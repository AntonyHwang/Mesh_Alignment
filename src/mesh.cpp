#include "mesh.h"
#include "new_nanoflann/nanoflann.hpp"
#include <random>

using namespace Eigen;
using namespace std;
using namespace nanoflann;


tuple<MatrixXd, MatrixXd, double> nearest_neighbor(MatrixXd Pv, MatrixXd Qv, int sample) {
    const size_t dim = 3;
    const size_t N_q = Qv.rows();
    const size_t N_p = Pv.rows();
    double dist = 0;

    Matrix<double, Dynamic, Dynamic>  mat(N_q, dim);

    for (size_t i = 0; i < N_q; i++)
        for (size_t d = 0; d < dim; d++)
            mat(i,d) = Qv(i,d);

//	cout << mat << endl;
    const size_t num_results = 1;
    typedef KDTreeEigenMatrixAdaptor<Matrix<double, Dynamic, Dynamic> > my_kd_tree_t;

    my_kd_tree_t mat_index(mat, 10 /* max leaf */ );
    mat_index.index->buildIndex();

    // Query point:
    vector<double> query_pt(dim);
    //matching results
    Matrix<double, Dynamic, Dynamic>  P(N_p, dim);
    Matrix<double, Dynamic, Dynamic>  Q(N_p, dim);
    int count = 0;
    for (size_t i = 0; i < Pv.rows(); i+=sample) {
        for (size_t d = 0; d < dim; d++)
            query_pt[d] = Pv(i, d);
        // do a knn search
        vector<size_t> ret_indexes(num_results);
        vector<double> out_dists_sqr(num_results);

        nanoflann::KNNResultSet<double> resultSet(num_results);

        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

//        std::cout << "knnSearch(nn=" << num_results << "): \n";
        if (out_dists_sqr[0] < 0.01) {
            int idx = ret_indexes[0];
            P.row(count) = Pv.row(i);
            Q.row(count) = Qv.row(idx);
            count ++;
            dist += out_dists_sqr[0];
        }
//        for (size_t i = 0; i < num_results; i++)
//            std::cout << "ret_index[" << i << "]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << endl;
    }
    P = P.topRows(count);
    dist /= P.rows();
    Q = Q.topRows(count);
    return {P, Q, dist};
}

tuple<Matrix3d, Vector3d, double> mesh::ICP(MatrixXd Pv, MatrixXd Qv, int step_size) {
    Vector3d t;
    MatrixXd P, Q, U, V, R, A;
    MatrixXd P_mean, Q_mean;
    double dist;
    tie(P, Q, dist) = nearest_neighbor(Pv, Qv, step_size);
    //compute the mean
    P_mean = P.colwise().sum()/P.rows();
    Q_mean = Q.colwise().sum()/Q.rows();
    //compute matrix A for SVD
    A = (Q - Q_mean.replicate(Q.rows(), 1)).transpose() * (P - P_mean.replicate(P.rows(), 1));
    //compute SVD
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    U = svd.matrixU();
    V = svd.matrixV();
    //compute rotation and translation
    R = V * U.transpose();
    t = P_mean.transpose() - (R * Q_mean.transpose());
    return {R, t, dist};
}
//function to add noise
MatrixXd mesh::Add_noise(MatrixXd m, double noise_val) {
    MatrixXd n_m = MatrixXd::Zero(m.rows(), m.cols());
    double noise = 0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    //distribution with zero mean
    normal_distribution<double> dist(0.0,noise_val);

    for (int i = 0; i < m.rows(); i++) {
        for (int d = 0; d < 3; d++) {
            noise = dist(generator);
            n_m(i, d) = m(i, d) + noise;
        }
    }
    return n_m;
}
//rotate the mesh vertices
pair<MatrixXd, MatrixXd> mesh::rotate(MatrixXd m1_V, double x, double y, double z) {
    Matrix3d r_x, r_y, r_z, R;
    MatrixXd V_mean, m2_V;
    //R matrix from axis angles
    r_x << 1, 0, 0,
            0, cos(x), -sin(x),
            0, sin(x), cos(x);
    r_y << cos(y), 0, sin(y),
            0, 1, 0,
            -sin(y), 0, cos(y);
    r_z << cos(z), -sin(z), 0,
            sin(z), cos(z), 0,
            0, 0, 1;

    R = r_z * r_y * r_x;

    V_mean = (m1_V.colwise().sum()/m1_V.rows()).replicate(m1_V.rows(),1);
    m1_V = m1_V - V_mean;
    m2_V = (R * m1_V.transpose()).transpose();

    return {m1_V, m2_V};
}
