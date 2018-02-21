#ifndef mesh_H
#define mesh_H

#include <igl/readOFF.h>

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Eigen/Dense>
using namespace Eigen;
using namespace std;


class mesh {
    public:
    tuple<Matrix3d, Vector3d, double> ICP(MatrixXd Pv, MatrixXd Qv, int step_size);
    MatrixXd Add_noise(MatrixXd m, double noise_val);
    pair<MatrixXd, MatrixXd> rotate(MatrixXd, double x, double y, double z);
};
#endif

