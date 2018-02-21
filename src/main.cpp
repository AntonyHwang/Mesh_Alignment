#include "acq/normalEstimation.h"
#include "acq/decoratedCloud.h"
#include "acq/cloudManager.h"

#include "nanogui/formhelper.h"
#include "nanogui/screen.h"

#include "igl/readOFF.h"
#include "igl/viewer/Viewer.h"
#include "mesh.h"
#include "Eigen/Dense"

#include <iostream>
#include <string>
#include <cmath>

namespace acq {

/** \brief                      Re-estimate normals of cloud \p V fitting planes
 *                              to the \p kNeighbours nearest neighbours of each point.
 * \param[in ] kNeighbours      How many neighbours to use (Typiclaly: 5..15)
 * \param[in ] vertices         Input pointcloud. Nx3, where N is the number of points.
 * \param[in ] maxNeighbourDist Maximum distance between vertex and neighbour.
 * \param[out] viewer           The viewer to show the normals at.
 * \return                      The estimated normals, Nx3.
 */
    NormalsT
    recalcNormals(
            int                 const  kNeighbours,
            CloudT              const& vertices,
            float               const  maxNeighbourDist
    ) {
        NeighboursT const neighbours =
                calculateCloudNeighbours(
                        /* [in]        cloud: */ vertices,
                        /* [in] k-neighbours: */ kNeighbours,
                        /* [in]      maxDist: */ maxNeighbourDist
                );

        // Estimate normals for points in cloud vertices
        NormalsT normals =
                calculateCloudNormals(
                        /* [in]               Cloud: */ vertices,
                        /* [in] Lists of neighbours: */ neighbours
                );

        return normals;
    } //...recalcNormals()

    void setViewerNormals(
            igl::viewer::Viewer      & viewer,
            CloudT              const& vertices,
            NormalsT            const& normals
    ) {
        // [Optional] Set viewer face normals for shading
        //viewer.data.set_normals(normals);

        // Clear visualized lines (see Viewer.clear())
        viewer.data.lines = Eigen::MatrixXd(0, 9);

        // Add normals to viewer
        viewer.data.add_edges(
                /* [in] Edge starting points: */ vertices,
                /* [in]       Edge endpoints: */ vertices + normals * 0.01, // scale normals to 1% length
                /* [in]               Colors: */ Eigen::Vector3d::Zero()
        );
    }

} //...ns acq

int main(int argc, char *argv[]) {

    mesh msh;

    MatrixXd V_1, V_2, V_3, V_4, V_5;
    MatrixXi F_1, F_2, F_3, F_4, F_5;
    //file path
    std::string m1Path = "../off_files/bun000.off";
    std::string m2Path = "../off_files/bun045.off";
    std::string m3Path = "../off_files/bun090.off";
    std::string m4Path = "../off_files/bun180.off";
    std::string m5Path = "../off_files/top2.off";

    igl::readOFF(m1Path, V_1, F_1);
    igl::readOFF(m2Path, V_2, F_2);
    igl::readOFF(m3Path, V_3, F_3);
    igl::readOFF(m4Path, V_4, F_4);
    igl::readOFF(m5Path, V_5, F_5);

    // How many neighbours to use for normal estimation, shown on GUI.
    int kNeighbours = 10;
    // Maximum distance between vertices to be considered neighbours (FLANN mode)
    float maxNeighbourDist = 0.15; //TODO: set to average vertex distance upon read

    double rot_x = 0;
    double rot_y = 0;
    double rot_z = 0;
    double noise_val = 0.0005;
    int max_iteration = 250;
    int step_size = 1;

    Eigen::Vector3d T;


    // Visualize the mesh in a viewer
    igl::viewer::Viewer viewer;
    {
        // Don't show face edges
        viewer.core.show_lines = false;
    }

    // Store cloud so we can store normals later
    acq::CloudManager cloudManager;
    // Read mesh from meshPath
    {
        int total_V = V_1.rows() + V_2.rows();
        int total_F = F_1.rows() + F_2.rows();
        int dim = 3;
        cout << V_2.rows();

        Eigen::MatrixXd disp_V;
        disp_V.resize(total_V, dim);
        disp_V << V_1, V_2;

        Eigen::MatrixXi disp_F;
        disp_F.resize(total_F,F_1.cols());
        disp_F << F_1,(F_2.array() + V_1.rows());

        Eigen::MatrixXd Color(disp_F.rows(),dim);

        RowVector3d m1_color(0, 0, 1);
        RowVector3d m2_color(1, 0, 0);
        Color << m1_color.replicate(F_1.rows(),1),
                 m2_color.replicate(F_2.rows(),1);


        // Store read vertices and faces
        cloudManager.addCloud(acq::DecoratedCloud(disp_V, disp_F));

        cloudManager.addCloud(acq::DecoratedCloud(V_1, F_1));
        cloudManager.setCloud(acq::DecoratedCloud(V_1, F_1),1);

        cloudManager.addCloud(acq::DecoratedCloud(V_2, F_2));
        cloudManager.setCloud(acq::DecoratedCloud(V_2, F_2),2);

        cloudManager.addCloud(acq::DecoratedCloud(V_3, F_3));
        cloudManager.setCloud(acq::DecoratedCloud(V_3, F_3),3);

        cloudManager.addCloud(acq::DecoratedCloud(V_4, F_4));
        cloudManager.setCloud(acq::DecoratedCloud(V_4, F_4),4);

        cloudManager.addCloud(acq::DecoratedCloud(V_5, F_5));
        cloudManager.setCloud(acq::DecoratedCloud(V_5, F_5),5);

        cloudManager.addCloud(acq::DecoratedCloud(V_1, F_1));
        cloudManager.setCloud(acq::DecoratedCloud(V_1, F_1),6);

        cloudManager.addCloud(acq::DecoratedCloud(V_2, F_2));
        cloudManager.setCloud(acq::DecoratedCloud(V_2, F_2),7);
        viewer.data.clear();

        // Show mesh
        viewer.data.set_mesh(
                cloudManager.getCloud(0).getVertices(),
                cloudManager.getCloud(0).getFaces()
        );
        viewer.data.set_colors(Color);
//
//        // Calculate normals on launch
//        cloudManager.getCloud(0).setNormals(
//                acq::recalcNormals(
//                        /* [in]      K-neighbours for FLANN: */ kNeighbours,
//                        /* [in]             Vertices matrix: */ cloudManager.getCloud(0).getVertices(),
//                        /* [in]      max neighbour distance: */ maxNeighbourDist
//                )
//        );
//        // Update viewer
//        acq::setViewerNormals(
//                viewer,
//                cloudManager.getCloud(0).getVertices(),
//                cloudManager.getCloud(0).getNormals()
//        );
    } //...read mesh







    // Extend viewer menu using a lambda function
    viewer.callback_init =
            [
                    &cloudManager, &kNeighbours, &maxNeighbourDist, &V_1, &V_2, &V_3, &V_4, &V_5, &F_1, &F_2, &F_3, &F_4, &F_5, &step_size, &max_iteration, &noise_val, &msh, &rot_x, &rot_y, &rot_z
            ] (igl::viewer::Viewer& viewer)
            {
                // Add an additional menu window
                viewer.ngui->addWindow(Eigen::Vector2i(900,10), "ICP Control");

                // Add new group
                viewer.ngui->addGroup("ICP:");

                // Add k-neighbours variable to GUI
                viewer.ngui->addVariable<int>(
                        /* Displayed name: */ "Step Size",

                        /*  Setter lambda: */ [&] (int val) { step_size = val; },

                        /*  Getter lambda: */ [&]() { return step_size; }
                );
                viewer.ngui->addVariable<int>(
                        /* Displayed name: */ "Num Iter",

                        /*  Setter lambda: */ [&] (int val) { max_iteration = val; },

                        /*  Getter lambda: */ [&]() { return max_iteration; }
                );
                // Add an additional menu window

                //ICP Button

                viewer.ngui->addButton(
                        /* displayed label: */  "Point To Point ICP",
                        [&](){
                            Eigen::MatrixXd R, Pv, Qv;
                            Eigen::MatrixXi Pf, Qf;
                            Eigen::Vector3d t;
                            double d_diff = 1000;
                            int count = 0;

                            //Get face and vertex matrix
                            Pv = cloudManager.getCloud(6).getVertices();
                            Pf = cloudManager.getCloud(6).getFaces();

                            Qv = cloudManager.getCloud(7).getVertices();
                            Qf = cloudManager.getCloud(7).getFaces();
                            //ICP algorithm
                            //Call function in class
                            int dim = 3;
                            double pre_diff = 1000;
                            double t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if ((d_diff > pre_diff && d_diff < 0.0001) || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            double time = (std::clock() - t_start)*1.0/CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;
                            //store mesh
                            cloudManager.setCloud(acq::DecoratedCloud(Pv, Pf),6);
                            cloudManager.setCloud(acq::DecoratedCloud(Qv, Qf),7);
                            MatrixXd result_V;
                            result_V.resize(Pv.rows() + Qv.rows(), dim);
                            MatrixXi result_F(Pf.rows() + Qf.rows(), Pf.cols());
                            result_V << Qv, Pv;
                            result_F << Qf, (Pf.array() + Qv.rows());
                            //set mesh color
                            RowVector3d m1_color(0, 0, 1);
                            RowVector3d m2_color(1, 0, 0);

                            Eigen::MatrixXd Color(result_F.rows(), dim);

                            Color << m1_color.replicate(F_1.rows(),1),
                                     m2_color.replicate(F_2.rows(),1);

                            cloudManager.setCloud(acq::DecoratedCloud(result_V, result_F),0);
                            viewer.data.clear();
                            viewer.data.set_mesh(result_V, result_F);
                            viewer.data.set_colors(Color);
                        }
                );

                viewer.ngui->addVariable<double>(
                        /* Displayed name: */ "Noise (0-0.01):",

                        /*  Setter lambda: */ [&] (double val) { noise_val = val; },

                        /*  Getter lambda: */ [&]() { return noise_val; }
                );
                viewer.ngui->addButton(
                        "Add Noise to M2",
                        [&](){
                            MatrixXd n_V_2;
                            MatrixXi n_F_2;

                            int total_V = V_1.rows() + V_2.rows();
                            int total_F = F_1.rows() + F_2.rows();

                            cloudManager.setCloud(acq::DecoratedCloud(V_2, F_2),2);
                            //call class function to add noise
                            n_V_2 = msh.Add_noise(cloudManager.getCloud(7).getVertices(), noise_val);
                            n_F_2 = cloudManager.getCloud(7).getFaces();

                            Eigen::MatrixXd disp_V;
                            disp_V.resize(total_V, 3);
                            disp_V << V_1, n_V_2;

                            Eigen::MatrixXi disp_F;
                            disp_F.resize(total_F,F_1.cols());
                            disp_F << F_1,(n_F_2.array() + V_1.rows());
                            //store current scene mesh
                            cloudManager.setCloud(acq::DecoratedCloud(V_1, F_1),6);
                            cloudManager.setCloud(acq::DecoratedCloud(n_V_2, n_F_2),7);

                            viewer.data.clear();

                            viewer.data.set_mesh(disp_V, disp_F);

                            Eigen::MatrixXd Color(disp_F.rows(), 3);

                            RowVector3d m1_color(0, 0, 1);
                            RowVector3d m2_color(1, 0, 0);

                            Color << m1_color.replicate(F_1.rows(),1),
                                     m2_color.replicate(n_F_2.rows(),1);

                            cloudManager.setCloud(acq::DecoratedCloud(disp_V, disp_F),0);
                            viewer.data.set_colors(Color);

                        }
                );
                viewer.ngui->addVariable<double>(
                        /* Displayed name: */ "X-rotation:",

                        /*  Setter lambda: */ [&] (double r_x) { rot_x = r_x; },

                        /*  Getter lambda: */ [&]() { return rot_x; }
                );
                viewer.ngui->addVariable<double>(
                        /* Displayed name: */ "Y-rotation:",

                        /*  Setter lambda: */ [&] (double r_y) { rot_y = r_y; },

                        /*  Getter lambda: */ [&]() { return rot_z; }
                );
                viewer.ngui->addVariable<double>(
                        /* Displayed name: */ "Z-rotation:",

                        /*  Setter lambda: */ [&] (double r_z) { rot_z = r_z; },

                        /*  Getter lambda: */ [&]() { return rot_z; }
                );
                viewer.ngui->addButton(
                        "Get Rotated M1",
                        [&](){
                            MatrixXd m1_V = V_1;
                            MatrixXd m2_V(m1_V.rows(), m1_V.cols());

                            int dim = 3;

                            double deg_to_rad = M_PI / 180;
                            double x = rot_x * deg_to_rad;
                            double y = rot_y * deg_to_rad;
                            double z = rot_z * deg_to_rad;

                            tie(m1_V, m2_V) = msh.rotate(m1_V, x, y, z);

                            MatrixXd result_V(m1_V.rows() * 2, dim);
                            MatrixXi result_F(F_1.rows() * 2, F_1.cols());
                            result_V << m1_V, m2_V;
                            result_F << F_1, (F_1.array() + V_1.rows());

                            RowVector3d m1_color(0, 0, 1);
                            RowVector3d m2_color(1, 0, 0);


                            Eigen::MatrixXd Color(result_F.rows(), dim);


                            Color << m1_color.replicate(F_1.rows(),1),
                                     m2_color.replicate(F_2.rows(),1);

                            cloudManager.setCloud(acq::DecoratedCloud(m1_V, F_1),6);
                            cloudManager.setCloud(acq::DecoratedCloud(m2_V, F_1),7);
                            cloudManager.setCloud(acq::DecoratedCloud(result_V, result_F),0);
                            viewer.data.clear();
                            viewer.data.set_mesh(result_V,result_F);
                            viewer.data.set_colors(Color);

                        }
                );
                viewer.ngui->addButton(
                        "Multi-Scan",
                        [&](){
                            //ICP from 0 degree to 90 degree
                            //from 90 to 180
                            //from 180 to 270
                            //from 315 to 270
                            cloudManager.setCloud(acq::DecoratedCloud(V_1, F_1),1);
                            cloudManager.setCloud(acq::DecoratedCloud(V_2, F_2),2);
                            cloudManager.setCloud(acq::DecoratedCloud(V_3, F_3),3);
                            cloudManager.setCloud(acq::DecoratedCloud(V_4, F_4),4);
                            cloudManager.setCloud(acq::DecoratedCloud(V_5, F_5),5);
                            viewer.data.clear();
                            int total_V = V_1.rows() + V_2.rows() + V_3.rows() + V_4.rows() + V_5.rows();
                            int total_F = F_1.rows() + F_2.rows() + F_3.rows() + F_4.rows() + F_5.rows();
                            int dim = 3;

                            Eigen::MatrixXd R, Pv, Qv;
                            Eigen::MatrixXi Pf, Qf;
                            Eigen::Vector3d t;
                            double d_diff = 1000;
                            double pre_diff = 1000;
                            double t_start;
                            double time;
                            int count = 0;

                            //ICP mesh 1-2
                            Pv = cloudManager.getCloud(2).getVertices();
                            Pf = cloudManager.getCloud(2).getFaces();

                            Qv = cloudManager.getCloud(1).getVertices();
                            Qf = cloudManager.getCloud(1).getFaces();

                            t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if (d_diff > pre_diff || d_diff < 0.000018 || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            time = (clock() - t_start) * 1.0 / CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;

                            MatrixXd V1_result;
                            MatrixXi F1_result;
                            V1_result = Qv;
                            F1_result = F_1;


                            //ICP mesh 2-3
                            Pv = cloudManager.getCloud(3).getVertices();
                            Pf = cloudManager.getCloud(3).getFaces();

                            Qv = cloudManager.getCloud(2).getVertices();
                            Qf = cloudManager.getCloud(2).getFaces();

                            d_diff = 1000;
                            pre_diff = 1000;
                            count = 0;

                            t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if (d_diff > pre_diff || d_diff < 0.000018 || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            time = (clock() - t_start) * 1.0 / CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;

                            MatrixXd V2_result;
                            MatrixXi F2_result;
                            V2_result = Qv;
                            F2_result = F_2;

                            //ICP mesh 123-4
                            Pv = cloudManager.getCloud(3).getVertices();
                            Pf = cloudManager.getCloud(3).getFaces();

                            Qv = cloudManager.getCloud(4).getVertices();
                            Qf = cloudManager.getCloud(4).getFaces();

                            d_diff = 1000;
                            pre_diff = 1000;
                            count = 0;

                            t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if (d_diff > pre_diff || d_diff < 0.000018 || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            time = (clock() - t_start) * 1.0 / CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;

                            MatrixXd V4_result;
                            MatrixXi F4_result;
                            V4_result = Qv;
                            F4_result = F_4;

                            //ICP mesh 5-4
                            Pv = cloudManager.getCloud(4).getVertices();
                            Pf = cloudManager.getCloud(4).getFaces();

                            Qv = cloudManager.getCloud(5).getVertices();
                            Qf = cloudManager.getCloud(5).getFaces();

                            d_diff = 1000;
                            pre_diff = 1000;
                            count = 0;

                            t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if (d_diff > pre_diff || d_diff < 0.000018 || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            time = (clock() - t_start) * 1.0 / CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;

                            MatrixXd V5_result;
                            MatrixXi F5_result;
                            V5_result = Qv;
                            F5_result = F_5;

                            MatrixXd result_V(total_V, dim);
                            MatrixXi result_F(total_F, F_1.cols());
                            result_V << V1_result, V2_result, V_3, V4_result, V5_result;

                            result_F <<  F_2,
                                        (F_2.array() + V_1.rows()),
                                        (F_3.array() + V_1.rows() + V_2.rows()),
                                        (F_4.array() + V_1.rows() + V_2.rows() + V_3.rows()),
                                        (F_5.array() + V_1.rows() + V_2.rows() + V_3.rows() + V_4.rows());

                            RowVector3d m1_color(0, 0, 1);
                            RowVector3d m2_color(0, 1, 0);
                            RowVector3d m3_color(1, 0, 0);
                            RowVector3d m4_color(0.5, 0.5, 0.5);
                            RowVector3d m5_color(0, 0.5, 0.5);

                            Eigen::MatrixXd Color(result_F.rows(), dim);


                            Color << m1_color.replicate(F_1.rows(),1),
                                     m2_color.replicate(F_2.rows(),1),
                                     m3_color.replicate(F_3.rows(),1),
                                     m4_color.replicate(F_4.rows(),1),
                                     m5_color.replicate(F_5.rows(),1);

                            viewer.data.clear();
                            cloudManager.setCloud(acq::DecoratedCloud(result_V, result_F),0);

                            viewer.data.set_mesh(result_V,result_F);
                            viewer.data.set_colors(Color);

                        }
                );
                viewer.ngui->addButton(
                        "Multi-Scan 2",
                        [&](){
                            cloudManager.setCloud(acq::DecoratedCloud(V_1, F_1),1);
                            cloudManager.setCloud(acq::DecoratedCloud(V_2, F_2),2);
                            cloudManager.setCloud(acq::DecoratedCloud(V_3, F_3),3);
                            cloudManager.setCloud(acq::DecoratedCloud(V_4, F_4),4);
                            cloudManager.setCloud(acq::DecoratedCloud(V_5, F_5),5);
                            viewer.data.clear();
                            int total_V = V_1.rows() + V_2.rows() + V_3.rows() + V_4.rows() + V_5.rows();
                            int total_F = F_1.rows() + F_2.rows() + F_3.rows() + F_4.rows() + F_5.rows();
                            int dim = 3;

                            Eigen::MatrixXd R, Pv, Qv;
                            Eigen::MatrixXi Pf, Qf;
                            Eigen::Vector3d t;
                            double d_diff = 1000;
                            double pre_diff = 1000;
                            double t_start;
                            double time;
                            int count = 0;
                            MatrixXd no;

                            //ICP mesh 1-2
                            tie(no, Pv) = msh.rotate(cloudManager.getCloud(3).getVertices(), 0, 0, 0);
                            Pf = cloudManager.getCloud(3).getFaces();

                            tie(no, Qv) = msh.rotate(cloudManager.getCloud(1).getVertices(), 0, -90, 0);
                            Qf = cloudManager.getCloud(1).getFaces();

                            t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if ((d_diff > pre_diff && d_diff < 0.0001) || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            time = (clock() - t_start) * 1.0 / CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;

                            MatrixXd result_V12(Qv.rows() + Qv.rows(), dim);
                            MatrixXi result_F12(Pf.rows() + Qf.rows(), Pf.cols());
                            result_V12 << Qv, Pv;
                            result_F12 << Qf, (Pf.array() + Qv.rows());

                            //ICP mesh 3-12

                            Pv = result_V12;
                            Pf = result_F12;

                            tie(no, Qv) = msh.rotate(cloudManager.getCloud(2).getVertices(), 0, -45, 0);
                            Qf = cloudManager.getCloud(2).getFaces();

                            d_diff = 1000;
                            pre_diff = 1000;
                            count = 0;

                            t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if ((d_diff > pre_diff && d_diff < 0.0001) || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            time = (clock() - t_start) * 1.0 / CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;

                            MatrixXd result_V32;
                            result_V32.resize(Pv.rows() + Qv.rows(), dim);
                            MatrixXi result_F32(Pf.rows() + Qf.rows(), Pf.cols());
                            result_V32 << Qv, Pv;
                            result_F32 << Qf, (Pf.array() + Qv.rows());


                            //ICP mesh 5-4
                            Pv = result_V32;
                            Pf = result_F32;

                            tie(no, Qv) = msh.rotate(cloudManager.getCloud(4).getVertices(), 0, 90, 0);
                            Qf = cloudManager.getCloud(4).getFaces();

                            d_diff = 1000;
                            pre_diff = 1000;
                            count = 0;

                            t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if ((d_diff > pre_diff && d_diff < 0.0001) || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            time = (clock() - t_start) * 1.0 / CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;

                            MatrixXd result_V54;
                            result_V54.resize(Pv.rows() + Qv.rows(), dim);
                            MatrixXi result_F54(Pf.rows() + Qf.rows(), Pf.cols());
                            result_V54 << Qv, Pv;
                            result_F54 << Qf, (Pf.array() + Qv.rows());

                            //ICP mesh 4_123
                            Pv = result_V54;
                            Pf = result_F54;

                            tie(no, Qv) = msh.rotate(cloudManager.getCloud(5).getVertices(), -90, 0, 0);
                            Qf = cloudManager.getCloud(5).getFaces();

                            d_diff = 1000;
                            pre_diff = 1000;
                            count = 0;

                            t_start = clock();
                            for (int i = 0; i < max_iteration; i++){
                                count++;
                                tie(R, t, d_diff) = msh.ICP(Pv, Qv, step_size);
                                d_diff = abs(d_diff);
                                Qv = (R * Qv.transpose()).transpose() + t.replicate(1, Qv.rows()).transpose();
                                if ((d_diff > pre_diff && d_diff < 0.0001) || i > max_iteration) {
                                    cout << "\nEnd distance: " << d_diff << "\n";
                                    break;
                                }
                                pre_diff = d_diff;
                            }
                            cout<< "Iterations: "<<count << "\n";
                            time = (clock() - t_start) * 1.0 / CLOCKS_PER_SEC;
                            cout << "Processing Time: " << time << " s" << endl;

                            MatrixXd result_V;
                            result_V.resize(Pv.rows() + Qv.rows(), dim);
                            MatrixXi result_F(Pf.rows() + Qf.rows(), Pf.cols());
                            result_V << Qv, Pv;
                            result_F << Qf, (Pf.array() + Qv.rows());

//                            MatrixXd result_V(total_V, dim);
//                            MatrixXi result_F(total_F, F_1.cols());
//                            result_V << V1_result, V2_result, V_3, V4_result, V5_result;
//                              result_F << F_1, (F_2.array() + V_1.rows()), (F_3.array() + V_2.rows() + V_1.rows()), (F_4.array() + V_3.rows() + V_2.rows() + V_1.rows()), (F_5.array() + V_4.rows() + V_3.rows() + V_2.rows() + V_1.rows());
//                            result_F <<  F_2,
//                                    (F_2.array() + V_1.rows()),
//                                    (F_3.array() + V_1.rows() + V_2.rows()),
//                                    (F_4.array() + V_1.rows() + V_2.rows() + V_3.rows()),
//                                    (F_5.array() + V_1.rows() + V_2.rows() + V_3.rows() + V_4.rows());

                            RowVector3d m1_color(0, 0, 1);
                            RowVector3d m2_color(0, 1, 0);
                            RowVector3d m3_color(1, 0, 0);
                            RowVector3d m4_color(0.5, 0.5, 0.5);
                            RowVector3d m5_color(0, 0.5, 0.5);
//
//
                            Eigen::MatrixXd Color(result_F.rows(), dim);
//
//
                            Color << m1_color.replicate(F_1.rows(),1),
                                    m2_color.replicate(F_2.rows(),1),
                                    m3_color.replicate(F_3.rows(),1),
                                    m4_color.replicate(F_4.rows(),1),
                                    m5_color.replicate(F_5.rows(),1);

                            viewer.data.clear();
                            cloudManager.setCloud(acq::DecoratedCloud(result_V, result_F),0);

                            viewer.data.set_mesh(result_V,result_F);
                            viewer.data.set_colors(Color);

                        }
                );
                viewer.ngui->addButton(
                        "Show M1 and M2",
                        [&]() {

                            int total_V = V_1.rows() + V_2.rows();
                            int total_F = F_1.rows() + F_2.rows();
                            int dim = 3;

                            Eigen::MatrixXd disp_V;
                            disp_V.resize(total_V, dim);
                            disp_V << V_1, V_2;

                            Eigen::MatrixXi disp_F;
                            disp_F.resize(total_F,F_1.cols());
                            disp_F << F_1,(F_2.array() + V_1.rows());

                            Eigen::MatrixXd Color(disp_F.rows(),dim);

                            RowVector3d m1_color(0, 0, 1);
                            RowVector3d m2_color(1, 0, 0);
                            Color << m1_color.replicate(F_1.rows(),1),
                                     m2_color.replicate(F_2.rows(),1);

                            cloudManager.setCloud(acq::DecoratedCloud(disp_V, disp_F),0);
                            cloudManager.setCloud(acq::DecoratedCloud(V_1, F_1),6);
                            cloudManager.setCloud(acq::DecoratedCloud(V_2, F_2),7);
                            // Show mesh
                            viewer.data.clear();
                            viewer.data.set_mesh(
                                    cloudManager.getCloud(0).getVertices(),
                                    cloudManager.getCloud(0).getFaces()
                            );
                            viewer.data.set_colors(Color);
                        }
                );

                // Generate menu
                viewer.screen->performLayout();

                return false;
            }; //...viewer menu


    // Start viewer
    viewer.launch();

    return 0;
} //...main()

