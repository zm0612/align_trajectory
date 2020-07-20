#include "compute_sim3.h"
#include <pangolin/pangolin.h>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    pangolin::CreateWindowAndBind("Main", 1024, 768);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    pangolin::OpenGlRenderState s_cam(
        pangolin::ProjectionMatrix(1024, 768, 2000, 2000, 512, 384 , 0.1, 1e6),
        pangolin::ModelViewLookAt(0, 400, 400, 0, 300, 300, pangolin::AxisY)
    );

    pangolin::Handler3D handler(s_cam);
    pangolin::View& d_cam = pangolin::CreateDisplay()
        .SetBounds(0.0, 1.0, 0.0, 1.0, -1024/768.0f)
        .SetHandler(&handler);

    pangolin::CreatePanel("ui").SetBounds(0.0, 1.0, 0.0, pangolin::Attach::Pix(300));

    ComputeSim3<double> computeSim3(0.01);

    std::string strTraj_1 = "./data/frame_trajectory.txt";
    std::string strTraj_2 = "./data/lidar-trajectory.txt";
    std::string savedTrajFile = "./";

    computeSim3.LoadTraj(strTraj_1, strTraj_2);

    std::vector<TrajPoint<double>> syncedTraj_1;
    std::vector<TrajPoint<double>> syncedTraj_2;
    computeSim3.SyncTraj(syncedTraj_1, syncedTraj_2);

    Eigen::Matrix4d sim3 = computeSim3.GetSim3();

    std::vector<Eigen::Vector3d> points_1;
    std::vector<Eigen::Vector3d> points_2;

    computeSim3.GetSyncedTrajPoints(points_1, points_2);

    pangolin::Var<bool> original_button("ui.Origin", false, true);
    pangolin::Var<bool> trans_button("ui.Trans", false, true);
    pangolin::Var<int> match_int("ui.Match", -1, 0, points_1.size()-1);
    pangolin::Var<double> error("ui.Error", 0);
    while (!pangolin::ShouldQuit()){
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        d_cam.Activate(s_cam);


        if (pangolin::Pushed(trans_button)){
            glColor3f(0.0, 1.0, 0.0);
            glPointSize(3.0);
            glBegin(GL_POINTS);
            for (int i = 0; i < points_1.size(); ++i) {
                Eigen::Vector3d points;
                points = sim3.block(0,0,3,3)*points_1[i] + sim3.block(0,3,3,1);
                glVertex3d(points(0), points(1), points(2));
            }

            glColor3f(1.0, 1.0, 1.0);
            for (int i = 0; i < points_2.size(); ++i) {
                glVertex3d(points_2[i](0), points_2[i](1), points_2[i](2));
            }
            glEnd();

            original_button = false;
            trans_button = true;
            match_int = -1;
        }

        if (pangolin::Pushed(original_button)){
            glColor3f(0.0, 1.0, 0.0);
            glPointSize(3.0);
            glBegin(GL_POINTS);
            for (int i = 0; i < points_1.size(); ++i) {
                glVertex3d(points_1[i](0), points_1[i](1), points_1[i](2));
            }

            glColor3f(1.0, 1.0, 1.0);
            for (int i = 0; i < points_2.size(); ++i) {
                glVertex3d(points_2[i](0), points_2[i](1), points_2[i](2));
            }
            glEnd();

            original_button = true;
            trans_button = false;
            match_int = -1;
        }

        if (match_int >= 0){
            glPointSize(1.0);
            glBegin(GL_POINTS);
            for (int i = 0; i < points_1.size(); ++i) {
                glColor3f(1.0, 1.0, 1.0);
                glVertex3d(points_1[i](0), points_1[i](1), points_1[i](2));
            }

            for (int i = 0; i < points_2.size(); ++i) {
                glVertex3d(points_2[i](0), points_2[i](1), points_2[i](2));
            }
            glEnd();

            glPointSize(10);
            glColor3f(1.0, 0.0, 0.0);
            glBegin(GL_POINTS);
            Eigen::Vector3d points;
            points = sim3.block(0,0,3,3)*points_1[match_int] + sim3.block(0,3,3,1);
            glVertex3d(points(0), points(1), points(2));
            glEnd();

            glPointSize(10);
            glColor3f(0.0, 1.0, 0.0);
            glBegin(GL_POINTS);
            glVertex3d(points_2[match_int](0), points_2[match_int](1), points_2[match_int](2));
            glEnd();

            error = (points-points_2[match_int]).norm();
        }

        pangolin::FinishFrame();
    }

    return 0;
}
