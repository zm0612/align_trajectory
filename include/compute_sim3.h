//
// Created by meng on 2020/7/17.
//
#ifndef COMPUTE_SIM3_COMPUTE_SIM3_H
#define COMPUTE_SIM3_COMPUTE_SIM3_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <pangolin/pangolin.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>


template<typename Type>
struct TrajPoint{
    Type timestamp;
    Eigen::Matrix<Type, 3, 1> translation;
    Eigen::Matrix<Type, 4, 1> rotation;
};

template<typename Type>
class ComputeSim3 {
public:
    using _tTraj = std::vector<TrajPoint<Type>>;
    using _tTrajPoints = std::vector<Eigen::Matrix<Type, 3, 1>>;
    ComputeSim3(Type syncThreshold = 0.01):_syncThreshold(syncThreshold){};
    ~ComputeSim3(){};

    void LoadTraj(const std::string& strTraj1, const std::string& strTraj2);

    void SyncTraj(_tTraj& syncedTraj_1, _tTraj& syncedTraj_2);

    void GetTraj(_tTraj& traj_1, _tTraj& traj_2) const ;

    void GetSyncedTrajPoints(_tTrajPoints& points_1, _tTrajPoints& points_2) const;

    Eigen::Matrix<Type, 4, 4> GetSim3();

    void GetSyncedTraj(_tTraj& traj_1, _tTraj& traj_2);

    void SaveSyncedTraj(const std::string& savedPath);

private:
    _tTraj _vTraj_1;
    _tTraj _vTraj_2;
    _tTraj _vSyncedTraj_1;
    _tTraj _vSyncedTraj_2;
    Type _syncThreshold;
};

template<typename Type>
void ComputeSim3<Type>::GetSyncedTrajPoints(_tTrajPoints &points_1, _tTrajPoints &points_2) const {
    points_1.reserve(_vSyncedTraj_1.size());
    points_2.reserve(_vSyncedTraj_2.size());

    for (int i = 0; i < _vSyncedTraj_1.size(); ++i) {
        points_1.push_back(_vSyncedTraj_1[i].translation);
        points_2.push_back(_vSyncedTraj_2[i].translation);
    }
}

template<typename Type>
Eigen::Matrix<Type, 4,4> ComputeSim3<Type>::GetSim3() {
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> trajPoints_1;
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> trajPoints_2;

    const int n = _vSyncedTraj_1.size();
    trajPoints_1.resize(3, n);
    trajPoints_2.resize(3, n);

    for (int i = 0; i < _vSyncedTraj_1.size(); ++i) {
        trajPoints_1.block(0, i, 3, 1) = _vSyncedTraj_1[i].translation;
        trajPoints_2.block(0, i, 3, 1) = _vSyncedTraj_2[i].translation;
    }

    Eigen::Matrix<Type, 3, 1> means_1;
    Eigen::Matrix<Type, 3, 1> means_2;

    Type one_over_n = 1 / static_cast<Type>(n);
    means_1 = trajPoints_1.rowwise().sum() * one_over_n;
    means_2 = trajPoints_2.rowwise().sum() * one_over_n;

    Eigen::Matrix<Type, 3, Eigen::Dynamic> demeans_1;
    Eigen::Matrix<Type, 3, Eigen::Dynamic> demeans_2;
    demeans_1 = trajPoints_1.colwise() - means_1;
    demeans_2 = trajPoints_2.colwise() - means_2;

    Type var_1 = demeans_1.rowwise().squaredNorm().sum() * one_over_n;

    Eigen::Matrix<Type, 3, 3> sigma;
    sigma = demeans_2 * demeans_1.transpose() * one_over_n;

    Eigen::JacobiSVD<Eigen::Matrix<Type, 3, 3>> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Matrix<Type, 3, 1> S;
    S = Eigen::Matrix<Type, 3, 1>::Ones();

    if (svd.matrixU().determinant() * svd.matrixV().determinant() < 0)
        S(2) = -1;

    Eigen::Matrix<Type, 4,4> sim3 = Eigen::Matrix<Type, 4,4>::Identity();

    sim3.block(0,0,3,3).noalias() = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();

    Type c = 1 / var_1 * svd.singularValues().dot(S);

    sim3.block(0,3,3,1) = means_2 - c * sim3.block(0,0,3,3)*means_1;
    sim3.block(0,0,3,3) = sim3.block(0,0,3,3)*c;

    return sim3;
}

template<typename Type>
void ComputeSim3<Type>::SaveSyncedTraj(const std::string& savedPath) {
    std::ofstream f_1;
    std::ofstream f_2;

    f_1.open(savedPath+"/"+"traj_1.txt");
    f_2.open(savedPath+"/"+"traj_2.txt");
    f_1 << std::fixed;
    f_2 << std::fixed;

    if (_vSyncedTraj_1.empty() || _vSyncedTraj_2.empty()){
        std::cerr << "_vSyncedTraj_1.empty() || _vSyncedTraj_2.empty()" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (_vSyncedTraj_1.size() != _vSyncedTraj_2.size()){
        std::cerr << "_vSyncedTraj_1.size() != _vSyncedTraj_2.size()" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < _vSyncedTraj_1.size(); ++i) {
        f_1 << std::setprecision(6) << _vSyncedTraj_1[i].timestamp << " " << std::setprecision(7) << _vSyncedTraj_1[i].translation(0)
            << " " << _vSyncedTraj_1[i].translation(1) << " " << _vSyncedTraj_1[i].translation(2) << " " << _vSyncedTraj_1[i].rotation(0)
            << " " << _vSyncedTraj_1[i].rotation(1) << " " << _vSyncedTraj_1[i].rotation(2)<< " " << _vSyncedTraj_1[i].rotation(3) << std::endl;

        f_2 << std::setprecision(6) << _vSyncedTraj_2[i].timestamp << " " << std::setprecision(7) << _vSyncedTraj_2[i].translation(0)
            << " " << _vSyncedTraj_2[i].translation(1) << " " << _vSyncedTraj_2[i].translation(2) << " " << _vSyncedTraj_2[i].rotation(0)
            << " " << _vSyncedTraj_2[i].rotation(1) << " " << _vSyncedTraj_2[i].rotation(2)<< " " << _vSyncedTraj_2[i].rotation(3) << std::endl;
    }

    f_1.close();
    f_2.close();
}

template<typename Type>
void ComputeSim3<Type>::GetSyncedTraj(_tTraj &traj_1, _tTraj &traj_2) {
    traj_1 = _vSyncedTraj_1;
    traj_2 = _vSyncedTraj_2;
}

template<typename Type>
void ComputeSim3<Type>::SyncTraj(_tTraj &syncedTraj_1, _tTraj &syncedTraj_2) {
    int size_1 = _vTraj_1.size();
    int size_2 = _vTraj_2.size();
    int currentFirst = 0;

    if (size_1 > size_2){
        _vSyncedTraj_1.reserve(size_2);
        _vSyncedTraj_2.reserve(size_2);
        for (int i = 0; i < size_2; ++i) {
            for (int j = currentFirst; j < size_1; ++j) {

                if (std::abs(_vTraj_1[j].timestamp - _vTraj_2[i].timestamp) > _syncThreshold &&
                             _vTraj_1[j].timestamp > _vTraj_2[i].timestamp){
                    break;
                }


                if (std::abs(_vTraj_1[j].timestamp - _vTraj_2[i].timestamp) <= _syncThreshold){
                    _vSyncedTraj_1.push_back(_vTraj_1[j]);
                    _vSyncedTraj_2.push_back(_vTraj_2[i]);
                    currentFirst++;
                    break;
                }
            }
        }
    } else{// size_1 <= size_2
        for (int i = 0; i < size_1; ++i) {
            for (int j = currentFirst; j < size_2; ++j) {
                if (std::abs(_vTraj_1[i].timestamp - _vTraj_2[j].timestamp) > _syncThreshold &&
                    _vTraj_2[j].timestamp > _vTraj_1[i].timestamp){
                    break;
                }

                if (std::abs(_vTraj_1[i].timestamp - _vTraj_2[j].timestamp) <= _syncThreshold){
                    _vSyncedTraj_1.push_back(_vTraj_1[i]);
                    _vSyncedTraj_2.push_back(_vTraj_2[j]);
                    currentFirst++;
                    break;
                }
            }
        }

    }
    syncedTraj_1 = _vSyncedTraj_1;
    syncedTraj_2 = _vSyncedTraj_2;
}

template<typename Type>
void ComputeSim3<Type>::GetTraj(_tTraj &traj_1, _tTraj &traj_2) const {
    traj_1 = _vTraj_1;
    traj_2 = _vTraj_2;
}

template<typename Type>
void ComputeSim3<Type>::LoadTraj(const std::string &strTraj1, const std::string& strTraj2) {
    std::ifstream fTraj_1;
    std::ifstream fTraj_2;

    fTraj_1.open(strTraj1.c_str());
    if (!fTraj_1.is_open()){
        std::cerr << "Fail to open traj_1" << std::endl;
        exit(EXIT_FAILURE);
    }

    while (!fTraj_1.eof()){
        std::string s;
        std::getline(fTraj_1, s);
        if (!s.empty()){
            std::stringstream ss;
            ss.setf(std::ios::fixed);
            ss << s;
            TrajPoint<Type> trajPoint;
            ss >> trajPoint.timestamp;
            ss >> trajPoint.translation(0) >> trajPoint.translation(1) >> trajPoint.translation(2);
            ss >> trajPoint.rotation(0) >> trajPoint.rotation(1) >> trajPoint.rotation(2) >> trajPoint.rotation(3);
            _vTraj_1.emplace_back(trajPoint);
        }
    }

    fTraj_2.open(strTraj2.c_str());
    if (!fTraj_2.is_open()){
        std::cerr << "Fail to open traj_2" << std::endl;
        exit(EXIT_FAILURE);
    }

    while (!fTraj_2.eof()){
        std::string  s;
        std::getline(fTraj_2, s);
        if (!s.empty()){
            std::stringstream ss;
            ss.setf(std::ios::fixed);
            ss << s;
            TrajPoint<Type> trajPoint;
            ss >> trajPoint.timestamp;
            ss >> trajPoint.translation(0) >> trajPoint.translation(1) >> trajPoint.translation(2);
            ss >> trajPoint.rotation(0) >> trajPoint.rotation(1) >> trajPoint.rotation(2) >> trajPoint.rotation(3);
            _vTraj_2.emplace_back(trajPoint);
        }
    }
    fTraj_1.close();
    fTraj_2.close();
}
#endif //COMPUTE_SIM3_COMPUTE_SIM3_H
