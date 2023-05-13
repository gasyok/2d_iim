#pragma once
#include <preprocess/preprocess.h>
#include <unordered_map>

using namespace Eigen;

class System : public PreProcess {
private:
    // vector<double> coord_x, coord_y;
    vector<vector<Vector3d>> u;
    vector<vector<Vector3d>> u_next;
    vector<double> boundary;
public:
    System(InitValues init);
    Vector3d GetValue(int i, int j);
    Vector3d equation(int i, int j);
    Vector3d irrEquation(int i, int j);
    void solve(double t);
    void shift();
};
