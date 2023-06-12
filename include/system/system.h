#pragma once
#include <preprocess/preprocess.h>
#include <unordered_map>

using namespace Eigen;

class System : public PreProcess {
private:
    double total_time;
    int total_steps;
public:
    System(int _M, double _x0, double _y0, double _A, double _omega, double _alpha);
    Vector3d GetValue(int i, int j);
    Vector3d equation(int i, int j);
    Vector3d irrEquation(int i, int j);
    void solve(double t);
    void sample();
};
