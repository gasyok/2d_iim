#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
using namespace Eigen;
using std::vector;

class InitValues {
private:
    // vector<double> coord_x, coord_y, velocity_x, velocity_y;
    double A, omega, alpha, x0, y0;
public:
    vector<vector<Vector3d>> u;
    Matrix3d A_minus, A_plus, B_minus, B_plus;
    double rho_minus, rho_plus, c_minus, c_plus, k_minus, k_plus;
    double tau, h;
    int Mx, My;

    InitValues();
    InitValues(double _tau, double _h, int _Mx, int _My, double _x0, double _y0, double _A, double _omega, double _alpha);
    void SetInitU(double x0, double y0, double A, double omega, double alpha);
    void SetInitRadU(double x0, double y0, double omega);
    void PrintInit();
};
