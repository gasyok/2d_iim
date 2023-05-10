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
    vector<double> coord_x, coord_y;
    vector<vector<double>> pressure, velocity_x, velocity_y;
    double tau, h, c, rho, min_x, max_x, min_y, max_y;
    int size_x, size_y;
public:
    Matrix3d A_minus, A_plus, B_minus, B_plus;
    double rho_minus, rho_plus, c_minus, c_plus, k_minus, k_plus;
    InitValues(double tau, double h, double c, double rho, double min_x,
               double max_x, double min_y, double max_y, double x0, double y0, double A, double sigma, double a, double b);
    InitValues();
    void SetPressure(double x0, double y0, double A, double sigma, double a, double b);
    void SetVelocity(double x0, double y0, double a, double b);

    double GetTau();
    double GetH();
    double GetRho();
    double GetSpeed();
    int GetSizeX();
    int GetSizeY();

    vector<double> GetCoordX();
    vector<double> GetCoordY();
    vector<vector<double>> GetVelocityX();
    vector<vector<double>> GetVelocityY();
    vector<vector<double>> GetPressure();
};
