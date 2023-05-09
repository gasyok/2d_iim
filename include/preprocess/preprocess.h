#pragma once
#include <curve/curve.h>

class PreProcess : public Curve {
private:
    pair<int, int> GetPoint(int l, pair<int, int> point);
    Vector2d bisection(double x1, double x2, double y1, double y2, double tol, int max_iter);
    Vector2d GetOrigin(pair<int, int> point);
    Eigen::Vector2d normalize(const Eigen::Vector2d& vec);
    Vector2d GetRotationalCoord(int l, pair<int, int>point);
    Matrix3d RotateMatrix(pair<int, int> point, Matrix3d matrix);
    vector<Matrix3d> GetDefaultQ (pair<int, int> point);
    Matrix3d OppositeQ(int i, int l, pair<int, int> point);
    Matrix3d BesideQ(int i, int l, pair<int, int> point);
    Matrix3d GetQmatrix(int i, int l, pair<int, int> point);
    Matrix3d GetFmatrix(int i, pair<int, int> point);
    vector<Matrix3d> CalcGammaMatrices(pair<int, int> point);
public:
    std::unordered_map<std::pair<int, int>, std::vector<Eigen::Matrix3d>, pair_hash> gamma_matrices;
    Matrix3d A_minus, A_plus, B_minus, B_plus;
    double rho_minus, rho_plus, c_minus, c_plus, k_minus, k_plus;
    PreProcess(InitValues init);
    void Solve();
};
