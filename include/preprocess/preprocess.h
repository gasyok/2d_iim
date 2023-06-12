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
    void Solve();
public:
    std::unordered_map<std::pair<int, int>, std::vector<Eigen::Matrix3d>, pair_hash> gamma_matrices;
    PreProcess(int _M, double _x0, double _y0, double _A, double _omega, double _alpha);
};
