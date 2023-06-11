#include <system/system.h>
#include <fstream>

Vector3d System::GetValue(int i, int j) {
    return u[i][j];
}
Vector3d System::equation(int i, int j) {
    int i1, i2, j1, j2;
    double K = tau / h;

    double c = c_minus;
    Matrix3d A = A_minus;
    Matrix3d B = B_minus;

    if (func(h * i, h * j) > 0.0) {
        A = A_plus;
        B = B_plus;
        c = c_plus;
    }
    
    i1 = (Mx + i + 1) % Mx;
    i2 = (Mx + i - 1) % Mx;
    j1 = (My + j + 1) % My;
    j2 = (My + j - 1) % My;

    return u[i][j] - 0.5 * tau / h * (A * (u[i1][j] - u[i2][j]) + B * (u[i][j1] - u[i][j2])) + \
            0.5 * (tau / h) * (tau / h) * c * c * (u[i1][j] + u[i2][j] + u[i][j1] + u[i][j2] - 4 * u[i][j]);
}
Vector3d System::irrEquation(int i, int j) {
    Vector3d res (0.0, 0.0, 0.0);

    int i1 = (Mx + (i - 1)) % Mx;
    int i2 = (Mx + (i + 1)) % Mx;

    int j1 = (My + (j - 1)) % My;
    int j2 = (My + (j + 1)) % My;

    const vector<pair<int, int>> offsets = {{i1, j}, {i, j}, {i2, j}, {i, j2}, {i, j1}, {i2, j1}};

    pair<int, int> point = std::make_pair(i, j);
    for (int l = 0; l < 6; ++l) {
        int new_i = offsets[l].first;
        int new_j = offsets[l].second;

        res += gamma_matrices[point][l] * u[new_i][new_j];
    }
    return u[i][j] + (tau / h) * res; 
}
void System::solve(double t) {
    vector<vector<Vector3d>> new_u = u;
    for (auto c : irregular_points) {
        new_u[c.first][c.second] = irrEquation(c.first, c.second);
    }
    for (auto c : regular_points) {
        new_u[c.first][c.second] = equation(c.first, c.second);
    }
    u = new_u;
}
System::System(double _tau, double _h, int _Mx, int _My, double _x0, double _y0, double _A, double _omega, double _alpha)
    : PreProcess(_tau, _h, _Mx, _My, _x0, _y0, _A, _omega, _alpha) {
}
