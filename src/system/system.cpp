#include <system/system.h>
#include <fstream>

Vector3d System::GetValue(int i, int j) {
    return u[i][j];
}
Vector3d System::equation(int i, int j) {
    int i1, i2, j1, j2;
    double K = k / h;

    double c = c_minus;
    Matrix3d A = A_minus;
    Matrix3d B = B_minus;

    if (func(get_x(i), get_y(j)) >= 0.0) {
        A = A_plus;
        B = B_plus;
        c = c_plus;
    }
    i1 = (x_size + i + 1) % x_size;
    i2 = (x_size + i - 1) % x_size;
    j1 = (y_size + j + 1) % y_size;
    j2 = (y_size + j - 1) % y_size;

    return u_next[i][j] - 0.5 * k / h * (A * (u_next[i1][j] - u_next[i2][j]) + B * (u_next[i][j1] - u_next[i][j2])) + \
            0.5 * (k / h) * (k / h) * c * c * (u_next[i1][j] + u_next[i2][j] + u_next[i][j1] + u_next[i][j2] - 4 * u_next[i][j]);
    // return u[i][j] -  0.5 * K * (A * (u[i1][j] - u[i2][j]) + B * (u[i][j1] - u[i][j2])) + 0.5 * K * K * c * c * (u[i1][j] + u[i2][j] + u[i][j1] + u[i][j2] - 4 * u[i][j]);
}
Vector3d System::irrEquation(int i, int j) {
    Vector3d res (0.0, 0.0, 0.0);

    int i1 = (x_size + (i - 1)) % x_size;
    int i2 = (x_size + (i + 1)) % x_size;

    int j1 = (y_size + (j - 1)) % y_size;
    int j2 = (y_size + (j + 1)) % y_size;
    const vector<pair<int, int>> offsets = {{i1, j}, {i, j}, {i2, j}, {i, j2}, {i, j1}, {i2, j1}};
    // const vector<pair<int, int>> offsets = {{i, j1}, {i, j}, {i, j2}, {i2, j}, {i1, j}, {i1, j2}};

    pair<int, int> point = std::make_pair(i, j);
    for (int l = 0; l < 6; ++l) {
        int new_i = offsets[l].first;
        int new_j = offsets[l].second;

        res += gamma_matrices[point][l] * u[new_i][new_j];
    }
    return u[i][j] + k / h * res; 
}
void System::solve(double t) {
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size ; ++j) {
            if (i == 0 || i == x_size - 1 || j == 0 || j == y_size - 1) {
                double x0 = 0;
                double y0 = -0.4;
                double dx = 1;
                double dy = 10;
                double d_norm = sqrt(dx * dx + dy * dy);
                double d_norm_x = dx / d_norm;
                double d_norm_y = dy / d_norm;
                double A = 1;
                double omega = 5;
                double _c, _rho, pressure;
                if (func(get_x(i), get_y(j)) <= 0) {
                    _c = c_minus;
                    _rho = rho_minus;
                }
                else {
                    _c = c_plus;
                    _rho = rho_plus;
                }
                double theta = d_norm_x * (get_x(i) - x0) + d_norm_y * (get_y(j) - y0) - _c * t;
                if (theta >= 0 && theta <= (1 / omega)) {
                    pressure = (0.5 * A * (1 - cos(2 * M_PI * omega * theta)));
                }
                else {
                    pressure = 0.0;
                }
                u[i][j](0) = d_norm_x / (_rho * _c) * pressure;
                u[i][j](1) = d_norm_y / (_rho * _c) * pressure;
                u[i][j](2) = pressure;
                continue;
            }
            pair<int, int> point = std::make_pair(i, j);
            if (IsInIrregular(point)) {
                u[i][j] = irrEquation(i, j);
            }
            else {
                u[i][j] = equation(i, j);
            }
            // u[i][j] = equation(i, j);
        }
    }
}
void System::shift() {
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            Vector3d tmp (u[i][j](0), u[i][j](1), u[i][j](2));
            u_next[i][j] = tmp;
        }
    }
}
System::System(InitValues init) : PreProcess(init) {
    vector<vector<double>> pressure = init.GetPressure();
    vector<vector<double>> velocity_x = init.GetVelocityX();
    vector<vector<double>> velocity_y = init.GetVelocityY();

    for (int i = 0; i < x_size; ++i) {
        u.push_back(vector<Vector3d>());
        u_next.push_back(vector<Vector3d>());
        for (int j = 0; j < y_size; ++j) {
            Vector3d tmp(velocity_x[i][j], velocity_y[i][j], pressure[i][j]);
            u[i].push_back(tmp);
            u_next[i].push_back(tmp);
        }
    }
}
