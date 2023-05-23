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

    if (func(get_x(i), get_y(j)) > 0.0) {
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
    double x0 = boundary[0];
    double y0 = boundary[1];
    double A = boundary[2];
    double omega = boundary[3];
    double alpha = boundary[4];
    Vector2d a_initial (cos(alpha), sin(alpha));
    Vector2d a_reflected (cos(alpha), -sin(alpha));
    Vector2d a_transmited (sqrt(1 - (c_plus / c_minus) * (c_plus / c_minus) * sin(alpha) * sin(alpha)), c_plus / c_minus * sin(alpha));
    double phi = 0.4;
    Matrix2d R;
    double cos_alpha_transm = sqrt(1 - (c_plus / c_minus) * (c_plus / c_minus) * sin(alpha) * sin(alpha));
    R << cos(phi), -sin(phi),
          sin(phi), cos(phi);
    double A_r = (rho_plus * c_plus * cos(alpha) - rho_minus * c_minus * cos_alpha_transm) / (rho_plus * c_plus * cos(alpha) + rho_minus * c_minus * cos_alpha_transm);
    double A_t = 1 + A_r;
    // double A_t = (2 * rho_plus * c_plus * cos(alpha)) / (rho_plus * c_plus * cos(alpha) + rho_minus * c_minus * cos_alpha_transm);

    a_initial = R * a_initial;
    a_reflected = R * a_reflected;
    a_transmited = R * a_transmited;
    double omega_t = omega * cos(alpha) / cos_alpha_transm;

    // for (int i = 1; i < x_size - 1; ++i) {
    //     for (int j = 1; j < y_size - 1; ++j) {
    //         pair<int, int> point = std::make_pair(i, j);
    //         if (IsInIrregular(point)) {
    //             u[i][j] = irrEquation(i, j);
    //         }
    //         else {
    //             u[i][j] = equation(i, j);
    //         }
    //     }
    // }
    for (auto c : irregular_points) {
        u[c.first][c.second] = irrEquation(c.first, c.second);
    }
    for (auto c : regular_points) {
        u[c.first][c.second] = equation(c.first, c.second);
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

    boundary = init.GetInitU();

    for (int i = 0; i < x_size; ++i) {
        u.push_back(vector<Vector3d>());
        u_next.push_back(vector<Vector3d>());
        for (int j = 0; j < y_size; ++j) {
            Vector3d tmp(velocity_x[i][j], velocity_y[i][j], pressure[i][j]);
            u[i].push_back(tmp);
            u_next[i].push_back(tmp);
        }
    }
    irregular_points = GetIrregularPoints();
    regular_points = GetRegularPoints();
}
