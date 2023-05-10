#include <preprocess/preprocess.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

pair<int, int> PreProcess::GetPoint(int l, pair<int, int> point) {
    int i = point.first;
    int j = point.second;
    int i1 = (x_size + (i - 1)) % x_size;
    int i2 = (x_size + (i + 1)) % x_size;

    int j1 = (y_size + (j - 1)) % y_size;
    int j2 = (y_size + (j + 1)) % y_size;

    const vector<pair<int, int>> offsets = {{i, j1}, {i, j}, {i, j2}, {i2, j}, {i1, j}, {i1, j2}};
    return offsets[l];
}
Vector2d PreProcess::bisection(double x1, double x2, double y1, double y2, double tol=1e-8, int max_iter = 1000) {
    double x0 = (x1 + x2) / 2.0;
    double y0 = (y1 + y2) / 2.0;

    int iter = 0;
    double prev_x0, prev_y0;
    const double epsilon = tol * 10;

    while ((std::abs(func(x0, y0)) > tol || std::abs(prev_x0 - x0) > tol || std::abs(prev_y0 - y0) > tol) && iter < max_iter) {
        prev_x0 = x0;
        prev_y0 = y0;

        // Find y0
        if (func(x0, y1) * func(x0, y0) <= 0.0) {
            y2 = y0;
        } else {
            y1 = y0;
        }
        y0 = (y1 + y2) / 2.0;

        // Find x0
        if (func(x1, y0) * func(x0, y0) <= 0.0) {
            x2 = x0;
        } else {
            x1 = x0;
        }
        x0 = (x1 + x2) / 2.0;

        iter++;
    }
    return Eigen::Vector2d(x0, y0);
}
Vector2d PreProcess::GetOrigin(pair<int, int> point) {
    int i = point.first;
    int j = point.second;
    int i1 = (x_size + (i - 1)) % x_size;
    int i2 = (x_size + (i + 1)) % x_size;

    int j1 = (y_size + (j - 1)) % y_size;
    int j2 = (y_size + (j + 1)) % y_size;


    const vector<pair<int, int>> offsets = {{i, j1}, {i, j}, {i, j2}, {i2, j}, {i1, j}};
    Vector2d res;

    for (const auto& offset : offsets) {
        int new_i = offset.first;
        int new_j = offset.second;

        if(is_opposite(i, j, new_i, new_j)) {
            res = bisection(get_x(i), get_x(new_i), get_y(j), get_y(new_j));
            break;
        }
    }
    return res;
}
Eigen::Vector2d PreProcess::normalize(const Eigen::Vector2d& vec) {
    double norm = vec.norm();
    if (norm > 0) {
        return vec / norm;
    } else {
        return Eigen::Vector2d::Zero();
    }
}
Vector2d PreProcess::GetRotationalCoord(int l, pair<int, int>point) {
    // Получение координат x_l, y_l
    std::pair<int, int> point_l = GetPoint(l, point);
    double x_l = get_x(point_l.first);
    double y_l = get_y(point_l.second);
    Vector2d vec_point (x_l, y_l);
    Vector2d origin = GetOrigin(point);

    // Создание матрицы обратного поворота R(-β)
    Eigen::Matrix2d R_inv;
    // Vector2d Deriv (DerivationX(origin(0), origin(1)), DerivationY(origin(0), origin(1)));
    // Vector2d Deriv (diff(func, origin(0), origin(1), true), diff(func, origin(0), origin(1), false));
    Vector2d Deriv(diff([this](double x, double y) { return func(x, y); }, origin(0), origin(1), true),
                   diff([this](double x, double y) { return func(x, y); }, origin(0), origin(1), false));

    Vector2d normalized = normalize(Deriv);

    R_inv << normalized(0), -normalized(1),
        normalized(1),  normalized(0);

    // Создание вектора (x_l - x₀, y_l - y₀)

    // Умножение матрицы обратного поворота на вектор (x_l - x₀, y_l - y₀)
    Eigen::Vector2d result = R_inv * (vec_point - origin);

    return result;
}
Matrix3d PreProcess::RotateMatrix(pair<int, int> point, Matrix3d matrix) {
    Matrix3d Q_zero;
    Vector2d origin = GetOrigin(point);

    Vector2d Deriv(diff([this](double x, double y) { return func(x, y); }, origin(0), origin(1), true),
                   diff([this](double x, double y) { return func(x, y); }, origin(0), origin(1), false));

    Vector2d normalized = normalize(Deriv);
    Q_zero << normalized(0), -normalized(1), 0,
        normalized(1), normalized(0), 0,
        0, 0, 1;
    return Q_zero * matrix * Q_zero.inverse();
}
vector<Matrix3d> PreProcess::GetDefaultQ (pair<int, int> point) {
    double rho_temp, k_temp, c_temp;
    const double eps = 1e-12;
    rho_temp = rho_minus / rho_plus;
    k_temp = k_minus / k_plus;
    c_temp = c_minus / c_plus;
    if (func(get_x(point.first), get_y(point.second)) > 0.0) {
        rho_temp = 1 / rho_temp;
        k_temp = 1 / k_temp;
        c_temp = 1 / c_temp;
    }
    Vector3d q1_diag (1, rho_temp, 1);
    Vector3d q2_diag (k_temp, 1,  1 / rho_temp);

    Matrix3d Zero1, Zero2;
    Zero1 << 0, 1, 0,
        0, 0, 0,
        0, 0, 0;

    Zero2 << 0, 1, 0,
        -1, 0, 0,
        0, 0, 0;

    Matrix3d q1 = q1_diag.asDiagonal();
    Matrix3d q2 = q2_diag.asDiagonal();
    Matrix3d q3 = rho_temp * (c_temp * c_temp - 1) * Zero1;
    Matrix3d q4 = Zero2;

    vector<Matrix3d> result {q1, q2, q3, q4};
    return result;
}
Matrix3d PreProcess::OppositeQ(int i, int l, pair<int, int> point) {
    Vector2d new_coord = GetRotationalCoord(l, point);
    double xi = new_coord(0);
    const double eps = 1e-12;
    double eta = new_coord(1);
    double c_temp = c_minus / c_plus;
    if (func(get_x(point.first), get_y(point.second)) > 0.0) {
        c_temp = 1 / c_temp;
    }
    Vector2d origin = GetOrigin(point);
    double curvature_value = curvature(origin(0), origin(1));

    vector<Matrix3d> matrices_default = GetDefaultQ(point);
    Matrix3d q1 = matrices_default[0];
    Matrix3d q2 = (xi / h) * matrices_default[1] + 0.5 * (1 / (h * h)) * curvature_value * (xi * xi * Matrix3d::Identity() - eta * eta * Matrix3d::Identity() + 2 * xi * eta * matrices_default[3]) * (matrices_default[1] - matrices_default[0]);
    Matrix3d q3 = (eta / h) * matrices_default[0] + (xi / h) * matrices_default[2] + 0.5 / (h * h) * curvature_value * (xi * xi * Matrix3d::Identity() - eta * eta * Matrix3d::Identity() + 2 * xi * eta * matrices_default[3]) * matrices_default[2];
    Matrix3d q4 = (xi / h) * (xi / h) * c_temp * c_temp * matrices_default[0];
    Matrix3d q5 = 2 * xi * eta / (h * h) * matrices_default[1];
    Matrix3d q6 = (eta / h) * (eta / h) * matrices_default[0] + 2 * xi * eta / (h * h) * matrices_default[2] + (xi / h) * (xi / h) * (c_temp * c_temp - 1) * matrices_default[0];

    vector<Matrix3d> matrices = {q1, q2, q3, q4, q5, q6};
    return matrices[i];
}
Matrix3d PreProcess::BesideQ(int i, int l, pair<int, int> point) {
    Vector2d new_coord = GetRotationalCoord(l, point);
    double xi = new_coord(0);
    double eta = new_coord(1);
    Matrix3d q1 = Matrix3d::Identity();
    Matrix3d q2 = xi * Matrix3d::Identity();
    Matrix3d q3 = eta * Matrix3d::Identity();
    Matrix3d q4 = xi * xi * Matrix3d::Identity();
    Matrix3d q5 = 2 * xi * eta * Matrix3d::Identity();
    Matrix3d q6 = eta * eta * Matrix3d::Identity();

    vector<Matrix3d> matrices = {q1, q2, q3, q4, q5, q6};
    return matrices[i];
}
Matrix3d PreProcess::GetQmatrix(int i, int l, pair<int, int> point) {
    pair<int, int> new_point = GetPoint(l, point);
    if (is_opposite(point.first, point.second, new_point.first, new_point.second)) {
        return OppositeQ(i, l, point);
    }
    else {
        return BesideQ(i, l, point);
    }
}
Matrix3d PreProcess::GetFmatrix(int i, pair<int, int> point) {
    Matrix3d A = A_minus, B = B_minus;
    double c = c_minus;

    Vector2d new_coord = GetRotationalCoord(1, point);
    double xi = new_coord(0);
    double eta = new_coord(1);
    const double eps = 1e-12;

    if (func(get_x(point.first), get_y(point.second)) > 0.0) {
        A = A_plus;
        B = B_plus;
        c = c_plus;
    }
    Matrix3d f1 = Matrix3d::Zero();
    Matrix3d f2 = -h * A;
    Matrix3d f3 = -h * B;
    Matrix3d f4 = h * k * c * c * Matrix3d::Identity() - 2 * h * xi * A;
    Matrix3d f5 = -2 * h * (eta * A + xi * B);
    Matrix3d f6 = h * k * c * c * Matrix3d::Identity() - 2 * h * eta * B;

    vector<Matrix3d> matrices = {f1, f2, f3, f4, f5, f6};
    return matrices[i];
}
vector<Matrix3d> PreProcess::CalcGammaMatrices(pair<int, int> point) {
    Matrix<double, 18, 18> Q;
    Eigen::Matrix<double, 18, 3> F;
    Q.setZero();
    F.setZero();

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            Q.block<3, 3>(3 * i, 3 * j) = GetQmatrix(i, j, point).transpose();
        }
    }
    for (int i = 0; i < 6; ++i) {
        F.block<3, 3>(3 * i, 0) = GetFmatrix(i, point).transpose();
    }
    vector<Matrix3d> block_matrices;
    Matrix<double, 18, 3> G = Q.colPivHouseholderQr().solve(F);

    for (int i = 0; i < 18; i += 3) {
        Matrix3d block = G.block<3, 3>(i, 0);
        Matrix3d gamma = RotateMatrix(point, block.transpose());
        block_matrices.push_back(gamma);
        // block_matrices.push_back(block.transpose());
    }
    return block_matrices;
}
PreProcess::PreProcess(InitValues init) : Curve(init) {
    A_minus = init.A_minus;
    A_plus = init.A_plus;
    B_minus = init.B_minus;
    B_plus = init.B_plus;
    rho_minus = init.rho_minus;
    rho_plus = init.rho_plus;
    c_minus = init.c_minus;
    c_plus = init.c_plus;
    k_minus = init.k_minus;
    k_plus = init.k_plus;
    Solve();
}
void PreProcess::Solve() {
    for (const auto& point: GetIrregularPoints()) {
        vector<Matrix3d> gammas = CalcGammaMatrices(point);
        gamma_matrices[point] = gammas;
    }
}
