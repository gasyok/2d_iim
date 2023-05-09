#include <init/init.h>
#include <cmath>

void InitValues::SetTau(double _tau) {
    tau = _tau;
}
void InitValues::SetH(double _h) {
    h = _h;
}

void InitValues::SetRho(double _rho) {
    rho = _rho;
}
void InitValues::SetSpeed(double _c) {
    c = _c;
}
void InitValues::SetPressure(double x0, double y0, double A, double omega, double a, double b) {
    double dx = a;
    double dy = b;
    double d_norm = sqrt(dx * dx + dy * dy);
    double d_norm_x = dx / d_norm;
    double d_norm_y = dy / d_norm;

    for (int i = 0; i < size_x; ++i) {
        pressure.push_back(vector<double>());
        for (int j = 0; j < size_y; ++j) {
            // double theta = d_norm_x * (coord_x[i] - x0) + d_norm_y * (coord_y[j] - y0);
            // if (theta <= 1 / omega && theta >= 0) {
            //     pressure[i].push_back(0.5 * 3 * (1 - cos(2 * M_PI * omega * theta)));
            // }
            // else {
            //     pressure[i].push_back(0);
            // }
            pressure[i].push_back(A * exp(-(pow(((coord_x[i] - x0) * d_norm_x + (coord_y[j] - y0) * d_norm_y), 2) / (2 * omega * omega))));
            // pressure[i].push_back(A * exp(-(((coord_x[i] - x0) * (coord_x[i] - x0) + (coord_y[j] - y0) * (coord_y[j] - y0)) / (2 * sigma * sigma))));
            // pressure[i].push_back(exp(-626 * ((coord_x[i] - x0) * (coord_x[i] - x0) + (coord_y[j] - y0) * (coord_y[j] - y0))));
        }
    }
}
void InitValues::SetVelocity(double x0, double y0, double a, double b) {
    double dx = a;
    double dy = b;
    double d_norm = sqrt(dx * dx + dy * dy);
    double d_norm_x = dx / d_norm;
    double d_norm_y = dy / d_norm;
    for (int i = 0; i < size_x; ++i) {
        velocity_x.push_back(vector<double>());
        velocity_y.push_back(vector<double>());
        for (int j = 0; j < size_y; ++j) {
            velocity_x[i].push_back(pressure[i][j] * d_norm_x);
            velocity_y[i].push_back(pressure[i][j] * d_norm_y);
            // velocity_x[i].push_back(0);
            // velocity_y[i].push_back(0);
        }
    }
}
vector<double> InitValues::GetCoordX() {
    return coord_x;
}
vector<double> InitValues::GetCoordY() {
    return coord_y;
}
vector<vector<double>> InitValues::GetVelocityX() {
    return velocity_x;
}
vector<vector<double>> InitValues::GetVelocityY() {
    return velocity_y;
}
vector<vector<double>> InitValues::GetPressure() {
    return pressure;
}
int InitValues::GetSizeX() {
    return size_x;
}
int InitValues::GetSizeY() {
    return size_y;
}
double InitValues::GetRho() {
    return rho;
}
double InitValues::GetSpeed() {
    return c;
}
double InitValues::GetTau() {
    return tau;
}
double InitValues::GetH() {
    return h;
}
void InitValues::SetSizeX(double a, double b) {
    size_x = (b - a) / h;
}
void InitValues::SetSizeY(double a, double b) {
    size_y = (b - a) / h;
}
InitValues::InitValues() {}
InitValues::InitValues(double tau, double h, double c, double rho, double min_x,
                       double max_x, double min_y, double max_y, double x0, double y0, double A, double sigma, double a, double b) {
    SetTau(tau);
    SetH(h);
    SetSpeed(c);
    SetRho(rho);
    SetSizeX(min_x, max_x);
    SetSizeY(min_y, max_y);

    // Define x coord
    coord_x.clear();
    for (int i = 0; i < size_x; ++i) {
        coord_x.push_back(min_x + i * h);
    }

    // Define y coord
    for (int i = 0; i < size_y; ++i) {
        coord_y.push_back(min_y + i * h);
    }
    SetPressure(x0, y0, A, sigma, a, b);
    SetVelocity(x0, y0, a, b);

    rho_minus = 1;
    rho_plus = 2;
    c_minus = 3;
    c_plus = 1;
    k_minus = c_minus * c_minus * rho_minus;
    k_plus = c_plus * c_plus * rho_plus;

    A_minus << 0, 0, 1 / rho_minus,
               0, 0, 0,
               k_minus, 0, 0;
    A_plus << 0, 0, 1 / rho_plus,
               0, 0, 0,
               k_plus, 0, 0;
    B_minus << 0, 0, 0,
               0, 0, 1 / rho_minus,
               0, k_minus, 0;
    B_plus << 0, 0, 0,
               0, 0, 1 / rho_plus,
               0, k_plus, 0;
}
