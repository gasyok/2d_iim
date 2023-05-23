#include <curve/curve.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::make_pair;

bool Curve::IsInIrregular(pair<int, int> point) {
    if (irregular_points.find(point) != irregular_points.end()) return true;
    return false;
}
double Curve::func(double x, double y) {

    // return y - 0.3 - tan(0.4) * x;
    // return (y - 0.5) - h / 2;
    // return x - 0.4 - h / 2;
    return (y - 0.5) * (y - 0.5) + (x - 0.3) * (x - 0.3) - (0.2 + h / 2) * (0.2 + h / 2);
    // return y + 0.25 - h / 2;
    // return 0.3 + tan(0.4) * y - x;
    // return (x - h / 2) * (x - h / 2) + (y - h / 2) * (y - h / 2);
}
double Curve::curvature(double x, double y) {
    double f_x = diff([this](double x, double y) { return func(x, y); }, x, y, true);
    double f_y = diff([this](double x, double y) { return func(x, y); }, x, y, false);

    double f_xx = diff([this, &f_x](double x, double y) { return diff([this](double x, double y) { return func(x, y); }, x, y, true);}, x, y, true);
    double f_xy = diff([this, &f_y](double x, double y) { return diff([this](double x, double y) { return func(x, y); }, x, y, true);}, x, y, false);
    double f_yx = diff([this, &f_x](double x, double y) { return diff([this](double x, double y) { return func(x, y); }, x, y, false);}, x, y, true);
    double f_yy = diff([this, &f_y](double x, double y) { return diff([this](double x, double y) { return func(x, y); }, x, y, false);}, x, y, false);

    double curvature = (f_xx * f_yy - f_xy * f_yx) / pow((f_x * f_x + f_y * f_y), 1.5);

    return curvature;
}
double Curve::get_x(int i) {
    return x[i];
}
double Curve::get_y(int j) {
    return y[j];
}
void Curve::FPrint() {
    std::ostringstream filename;
    filename << "../bin/animation/curve_points.bin";
    std::ofstream file_curve(filename.str(), std::ios::binary);
    for (auto c : irregular_points) {
        double tmp = 0.3 + h / 2;
        float _x = x[c.first];
        float _y = y[c.second];
        file_curve.write(reinterpret_cast<char*>(&_x), sizeof(_x));
        file_curve.write(reinterpret_cast<char*>(&_y), sizeof(_x));
    }
    file_curve.close();
}
bool Curve::is_opposite(int i, int j, int new_i, int new_j) {
    if (func(x[i], y[j]) * func(x[new_i], y[new_j]) < 0) return true;
    return false;
}
bool Curve::IsIrregular(int i, int j) {
    const std::vector<std::pair<int, int>> offsets = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    for (const auto& offset : offsets) {
        int new_i = i + offset.first;
        int new_j = j + offset.second;
        // int new_i = (x_size + i + offset.first) % x_size;
        // int new_j = (y_size + j + offset.second) % y_size;

        if (is_opposite(i, j, new_i, new_j)) return true;
    }
    return false;
}
void Curve::SetIrregularPoints() {
    for (int i = 1; i < x_size - 1; ++i) {
        for (int j = 1; j < y_size - 1; ++j) {
            pair<int, int> current_point = make_pair(i, j);
            // if (i == 0 || i == x_size - 1 || j == 0 || j == y_size - 1) {
            //     regular_points.insert(current_point);
            //     continue;
            // }
            if (IsIrregular(i, j)) {
                irregular_points.insert(current_point);
            }
            else {
                regular_points.insert(current_point);
            }
        }
    }
}
unordered_set<pair<int, int>, pair_hash> Curve::GetIrregularPoints() {
    return irregular_points;
}
unordered_set<pair<int, int>, pair_hash> Curve::GetRegularPoints() {
    return regular_points;
}
Curve::Curve() : Curve(InitValues()) {}
Curve::Curve(InitValues init) : h(init.GetH()), k(init.GetTau()), x(init.GetCoordX()), y(init.GetCoordY()), x_size(init.GetSizeX()), y_size(init.GetSizeY()) {
    SetIrregularPoints();
    FPrint();
}
