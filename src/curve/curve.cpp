#include <curve/curve.h>
#include <cmath>
using std::make_pair;

bool Curve::IsInIrregular(pair<int, int> point) {
    if (irregular_points.find(point) != irregular_points.end()) return true;
    return false;
}
double Curve::func(double x, double y) {

    return y - h / 2;
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
bool Curve::IsIrregular(int i, int j) {
    const std::vector<std::pair<int, int>> offsets = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    for (const auto& offset : offsets) {
        int new_i = (x_size + i + offset.first) % x_size;
        int new_j = (y_size + j + offset.second) % y_size;

        if (func(x[i], y[j]) * func(x[new_i], y[new_j]) < 0) {
            return true;
        }
    }
    return false;
}
void Curve::SetIrregularPoints() {
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            if (IsIrregular(i, j)) {
                pair<int, int> current_point = make_pair(i, j);
                irregular_points.insert(current_point);
            }
        }
    }
}
unordered_set<pair<int, int>, pair_hash> Curve::GetIrregularPoints() {
    return irregular_points;
}
Curve::Curve() {}
Curve::Curve(InitValues init) {
    h = init.GetH();
    k = init.GetTau();
    x = init.GetCoordX();
    y = init.GetCoordY();
    x_size = init.GetSizeX();
    y_size = init.GetSizeY();
    SetIrregularPoints();
}
