#include <curve/curve.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::make_pair;

double Curve::func(double x, double y) {

    // return y - 0.3 - tan(0.4) * x;
    // return (y - 0.5) - h / 2;
    return x - 0.4 - h / 2;
    // return (y - 0.5) * (y - 0.5) + (x - 0.3) * (x - 0.3) - (0.2 + h / 2) * (0.2 + h / 2);
    // return y + 0.25 - h / 2;
    // return 0.3 + tan(0.4) * y - x;
    // return (x - h / 2) * (x - h / 2) + (y - h / 2) * (y - h / 2);
}
void Curve::FPrint() {
    std::ostringstream filename;
    filename << "../bin/animation/curve_points.bin";
    std::ofstream file_curve(filename.str(), std::ios::binary);
    for (auto c : irregular_points) {
        float _x = c.first * h;
        float _y = c.second * h;
        file_curve.write(reinterpret_cast<char*>(&_x), sizeof(_x));
        file_curve.write(reinterpret_cast<char*>(&_y), sizeof(_x));
    }
    file_curve.close();
}
bool Curve::is_opposite(int i, int j, int new_i, int new_j) {
    if (func(h * i, h * j) * func(h * new_i, h * new_j) < 0) return true;
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
    for (int i = 1; i < Mx - 1; ++i) {
        for (int j = 1; j < My - 1; ++j) {
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
Curve::Curve(double _tau, double _h, int _Mx, int _My, double _x0, double _y0, double _A, double _omega, double _alpha)
: InitValues(_tau, _h, _Mx, _My, _x0, _y0, _A, _omega, _alpha) {
    SetIrregularPoints();
    FPrint();
}
