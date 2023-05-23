#pragma once
#include <init/init.h>
#include <unordered_set>
#include <unordered_map>
#include <utility>

using std::pair;
using std::unordered_set;

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

class Curve {
private:
    unordered_set<pair<int, int>, pair_hash> irregular_points;
    unordered_set<pair<int, int>, pair_hash> regular_points;
    vector<double> x;
    vector<double> y;
    bool IsIrregular(int i, int j);

public:
    bool is_opposite(int i, int j, int new_i, int new_j);
    double diff(std::function<double(double, double)> func, double x, double y, bool is_x_derivative) {
        const double epsilon = 1e-6;
        if (is_x_derivative) {
            return (func(x + epsilon, y) - func(x - epsilon, y)) / (2 * epsilon);
        } else {
            return (func(x, y + epsilon) - func(x, y - epsilon)) / (2 * epsilon);
        }
    }
    void FPrint();
    double func(double x, double y);
    int x_size, y_size;
    double h, k;
    Curve();
    Curve(InitValues init);
    double curvature(double x, double y);
    void SetIrregularPoints();
    unordered_set<pair<int, int>, pair_hash> GetIrregularPoints();
    unordered_set<pair<int, int>, pair_hash> GetRegularPoints();
    double get_x(int i);
    double get_y(int j);
    bool IsInIrregular(pair<int, int> point);
};
