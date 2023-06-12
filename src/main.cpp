#include <system/system.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::ofstream;

int main() {
    int M = 200;
    double x0 = 0.3;
    double y0 = 1;
    double A = 1;
    double omega = 0.1;
    double alpha = M_PI / 2 - 0.4;
    System mesh(M, x0, y0, A, omega, alpha);
    mesh.sample();
    return 0;
}
