#include <system/system.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::ofstream;

int main() {
    System mesh(0.0004, 0.01, 100, 100, 0.6, 0.5, 1, 10, M_PI / 2 - 0.4);
    int N = 800;
    int Mx = mesh.Mx;
    int My = mesh.My;
    for (int n = 0; n < N; ++n) {
        mesh.solve(n * mesh.tau);
        // Формирование имени файла с индексом временного шага
        std::ostringstream filename;
        filename << "../bin/animation/velocity_out_" << std::setfill('0') << std::setw(5) << n << ".bin";
        std::ofstream file_velocity(filename.str(), std::ios::binary);
        // Запись количества точек данных для текущего временного шага (4 байта, little-endian)
        for (int i = 0; i < Mx; ++i) {
            for (int j = 0; j < My; ++j) {
                float x_coord = mesh.h * i;
                float y_coord = mesh.h * j;
                float pressure_value = mesh.GetValue(i, j)(2);
                // Запись x, y и p(x, y) (каждый параметр - 4 байта, little-endian)
                file_velocity.write(reinterpret_cast<char*>(&x_coord), sizeof(x_coord));
                file_velocity.write(reinterpret_cast<char*>(&y_coord), sizeof(y_coord));
                file_velocity.write(reinterpret_cast<char*>(&pressure_value), sizeof(pressure_value));
            }
        }
        file_velocity.close();
    }
    return 0;
}
