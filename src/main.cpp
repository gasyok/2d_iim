#include <system/system.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::ofstream;

int main() {
    InitValues init;
    System mesh(init);
    vector<double> coord_x = init.GetCoordX();
    vector<double> coord_y = init.GetCoordY();
    int N = 800;
    int Mx = init.GetSizeX();
    int My = init.GetSizeY();
    for (int n = 0; n < N; ++n) {
        mesh.solve(n * init.GetTau());
        mesh.shift();
        // Формирование имени файла с индексом временного шага
        std::ostringstream filename;
        filename << "../bin/animation/velocity_out_" << std::setfill('0') << std::setw(5) << n << ".bin";
        std::ofstream file_velocity(filename.str(), std::ios::binary);
        // Запись количества точек данных для текущего временного шага (4 байта, little-endian)
        for (int i = 0; i < Mx; ++i) {
            for (int j = 0; j < My; ++j) {
                float x_coord = coord_x[i];
                float y_coord = coord_y[j];
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
