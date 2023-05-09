#include <system/system.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::ofstream;

int main() {
    // InitValues(double tau, double h, double c, double rho, double min_x,
    //            double max_x, double min_y, double max_y, double x0, double y0, A, sigma);
    // InitValues init(0.0002, 0.01, 3, 2, -0.5, 0.5, -0.5, 0.5, 0, 0, 1, 0.3, -1 / sqrt(2), -1 / sqrt(2));
    InitValues init(0.0002, 0.01, 3, 2, -0.5, 0.5, -0.5, 0.5, 0.25, 0.25, 1, 0.05, -1, -3);
    vector<double> coord_x = init.GetCoordX();
    vector<double> coord_y = init.GetCoordY();
    int N = 2;
    int Mx = init.GetSizeX();
    int My = init.GetSizeY();
    System mesh(init);
    for (int n = 0; n < N; ++n) {
        mesh.solve();
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

    // int result = std::system("python3 message_post.py");
    // if (result == 0) {
    //     std::cout << "Python script executed successfully!" << std::endl;
    // } else {
    //     std::cerr << "Failed to execute Python script. Error code: " << result << std::endl;
    // }
    return 0;
}
