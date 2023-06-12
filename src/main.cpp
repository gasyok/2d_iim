#include <system/system.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::ofstream;

int main() {
    // int N = 800;
    // int Mx = mesh.Mx;
    // int My = mesh.My;
    int M = 400;
    double x0 = 1;
    double y0 = 1;
    double A = 1;
    double omega = 0.3;
    double alpha = M_PI / 2 - 0.4;

    System mesh(M, x0, y0, A, omega, alpha);
    mesh.sample();
    //
    // std::ofstream file_list("../bin/animation/velocity_out.pvd");
    // file_list << "<?xml version=\"1.0\"?>\n";
    // file_list << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    // file_list << "  <Collection>\n";
    //
    // for (int n = 0; n < N; ++n) {
    //     mesh.solve(n * mesh.tau);
        // std::ostringstream filename;
        // filename << "velocity_out_" << std::setfill('0') << std::setw(5) << n << ".vtk";
        // std::ofstream file_velocity("../bin/animation/" + filename.str());
        // 
        // // Запись данных в формате VTK, как в предыдущем примере...
        // file_velocity << "# vtk DataFile Version 3.0\n";
        // file_velocity << "Pressure data\n";
        // file_velocity << "ASCII\n";
        // file_velocity << "DATASET STRUCTURED_POINTS\n";
        // file_velocity << "DIMENSIONS " << Mx << " " << My << " 1\n";
        // file_velocity << "ORIGIN 0 0 0\n";
        // file_velocity << "SPACING " << mesh.h << " " << mesh.h << " 1\n";
        // file_velocity << "POINT_DATA " << (Mx * My) << "\n";
        // file_velocity << "SCALARS pressure float 1\n";
        // file_velocity << "LOOKUP_TABLE default\n";
        //
        // for (int i = 0; i < Mx; ++i) {
        //     for (int j = 0; j < My; ++j) {
        //         float pressure_value = mesh.GetValue(i, j)(2);
        //         file_velocity << pressure_value << "\n";
        //     }
        // }
        // file_velocity.close();
        // file_list << "    <DataSet timestep=\"" << n << "\" part=\"0\" file=\"" << filename.str() << "\"/>\n";
    // }

    // file_list << "  </Collection>\n";
    // file_list << "</VTKFile>\n";
    // file_list.close();
    return 0;
}
