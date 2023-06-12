#include <system/system.h>
#include <fstream>

Vector3d System::GetValue(int i, int j) {
    return u[i][j];
}
Vector3d System::equation(int i, int j) {
    int i1, i2, j1, j2;
    double K = tau / h;

    double c = c_minus;
    Matrix3d A = A_minus;
    Matrix3d B = B_minus;

    if (func(h * i, h * j) > 0.0) {
        A = A_plus;
        B = B_plus;
        c = c_plus;
    }
    
    i1 = (M + i + 1) % M;
    i2 = (M + i - 1) % M;
    j1 = (M + j + 1) % M;
    j2 = (M + j - 1) % M;

    return u[i][j] - 0.5 * tau / h * (A * (u[i1][j] - u[i2][j]) + B * (u[i][j1] - u[i][j2])) + \
            0.5 * (tau / h) * (tau / h) * c * c * (u[i1][j] + u[i2][j] + u[i][j1] + u[i][j2] - 4 * u[i][j]);
}
Vector3d System::irrEquation(int i, int j) {
    Vector3d res (0.0, 0.0, 0.0);

    int i1 = (M + (i - 1)) % M;
    int i2 = (M + (i + 1)) % M;

    int j1 = (M + (j - 1)) % M;
    int j2 = (M + (j + 1)) % M;

    const vector<pair<int, int>> offsets = {{i1, j}, {i, j}, {i2, j}, {i, j2}, {i, j1}, {i2, j1}};

    pair<int, int> point = std::make_pair(i, j);
    for (int l = 0; l < 6; ++l) {
        int new_i = offsets[l].first;
        int new_j = offsets[l].second;

        res += gamma_matrices[point][l] * u[new_i][new_j];
    }
    return u[i][j] + (tau / h) * res; 
}
void System::solve(double t) {
    vector<vector<Vector3d>> new_u = u;
    for (auto c : irregular_points) {
        new_u[c.first][c.second] = irrEquation(c.first, c.second);
    }
    for (auto c : regular_points) {
        new_u[c.first][c.second] = equation(c.first, c.second);
    }
    // for (int i = 1; i < M-1; ++i) {
    //     for (int j = 1; j < M-1; ++j) {
    //         new_u[i][j] = equation(i, j);
    //     }
    // }
    u = new_u;
}
void System::sample() {
    int n = 0;
    std::ofstream file_list("../bin/animation/velocity_out.pvd");
    file_list << "<?xml version=\"1.0\"?>\n";
    file_list << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file_list << "  <Collection>\n";
    while (n < total_steps) {
        std::ostringstream filename;
        filename << "velocity_out_" << std::setfill('0') << std::setw(5) << n << ".vtk";
        std::ofstream file_velocity("../bin/animation/" + filename.str());

        // Запись данных в формате VTK, как в предыдущем примере...
        file_velocity << "# vtk DataFile Version 3.0\n";
        file_velocity << "Pressure data\n";
        file_velocity << "ASCII\n";
        file_velocity << "DATASET STRUCTURED_POINTS\n";
        file_velocity << "DIMENSIONS " << M << " " << M << " 1\n";
        file_velocity << "ORIGIN 0 0 0\n";
        file_velocity << "SPACING " << h << " " << h << " 1\n";
        file_velocity << "POINT_DATA " << (M * M) << "\n";
        file_velocity << "SCALARS pressure double 1\n";
        file_velocity << "LOOKUP_TABLE default\n";

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                float pressure_value = u[i][j](2);
                file_velocity << pressure_value << "\n";
            }
        }
        file_velocity.close();
        file_list << "    <DataSet timestep=\"" << n << "\" part=\"0\" file=\"" << filename.str() << "\"/>\n";

        n++;
        solve(n * tau);
    }
    file_list << "  </Collection>\n";
    file_list << "</VTKFile>\n";
    file_list.close();
}

System::System(int _M, double _x0, double _y0, double _A, double _omega, double _alpha)
    : PreProcess(_M, _x0, _y0, _A, _omega, _alpha) {
    total_steps = 400;
}
