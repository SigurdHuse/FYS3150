#include "solver.hpp"

int main(int argc, const char *argv[])
{
    // assert(argc == 11);
    int side_length, time_steps;
    double T, x_c, sigma_x, p_x, sigma_y, y_c, p_y, v0;
    std::string name;
    side_length = atoi(argv[1]);
    T = atof(argv[2]);
    time_steps = atoi(argv[3]);

    sigma_x = atof(argv[4]);
    x_c = atof(argv[5]);
    p_x = atof(argv[6]);

    sigma_y = atof(argv[7]);
    y_c = atof(argv[8]);
    p_y = atof(argv[9]);

    v0 = std::stod(argv[10]);
    name = argv[11];

    std::cout << v0 + 1 << "\n";
    Solver solver(side_length, T, time_steps, v0, name);
    solver.set_initial_state(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    solver.initialise_V();
    solver.simulate();
}