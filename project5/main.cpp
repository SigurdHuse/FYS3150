#include "solver.hpp"

#include <sys/stat.h>

/**
 * @brief Simulates the time evolution of the probability of a wave packet using the Crank-Nicolson method to
 * approximate the Schrodinger equation.
 *
 * @param side_length Integer, number of points along x- and y-axis with length 1
 * @param T Double, total time of simulation
 * @param time_steps Integer, number of time steps in simulation
 * @param sigma_x Double, sigma_x parameter in simulation
 * @param x_c Double, x-coordinate of starting position of wave packet
 * @param p_x Double, p_x in simulation
 * @param sigma_y Double, sigma_y in simulation
 * @param y_c Double, y-coordinate of starting position of wave packet
 * @param p_y Double, p_y in simulation
 * @param v0 Double, value of potential in the wall in simulation
 * @param name String, name of simulation
 */
void simulate(int side_length, double T, int time_steps, double sigma_x, double x_c, double p_x,
              double sigma_y, double y_c, double p_y, double v0, std::string name)
{
    Solver solver(side_length, T, time_steps, v0, name);
    solver.set_initial_state(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    solver.initialise_V("config.txt");
    solver.simulate();
}

int main(int argc, const char *argv[])
{
    std::string folder_name = "data";
    mkdir(folder_name.c_str(), 0777);
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

    simulate(side_length, T, time_steps, sigma_x, x_c, p_x, sigma_y, y_c, p_y, v0, name);
}