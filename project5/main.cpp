#include "solver.hpp"

#include <sys/stat.h>

/**
 * @brief Simulates the time evolution of the probability of a wave packet using the Crank-Nicolson method to
 * approximate the Schrodinger equation.
 *
 * @param side_length Integer, number of points along x- and y-axis with length 1
 * @param T Double, total time of simulation
 * @param time_steps Integer, number of time steps in simulation
 * @param sigma_x Double, initial width of wave packet in x-direction
 * @param x_c Double, x-coordinate of starting position of wave packet
 * @param p_x Double, wave packet momenta in x-direction
 * @param sigma_y Double, initial width of wave packet in y-direction
 * @param y_c Double, y-coordinate of starting position of wave packet
 * @param p_y Double, wave packet momenta in y-direction
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

void generate_all_data()
{
    int side_length = 200, time_steps = 320;
    double T = 0.008, x_c = 0.25, sigma_x = 0.05, p_x = 200;
    double sigma_y = 0.05, y_c = 0.5, p_y = 0, v0 = 0;

    Solver solver1(side_length, T, time_steps, v0, "No_slit_sigma_y_005");
    solver1.set_initial_state(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    solver1.initialise_V("configs/config1.txt");

    v0 = 1e10;
    sigma_y = 0.1;
    Solver solver2(side_length, T, time_steps, v0, "Two_slits_sigma_y_0.1");
    solver2.set_initial_state(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    solver2.initialise_V("configs/config2.txt");

    T = 0.002;
    time_steps = 80;
    sigma_y = 0.2;

    Solver solver3(side_length, T, time_steps, v0, "Two_slits_sigma_y_0.2");
    solver3.set_initial_state(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    solver3.initialise_V("configs/config3.txt");

    Solver solver4(side_length, T, time_steps, v0, "One_slit_sigma_y_0.2");
    solver4.set_initial_state(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    solver4.initialise_V("configs/config4.txt");

    Solver solver5(side_length, T, time_steps, v0, "Three_slits_sigma_y_0.2");
    solver5.set_initial_state(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    solver5.initialise_V("configs/config5.txt");

    T = 0.2;
    time_steps = 8000;

    Solver solver6(side_length, T, time_steps, v0, "Seven_slits_sigma_y_0.2");
    solver6.set_initial_state(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
    solver6.initialise_V("configs/config6.txt");

    // solver1.simulate();
    // solver2.simulate();
    // solver3.simulate();
    // solver4.simulate();
    // solver5.simulate();
    solver6.simulate();
}

int main(int argc, const char *argv[])
{
    std::string folder_name = "data";
    mkdir(folder_name.c_str(), 0777);

    if (argc == 1)
    {
        std::cout << "Generating all data from report!\n";
        generate_all_data();
    }
    else
    {
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
}