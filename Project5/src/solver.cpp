#include "solver.hpp"
#include "grid.hpp"

/*
NOTE (if you have read the report)

If reviewing the code after reading the report, note that where as the coordinate of u_{ij}^n in the paper is (x_i,y_i).
The opposite is true in this code as the elements (x_i, y_j) corresponds to U(j,i). This could be changed in the future,
but as the code might seems weird without this information it is reported on. Hopefully it is not too confusing, as it is simply
the way I think about matricies when coding them.
*/

/*
NOTE (if you have not read the report)

You should read the report! All definitons come from there, and a more in depth explantation of what is simulated is presented.
*/

/**
* @brief Constructor for Solver class.

* @param side_length Integer, determines the number points along the x-and y-axis in the discretised box.
* @param time Double, sets the time for how long the system is to be simulated for.
* @param time_delta Double. sets how big the time step is.
* @param v0 Double, value of walls in potentital.
* @param file_name String, name of file to store data from simulation too.
*/
Solver::Solver(int side_length, double time, int time_delta, double v0, std::string file_name)
{
    // Initialises A and B matrix as defined in the report
    A_matrix.set_side_length(side_length);
    A_matrix.set_time(time);
    A_matrix.set_time_step(time_delta);

    B_matrix.set_side_length(side_length);
    B_matrix.set_time(time);
    B_matrix.set_time_step(time_delta);

    current_state = arma::cx_mat(side_length, side_length);
    current_state_vec = arma::cx_colvec((side_length - 2) * (side_length - 2));
    states = arma::cx_cube(side_length, side_length, time_delta + 1);

    M = side_length;
    M_squared = (side_length - 2) * (side_length - 2);
    h = 1.0 / (side_length - 1);
    time_steps = time_delta;
    v_0 = v0;
    name = file_name;

    filename = "data/" + file_name + "_M_" + std::to_string(side_length) + "_dt_" + std::to_string(time_delta);
}

/**
 * @brief Reads values from specified file.
 *
 * @param thickness Reference to integer determining the thickness of the wall.
 * @param x_pos Reference to integer determning x-coordinate to center wall around.
 * @param seperation Reference to integer determining length of slits.
 * @param aperture Reference to integer determining seperation between slits.
 * @param slits Reference to integer determining number of slits.
 * @param filename String, file to read data from.
 */
void get_values_from_file(double &thickness, double &x_pos, double &seperation, double &aperture, int &slits, std::string filename)
{
    std::string line;
    std::ifstream infile(filename);
    std::getline(infile, line);
    std::getline(infile, line);
    thickness = std::atof(line.c_str());
    std::getline(infile, line);
    std::getline(infile, line);
    x_pos = std::atof(line.c_str());
    std::getline(infile, line);
    std::getline(infile, line);
    seperation = std::atof(line.c_str());
    std::getline(infile, line);
    std::getline(infile, line);
    aperture = std::atof(line.c_str());
    std::getline(infile, line);
    std::getline(infile, line);
    slits = std::atoi(line.c_str());
    infile.close();
}

/**
 * @brief Initialises the potenital matrix V, and has to be called as a class function for the system to be able to simulate.
 *
 * @param name String, file name for the file from which the configurations for potenital are read from using the get_values_from_file() function.
 */
void Solver::initialise_V(std::string name)
{
    double thickness, x_pos, seperation, aperture;
    int slits;
    V = arma::mat(M - 2, M - 2);

    get_values_from_file(thickness, x_pos, seperation, aperture, slits, name);
    filename += "_slits_" + std::to_string((int)slits) + ".bin";

    int start_x = x_pos / h;
    int wall_thickness = thickness / h;
    start_x -= wall_thickness / 2;
    int length_of_wall = seperation / h;
    int opening = aperture / h;

    // Subtract two as V is of size (M-1)x(M-1)
    int start_y = ypos / h - 2;

    if (slits & 1)
    {
        // Goes half a space down as centre is a space
        start_y -= opening / 2;

        // number of walls and slits below center
        start_y -= (slits - 1) / 2 * length_of_wall;
        start_y -= (slits - 1) / 2 * opening;
    }
    else
    {
        // Goes half a wall down as centre is wall
        start_y -= length_of_wall / 2;

        // Number of slits and walls below center
        start_y -= slits / 2 * opening;
        start_y -= (slits - 2) / 2 * length_of_wall;
    }

    // Initialises wall
    for (int y = 0; y < M - 2; ++y)
    {
        for (int x = start_x; x <= start_x + wall_thickness; ++x)
        {
            V(y, x) = v_0;
        }
    }

    // Makes slits in the wall
    for (int slit = 1; slit <= slits; ++slit)
    {
        for (int y = start_y; y <= start_y + opening; ++y)
        {
            for (int x = start_x; x <= start_x + wall_thickness; ++x)
            {
                V(y, x) = 0;
            }
        }
        start_y += length_of_wall + opening + 1;
    }

    // Used for debugging
    // V.save(std::to_string(slits) + "slit.txt", arma::raw_ascii);

    // Fills A and B matrix using V
    fill_matrices();
}

// Fills A and B matrix
void Solver::fill_matrices()
{
    A_matrix.fill_matrix(V, 1);
    B_matrix.fill_matrix(V, 0);
}

/**
 * @brief Sets inital state of system with Dirichlet boundary conditions and normalises all grid squares to have absolute value 1.
 *
 * @param x_c Double determining initial x-coordinate of the wave packet.
 * @param y_c Double determining initial y-coordinate of the wave packet.
 * @param sigma_x Double determining initial width of the of the wave packet in x-direction.
 * @param sigma_y Double determining initial width of the of the wave packet in y-direction.
 * @param p_x Double determining wave packet momenta in x-direction.
 * @param p_y Double determining wave packet momenta in y-direction.
 */
void Solver::set_initial_state(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y)
{
    double real_comp, imag_comp;
    double sigma_x_squared_2 = 2 * sigma_x * sigma_x;
    double sigma_y_squared_2 = 2 * sigma_y * sigma_y;
    ypos = y_c;

    long double normalization_factor = 0;
    for (int y = 1; y < M - 1; ++y)
    {
        for (int x = 1; x < M - 1; ++x)
        {
            real_comp = -pow(h * x - x_c, 2) / sigma_x_squared_2 - pow(h * y - y_c, 2) / sigma_y_squared_2;
            imag_comp = p_x * (h * x - x_c) + p_y * (h * y - y_c);
            std::complex<double> exponent(real_comp, imag_comp);
            current_state(y, x) = exp(exponent);
            normalization_factor += norm(current_state(y, x));
        }
    }
    // We take the sqrt as we are taking the sum of each value squared
    current_state /= sqrt(normalization_factor);
}

// Converts current state to vector using class function indicies_to_index().
void Solver::convert_current_state_to_vector()
{
    for (int y = 1; y < M - 1; ++y)
    {
        for (int x = 1; x < M - 1; ++x)
        {
            current_state_vec[indicies_to_index(y, x)] = current_state(y, x);
        }
    }
}

//

/**
 * @brief Updates current state of system from vector, given that the vector contains all interal points of the grid. The vector also
 * needs to contain the columns of the matrix one after the other.
 *
 * @param v Complex arma vector, representing the system one colummn after another.
 */
void Solver::update_current_state(arma::cx_vec &v)
{
    // v.raw_print();
    for (int y = 1; y < M - 1; ++y)
    {
        for (int x = 1; x < M - 1; ++x)
        {
            current_state(y, x) = v(indicies_to_index(y, x));
        }
    }
    // print_current();
}

//

/**
 * @brief Translates pair of indicies to single index, to be able to go from matrix representation to vector representation.
 * Might seem weird after reading the report, for clarification see NOTE on top of file.
 *
 * @param y Integer, y-coordinate in matrix
 * @param x Integer, x-coordinate in matrix
 * @return Integer, the index in the vector representation of the matrix.
 */
int Solver::indicies_to_index(int y, int x)
{
    return y - 1 + (M - 2) * (x - 1);
}

// Computes b in Crank-Nicolson method
arma::cx_vec Solver::compute_b()
{
    convert_current_state_to_vector();
    return B_matrix.multiply_matrix_with_vector(current_state_vec);
}

// Finds next state using Crank-Nicolson
void Solver::find_next_state()
{
    arma::cx_vec b = compute_b();

    arma::cx_vec ans;

    // Solves the matrix equation
    arma::spsolve(ans, A_matrix.get_matrix(), b, "superlu");
    update_current_state(ans);
}

// Prints the current state of the system
void Solver::print_current()
{
    current_state.raw_print();
}

// Prints current state of system as vector
void Solver::print_current_state_vector()
{
    convert_current_state_to_vector();
    current_state_vec.raw_print();
}

// Simulates the system
void Solver::simulate()
{
    // Loops over all time steps and saves the current state in a cube
    for (int i = 0; i < time_steps; ++i)
    {
        states.slice(i) = current_state;
        find_next_state();
    }
    states.slice(time_steps) = current_state;
    write_states_to_file();
}

// Writes all computed states to file
void Solver::write_states_to_file()
{
    states.save(filename);
}

// Computes the sum of probabilities in the system
double Solver::compute_sum_of_probabilities()
{
    double ans = 0;
    for (int y = 1; y < M - 1; ++y)
    {
        for (int x = 1; x < M - 1; ++x)
        {
            ans += norm(current_state(y, x));
        }
    }
    return ans;
}