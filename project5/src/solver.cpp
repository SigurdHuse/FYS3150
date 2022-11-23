#include "solver.hpp"
#include "grid.hpp"

// Constructor
Solver::Solver(int side_length, double time, int time_delta)
{
    A_matrix.set_side_length(side_length);
    A_matrix.set_time(time);
    A_matrix.set_time_step(time_delta);

    B_matrix.set_side_length(side_length);
    B_matrix.set_time(time);
    B_matrix.set_time_step(time_delta);

    current_state = arma::cx_mat(side_length, side_length);
    states = arma::cx_cube(side_length, side_length, time_delta + 1);

    M = side_length;
    h = 1.0 / side_length;
    time_steps = time_delta;
    V = arma::cx_mat(side_length - 2, side_length - 2);

    filename = "Data_M_" + std::to_string(side_length) + "_dt_" + std::to_string(time_delta) + ".bin";
    fill_matrices();
}

void Solver::fill_matrices()
{
    A_matrix.fill_matrix(V, 1);
    B_matrix.fill_matrix(V, 0);
}

// Sets inital state of system with Dirichlet boundary conditions
// and normalises all grid squares to have absolute value 1
void Solver::set_initial_state(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y)
{
    double real_comp, imag_comp;
    double sigma_x_squared_2 = 2 * sigma_x * sigma_x;
    double sigma_y_squared_2 = 2 * sigma_y * sigma_y;

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
    // std::cout << "NOrmal:" << normalization_factor << "\n";
    current_state /= sqrt(normalization_factor);
}

// Converts current state to vector
arma::cx_vec Solver::convert_current_state_to_vector()
{
    arma::cx_vec ans = arma::cx_vec((M - 2) * (M - 2));

    for (int y = 1; y < M - 1; ++y)
    {
        for (int x = 1; x < M - 1; ++x)
        {
            ans[indicies_to_index(y, x)] = current_state(y, x);
        }
    }
    return ans;
}

void Solver::update_current_state(arma::cx_vec v)
{
    for (int y = 1; y < M - 1; ++y)
    {
        for (int x = 1; x < M - 1; ++x)
        {
            current_state(y, x) = v[indicies_to_index(y, x)];
        }
    }
}

// Translates pair of indicies to single index
int Solver::indicies_to_index(int y, int x)
{
    return y - 1 + (M - 2) * (x - 1);
}

arma::cx_vec Solver::compute_b()
{
    return B_matrix.multiply_matrix_with_vector(convert_current_state_to_vector());
}

// Finds next state using Crank-Nicolson
void Solver::find_next_state()
{

    arma::cx_vec b = compute_b();

    arma::cx_vec ans;
    arma::spsolve(ans, A_matrix.get_matrix(), b, "superlu");
    update_current_state(ans);
}

// Prints the current state of the system
void Solver::print_current()
{
    current_state.raw_print();
}

// Simulates the system
void Solver::simulate()
{
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