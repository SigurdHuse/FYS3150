#ifndef __grid_hpp__
#define __grid_hpp__

#include <armadillo>
#include <vector>
#include <complex>

class Grid
{
private:
    int M, M_squared;
    double T, h, delta_t;

    // Matrix optimized to contain mainly zeros
    arma::sp_cx_mat matrix;

public:
    // Constructor
    Grid(int side_length, double time, int time_steps);

    // Empty constructor to use in solver.cpp class
    Grid();

    // Translates pair of indicies to single index
    int indicies_to_index(int i, int j);

    // Fills the matrix with the constant r
    void fill_matrix_with_r(std::complex<double> r);

    // Fills matrix with values from a vector
    void fill_matrix_from_vector(arma::cx_vec v);

    // Fills matrix with values to do the Crank-Nicolson method
    void fill_matrix(arma::mat V, bool A_matrix);

    // Prints matrix
    void print_matrix();

    // Multiplies matrix with a vector v
    arma::cx_vec multiply_matrix_with_vector(arma::cx_colvec v);

    // Setters meant to be used with empty constructor in grid.cpp
    void set_side_length(int side_length);

    void set_time(double time);

    void set_time_step(int time_steps);

    // Returns the sp_cx_mat matrix
    arma::sp_cx_mat get_matrix();
};

#endif