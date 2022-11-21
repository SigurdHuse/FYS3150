#ifndef __grid_hpp__
#define __grid_hpp__

#include <armadillo>
#include <vector>
#include <complex>

class Grid
{
private:
    int M, dim;
    double T, h, delta_t;
    // Matrix optimised to contain mainly zeros
    // Should maybe be sp_cx_mat if larger than 100 x 100
    arma::cx_mat matrix;

public:
    // Constructor
    Grid(int side_length, int dimension, double time, int time_steps);

    // Translates pair of indicies to single index
    int indicies_to_index(int i, int j);

    // Fills the matrix with the constant r
    void fill_matrix_with_r(std::complex<double> r);

    // Fills matrix with values from a vector
    void fill_matrix_from_vector(std::vector<std::complex<double>> v);

    // Fills matrix with values to compute Crank-Nicolson
    void fill_matrix(arma::cx_mat V, bool A_matrix);

    // Prints matrix
    void print_matrix();

    // Solves system of equations and gets next time step
    arma::cx_colvec get_next_step(arma::cx_colvec cur);
};

#endif