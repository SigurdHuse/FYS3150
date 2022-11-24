#ifndef __solver_hpp__
#define __solver_hpp__

#include "grid.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

class Solver
{
private:
    Grid A_matrix;
    Grid B_matrix;

    std::string filename;
    int M, M_squared, time_steps;
    double h;
    long double v_0;
    arma::cx_mat current_state;
    arma::mat V;
    arma::cx_cube states;
    arma::cx_colvec current_state_vec;

public:
    // Constructor
    Solver(int side_length, double time, int time_delta, long double v0);

    // Computes b component of Crank-Nicolson method
    arma::cx_vec compute_b();

    // Sets inital state of system with Dirichlet boundary conditions
    void set_initial_state(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);

    // Converts current state to vector
    void convert_current_state_to_vector();

    // Translates pair of indicies to single index
    int indicies_to_index(int y, int x);

    // Print current state
    void print_current();

    // Finds next state using Crank-Nicolson
    void find_next_state();

    // Updates current state from vector
    void update_current_state(arma::cx_vec &v);

    // Simulates system using Crank-Nicolson
    void simulate();

    double compute_sum_of_probabilities();

    // Fills A_matrix and B_matrix
    void fill_matrices();

    // Writes current state to file
    void write_states_to_file();

    // Prints current state of system as vector
    void print_current_state_vector();

    // Initialises the potential V
    void initialise_V();
};

#endif