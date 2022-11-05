#ifndef __system_hpp__
#define __system_hpp__

#include <armadillo>
#include <utility>
#include <random>
#include <string>

class System
{
    /*
    Class representing the lxl grid in a Ising model
    */
private:
    // Side length of grid
    int L;
    // Area of grid
    int N;

    // Temprature of system
    double T;
    double beta;

    // Values of system
    int energy, magnetism;

public:
    // Random number generator;
    std::mt19937 engine;

    // Distributions
    std::uniform_int_distribution<int> gen_0_l;
    std::uniform_int_distribution<int> gen_int_0_1;
    std::uniform_real_distribution<double> uniform_dist;

    // The grid
    arma::Mat<int> grid;

    // Ratio between p-values in our MC simulation
    std::vector<double> delta_E_values;

    // Neighbors
    std::vector<std::vector<std::pair<int, int>>> neig;

    // Constructor
    System(int length, double T, unsigned int seed, bool random);

    // Computes the energy of the system
    void set_energy();

    // Computes the magnetisation of the system
    void set_magnetisation();

    // Computes the specific heat capacity of the system
    double compute_specific_heat_capacity();

    // Computes the suspectibility of the system
    double compute_suspectibility();

    // Fills the system with randomly generated up and down spins
    void fill_with_random_spins();

    // Perform one MC cycle
    void one_MC_cycle();

    // Fills the grid with only positive spins
    void fill_with_positive();

    // Fills in the neigbhor matrix
    void fill_neig();

    // Returns the magnetism of the system
    int get_magnetism();

    // Returns the energy of the system
    int get_energy();
};

#endif