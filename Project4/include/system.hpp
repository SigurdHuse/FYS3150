#ifndef __system_hpp__
#define __system_hpp__

#include <armadillo>

class System
{
    /*
    Class representing the lxl grid in a Ising model
    */
private:
    // Side length of grid
    int l;

    // Temprature of system
    double T;

public:
    // The grid
    arma::Mat<int> grid;

    // Constructor
    System(int l, double T);

    // Computes the energy of the system
    int compute_engergy();

    // Computes the magnetisation of the system
    int compute_magnetisation();

    // Computes the specific heat capacity of the system
    double compute_specific_heat_capacity();

    // Computes the suspectibility of the system
    double compute_suspectibility();
};

#endif