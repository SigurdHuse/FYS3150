#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <Particle.hpp>
#include <bits/stdc++.h>

class PenningTrap
{
private:
    double B0_, V0_, d_, f_, omega_v_;
    bool supposed_to_be_time_dependent;
    long double t_;
    std::vector<Particle> particles_;

public:
    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in);

    // Constructor for time ans positon dependent V_0
    PenningTrap(double B0_in, double V0_in, double d_in, double f, double omega_v);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);

    // External electric field at point r=(x,y,z) with positon dependency
    arma::vec external_E_field(arma::vec r, bool use_distance);

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // External magnetic field at point r=(x,y,z) with positon dependency
    arma::vec external_B_field(arma::vec r, bool use_distance);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    // With or without interactions between particles and with or without time_dependence
    void evolve_RK4(long double dt, bool interaction);
    void evolve_RK4(long double dt, bool interaction, bool time_dependent);

    // Evolve the system one time step (dt) using Forward Euler
    // With or without interactions between particles and with or without time_dependence
    void evolve_forward_Euler(long double dt, bool interaction);
    void evolve_forward_Euler(long double dt, bool interaction, bool time_dependent);

    // Writes current positon of particles to file
    void write_positions_to_file(std::ofstream &outfile, int width, int prec);

    // Writes current velocities of particles to file
    void write_velocities_to_file(std::ofstream &outfile, int width, int prec);

    // Right side of ODE
    std::vector<arma::vec> f(double omega0, double omegaz, const arma::vec &deriv, const arma::vec &pos, const arma::vec &forces);

    // Compute values in ODE when it's time dependent
    double get_V0(double offset);
    double get_omega0(double q, double m, arma::vec pos);
    double get_omegaz(double offset, double q, double m, arma::vec pos);

    // Gets number of particles in the trap
    int get_number_of_particles_in_trap();

    // Fill trap with n randomly generated particles with charge q and mass m
    void fill_trap(double q, double m, int n);
};

#endif