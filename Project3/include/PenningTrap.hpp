#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <Particle.hpp>
#include <bits/stdc++.h>

class PenningTrap
{
private:
    double B0_, V0_, d_;
    std::vector<Particle> particles_;

public:
    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    // With or without interactions between particles
    void evolve_RK4(double dt, bool interaction);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt, bool interaction);

    // Writes current positon of particles to file
    void write_positions_to_file(std::ofstream &outfile, int width, int prec);

    // Writes current velocities of particles to file
    void write_velocities_to_file(std::ofstream &outfile, int width, int prec);

    // Right side of ODE
    std::vector<arma::vec> f(double omega0, double omegaz, const arma::vec &deriv, const arma::vec &pos, const arma::vec &forces);
};

#endif