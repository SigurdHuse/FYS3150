#include <PenningTrap.hpp>

const double ke = 1.38935333 * 1e5;

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
}

void PenningTrap::add_particle(Particle p_in)
{
    particles_.push_back(p_in);
}

// Loop trough particles and compute E(r)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    arma::vec E = arma::vec(3).fill(0.);
    for (auto p : particles_)
    {
        E += (r - p.r_) / pow(arma::abs(r - p.r_), 3) * p.q_;
    }
    return ke * E;
}