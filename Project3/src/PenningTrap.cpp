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

// E = -\delta V
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    arma::vec ans = arma::vec(3);
    ans(0) = -r(0);
    ans(1) = -r(1);
    ans(2) = r(2) * 2;
    return V0_ * ans / d_ / d_;
    /* Will use later
    arma::vec E = arma::vec(3).fill(0.);
    for (auto p : particles_)
    {
        E += (r - p.r_) / pow(arma::abs(r - p.r_), 3) * p.q_;
    }
    return ke * E;
    */
}

arma::vec PenningTrap::external_B_field(arma::vec r)
{
    arma::vec field = {0, 0, B0_};
    return r * field;
}