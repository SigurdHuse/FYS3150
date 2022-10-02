#include <PenningTrap.hpp>

const double ke = 1.38935333 * 1e5;

// Constructor for Penning trap
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
}

// Adds particle to vector which stores particles in the Penning trap
void PenningTrap::add_particle(Particle p_in)
{
    particles_.push_back(p_in);
}

// Compute external E field from formula, E = -\delta V
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    arma::vec ans = arma::vec(3);
    ans(0) = r(0);
    ans(1) = r(1);
    ans(2) = -r(2) * 2;
    return V0_ * ans / d_ / d_;
}

// Computes external magnetic field B from definition
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    arma::vec field = {0, 0, B0_};
    return field;
}

// Computes force on particle i from particle j
arma::vec PenningTrap::force_particle(int i, int j)
{
    Particle p1 = particles_[i], p2 = particles_[j];
    return ke * p2.q_ * (p1.r_ - p2.r_) / pow(arma::norm(p1.r_ - p2.r_), 3);
}

// Compute the Lorentz force
arma::vec PenningTrap::total_force_external(int i)
{
    Particle p = particles_[i];
    return p.q_ * external_E_field(p.r_) + p.q_ * arma::cross(p.v_, external_B_field(p.r_));
}

// Compute the total force from the other n particles on particle i
// by using force_particle(i,j) for j = 0,1,...,i-1,i+1,...,n.
arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec ans = arma::vec(3).fill(0.);
    int n = particles_.size();
    for (int j = 0; j < i; ++j)
    {
        ans += force_particle(i, j);
    }
    for (int j = i + 1; j < n; ++j)
    {
        ans += force_particle(i, j);
    }
    return ans;
}

// Computes sum of force from fields and other particles
arma::vec PenningTrap::total_force(int i)
{
    return total_force_particles(i) + total_force_external(i);
}