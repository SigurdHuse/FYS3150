#include <PenningTrap.hpp>

const double ke = 1.38935333 * 1e5;

// Constructor for Penning trap
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
    std::vector<arma::vec> particles_;
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
    (p.q_ * external_E_field(p.r_)).print();
    std::cout << "electric\n";
    (p.q_ * arma::cross(p.v_, external_B_field(p.r_))).print();
    std::cout << "B\n";
    return p.q_ * external_E_field(p.r_) + p.q_ * arma::cross(p.v_, external_B_field(p.r_));
}

// Compute the total force from the other n particles on particle i
// by using force_particle(i,j) for j = 0,1,...,i-1,i+1,...,n.
arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec ans = arma::vec(3).fill(0.);
    Particle p = particles_[i];
    int n = particles_.size();
    for (int j = 0; j < i; ++j)
    {
        ans += force_particle(i, j);
    }
    for (int j = i + 1; j < n; ++j)
    {
        ans += force_particle(i, j);
    }
    return ans / p.m_ * p.q_;
}

// Computes sum of force from fields and other particles
arma::vec PenningTrap::total_force(int i)
{
    total_force_external(i).print();
    std::cout << "\n";
    return total_force_particles(i) + total_force_external(i);
}

// Right side of ODE without force applied
arma::vec PenningTrap::f(double omega0, double omegaz, arma::vec deriv, arma::vec pos)
{
    arma::vec ans = {omega0 * deriv(1) + omegaz / 2 * pos(0),
                     -omega0 * deriv(0) + omegaz / 2 * pos(1),
                     -omegaz * pos(2)};
    return ans;
}

// Go one step forward in simulation with forward Euler
// The position of each particle is r, and it's derivative is v
void PenningTrap::evolve_forward_Euler(double dt, bool interaction)
{
    int n = particles_.size();
    std::vector<arma::vec> forces(n), vs(n);
    for (int i = 0; i < n; ++i)
    {
        if (interaction)
        {
            forces[i] = total_force_particles(i);
        }
        else
        {
            forces[i] = {0, 0, 0};
        }
        vs[i] = particles_[i].v_;
    }
    double omega0, omegaz;
    const double dd = 1. / d_ / d_, V02 = V0_ * 2;
    for (int i = 0; i < n; ++i)
    {
        omega0 = particles_[i].q_ * B0_ / particles_[i].m_;
        omegaz = V02 * particles_[i].q_ / particles_[i].m_ * dd;
        arma::vec tmp = f(omega0, omegaz, particles_[i].v_, particles_[i].r_);
        tmp += forces[i];
        particles_[i].v_ += dt * tmp;
    }

    for (int i = 0; i < n; ++i)
    {
        particles_[i].r_ += dt * vs[i];
    }
}

void PenningTrap::evolve_RK4(double dt, bool interaction)
{
    int n = particles_.size();
    std::vector<arma::vec> forces(n), vs(n);
    double omega0, omegaz;
    const double dd = 1. / d_ / d_, V02 = V0_ * 2;
    for (int i = 0; i < n; ++i)
    {
        if (interaction)
        {
            forces[i] = total_force_particles(i);
        }
        else
        {
            forces[i] = {0, 0, 0};
        }
        vs[i] = particles_[i].v_;
    }
    for (int i = 0; i < n; ++i)
    {
        omega0 = particles_[i].q_ * B0_ / particles_[i].m_;
        omegaz = V02 * particles_[i].q_ / particles_[i].m_ * dd;
        arma::vec k1_r = (f(omega0, omegaz, particles_[i].v_, particles_[i].r_) + forces[i]) * dt;
        arma::vec k1_v = vs[i] * dt;

        arma::vec k2_r = (f(omega0, omegaz, particles_[i].v_ + k1_r / 2, particles_[i].r_ + k1_v / 2) + forces[i]) * dt;
        arma::vec k2_v = (vs[i] + k1_v / 2) * dt;

        arma::vec k3_r = (f(omega0, omegaz, particles_[i].v_ + k2_r / 2, particles_[i].r_ + k2_v / 2) + forces[i]) * dt;
        arma::vec k3_v = (vs[i] + k2_v / 2) * dt;

        arma::vec k4_r = (f(omega0, omegaz, particles_[i].v_ + k3_r, particles_[i].r_ + k3_v) + forces[i]) * dt;
        arma::vec k4_v = (vs[i] + k3_v) * dt;

        particles_[i].r_ += (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6.;
        particles_[i].v_ += (k1_r + 2 * k2_r + 2 * k3_r + k4_r) / 6.;
    }
}

// Writes all x-positions to first row, all y-positions to second row and all z-positons to third row
void PenningTrap::write_positions_to_file(std::ofstream &outfile, int width, int prec)
{
    int n = particles_.size();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            outfile << std::setw(width) << std::setprecision(prec) << std::scientific << particles_[j].r_(i) << ", ";
        }
        outfile << std::endl;
    }
}

// Writes all x-velocities to first row, all y-velocities to second row and all z-velocities to third row
void PenningTrap::write_velocities_to_file(std::ofstream &outfile, int width, int prec)
{
    int n = particles_.size();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            outfile << std::setw(width) << std::setprecision(prec) << std::scientific << particles_[j].v_(i) << ", ";
        }
        outfile << std::endl;
    }
}