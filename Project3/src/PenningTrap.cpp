#include <PenningTrap.hpp>

const double ke = 1.38935333 * 1e5;

// Constructor for Penning trap
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
    // Time starts at t = 0
    t_ = 0;
    std::vector<arma::vec> particles_;
}

// Constructor for time and positon dependent V_0
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double f, double omega_v)
{
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
    // Time starts at t = 0
    t_ = 0;
    f_ = f;
    omega_v_ = omega_v;
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

// Compute external E field from formula, E = -\delta V with positon dependency
arma::vec PenningTrap::external_E_field(arma::vec r, bool use_distance)
{
    if (use_distance == 0 || arma::norm(r) > d_)
    {
        return external_E_field(r);
    }
    arma::vec ans = {0, 0, 0};
    return V0_ * ans / d_ / d_;
}

// Computes external magnetic field B from definition
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    arma::vec field = {0, 0, B0_};
    return field;
}

// Computes external magnetic field B from definition with positon dependency
arma::vec PenningTrap::external_B_field(arma::vec r, bool use_distance)
{
    if (use_distance == 0 || arma::norm(r) > d_)
    {
        return external_B_field(r);
    }
    arma::vec ans = {0, 0, 0};
    return ans;
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
    return total_force_particles(i) + total_force_external(i);
}

// Right side of ODE
std::vector<arma::vec> PenningTrap::f(double omega0, double omegaz, const arma::vec &deriv, const arma::vec &pos, const arma::vec &forces)
{
    arma::vec double_deriv = {omega0 * deriv(1) + omegaz / 2 * pos(0),
                              -omega0 * deriv(0) + omegaz / 2 * pos(1),
                              -omegaz * pos(2)};
    double_deriv += forces;
    arma::vec vel = deriv;
    return {double_deriv, vel};
}

// Go one step forward in simulation with forward Euler
// The position of each particle is r, and it's derivative is v
void PenningTrap::evolve_forward_Euler(long double dt, bool interaction)
{
    int n = particles_.size();
    std::vector<arma::vec> forces(n);
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
    }
    double omega0, omegaz;
    const double dd = 1. / d_ / d_, V02 = V0_ * 2;
    for (int i = 0; i < n; ++i)
    {
        omega0 = particles_[i].q_ * B0_ / particles_[i].m_;
        omegaz = V02 * particles_[i].q_ / particles_[i].m_ * dd;
        std::vector<arma::vec> tmp = f(omega0, omegaz, particles_[i].v_, particles_[i].r_, forces[i]);
        particles_[i].v_ += dt * tmp[0];
        particles_[i].r_ += dt * tmp[1];
    }
}

// Go one step forward in simulation with forward Euler with time and position dependecy
void PenningTrap::evolve_forward_Euler(long double dt, bool interaction, bool time_interaction)
{
    if (time_interaction == 0)
    {
        evolve_forward_Euler(dt, interaction);
        return;
    }
    int n = particles_.size();
    std::vector<arma::vec> forces(n);
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
    }
    double omega0, omegaz;
    const double dd = 1. / d_ / d_, V02 = V0_ * 2;
    for (int i = 0; i < n; ++i)
    {
        omega0 = get_omega0(particles_[i].q_, particles_[i].m_, particles_[i].r_);
        omegaz = get_omegaz(0., particles_[i].q_, particles_[i].m_, particles_[i].r_);
        std::vector<arma::vec> tmp = f(omega0, omegaz, particles_[i].v_, particles_[i].r_, forces[i]);
        particles_[i].v_ += dt * tmp[0];
        particles_[i].r_ += dt * tmp[1];
    }
    t_ += dt;
}

// Go one step forward in simulation with RK4
void PenningTrap::evolve_RK4(long double dt, bool interaction)
{
    int n = particles_.size();
    std::vector<arma::vec> forces(n);
    double omega0, omegaz;
    const double dd = 1. / d_ / d_;
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
    }
    for (int i = 0; i < n; ++i)
    {
        omega0 = particles_[i].q_ * B0_ / particles_[i].m_;
        omegaz = V0_ * particles_[i].q_ / particles_[i].m_ * dd;
        std::vector<arma::vec> k1 = f(omega0, omegaz, particles_[i].v_, particles_[i].r_, forces[i]);
        k1[0] *= dt;
        k1[1] *= dt;

        std::vector<arma::vec> k2 = f(omega0, omegaz, particles_[i].v_ + k1[0] / 2, particles_[i].r_ + k1[1] / 2, forces[i]);
        k2[0] *= dt;
        k2[1] *= dt;

        std::vector<arma::vec> k3 = f(omega0, omegaz, particles_[i].v_ + k2[0] / 2, particles_[i].r_ + k2[1] / 2, forces[i]);
        k3[0] *= dt;
        k3[1] *= dt;

        std::vector<arma::vec> k4 = f(omega0, omegaz, particles_[i].v_ + k3[0], particles_[i].r_ + k3[1], forces[i]);
        k4[0] *= dt;
        k4[1] *= dt;

        particles_[i].r_ += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.;
        particles_[i].v_ += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.;
    }
}

// Go one step forward in simulation with RK4 with time and position dependecy
void PenningTrap::evolve_RK4(long double dt, bool interaction, bool time_dependent)
{
    if (time_dependent == 0)
    {
        evolve_RK4(dt, interaction);
        return;
    }
    int n = particles_.size();
    std::vector<arma::vec> forces(n);
    double omega0, omegaz;
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
    }
    for (int i = 0; i < n; ++i)
    {
        omega0 = get_omega0(particles_[i].q_, particles_[i].m_, particles_[i].r_);
        omegaz = get_omegaz(0., particles_[i].q_, particles_[i].m_, particles_[i].r_);
        std::vector<arma::vec> k1 = f(omega0, omegaz, particles_[i].v_, particles_[i].r_, forces[i]);
        k1[0] *= dt;
        k1[1] *= dt;

        omega0 = get_omega0(particles_[i].q_, particles_[i].m_, particles_[i].r_ + k1[1] / 2);
        omegaz = get_omegaz(dt / 4., particles_[i].q_, particles_[i].m_, particles_[i].r_ + k1[1] / 2);
        std::vector<arma::vec> k2 = f(omega0, omegaz, particles_[i].v_ + k1[0] / 2, particles_[i].r_ + k1[1] / 2, forces[i]);
        k2[0] *= dt;
        k2[1] *= dt;

        omega0 = get_omega0(particles_[i].q_, particles_[i].m_, particles_[i].r_ + k2[1] / 2);
        omegaz = get_omegaz(dt / 2., particles_[i].q_, particles_[i].m_, particles_[i].r_ + k2[1] / 2);
        std::vector<arma::vec> k3 = f(omega0, omegaz, particles_[i].v_ + k2[0] / 2, particles_[i].r_ + k2[1] / 2, forces[i]);
        k3[0] *= dt;
        k3[1] *= dt;

        omega0 = get_omega0(particles_[i].q_, particles_[i].m_, particles_[i].r_ + k3[1]);
        omegaz = get_omegaz(dt * 3. / 4., particles_[i].q_, particles_[i].m_, particles_[i].r_ + k3[1]);
        std::vector<arma::vec> k4 = f(omega0, omegaz, particles_[i].v_ + k3[0], particles_[i].r_ + k3[1], forces[i]);
        k4[0] *= dt;
        k4[1] *= dt;

        particles_[i].r_ += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.;
        particles_[i].v_ += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.;
    }
    t_ += dt;
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

// Computes V0 when it's time dependent
double PenningTrap::get_V0(double offset)
{
    return V0_ * (1. + f_ * std::cos(omega_v_ * (t_ + offset)));
}

// Computes omega0 when it's time and position dependent
double PenningTrap::get_omega0(double q, double m, arma::vec pos)
{
    if (arma::norm(pos) > d_)
    {
        return 0;
    }
    return q * B0_ / m;
}

// Computes omegaz when it's time and position dependent
double PenningTrap::get_omegaz(double offset, double q, double m, arma::vec pos)
{
    if (arma::norm(pos) > d_)
    {
        return 0;
    }
    return get_V0(offset) * q / m / d_ / d_;
}

// Counts number of particles in the trap (particles which satisfies |r| < d)
int PenningTrap::get_number_of_particles_in_trap()
{
    int ans = 0;
    for (auto p : particles_)
    {
        ans += (arma::norm(p.r_) < d_);
    }
    return ans;
}

// Fill trap with n randomly generated particles with charge q and mass m
void PenningTrap::fill_trap(double q, double m, int n)
{
    for (int i = 0; i < n; ++i)
    {
        arma::vec r = arma::vec(3).randn() * 0.1 * d_; // random initial position
        arma::vec v = arma::vec(3).randn() * 0.1 * d_; // random initial velocity
        Particle p(q, m, r, v);
        particles_.push_back(p);
    }
}