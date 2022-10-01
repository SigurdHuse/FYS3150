#include <PenningTrap.hpp>

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
