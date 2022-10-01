#include "include/PenningTrap.hpp"

using namespace std;

const double ke = 1.38935333 * 1e5, T = 9.64852559 * 1e1;
const double B0 = 9.65 * 1e1, V0 = 9.65 * 1e8, d = 1e4;

// Need to find weight of singly-charged Calcium ion
const double m = 46;

void test_add_particle()
{
    PenningTrap test(B0, V0, d);
    arma::vec r = {1, 0, 0}, v = {1, 1, 1};
    Particle p(ke, m, r, v);
    test.add_particle(p);
}

int main()
{
    test_add_particle();
}
