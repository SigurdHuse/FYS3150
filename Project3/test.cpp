#include "PenningTrap.hpp"

using namespace std;

const double ke = 1.38935333 * 1e5, T = 9.64852559 * 1e1;
const double B0 = 9.65 * 1e1, V0 = 9.65 * 1e8, d = 1e4;
const double tol = 1e8;

// Need to find weight of singly-charged Calcium ion
const double m = 46;

void test_add_particle()
{
    PenningTrap test(B0, V0, d);
    arma::vec r = {1, 0, 0}, v = {1, 1, 1};
    Particle p(ke, m, r, v);
    test.add_particle(p);
}

void test_E_field()
{
    PenningTrap test(4, V0, d);
    arma::vec expected = {0., 0., 8.};
    expected *= ke;
    arma::vec r = {1.0, 1.0, 2.0};
    arma::vec computed = test.external_E_field(r);
    for (int i = 0; i < 3; ++i)
    {
        assert((expected(i) - computed(i)) < tol);
    }
}

void test_B_field()
{
    PenningTrap test(B0, V0, d);
    arma::vec expected = {-1., -1., 4.};
    expected *= ke;
    arma::vec r = {1.0, 1.0, 2.0};
    arma::vec computed = test.external_E_field(r);
    for (int i = 0; i < 3; ++i)
    {
        assert((expected(i) - computed(i)) < tol);
    }
}

int main()
{
    test_add_particle();
    test_E_field();
    test_B_field();
}
