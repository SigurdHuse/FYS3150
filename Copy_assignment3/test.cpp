#include "PenningTrap.hpp"

using namespace std;

const double ke = 1.38935333 * 1e5, T = 9.64852559 * 1e1;
const double B0 = 9.65 * 1e1, V0 = 9.65 * 1e8, d = 1e4;
const double tol = 1e-8;

const double m = 40.0775;

void test_add_particle()
{
    PenningTrap test(B0, V0, d);
    arma::vec r = {1, 0, 0}, v = {1, 1, 1};
    Particle p(ke, m, r, v);
    test.add_particle(p);
    cout << "add_particle() test passed" << endl;
}

void test_E_field()
{
    PenningTrap test(4, 2, 1.);
    arma::vec expected = {2., 2., -8.};
    arma::vec r = {1.0, 1.0, 2.0};
    arma::vec computed = test.external_E_field(r);
    for (int i = 0; i < 3; ++i)
    {
        assert(abs(expected(i) - computed(i)) < tol);
    }
    cout << "external_E_field() test passed" << endl;
}

void test_B_field()
{
    PenningTrap test1(1., V0, d), test2(3.4, V0, d);
    arma::vec expected1 = {0., 0., 1.}, expected2 = {0., 0., 3.4};
    arma::vec r = {1.0, 1.0, 2.0};
    arma::vec computed1 = test1.external_B_field(r), computed2 = test2.external_B_field(r);
    for (int i = 0; i < 3; ++i)
    {
        assert((expected1(i) - computed1(i)) < tol);
        assert((expected2(i) - computed2(i)) < tol);
    }
    cout << "external_B_field() test passed" << endl;
}

void test_force_particle()
{
    PenningTrap test(B0, V0, d);
    arma::vec r1 = {1., 1., 1.}, r2 = {2., 2., 2.};
    Particle p1(1, 2, r1, r1), p2(4, 1, r2, r2);
    test.add_particle(p1);
    test.add_particle(p2);

    arma::vec expected1 = {-4., -4., -4.};
    expected1 *= ke / (3. * sqrt(3));
    arma::vec expected2 = {1., 1., 1.};
    expected2 *= ke / (3. * sqrt(3));

    arma::vec computed1 = test.force_particle(0, 1);
    arma::vec computed2 = test.force_particle(1, 0);

    for (int i = 0; i < 3; ++i)
    {
        assert(abs(computed1(i) - expected1(i)) < tol);
        assert(abs(computed2(i) - expected2(i)) < tol);
    }
    cout << "force_particle() test passed" << endl;
}

void test_total_force_external()
{
    PenningTrap test(3., 2., 1.);
    arma::vec r = {1., 1., 1.}, v = {-1., 2, 3};
    Particle p(1., 1., r, v);
    test.add_particle(p);
    arma::vec expected = {8., 5., -4.};
    arma::vec computed = test.total_force_external(0);
    for (int i = 0; i < 3; ++i)
    {
        assert(abs(expected(i) - computed(i)) < tol);
    }
    cout << "total_force_external() test passed" << endl;
}

void test_total_force_particles()
{
    PenningTrap test(1., 1., 1.);
    arma::vec r1 = {1., 1., 1.}, v = {-1., 2, 3};
    arma::vec r2 = {2., 2., 2.}, r3 = {-3., -3., 0.};
    Particle p1(1., 1., r1, v), p2(2., 2., r2, v), p3(3., 3., r3, v);

    test.add_particle(p1);
    test.add_particle(p2);
    test.add_particle(p3);

    arma::vec tmp1 = {-2, -2, -2}, tmp2 = {12, 12, 3};
    tmp1 /= 3 * sqrt(3);
    tmp2 /= 33 * sqrt(33);
    arma::vec expected1 = ke * (tmp1 + tmp2);
    tmp1 = {1., 1., 1.}, tmp2 = {15., 15., 6.};
    tmp1 /= 3 * sqrt(3);
    tmp2 /= 54 * sqrt(54);
    arma::vec expected2 = ke * (tmp1 + tmp2);
    tmp1 = {-4., -4., -1.}, tmp2 = {-10., -10., -4.};
    tmp1 /= 33 * sqrt(33);
    tmp2 /= 54 * sqrt(54);
    arma::vec expected3 = ke * (tmp1 + tmp2);

    arma::vec computed1 = test.total_force_particles(0);
    arma::vec computed2 = test.total_force_particles(1);
    arma::vec computed3 = test.total_force_particles(2);

    for (int i = 0; i < 3; ++i)
    {
        assert(abs(expected1(i) - computed1(i)) < tol);
        assert(abs(expected2(i) - computed2(i)) < tol);
        assert(abs(expected3(i) - computed3(i)) < tol);
    }
    cout << "total_force_particles() test passed" << endl;
}

void test_get_total_of_particles_in_trap()
{
    PenningTrap test(1., 1., 100.);
    arma::vec r1 = {1., 1., 1.}, r2 = {100, 100, 100}, v = {-1., 2, 3};
    Particle p1(1, 1, r1, v), p2(1, 1, r2, v);
    assert(test.get_number_of_particles_in_trap() == 0);
    test.add_particle(p1);
    assert(test.get_number_of_particles_in_trap() == 1);
    for (int i = 0; i < 9; ++i)
    {
        test.add_particle(p1);
    }
    assert(test.get_number_of_particles_in_trap() == 10);
    test.add_particle(p2);
    assert(test.get_number_of_particles_in_trap() == 10);
    cout << "get_total_of_particles_in_trap() test passed" << endl;
}

void test_fill_trap()
{
    PenningTrap test(1., 1., 100.);
    test.fill_trap(1, m, 100);
    assert(test.get_number_of_particles_in_trap() == 100);
    cout << "fill_trap() test passed" << endl;
}

int main()
{
    test_add_particle();
    test_E_field();
    test_B_field();
    test_force_particle();
    test_total_force_external();
    test_total_force_particles();
    test_get_total_of_particles_in_trap();
    test_fill_trap();
}
