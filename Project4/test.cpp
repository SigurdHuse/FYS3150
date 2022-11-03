#include "system.hpp"
#include <assert.h>
#include <vector>
#include <algorithm>
#include <armadillo>

// Test compute_magnetisation in System class works as expected
void test_compute_magnetisation()
{
    System test1(5, 1), test2(3, 10);
    test1.grid.fill(1);
    test2.grid.fill(-1);
    assert(test1.compute_magnetisation() == 25);
    assert(test2.compute_magnetisation() == -9);
}

// Test compute_engergy in System class works as expected
void test_compute_engergy()
{
    std::vector<std::vector<int>> pos = {{1, 1, 1, 1},
                                         {-1, 1, 1, 1},
                                         {-1, -1, 1, 1},
                                         {-1, -1, -1, 1},
                                         {-1, -1, -1, -1}};
    std::vector<int> expected = {-8, 0, 0, 0, 0, 0, 0, 8, 8, 0, 0, 0, 0, 0, 0, -8};

    System test(2, 3);
    int cnt = 0;
    for (auto v : pos)
    {
        // std::sort(v.begin(), v.end());
        do
        {
            test.grid(0, 0) = v[0];
            test.grid(0, 1) = v[1];
            test.grid(1, 0) = v[2];
            test.grid(1, 1) = v[3];

            // test.grid.print();
            // std::cout << test.compute_engergy() << " " << expected[cnt] << "\n";

            assert(expected[cnt] == test.compute_energy());
            cnt++;
        } while (std::next_permutation(v.begin(), v.end()));
    }
    // std::cout << expected.size() << "\n";
}

// Test random_coordinate and r in System class works as expected
void test_random_r()
{
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    long double r = 0;
    for (int i = 0; i < 1e6; ++i)
    {
        r += uniform_dist(engine);
    }
    std::cout << r / 1e6 << std::endl;
}

// Test delta_E_values in System class is computed as expected
void test_delta_E_values()
{
    System test1(20, 1), test2(20, 2);
    std::vector<double> expected1 = {exp(8),
                                     exp(4),
                                     exp(0),
                                     exp(-4),
                                     exp(-8)};
    std::vector<double> expected2 = {exp(4), exp(2), exp(0), exp(-2), exp(-4)};
    double tol = 1e-8;
    for (int i = -8, j = 0; i <= 8; i += 4, j++)
    {
        // std::cout << test.delta_E_values[i + 8] << "\n";
        // std::cout << expected[j] << "\n";
        assert(abs(test1.delta_E_values[i + 8] - expected1[j]) < tol);
        assert(abs(test2.delta_E_values[i + 8] - expected2[j]) < tol);
    }
}

int main()
{
    test_compute_magnetisation();
    test_compute_engergy();
    test_random_r();
    test_delta_E_values();
}