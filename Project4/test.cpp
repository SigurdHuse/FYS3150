#include "system.hpp"
#include <assert.h>
#include <vector>
#include <algorithm>

void test_compute_magnetisation()
{
    System test1(5, 10), test2(3, 10);
    test1.grid.fill(1);
    test2.grid.fill(-1);
    assert(test1.compute_magnetisation() == 25);
    assert(test2.compute_magnetisation() == -9);
}

void test_compute_engergy()
{
    std::vector<std::vector<int>> pos = {{1, 1, 1, 1},
                                         {-1, 1, 1, 1},
                                         {-1, -1, 1, 1},
                                         {-1, -1, -1, 1},
                                         {-1, -1, -1, -1}};
    System test(2, 3);
    for (auto v : pos)
    {
        // std::sort(v.begin(), v.end());
        do
        {
            test.grid(0, 0) = v[0];
            test.grid(0, 1) = v[1];
            test.grid(1, 0) = v[2];
            test.grid(1, 1) = v[3];
            test.grid.print();
            std::cout << test.compute_engergy() << "\n";
            /*
            for (int i : v)
                std::cout << i << " ";
            std::cout << "\n";
            */

        } while (std::next_permutation(v.begin(), v.end()));
    }
}

int main()
{
    test_compute_magnetisation();
    test_compute_engergy();
}