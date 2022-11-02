#include "system.hpp"

int main()
{
    System test(2, 1, 1);
    test.fill_with_random_spins();
    // std::ofstream outfile;
    // outfile.open("energy.txt");
    int runs = 1e6;
    long double ans = test.compute_energy();
    for (int i = 0; i < runs; ++i)
    {
        test.one_MC_cycle();
        // test.grid.raw_print();
        // std::cout << "\n\n";
        ans += test.compute_energy();
        // std::cout << ans << "\n";
    }
    std::cout << ans << "\n";
    std::cout << ans / runs << "\n";
}