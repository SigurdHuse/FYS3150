#include "PenningTrap.hpp"

const double ke = 1.38935333 * 1e5, T = 9.64852559 * 1e1;
const double B0 = 9.65 * 1e1, V0 = 9.65 * 1e8, d = 1e4;

int main()
{
    PenningTrap test(B0, V0, d);
    arma::vec r = {10, 10, 10}, v = {-3, 5, 2};
    Particle p(1, 100, r, v);
    test.add_particle(p);
    const int time = 100, n = 10000;
    const double dt = (double)time / n;
    std::ofstream outfile;
    outfile.open("one_particle_n_10000.txt");
    for (int i = 0; i < n; ++i)
    {
        outfile << dt * i << "\n";
        test.write_positions_to_file(outfile, 14, 6);
        test.evolve_forward_Euler(dt);
    }
    outfile << time << "\n";
    test.write_positions_to_file(outfile, 14, 6);
    outfile.close();
}