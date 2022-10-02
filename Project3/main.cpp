#include "PenningTrap.hpp"

const double ke = 1.38935333 * 1e5, T = 9.64852559 * 1e1;
const double B0 = 9.65 * 1e1, V0 = 9.65 * 1e8, d = 1e4;

const double m = 40.0775, q = 1.;

const int width = 14, prec = 7;

void one_particle_100_mus(std::string &filename)
{
    // Constructs objects
    PenningTrap test(B0, V0, d);
    arma::vec r = {10., 10., 10.}, v = {-3., 5., 2.};
    Particle p(q, m, r, v);
    test.add_particle(p);

    const int time = 100, n = 2;
    const double dt = (double)time / n;
    std::ofstream outfile;
    outfile.open(filename);

    // Simulates and writes to file
    for (int i = 0; i < n; ++i)
    {
        outfile << dt * i << "\n";
        test.write_positions_to_file(outfile, width, prec);
        test.evolve_RK4(dt, 1);
    }
    outfile << time << "\n";
    test.write_positions_to_file(outfile, width, prec);
    outfile.close();
}

void two_particles(std::string filename, bool interaction)
{
    // Constructs objects
    PenningTrap test(B0, V0, d);
    arma::vec r1 = {10., 10., 10.}, v1 = {-3., 5., 2.};
    arma::vec r2 = {-10., -10., -5.}, v2 = {-7., 10., 3.};
    Particle p1(q, m, r1, v1), p2(q, m, r2, v2);
    test.add_particle(p1);
    test.add_particle(p2);

    const int time = 100, n = 2;
    const double dt = (double)time / n;
    std::ofstream outfile, time_file;
    time_file.open("time_" + filename);
    outfile.open(filename);

    // Simulates and writes to file
    for (int i = 0; i < n; ++i)
    {
        time_file << dt * i << "\n";
        test.write_positions_to_file(outfile, width, prec);
        test.evolve_RK4(dt, interaction);
    }
    time_file << time << "\n";
    test.write_positions_to_file(outfile, width, prec);
    outfile.close();
    time_file.close();
}

int main()
{
    two_particles("two_particles_with_interaction.txt", 1);
    // two_particles("two_particles_without_interaction.txt", 0);
}