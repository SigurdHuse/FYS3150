#include "PenningTrap.hpp"

const double ke = 1.38935333 * 1e5, T = 9.64852559 * 1e1;
const double B0 = 9.65 * 1e1, V0 = 2.41 * 1e6, d = 500;

const double m = 40.0775, q = 1.;

const int width = 20, prec = 9;

void one_particle_50_mus(std::string filename)
{
    // Constructs objects
    PenningTrap test(B0, V0, d);
    arma::vec r = {20., 0, 20.}, v = {0, 25., 0};
    Particle p(q, m, r, v);
    test.add_particle(p);

    const int time = 50, n = 100000;
    const double dt = (double)time / n;
    std::ofstream outfile, time_file;
    time_file.open("time_" + filename);
    outfile.open(filename);

    // Simulates and writes to file
    for (int i = 0; i < n; ++i)
    {
        time_file << dt * i << "\n";
        test.write_positions_to_file(outfile, width, prec);
        test.evolve_RK4(dt, 1);
    }
    time_file << time << "\n";
    test.write_positions_to_file(outfile, width, prec);
    time_file.close();
    outfile.close();
}

void two_particles(std::string filename, bool interaction, bool vel)
{
    // Constructs objects
    PenningTrap test(B0, V0, d);
    arma::vec r1 = {20., 0., 20.}, v1 = {0., 25., 0.};
    arma::vec r2 = {25., 25., 0.}, v2 = {0., 40., 5.};
    Particle p1(q, m, r1, v1), p2(q, m, r2, v2);
    test.add_particle(p1);
    test.add_particle(p2);

    const int time = 30, n = 100000;
    const double dt = (double)time / n;
    std::ofstream outfile, time_file;
    time_file.open("time_" + filename);
    outfile.open(filename);

    // Simulates and writes to file
    for (int i = 0; i < n; ++i)
    {
        time_file << dt * i << "\n";
        if (vel)
        {
            test.write_velocities_to_file(outfile, width, prec);
        }
        else
        {
            test.write_positions_to_file(outfile, width, prec);
        }
        test.evolve_RK4(dt, interaction);
    }
    time_file << time << "\n";
    test.write_positions_to_file(outfile, width, prec);
    outfile.close();
    time_file.close();
}

void one_particle_different_h(std::string filename, std::vector<int> nvals, bool euler)
{
    // Constructs objects

    const int time = 50;
    for (int n : nvals)
    {
        PenningTrap test(B0, V0, d);
        arma::vec r = {10., 0, 10.}, v = {0, 5., 0};
        Particle p(q, m, r, v);
        test.add_particle(p);
        const double h = (double)time / n;
        std::ofstream outfile;
        std::string file = filename + "_n_" + std::to_string(n) + "_method_";

        if (euler)
        {
            file += "euler";
        }
        else
        {
            file += "RK4";
        }

        file += ".txt";
        outfile.open(file);

        // Simulates and writes to file
        for (int i = 0; i < n; ++i)
        {
            test.write_positions_to_file(outfile, width, prec);
            if (euler)
            {
                test.evolve_forward_Euler(h, 1);
            }
            else
            {
                test.evolve_RK4(h, 1);
            }
        }
        test.write_positions_to_file(outfile, width, prec);
        outfile.close();
    }
}

void one_particle_100_mus_test(std::string filename)
{
    // Constructs objects
    PenningTrap test(B0, V0, d, 0.5, 1000.);
    arma::vec r = {10., 0, 10.}, v = {0, 5., 0};
    Particle p(q, m, r, v);
    test.add_particle(p);

    const int time = 100, n = 100000;
    const long double dt = (long double)time / n;
    std::ofstream outfile, time_file;
    time_file.open("time_" + filename);
    outfile.open(filename);

    // Simulates and writes to file
    for (int i = 0; i < n; ++i)
    {
        time_file << dt * i << "\n";
        test.write_positions_to_file(outfile, width, prec);
        test.evolve_RK4(dt, 1, 1);
    }
    time_file << time << "\n";
    test.write_positions_to_file(outfile, width, prec);
    time_file.close();
    outfile.close();
}

void hundred_particles_time_dependent(double f, std::string filename)
{
    double length = 500, V = 241250, time = 500;
    int steps = 500000;
    const double M = 1e6, dt = time / steps;

    std::ofstream outfile;
    outfile.open(filename);

    for (double omega_V = 0.2; omega_V <= 2.5; omega_V += 0.02)
    {
        PenningTrap trap(B0, length, V, f, omega_V * M);
        trap.fill_trap(1, m, 10);
        outfile << std::setw(width) << std::setprecision(prec) << std::scientific << omega_V << std::endl;
        for (int i = 0; i < steps; ++i)
        {
            trap.evolve_RK4(dt, 0, 1);
        }
        outfile << std::setw(width) << std::setprecision(prec) << std::scientific << trap.get_number_of_particles_in_trap() << std::endl;
    }
}

int main()
{
    two_particles("two_particles_with_interaction.txt", 1, 0);
    two_particles("two_particles_without_interaction.txt", 0, 0);
    two_particles("two_particles_with_interaction_vel.txt", 1, 1);
    two_particles("two_particles_without_interaction_vel.txt", 0, 1);
    one_particle_50_mus("one_particle_n_10000.txt");
    std::vector<int> nvals = {4000, 8000, 16000, 32000};
    one_particle_different_h("one_particle", nvals, 0);
    one_particle_different_h("one_particle", nvals, 1);
    one_particle_100_mus_test("test.txt");
    // hundred_particles_time_dependent(0.1, "hundred_particles_f_0.1.txt");
    //  hundred_particles_time_dependent(0.4, "hundred_particles_f_0.4.txt");
    //  hundred_particles_time_dependent(0.7, "hundred_particles_f_0.7.txt");
}