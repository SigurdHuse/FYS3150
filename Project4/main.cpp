#include "system.hpp"

#include <iostream>
#include <fstream>
#include <stdio.h>

void do_a_MC_simulation_with_random_start(int l, int T, int runs, int walkers)
{
#pragma omp parallel for
    for (int i = 0; i < walkers; ++i)
    {
        std::ofstream magnetism_out, energy_out;
        energy_out.open("Energy_states_l_" + std::to_string(l) + "_" +
                        std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt");
        magnetism_out.open("Magnetism_states_l_" + std::to_string(l) + "_" +
                           std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt");
        System cur(l, T);
        cur.fill_with_random_spins();

        // energy_out << cur.compute_energy() << "\n";
        // magnetism_out << cur.compute_magnetisation() << "\n";
        for (int i = 0; i < runs; ++i)
        {
            cur.one_MC_cycle();
            energy_out << cur.compute_energy() << "\n";
            magnetism_out << cur.compute_magnetisation() << "\n";
        }
        energy_out.close();
        magnetism_out.close();
    }

    std::ofstream main_energy, main_magnetism;
    std::ifstream cur_energy, cur_magnetism;
    main_energy.open("Energy_states_l_" + std::to_string(l) + "_" + std::to_string(runs) + ".txt");
    main_magnetism.open("Magnetism_states_l_" + std::to_string(l) + "_" + std::to_string(runs) + ".txt");

    int cnt = 0;
    for (int i = 0; i < walkers; ++i)
    {
        std::string filename1 = "Energy_states_l_" + std::to_string(l) + "_" +
                                std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt";
        std::string filename2 = "Magnetism_states_l_" + std::to_string(l) + "_" +
                                std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt";

        cur_energy.open(filename1);
        cur_magnetism.open(filename2);

        std::string line;
        while (std::getline(cur_energy, line))
        {
            cnt++;
            main_energy << line << '\n';
        }

        while (std::getline(cur_magnetism, line))
        {
            main_magnetism << line << '\n';
        }

        cur_energy.close();
        cur_magnetism.close();

        remove(filename1.c_str());
        remove(filename2.c_str());
    }
    std::cout << cnt << "\n";
    main_energy.close();
    main_magnetism.close();
}

int main()
{
    System test(20, 1);
    test.fill_with_random_spins();
    // std::ofstream outfile;
    // outfile.open("energy.txt");
    int T = 1, runs = 1e5, walkers = 10;
    do_a_MC_simulation_with_random_start(20, T, runs, walkers);
}