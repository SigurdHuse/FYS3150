#include "system.hpp"
#include "omp.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <chrono>
#include <assert.h>

/**
 * @brief Function performing a Monte Carlo simulation in Ising model with side length L.
 * Writing the energy and magnetisation of each Monte Carlo cycle to a file. Each cycle consists of
 * L x L iterations of the Metropolis rule. To speed up the simulation, the function has a specified
 * number of walkers which are run in parallel. Each of these walkers do runs number of cycles and
 * writes the energy and magnetisation of the state to their seperate file. In the end, each of these
 * documents are combined into one where length, temperature and number of runs per
 * walker is specified. To make the result reproducable the seed used for the simulation is also written to the file.
 *
 * @param L Integer specifieng the side length of our model.
 * @param T Double specifieng the temperatur of our system.
 * @param runs Integer number of runs for each walker.
 * @param walkers Integer number of walkers for the simulation.
 * @param seed Unsigned integer seed for the random number generation.
 * @param random_start Bool indiciating if the inital system state should be random or if all spins should be equal
 */

void do_a_MC_simulation(int L, double T, int runs, int walkers, unsigned int seed, bool random_start)
{
#pragma omp parallel for
    for (int i = 0; i < walkers; ++i)
    {
        std::ofstream magnetism_out, energy_out;
        energy_out.open("Energy_states_L_" + std::to_string(L) + "_T_" + std::to_string(T) + "_" +
                        std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt");

        magnetism_out.open("Magnetism_states_L_" + std::to_string(L) + "_T_" + std::to_string(T) + "_" +
                           std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt");

        // Each walkers seed is different to minimize dependence
        System cur(L, T, seed + i * 1e6, random_start);

        for (int i = 0; i < runs; ++i)
        {
            energy_out << cur.get_energy() << "\n";
            magnetism_out << cur.get_magnetism() << "\n";
            cur.one_MC_cycle();
        }
        energy_out.close();
        magnetism_out.close();
    }

    std::ofstream main_energy, main_magnetism;
    std::ifstream cur_energy, cur_magnetism;

    std::string filename_energy, filename_magnetism;

    filename_energy = "data/Energy_states_L_" + std::to_string(L) + "_T_" + std::to_string(T) + "_" + std::to_string(runs);
    filename_magnetism = "data/Magnetism_states_L_" + std::to_string(L) + "_T_" + std::to_string(T) + "_" + std::to_string(runs);

    if (random_start)
    {
        filename_energy += "_random";
        filename_magnetism += "_random";
    }
    else
    {
        filename_energy += "_positiv";
        filename_magnetism += "_positiv";
    }

    main_energy.open(filename_energy + ".txt");
    main_magnetism.open(filename_magnetism + ".txt");

    main_energy << "Base seed: " << seed << "\n";
    main_magnetism << "Base seed: " << seed << "\n";

    // Writes all the files from the walkers to a main file, and removes the files from the walkers
    for (int i = 0; i < walkers; ++i)
    {
        std::string filename1 = "Energy_states_L_" + std::to_string(L) + "_T_" + std::to_string(T) + "_" +
                                std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt";

        std::string filename2 = "Magnetism_states_L_" + std::to_string(L) + "_T_" + std::to_string(T) + "_" +
                                std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt";

        cur_energy.open(filename1);
        cur_magnetism.open(filename2);

        std::string line;
        while (std::getline(cur_energy, line))
        {
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
    main_energy.close();
    main_magnetism.close();
}

int main(int argc, const char *argv[])
{
    assert(argc == 6 || argc == 8);
    if (argc == 6)
    {
        int L = atoi(argv[1]), runs = atoi(argv[3]), walkers = atoi(argv[4]);
        double T = atof(argv[2]);
        bool random = atoi(argv[5]);
        // std::chrono::system_clock::now().time_since_epoch().count()
        unsigned int base_seed = 1e6;
        do_a_MC_simulation(L, T, runs, walkers, base_seed, random);
    }
    if (argc == 8)
    {
        int L = atoi(argv[1]), runs = atoi(argv[5]), walkers = atoi(argv[6]);
        double T_start = atof(argv[2]), T_end = atof(argv[3]), T_inc = atof(argv[4]);
        bool random = atoi(argv[7]);
        // std::chrono::system_clock::now().time_since_epoch().count()
        unsigned int base_seed = 1e6;
        for (double T = T_start; T <= T_end; T += T_inc)
        {
            do_a_MC_simulation(L, T, runs, walkers, base_seed, random);
        }
    }
}