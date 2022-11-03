#include "system.hpp"
#include "omp.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <chrono>
#include <assert.h>

/**
 * @brief Function performing a Mone Carlo simulation in Ising model with side length l. To speed up the
 * simulation, the function has a specified number of walkers which are run in parallel. Each of these walkers
 * do runs number of cycles and writes the energy and magnetisation of the state to their seperate file. In
 * the end, each of these documents are combined into one where length, temperature and number of runs per
 * walker is specified. To make the result reproducable the seed used for the simulation is also written to the file.
 *
 * @param l Integer specifieng the side length of our model.
 * @param T Double specifieng the temperatur of our system.
 * @param runs Integer number of runs for each walker.
 * @param walkers Integer number of walkers for the simulation.
 * @param seed Unsigned integer seed for the random number generation.
 * @param random_start Bool indiciating if the inital system state should be random or if all spins should be equal
 */

void do_a_MC_simulation(int l, double T, int runs, int walkers, unsigned int seed, bool random_start)
{
#pragma omp parallel for
    for (int i = 0; i < walkers; ++i)
    {
        std::ofstream magnetism_out, energy_out;
        energy_out.open("Energy_states_l_" + std::to_string(l) + "_T_" + std::to_string(T) + "_" +
                        std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt");

        magnetism_out.open("Magnetism_states_l_" + std::to_string(l) + "_T_" + std::to_string(T) + "_" +
                           std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt");
        System cur(l, T);
        cur.engine.seed(seed + i * 1e6);
        if (random_start)
        {
            cur.fill_with_random_spins();
        }
        else
        {
            cur.fill_with_positive();
        }
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

    std::string filename_energy, filename_magnetism;

    filename_energy = "Energy_states_l_" + std::to_string(l) + "_T_" + std::to_string(T) + "_" + std::to_string(runs);
    filename_magnetism = "Magnetism_states_l_" + std::to_string(l) + "_T_" + std::to_string(T) + "_" + std::to_string(runs);

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

    for (int i = 0; i < walkers; ++i)
    {
        std::string filename1 = "Energy_states_l_" + std::to_string(l) + "_T_" + std::to_string(T) + "_" +
                                std::to_string(runs) + "_walker_" + std::to_string(i) + ".txt";

        std::string filename2 = "Magnetism_states_l_" + std::to_string(l) + "_T_" + std::to_string(T) + "_" +
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
    assert(argc == 6);
    int l = atoi(argv[1]), runs = atoi(argv[3]), walkers = atoi(argv[4]);
    double T = atof(argv[2]);
    bool random = atoi(argv[5]);
    unsigned int base_seed = std::chrono::system_clock::now().time_since_epoch().count();
    do_a_MC_simulation(l, T, runs, walkers, base_seed, random);
}