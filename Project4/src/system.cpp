#include "system.hpp"

System::System(int length, double temp, unsigned int seed, bool random)
{
    l = length;
    N = length * length;

    T = temp;
    beta = 1. / temp;

    // Set seed of random number generator
    engine = std::mt19937(seed);

    // Sets parameters of distributions
    gen_0_l.param(std::uniform_int_distribution<int>::param_type(0, length - 1));
    gen_int_0_1.param(std::uniform_int_distribution<int>::param_type(0, 1));
    uniform_dist.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));

    // Grid containing the spins
    grid = arma::Mat<int>(length, length);

    // Neighbors of each spin, with peridodic boundary conditions
    neig = std::vector<std::vector<std::pair<int, int>>>(length, std::vector<std::pair<int, int>>(length));

    // Vector containing the different values of the ratio between probabilites of states
    delta_E_values = std::vector<double>(17);

    // Different values of p(x')/p(x_i)
    for (int i = -8; i <= 8; i += 4)
    {
        delta_E_values[i + 8] = exp(-beta * i);
    }

    // Fills neig matrix
    fill_neig();

    // Fills grid with values
    if (random)
    {
        fill_with_random_spins();
    }
    else
    {
        fill_with_positive();
    }

    // Computes energy of inital state
    set_energy();

    // Computes magnetisation of inital state
    set_magnetisation();
}

void System::fill_neig()
{
    // right neighbor is first, second is down neighbor
    // So first is x-coordinate, while second is y-coordinate
    for (int i = 0; i < l; ++i)
    {
        for (int j = 0; j < l; ++j)
        {
            neig[i][j].first = (j + 1) % l;
            neig[i][j].second = (i + 1) % l;
        }
    }
}

// Computes the magnetisation of the system by summing over all spins
void System::set_magnetisation()
{
    int ans = 0;
    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < l; ++i)
    {
        for (int j = 0; j < l; ++j)
        {
            // #pragma omp atomic
            ans += grid(i, j);
        }
    }
    magnetism = ans;
}

// Computes the energy of the system by taking the product and summing over all neighboring pairs
void System::set_energy()
{
    int ans = 0;
    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < l; ++i)
    {
        for (int j = 0; j < l; ++j)
        {
            // #pragma omp atomic
            ans += grid(i, j) * grid(i, neig[i][j].first);
            // #pragma omp atomic
            ans += grid(i, j) * grid(neig[i][j].second, j);
        }
    }

    energy = -ans;
}

// Fills the grid by generating 0's and 1's from a uniform distribution
void System::fill_with_random_spins()
{
    // Fill grid
    for (int i = 0; i < l; ++i)
    {
        for (int j = 0; j < l; ++j)
        {
            grid(i, j) = 2 * gen_int_0_1(engine) - 1;
        }
    }
}

// Perform a single MC cycle
void System::one_MC_cycle()
{
    for (int i = 0; i < N; ++i)
    {

        int x = gen_0_l(engine), y = gen_0_l(engine);

        int delta_E;
        delta_E = 2 * grid(y, x) * (grid((y - 1 + l) % l, x) + grid(y, (x - 1 + l) % l) + grid(y, neig[y][x].first) + grid(neig[y][x].second, x));

        double A = std::min(1., delta_E_values[delta_E + 8]);
        double r = uniform_dist(engine);

        if (r <= A)
        {
            grid(y, x) *= -1;
            energy += delta_E;
            magnetism += 2 * grid(y, x);
        }
    }
}

// Fills the grid with only positive spins
void System::fill_with_positive()
{
    for (int i = 0; i < l; ++i)
    {
        for (int j = 0; j < l; ++j)
        {
            grid(i, j) = 1;
        }
    }
}

// Returns the magnetism of the system
int System::get_magnetism()
{
    return magnetism;
}

// Returns the energy of the system
int System::get_energy()
{
    return energy;
}