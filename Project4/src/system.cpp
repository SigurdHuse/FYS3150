#include "system.hpp"

#include <chrono>

System::System(int length, double temp)
{
    l = length;
    N = length * length;

    T = temp;
    beta = 1. / T;

    grid = arma::Mat<int>(length, length);
    neig = std::vector<std::vector<std::pair<int, int>>>(length, std::vector<std::pair<int, int>>(length));
    delta_E_values = std::vector<double>(17);

    // right neighbor is first, down neighbor is second
    for (int i = 0; i < length - 1; ++i)
    {
        for (int j = 0; j < length - 1; ++j)
        {
            neig[i][j].first = i + 1;
            neig[i][j].second = j + 1;
        }
    }

    for (int i = 0; i < length - 1; ++i)
    {
        neig[length - 1][i].first = i + 1;
        neig[length - 1][i].second = 0;
    }

    for (int i = 0; i < length - 1; ++i)
    {
        neig[i][length - 1].first = 0;
        neig[i][length - 1].second = i + 1;
    }

    neig[length - 1][length - 1].first = 0;
    neig[length - 1][length - 1].second = 0;

    for (int i = -8; i <= 8; i += 4)
    {
        delta_E_values[i + 8] = exp(-beta * i);
        // std::cout << delta_E_values[i + 8] << "\n";
    }

    // seed = std::chrono::system_clock::now().time_since_epoch().count();
    // std::cout << std::chrono::system_clock::now().time_since_epoch().count() << "\n";
}

// Computes the magnetisation of the system by summing over all spins
int System::compute_magnetisation()
{
    int ans = 0;
#pragma omp parallel for collapse(2)
    for (int i = 0; i < l; ++i)
    {
        for (int j = 0; j < l; ++j)
        {
#pragma omp atomic
            ans += grid(i, j);
        }
    }
    return ans;
}

// Computes the energy of the system by taking the product and summing over all neighboring pairs
int System::compute_energy()
{
    int ans = 0;
#pragma omp parallel for collapse(2)
    for (int i = 0; i < l; ++i)
    {
        for (int j = 0; j < l; ++j)
        {
#pragma omp atomic
            ans += grid(i, j) * grid(i, neig[i][j].first);
#pragma omp atomic
            ans += grid(i, j) * grid(neig[i][j].second, j);
        }
    }

    return -ans;
}

// Fills the grid by generating 0's and 1's from a uniform distribution
void System::fill_with_random_spins()
{
    // Binary generator
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<int> gen(0, 1);

    // Fill grid
    for (int i = 0; i < l; ++i)
    {
        for (int j = 0; j < l; ++j)
        {
            grid(i, j) = 2 * gen(engine) - 1;
        }
    }
}

// Generates a random coordinate in the system
std::pair<int, int> System::generate_random_coordinate()
{
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<int> gen(0, l - 1);
    // Number generator
    std::pair<int, int> ans;
    ans.first = gen(engine);
    ans.second = gen(engine);
    // std::cout << ans.first << " " << ans.second << "\n";
    return ans;
}

// Perform a single MC cycle
void System::one_MC_cycle()
{
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    for (int i = 0; i < N; ++i)
    {
        auto cor = generate_random_coordinate();
        int x = cor.first, y = cor.second;
        int E_before, E_after;
        E_before = grid((y == 0) * l + y - 1, x) * grid(y, x) + grid(y, (x == 0) * l + x - 1) * grid(y, x);
        E_before += grid(y, neig[y][x].first) * grid(y, x) + grid(neig[y][x].second, x);
        grid(y, x) *= -1;
        E_after = grid((y == 0) * l + y - 1, x) * grid(y, x) + grid(y, (x == 0) * l + x - 1) * grid(y, x);
        E_after += grid(y, neig[y][x].first) * grid(y, x) + grid(neig[y][x].second, x) * grid(y, x);

        // std::cout << E_after - E_before << "\n";
        double A = std::min(1., delta_E_values[E_after - E_before + 8]);
        double r = uniform_dist(engine);
        // std::cout << r << "\n";
        // std::cout << A << " " << r << "\n";
        //   std::cout << "Number:" << -2 * (r > A) + 1 << "\n";
        grid(y, x) *= -2 * (r > A) + 1;
    }
}
