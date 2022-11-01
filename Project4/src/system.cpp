#include "system.hpp"

System::System(int length, double temp)
{
    l = length;
    grid = arma::Mat<int>(length, length);
    T = temp;
}

// Computes the magnetisation of the system by summing over all spins
int System::compute_magnetisation()
{
    int ans = 0;
    return arma::sum(arma::sum(grid));
}

// Computes the energy of the system by taking the product and summing over all neighboring pairs
int System::compute_engergy()
{
    int ans = 0;
    for (int i = 0; i < l - 1; ++i)
    {
        for (int j = 0; j < l - 1; ++j)
        {
            ans += grid(i, j) * grid(i + 1, j);
            ans += grid(i, j) * grid(i, j + 1);
        }
    }

    for (int i = 0; i < l - 1; ++i)
    {
        ans += grid(l - 1, i) * grid(l - 1, i + 1);
        ans += grid(l - 1, i) * grid(0, i);
    }

    for (int i = 0; i < l - 1; ++i)
    {
        ans += grid(i, l - 1) * grid(i + 1, l - 1);
        ans += grid(i, l - 1) * grid(i, 0);
    }

    ans += grid(l - 1, l - 1) * grid(0, l - 1);
    ans += grid(l - 1, l - 1) * grid(l - 1, 0);
    return ans;
}