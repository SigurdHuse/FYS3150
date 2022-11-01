#include "system.hpp"

System::System(int l)
{
    grid = arma::Mat<int>(l, l);
}

int System::compute_magnetisation()
{
    int ans = 0;
    return arma::sum(arma::sum(grid));
}