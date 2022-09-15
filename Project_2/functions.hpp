#include <bits/stdc++.h>
#include <armadillo>

double max_offdiag_symmetric(const arma::mat A, int &k, int &l)
{
    int n = A.n_cols;
    double mx = DBL_MIN;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            if (std::abs(A(i, j)) > mx)
            {
                std::cout << A(i, j) << std::endl;
                mx = std::abs(A(i, j));
                k = i;
                l = j;
            }
        }
    }
    return mx;
}

void test_max_offdiag_symmetric()
{
    arma::mat A = {{1, 0, 0, 0.5}, {0, 1, -0.7, 0}, {0, -0.7, 1, 0}, {0.5, 0, 0, 1}};
    int k, l;
    double ans;
    ans = max_offdiag_symmetric(A, k, l);
    assert(abs(ans - 0.7) < 1e-8);
    assert(k == 1);
    assert(l == 2);
}