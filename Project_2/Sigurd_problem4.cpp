#include "Sigurd_Jacobi.hpp"

int main()
{
    const int n = 6;
    arma::mat A = generate_A(n);
    arma::vec eig(n - 1);
    arma::mat eigvectors = arma::mat(n - 1, n - 1, arma::fill::eye);
    bool converged = 0;
    const int maxiter = 1000;
    int iterations = 0;
    jacobi_eigensolver(A, 1e-4, eig, eigvectors, maxiter, iterations, converged);
}