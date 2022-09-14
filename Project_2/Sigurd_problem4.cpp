#include "Sigurd_Jacobi.hpp"

int main(){
    const int n = 6;
    arma::mat A = generate_A(n);
    arma::vec eig;
    arma::mat eigvectors = arma::mat(n,n,arma::fill::eye);
    bool converged = 0;
    const int maxiter = 10;
    int iterations = 10;
    jacobi_eigensolver(A, 1e-4, eig, eigvectors, maxiter, iterations, converged);
}