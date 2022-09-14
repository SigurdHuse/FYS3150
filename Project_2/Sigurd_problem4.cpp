#include "Sigurd_Jacobi.hpp"

int main(){
    const int n = 6;
    arma::mat A = generate_A(n);
    arma::vec eig(n);
    arma::mat eigvectors = arma::mat(n,n,arma::fill::eye);
    bool converged = 0;
    const int maxiter = 1000;
    int iterations = 0;
    jacobi_eigensolver(A, 1e-8, eig, eigvectors, maxiter, iterations, converged);
}