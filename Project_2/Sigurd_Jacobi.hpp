#include <bits/stdc++.h>
#include <armadillo>

const double PI = 3.14159265359;

arma::mat generate_A(int n){
    arma::mat A = arma::mat(n,n);
    double h = 1./n;
    double a = -1./h/h, d = 2./h/h;
    A(0,0) = d; A(0,1) = a;
    A(n-1,n-1) = d, A(n-1,n-2) = a;
    for(int i = 1; i < n-1; ++i){
        A(i,i-1) = a; A(i, i+1) = a;
        A(i,i) = d;
    }
    return A;
}

double max_offdiag_symmetric(arma::mat A, int& k, int &l){
    int n = A.n_cols;
    double mx = DBL_MIN;
    //std::cout << mx << std::endl;
    for(int i = 0; i < n; ++i){
        for(int j = i + 1; j < n; ++j){
            if(abs(A(i,j)) > mx){
                mx = abs(A(i,j));
                k = i;
                l = j;
            }
        }
    }
    return mx;
}

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
    assert (k < l);
    double theta;
    if(A(k,k) == A(l,l)){
        if(A(k,l) > 0) theta = PI/4.;
        else theta = -PI/4.;
    }
    else theta = 0.5*atan2(2*A(k,l),A(k,k) - A(l,l));
    R(k,k) = cos(theta); R(l,l) = R(k,k);
    R(k,l) = sin(theta); R(l,k) = -R(k,l);
    A = R.t() * A * R;
    return;
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){
    int k, l;
    double mx = DBL_MAX;
    arma::mat R = arma::mat(A.n_cols,A.n_cols, arma::fill::eye);
    while(mx >= eps){
        mx = max_offdiag_symmetric(A,k,l);
        jacobi_rotate(A, R, k, l);
        eigenvectors = eigenvectors*R;
        R(k,k) = 0; R(l,l) = 0;
        R(k,l) = 0; R(l,k) = 0;
    }
    A.print();
}