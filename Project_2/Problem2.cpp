#include <bits/stdc++.h>
#include <armadillo>

using namespace std; 

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

arma::mat solve(arma::mat A){
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    int n = eigvec.n_cols;
    for(int i = 0; i < n; ++i){
        eigvec.col(i) = arma::normalise(eigvec.col(i));
    } 
    return eigvec;
}


int main()
{
	auto A = generate_A(6);
    auto ans = solve(A);
    cout << ans << endl;    
}