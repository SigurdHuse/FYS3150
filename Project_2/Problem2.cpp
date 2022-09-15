#include <bits/stdc++.h>
#include <armadillo>

using namespace std;

const double PI = 3.14159265359;

arma::mat generate_A(int n)
{
    arma::mat A = arma::mat(n - 1, n - 1);
    double h = 1. / n;
    double a = -1. / h / h, d = 2. / h / h;
    A(0, 0) = d;
    A(0, 1) = a;
    A(n - 2, n - 2) = d, A(n - 2, n - 3) = a;
    for (int i = 1; i < n - 2; ++i)
    {
        A(i, i - 1) = a;
        A(i, i + 1) = a;
        A(i, i) = d;
    }
    return A;
}

pair<arma::mat, arma::vec> solve(arma::mat A)
{
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    int n = eigvec.n_cols;

    for (int i = 0; i < n; ++i)
    {
        eigvec.col(i) = arma::normalise(eigvec.col(i));
    }

    return make_pair(eigvec, eigval);
}

arma::vec generate_v(int n, double i)
{
    arma::vec v(n - 1);
    for (int j = 1; j < n; ++j)
    {
        v[j - 1] = sin((j * i * PI) / (n));
    }
    return arma::normalise(v);
}

int main()
{
    const int n = 6;
    auto A = generate_A(n);
    auto ans = solve(A);
    arma::mat eigenvectors = ans.first;
    arma::vec eigenvalues = ans.second;
    arma::mat analytic_eigenvectors(n - 1, n - 1);
    arma::vec analytic_eigenvalues(n - 1);
    for (int i = 1; i < n; ++i)
    {
        analytic_eigenvectors.col(i - 1) = generate_v(n, i);
    }
    const double h = 1. / n;
    const double d = 2. / h / h, a = -1. / h / h;
    for (int i = 1; i < n; ++i)
    {
        analytic_eigenvalues(i - 1) = d + 2 * a * cos(i * PI / (n));
    }
    eigenvectors.print();
    cout << endl;
    analytic_eigenvectors.print();
    cout << endl;
    eigenvalues.raw_print();
    cout << endl;
    analytic_eigenvalues.raw_print();
    cout << endl;
}