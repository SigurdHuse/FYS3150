#include <bits/stdc++.h>
#include <armadillo>

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

double max_offdiag_symmetric(const arma::mat A, int &k, int &l)
{
    int n = A.n_cols;
    double mx = DBL_MIN;
    // Loops through each row and column from one off the diagnonal
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            // Checks if the current element is larger than the previous
            if (std::abs(A(i, j)) > mx)
            {
                mx = std::abs(A(i, j));
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
void jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l)
{
    // Compute the angle we need to rotate
    double tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
    double t, c, s;

    if (tau > 0)
        t = 1. / (tau + sqrt(1 + tau * tau));
    else
        t = -1. / (-tau + sqrt(1 + tau * tau));

    c = 1 / sqrt(1 + t * t);
    s = c * t;

    // Update A by performing A_(m+1) = S_m.T * A_(m) * R_(m)
    // We use the fact that A is symmetric to reduce computations
    double tmp = A(k, k);

    A(k, k) = A(k, k) * c * c - 2 * A(k, l) * c * s + A(l, l) * s * s;
    A(l, l) = A(l, l) * c * c + 2 * A(k, l) * c * s + tmp * s * s;
    A(k, l) = 0;
    A(l, k) = 0;

    for (int i = 0; i < (int)A.n_cols; ++i)
    {
        if (i == k || i == l)
        {
            continue;
        }
        tmp = A(i, k);
        A(i, k) = A(i, k) * c - A(i, l) * s;
        A(k, i) = A(i, k);

        A(i, l) = A(i, l) * c + tmp * s;
        A(l, i) = A(i, l);
    }

    // Update R by performing R_(m+1) = R_m * S_m

    for (int i = 0; i < (int)R.n_cols; ++i)
    {
        tmp = R(i, k);
        R(i, k) = R(i, k) * c - R(i, l) * s;
        R(i, l) = R(i, l) * c + tmp * s;
    }
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
void jacobi_eigensolver(arma::mat &A, double eps, arma::vec &eigenvalues, arma::mat &eigenvectors,
                        const int maxiter, int &iterations, bool &converged)
{
    int k, l;
    double mx = max_offdiag_symmetric(A, k, l);
    iterations = 0;
    while (mx >= eps && iterations < maxiter)
    {
        // Applies one Jacobi rotation
        jacobi_rotate(A, eigenvectors, k, l);
        iterations++;
        mx = max_offdiag_symmetric(A, k, l);
    }
    // Check if we converged
    if (iterations == maxiter)
        converged = 0;

    // Normalises the eigenvectors
    for (int i = 0; i < (int)A.n_cols; ++i)
    {
        eigenvalues(i) = A(i, i);
        eigenvectors.col(i) = arma::normalise(eigenvectors.col(i));
    }

    // Insertion sort to get the smallest eigenvalues first
    for (int i = 0; i < (int)A.n_cols; ++i)
    {

        int e = i;
        for (int j = i + 1; j < (int)A.n_cols; ++j)
        {
            if (eigenvalues(j) < eigenvalues(e))
            {
                e = j;
            }
        }
        if (i != e)
        {
            auto tmp = eigenvectors.col(i);
            eigenvectors.col(i) = eigenvectors.col(e);
            eigenvectors.col(e) = tmp;
            std::swap(eigenvalues(i), eigenvalues(e));
        }
    }
}
