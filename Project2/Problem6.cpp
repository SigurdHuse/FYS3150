#include "Sigurd_Jacobi.hpp"

using namespace std;

void write_to_file(ofstream &outfile, const int precision, const int width, vector<double> vals)
{
    int n = vals.size();
    for (int i = 0; i < n - 1; ++i)
    {
        outfile << setw(width) << setprecision(precision) << scientific << vals[i] << ", ";
    }
    outfile << setw(width) << setprecision(precision) << scientific << vals[n - 1] << endl;
}

int main()
{
    vector<int> n_values = {10, 100};
    const int width = 14, prec = 6;
    for (int n : n_values)
    {
        bool converged = 0;
        const int maxiter = 1000;
        int iterations = 0;
        double h = 1. / n;
        ofstream outfile;
        outfile.open("values_problem_6_n_" + to_string(n) + ".txt");
        arma::mat A = generate_A(n);
        arma::vec eigval(n - 1);
        arma::mat eigvec = arma::mat(n - 1, n - 1, arma::fill::eye);
        jacobi_eigensolver(A, 1e-5, eigval, eigvec, maxiter, iterations, converged);
        write_to_file(outfile, prec, width, {0, 0, 0, 0});
        for (int i = 0; i < n - 1; ++i)
        {
            write_to_file(outfile, prec, width, {(i + 1) * h, eigvec(i, 0), eigvec(i, 1), eigvec(i, 2)});
        }
        write_to_file(outfile, prec, width, {1, 0, 0, 0});
        outfile.close();
    }
}