#include "Sigurd_Jacobi.hpp"

using namespace std;

void write_to_file(ofstream &outfile, int precision, int width, double val1, double val2)
{
    outfile << setw(width) << setprecision(precision) << scientific << val1 << ", "
            << setw(width) << setprecision(precision) << scientific << val2 << endl;
}

int main()
{
    const int number_of_values = 6, runs = 10;
    vector<int> n_values(number_of_values, 4);
    vector<double> values1(number_of_values, 0), values2(number_of_values, 0);
    for (int i = 1; i < number_of_values; ++i)
    {
        n_values[i] = n_values[i - 1] * 2;
    }
    ofstream outfile;
    outfile.open("run_times_problem_5_not_dense.txt");
    const int width = 14, prec = 6;

    for (int _ = 0; _ < runs; _++)
    {
        for (int i = 0; i < number_of_values; ++i)
        {
            int n = n_values[i];
            arma::mat A = arma::mat(n - 1, n - 1).randn();
            A = arma::symmatu(A);
            arma::vec eigval(n - 1);
            arma::mat eigvectors = arma::mat(n - 1, n - 1, arma::fill::eye);
            bool converged = 0;
            const int maxiter = 1000;
            int iterations = 0;
            auto t1 = std::chrono::high_resolution_clock::now();
            jacobi_eigensolver(A, 1e-4, eigval, eigvectors, maxiter, iterations, converged);
            auto t2 = std::chrono::high_resolution_clock::now();
            double duration_seconds = std::chrono::duration<double>(t2 - t1).count();
            values1[i] += duration_seconds;
        }
    }

    for (int i = 0; i < number_of_values; i++)
    {
        write_to_file(outfile, prec, width, n_values[i], values1[i] / runs);
    }
    outfile.close();
    outfile.open("run_times_problem_5_dense.txt");

    for (int _ = 0; _ < runs; _++)
    {
        for (int i = 0; i < number_of_values; ++i)
        {
            int n = n_values[i];
            arma::mat A = generate_A(n);
            arma::vec eigval(n - 1);
            arma::mat eigvectors = arma::mat(n - 1, n - 1, arma::fill::eye);
            bool converged = 0;
            const int maxiter = 1000;
            int iterations = 0;
            auto t1 = std::chrono::high_resolution_clock::now();
            jacobi_eigensolver(A, 1e-4, eigval, eigvectors, maxiter, iterations, converged);
            auto t2 = std::chrono::high_resolution_clock::now();
            double duration_seconds = std::chrono::duration<double>(t2 - t1).count();
            values2[i] += duration_seconds;
        }
    }
    for (int i = 0; i < number_of_values; i++)
    {
        write_to_file(outfile, prec, width, n_values[i], values2[i] / runs);
    }
    outfile.close();
}