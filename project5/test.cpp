#include "grid.hpp"
#include "solver.hpp"

void test_fill_matrix_from_vector()
{
    int M = 5;
    Grid A(M, 1, 100), B(M, 1, 100);
    std::vector<std::complex<double>> a = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::complex<double> r(0, 1);

    A.fill_matrix_from_vector(a);
    B.fill_matrix_from_vector(a);

    A.fill_matrix_with_r(r);
    B.fill_matrix_with_r(r);

    std::cout << "A matrix: \n";
    A.print_matrix();

    std::cout << "\nB matrix: \n";
    B.print_matrix();
}

void test_fill_matrix()
{
    int M = 5;
    double r = 1;
    Grid A(M, 1, 100), B(M, 1, 100);

    double delta_t = 0.1, h = 0.1;
    arma::cx_mat V = arma::cx_mat(M - 2, M - 2);

    A.fill_matrix(V, 1);
    B.fill_matrix(V, 0);
    std::cout << "A matrix: \n";
    A.print_matrix();
    std::cout << "\nB matrix: \n";
    B.print_matrix();
}

void test_solver()
{
    Solver solver(200, 0.008, 320);
    solver.set_initial_state(0.25, 0.5, 0.05, 0.05, 200, 0);

    // solver.print_current();
    std::cout << solver.compute_sum_of_probabilities() << "\n";
    solver.simulate();
    // solver.print_current();
    std::cout << solver.compute_sum_of_probabilities() << "\n";

    // std::complex<double> test(2.2, 2.1);
    //  std::cout << abs(test) << "\n";
    //  std::cout << norm(test) << "\n";
}

int main()
{
    test_fill_matrix();
    // test_solver();
}