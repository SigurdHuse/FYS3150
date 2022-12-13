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

    std::cout << "A matrix filled from vector: \n";
    A.print_matrix();

    std::cout << "\nB matrix filled from vector: \n";
    B.print_matrix();
}

void test_fill_matrix()
{
    int M = 5;
    double r = 1;
    Grid A(M, 1, 100), B(M, 1, 100);

    double delta_t = 0.1, h = 0.1;
    arma::mat V = arma::mat(M - 2, M - 2);

    A.fill_matrix(V, 1);
    B.fill_matrix(V, 0);
    std::cout << "A matrix filled from V: \n";
    A.print_matrix();
    std::cout << "\nB matrix filled from V: \n";
    B.print_matrix();

    int side_length = 6, time_steps = 320;
    double T = 0.008, x_c = 0.25, sigma_x = 0.05, p_x = 200;
    double sigma_y = 0.05, y_c = 0.5, p_y = 0, v0 = 0;

    // arma::cx_mat tmp1(A.get_matrix());
    // tmp1.save("A3.txt", arma::raw_ascii);
    // arma::cx_mat tmp2(B.get_matrix());
    // tmp2.save("B3.txt", arma::raw_ascii);
}

int main()
{
    test_fill_matrix_from_vector();
    test_fill_matrix();
}