#include "grid.hpp"

// Constructor
Grid::Grid(int side_length, double time, int time_steps)
{
    M = side_length;
    M_squared = pow(side_length - 2, 2);
    h = 1.0 / side_length;
    delta_t = (double)time / time_steps;
    matrix = arma::sp_cx_mat(M_squared, M_squared);
}

Grid::Grid()
{
    ;
}

void Grid::set_side_length(int side_length)
{
    M = side_length;
    M_squared = pow(side_length - 2, 2);
    matrix = arma::sp_cx_mat((side_length - 2) * (side_length - 2), (side_length - 2) * (side_length - 2));
    h = 1.0 / side_length;
    // std::cout << "h : " << h << "\n";
}

void Grid::set_time(double time)
{
    T = time;
}

void Grid::set_time_step(int time_steps)
{
    delta_t = (double)T / time_steps;
    // std::cout << "delta_t:" << delta_t << "\n";
}

arma::sp_cx_mat Grid::get_matrix()
{
    return matrix;
}

// Fills the matrix with the constant r
void Grid::fill_matrix_with_r(std::complex<double> r)
{
    for (int i = 0; i < M_squared - 1; ++i)
    {
        if ((i + 1) % (M - 2) == 0)
        {
            continue;
        }
        matrix(i, i + 1) = r;
        matrix(i + 1, i) = r;
    }
    int offset = M - 2;
    for (int i = 0; i < M_squared - offset; ++i)
    {
        matrix(i, i + offset) = r;
        matrix(i + offset, i) = r;
    }
}

// Prints matrix
void Grid::print_matrix()
{
    arma::cx_mat tmp(matrix);
    tmp.print();
}

// Translates pair of indicies to single index
int Grid::indicies_to_index(int i, int j)
{
    return i + (M - 2) * j;
}

// Fills matrix with values from a vector
void Grid::fill_matrix_from_vector(arma::cx_vec v)
{
    for (int i = 0; i < M_squared; ++i)
    {
        matrix(i, i) = v[i];
    }
}

// r = i * delta_t / 2 / h / h by definition
void Grid::fill_matrix(arma::mat V, bool A_matrix)
{
    arma::cx_vec v(M_squared);
    std::complex<double> r(0, delta_t / 2 / h / h);
    std::complex<double> r4(0.0, 4 * r.imag());
    // V.raw_print();
    if (A_matrix)
    {
        for (int i = 0; i < M - 2; ++i)
        {
            for (int j = 0; j < M - 2; ++j)
            {
                std::complex<double> tmp(0.0, V(i, j) * delta_t / 2.0);
                std::complex<double> val(1.0 + tmp.real(), r4.imag() + tmp.imag());

                v[indicies_to_index(i, j)] = val;
            }
        }
        fill_matrix_with_r(-r);
    }
    else
    {
        for (int i = 0; i < M - 2; ++i)
        {
            for (int j = 0; j < M - 2; ++j)
            {
                std::complex<double> tmp(V(i, j) * delta_t / 2.0);
                std::complex<double> val(1.0 - tmp.real(), -r4.imag() - tmp.imag());

                v[indicies_to_index(i, j)] = val;
            }
        }
        fill_matrix_with_r(r);
    }
    fill_matrix_from_vector(v);
}

arma::cx_vec Grid::multiply_matrix_with_vector(arma::cx_colvec v)
{
    return matrix * v;
}