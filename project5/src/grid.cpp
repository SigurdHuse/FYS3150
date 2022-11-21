#include "grid.hpp"

Grid::Grid(int side_length, int dimension, double time, int time_steps)
{
    M = side_length;
    h = 1.0 / side_length;
    delta_t = (double)time / time_steps;
    dim = dimension;
    matrix = arma::cx_mat(dimension * (side_length - 2), dimension * (side_length - 2));
}

// Fills the matrix with the constant r
void Grid::fill_matrix_with_r(std::complex<double> r)
{
    for (int i = 0; i < dim * M - 7; ++i)
    {
        if ((i + 1) % 3 == 0)
        {
            continue;
        }
        matrix(i, i + 1) = r;
        matrix(i + 1, i) = r;
    }
    for (int i = 0; i < dim * M - 9; ++i)
    {
        matrix(i, i + 3) = r;
        matrix(i + 3, i) = r;
    }
}

// Prints matrix
void Grid::print_matrix()
{
    matrix.raw_print();
}

// Translates pair of indicies to single index
int Grid::indicies_to_index(int i, int j)
{
    return (M - 2) * i + j;
}

// Fills matrix with values from a vector
void Grid::fill_matrix_from_vector(std::vector<std::complex<double>> v)
{
    for (int i = 0; i < dim * (M - 2); ++i)
    {
        matrix(i, i) = v[i];
    }
}

// r = i * delta_t / 2 / h / h by definition
void Grid::fill_matrix(arma::cx_mat V, bool A_matrix)
{
    std::vector<std::complex<double>> v((M - 2) * (M - 2));
    std::complex<double> r(0, delta_t / 2 / h / h);
    std::complex<double> r4(0.0, 4 * r.imag());

    if (A_matrix)
    {
        for (int i = 0; i < M - 2; ++i)
        {
            for (int j = 0; j < M - 2; ++j)
            {
                std::complex<double> tmp(V(i, j) * delta_t / 2.0);
                std::complex<double> val(1.0 + tmp.real(), r4.imag() + tmp.imag());

                v[indicies_to_index(i, j)] = val;
            }
        }
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
    }
    fill_matrix_with_r(r);
    fill_matrix_from_vector(v);
}

// arma::cx_colvec Grid::get_next_step(arma::cx_colvec cur)