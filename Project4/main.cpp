#include <armadillo>
#include <random>

int main()
{
    auto gen = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
    arma::Mat<int> tmp = arma::Mat<int>(10, 10);
    for (int i = 0; i < 10; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            tmp(i, j) = 2 * (gen() == 1) - 1;
        }
    }

    tmp.print();
}