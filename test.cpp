#include <type_traits>
#include <algorithm>
#include <sstream>
#include <array>
#include <numeric>
#include <limits>
#include <fstream>
#include "tiny_test.hpp"
#include "matrix.h"

using testing::PrettyTest;
using testing::make_test;
using testing::TestGroup;
using namespace std::chrono_literals;

//#define LONG_TEST


template<typename T, typename U>
concept CanMultiply = requires(T first, U second) {
    {first * second};
};

template<typename T, typename U, typename R>
concept MultiplicationResult = requires(T first, U second) {
    {first * second} -> std::same_as<R>;
};

template<typename T, typename U>
concept CanAdd = requires(T first, U second) {
    {first + second};
    {first += second};
};

template<typename T, typename U>
concept CanSubtract = requires(T first, U second) {
    {first - second};
    {first -= second};
};

template<typename T, typename U, typename R>
concept SubtractionResult = CanSubtract<T, U> && requires(T first, U second) {
    {first - second} -> std::same_as<R>;
};

void for_each_index(size_t size, const auto& functor) {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            functor(i, j);
        }
    }
}

TestGroup all_tests[] = {
    {"Matrix",
        make_test<PrettyTest>("static checks", [](auto& test) {
            using M12 = Matrix<1, 2, size_t>;
            using M21 = Matrix<2, 1, size_t>;
            using M22 = Matrix<2, 2, size_t>;
            using M11 = Matrix<1, 1, size_t>;

            test.check(std::same_as<M11, SquareMatrix<1, size_t>>);
            test.check(std::same_as<Matrix<2, 2>, SquareMatrix<2>>);

            // Old way of doing such checks:
            // std::is_same_v<decltype(std::declval<M12>() * std::declval<M21>()), M11>;

            test.check(MultiplicationResult<M12, M21, M11>);
            test.check(MultiplicationResult<M21, M12, M22>);
            test.check(MultiplicationResult<M11, M11, M11>);
            test.check(MultiplicationResult<M22, M22, M22>);
            test.check(requires(M11 matrix) {
                {matrix *= matrix} -> std::same_as<M11&>;       
            });

            test.check(!CanMultiply<M12, M12>);
            test.check(!CanMultiply<M21, M21>);
            
            // In other form to show off requires expression
            test.check(requires(M12 m12, M21 m21, M11 m11, M22 m22) {
                {m12 + m12} -> std::same_as<M12>;
                {m12 += m12} -> std::same_as<M12&>;
                {m21 + m21} -> std::same_as<M21>;
                {m21 += m21} -> std::same_as<M21&>;
                {m11 + m11} -> std::same_as<M11>;
                {m11 += m11} -> std::same_as<M11&>;
                {m22 + m22} -> std::same_as<M22>;
                {m22 += m22} -> std::same_as<M22&>;
            });

            test.check(!CanAdd<M12, M11>);
            test.check(!CanAdd<M12, M21>);
            test.check(!CanAdd<M11, M22>);

            test.check(SubtractionResult<M12, M12, M12>);
            test.check(SubtractionResult<M21, M21, M21>);
            test.check(SubtractionResult<M11, M11, M11>);
            test.check(SubtractionResult<M22, M22, M22>);

            test.check(!CanSubtract<M12, M11>);
            test.check(!CanSubtract<M12, M21>);
            test.check(!CanSubtract<M11, M22>);
        }),

        make_test<PrettyTest>("operations", [](auto& test) {
            using field_t = Rational;
            SquareMatrix<3, field_t> matrix = {
                {1, 2, 3},
                {4, 5, 6},
                {7, 8, 9}
            };

            test.equals(matrix.det(), field_t(0));
            
            {
                auto transposed = matrix.transposed();
                for_each_index(3, [&](size_t i, size_t j) {
                    test.equals(transposed[i][j], matrix[j][i]);
                });
            }
            
            test.equals(matrix.rank(),  size_t(2));
            test.equals(matrix.trace(), field_t(15));
            
            auto&& row = matrix.getRow(1);
            field_t row_true[] = {4, 5, 6};
            test.check(std::equal(row.begin(), row.end(), row_true));

            auto&& col = matrix.getColumn(1);
            field_t col_true[] = {2, 5, 8};
            test.check(std::equal(col.begin(), col.end(), col_true));
        })
#ifdef LONG_TEST
        , testing::make_timed_test<testing::PrettyTest>("perfomance", [](auto& test){
            //TODO: extract matrix parsing from here
            constexpr size_t size = 20;
            SquareMatrix<size> test_matrix;
            SquareMatrix<size, double> correct_result;
            
            std::ifstream input("data/matr.txt");

            auto read_matrix = [&input](auto& matrix) {
                for_each_index(size, [&](size_t i, size_t j) {
                    input >> matrix[i][j];
                });
            };

            read_matrix(test_matrix);
            read_matrix(correct_result);
            auto inv = test_matrix.invert();

            for_each_index(size, [&](size_t i, size_t j) {
                auto diff = std::abs(double(inv[i][j]) - correct_result[i][j]);
                auto epsilon = 1e-5;
                test.check(diff < epsilon);
            });
        })
#endif // LONG_TEST
    }
};


int main() {
    bool success = true;
    for (auto& group : all_tests) {
        success &= group.run();
    }
    return success ? 0 : 1;
}

