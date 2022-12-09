#include <iostream>
#include "Complex.h"
#include "Matrix.h"
#include "Files.h"
#include <random>


int main() {
/*
    linalg::Matrix<int> check(1, 1);
    check(0, 0) = 2;
    std::cout << check(0, 0) << "\n";

    const linalg::Matrix<int> check2 = {{1, 2}, {2, 3}};
    std::cout << check2(0, 1) << "\n";
    // check2(0, 1) = 5;

    linalg::Matrix<int> m1 = {{1, 2}, {3, 4}};
    linalg::Matrix<int> m2 = {{3, 4}, {5, 6}};

    linalg::Matrix<int> mul;

    mul = m1 * m2;

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << mul(i, j) << " ";
        }
        std::cout << "\n";
    }

    std::cout << mul.trace() << "\n";

    linalg::Matrix<int> d = {{1, 5, 6}, {4, 5, 3}, {7, 9, 2}};

    std::cout << d.det() << "\n";

    linalg::Matrix<int> t = {{1, 2, 3}, {4, 5, 6}};

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << t(i, j) << " ";
        }
        std::cout << "\n";
    }

    linalg::Matrix<int> t2;
    t2 = transpose(t);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << t2(i, j) << " ";
        }
        std::cout << "\n";
    }

    std::cout << t.norm() << "\n";

    std::cout << t << "\n";

    linalg::Matrix<int> init = {{1, 6, 5}, {4, 2, 7}, {7, 4, 5}};
    linalg::Matrix<int> power = pow(init, 3);

    std::cout << power << "\n";
    std::cout << init * init * init << "\n";

    linalg::Matrix<int> mat = {{1, 6, 5}, {4, 2, 7}, {1, 6, 5}};

    std::cout << init.rank() << "\n" << mat.rank() << "\n";

    std::cout << inv(init) << "\n";

    linalg::Matrix<double> dd = {{1, 6, 5}, {4, 2, 7}, {7, 4, 5}};

    std::cout << dd * inv(init) << "\n";
*/
/*
    linalg::Matrix<Complex<int>> init = {{Complex<int>(1, 1), Complex<int>(6, 6), Complex<int>(5, 5)},
                         {Complex<int>(4, 4), Complex<int>(2, 2), Complex<int>(7, 7)},
                         {Complex<int>(7, 7), Complex<int>(4, 4), Complex<int>(5, 5)}};

    std::cout << init << "\n";
    std::cout << init.det() << "\n" << init.norm() << "\n" << init.trace() << "\n" << init.rank() << "\n" << inv(init) << "\n";
*/
    linalg::Matrix<int> init1 = {{1, 2, 3}, {2, 4, 1}, {6, 1, 5}};

    linalg::TextMode<int>::write("int_matrix.txt", init1);

    linalg::Matrix<int> new_m1 = linalg::TextMode<int>::read("int_matrix.txt");

    std::cout << new_m1 << "\n";

    linalg::BinaryMode<int>::write("int_matrix.bin", init1);

    linalg::Matrix<int> new_m2 = linalg::BinaryMode<int>::read("int_matrix.bin");

    std::cout << new_m2 << "\n";

    linalg::Matrix<Complex<int>> init1c = {{Complex<int>(1, 1), Complex<int>(2, 2), Complex<int>(3, 3)},
                                           {Complex<int>(2, 2), Complex<int>(4, 4), Complex<int>(1, 1)},
                                           {Complex<int>(6, 6), Complex<int>(1, 1), Complex<int>(5, 5)}};

    linalg::TextMode<Complex<int>>::write("intc_matrix.txt", init1c);

    linalg::Matrix<Complex<int>> new_mc1 = linalg::TextMode<Complex<int>>::read("intc_matrix.txt");

    std::cout << new_mc1 << "\n";

    linalg::BinaryMode<Complex<int>>::write("intc_matrix.bin", init1c);

    linalg::Matrix<Complex<int>> new_m3 = linalg::BinaryMode<Complex<int>>::read("intc_matrix.bin");

    std::cout << new_m3 << "\n";

    return 0;
}
