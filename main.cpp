#include <iostream>
#include "Complex.h"
#include "Matrix.h"
#include "Files.h"
#include <random>


int main() {
    linalg::Matrix<int> init1 = {{1, 2, 3}, {2, 4, 1}, {6, 1, 5}};

    std::cout << init1.det() << " " << init1.norm() << " " << init1.trace() << "\n" << inv(init1) << "\n";

    std::cout << pow(init1, 4) << "\n";

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

    linalg::Matrix<int> matr1 = {1, 2};
    linalg::Matrix<int> matr2 = {{1, 2}, {2, 3}};

    try {
        std::cout << matr1 + matr2 << "\n";
    }
    catch (std::logic_error& e) {
        std::cout << e.what() << "\n";
    }

    std::cout << "after throw\n";

    return 0;
}
