#include <iostream>
#include "Complex.h"
#include "Matrix.h"
#include <random>


int main() {
/*
    Matrix<int> m(2, 2);
    m[0][0] = std::rand() % 10;
    m[0][1] = std::rand() % 10;
    m[1][0] = std::rand() % 10;
    m[1][1] = std::rand() % 10;

    std::cout << m[0][0] << " " << m[0][1] << "\n" << m[1][0] << " " << m[1][1] << "\n";
    std::cout << m[1][1] << "\n";

    Matrix<int> m1(2, 2);
    m1[0][0] = std::rand() % 10;
    m1[0][1] = std::rand() % 10;
    m1[1][0] = std::rand() % 10;
    m1[1][1] = std::rand() % 10;

    Matrix<int> m2(2, 2);
    m2[0][0] = std::rand() % 10;
    m2[0][1] = std::rand() % 10;
    m2[1][0] = std::rand() % 10;
    m2[1][1] = std::rand() % 10;

    Matrix<int> m_mul = m1 * m2;

    std::cout << m_mul[0][0] << " " << m_mul[0][1] << "\n" << m_mul[1][0] << " " << m_mul[1][1] << "\n";

    Matrix<int> m3 = {{1, 2, 3, 4}, {5, 6, 7, 8}};
    Matrix<int> m4 = {1, 2, 3, 4};

    std::cout << m3[0][1] << "\n" << m4[1][0] << "\n";

    Matrix<int> power = pow(m, 3);

    std::cout << power[0][0] << " " << power[0][1] << "\n" << power[1][0] << " " << power[1][1] << "\n";
    */
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

    return 0;
}
