#pragma once
#include <cmath>
#include <iomanip>


template <class T = double>
class Complex {
    T re;
    T im;
public:
    Complex(T re = 0, T im = 0) {
        this->re = re;
        this->im = im;
    }

    T real() {
        return re;
    }

    void real(T r) {
        this->re = r;
    }

    T imag() {
        return im;
    }

    void imag(T i) {
        this->im = i;
    }

    Complex<T> & operator=(T other){
        this->re = other;
        this->im = 0;
        return *this;
    }

    Complex<T> & operator+=(T other) {
        this->re += other;
        return *this;
    }

    friend Complex<T> operator+(Complex<T> lhs, T rhs) {
        lhs += rhs;
        return lhs;
    }

    Complex<T> & operator-=(T other) {
        this->re -= other;
        return *this;
    }

    friend Complex<T> operator-(Complex<T> lhs, T rhs) {
        lhs -= rhs;
        return lhs;
    }

    Complex<T> & operator*=(T other) {
        this->re *= other;
        this->im *= other;
        return *this;
    }

    friend Complex<T> operator*(Complex<T> lhs, T rhs) {
        lhs *= rhs;
        return lhs;
    }


    Complex<T> & operator=(Complex<T> other){
        this->re = other.re;
        this->im = other.im;
        return *this;
    }

    Complex<T> & operator+=(Complex<T> other) {
        this->re += other.re;
        this->im += other.im;
        return *this;
    }

    friend Complex<T> operator+(Complex<T> lhs, Complex<T> rhs) {
        lhs += rhs;
        return lhs;
    }

    Complex<T> & operator-=(Complex<T> other) {
        this->re -= other.re;
        this->im -= other.im;
        return *this;
    }

    friend Complex<T> operator-(Complex<T> lhs, Complex<T> rhs) {
        lhs -= rhs;
        return lhs;
    }

    Complex<T> & operator*=(Complex<T> other) {
        T a = this->re;
        T b = this->im;
        T c = other.re;
        T d = other.im;

        // (a + bi)(c + di) = ac - bd + i(ad + bc);
        this->re = a * c - b * d;
        this->im = a * d + b * c;
        return *this;
    }

    friend Complex<T> operator*(Complex<T> lhs, Complex<T> rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend std::ostream & operator<<(std::ostream & out, const Complex<T> & c) {
        if (c.re >= 0) {
            if (c.im >= 0) {
                out << std::scientific << std::setprecision(4) << " " << c.re << " " << "+ " << c.im << " * i ";
            }
            else {
                out << std::scientific << std::setprecision(4) << " " << c.re << " " << "- " << -c.im << " * i ";
            }
        }
        else {
            if (c.im >= 0) {
                out << std::scientific << std::setprecision(4) << "-" << -c.re << " " << "+ " << c.im << " * i ";
            }
            else {
                out << std::scientific << std::setprecision(4) << "-" << -c.re << " " << "- " << -c.im << " * i ";
            }
        }

        return out;
    }

    double abs() {
        return sqrt(re * re + im * im);
    }
};
