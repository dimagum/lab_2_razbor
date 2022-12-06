#pragma once

template <class T = double>
class Matrix {
    T * m_ptr;
    unsigned m_rows;
    unsigned m_columns;
public:
    Matrix(unsigned r = 0, unsigned c = 1) {
        m_rows = r;
        m_columns = c;

        m_ptr = new T[r * c];

        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                m_ptr[i * m_columns + j] = 0;
            }
        }
    }
    Matrix(const Matrix<T> & other) {
        m_rows = other.m_rows;
        m_columns = other.m_columns;

        m_ptr = new T[other.m_rows * other.m_columns];

        for (int i = 0; i < other.m_rows; ++i) {
            for (int j = 0; j < other.m_columns; ++j) {
                m_ptr[i * m_columns + j] = other.m_ptr[i * m_columns + j];
            }
        }
    }
    Matrix(Matrix<T> && other) { // && - r-value ссылка;
        m_rows = other.m_rows;
        m_columns = other.m_columns;
        m_ptr = other.m_ptr;

        other.m_ptr = nullptr;
        other.m_rows = 0;
        other.m_columns = 0;
    }
    Matrix(const std::initializer_list<T> & lst) {
        m_rows = lst.size();
        m_columns = 1;
        m_ptr = new T[m_rows * m_columns];
        for (int i = 0; i < m_rows; ++i) {
            m_ptr[i * m_columns + 0] = *(lst.begin() + i);
        }
    }
    Matrix(const std::initializer_list<std::initializer_list<T>> & lst) {
        m_rows = lst.size();
        m_columns = lst.begin()->size();

        m_ptr = new T[m_rows * m_columns];
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_columns; ++j) {
                m_ptr[i * m_columns + j] = *((lst.begin() + i)->begin() + j);
            }
        }

    }
    ~Matrix() {
        delete [] m_ptr;
    }
/*
    class Proxy {
        T* ptr;
        unsigned cols;
    public:
        Proxy(T* other, unsigned c) {
            ptr = other;
            cols = c;
        }
        T & operator[](unsigned idx) {
            if (idx >= cols) {
                throw std::logic_error("index out of range\n");
            }
            return ptr[idx];
        }
    };

    Proxy operator[](unsigned idx) {
        if (idx >= m_rows) {
            throw std::logic_error("index out of range\n");
        }
        return Proxy(&m_ptr[idx], m_columns);
    }
    Proxy operator[](unsigned idx) const {
        if (idx >= m_rows) {
            throw std::logic_error("index out of range\n");
        }
        return Proxy(&m_ptr[idx], m_columns);
    }

    T * operator[](unsigned idx) {
        return m_ptr[idx];
    }
    // m[0][0] = 1;
    // a = m[0][0];

    T * operator[](unsigned idx) const {
        if (idx >= m_rows) {
            throw std::logic_error("index out of range\n");
        }
        return m_ptr[idx];
    }
*/
    // m(1, 1);
    T & operator()(int i, int j) {
        if (i >= m_rows || j >= m_columns) {
            throw std::logic_error("index out of range\n");
        }
        return m_ptr[i * m_columns + j];
    }

    const T & operator()(int i, int j) const {
        if (i >= m_rows || j >= m_columns) {
            throw std::logic_error("index out of range\n");
        }
        return m_ptr[i * m_columns + j];
    }

    Matrix<T> & operator=(const Matrix<T> & rhs) {
        if (this == &rhs) {
            return *this;
        }

        delete [] m_ptr;

        m_rows = rhs.m_rows;
        m_columns = rhs.m_columns;

        m_ptr = new T[m_rows * m_columns];
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_columns; ++j) {
                m_ptr[i * m_columns + j] = rhs.m_ptr[i * m_columns + j];
            }
        }

        return *this;
    }

    Matrix<T> & operator=(Matrix<T> && rhs)  noexcept {
        if (this == &rhs) {
            return *this;
        }

        delete [] m_ptr;

        m_rows = rhs.m_rows;
        m_columns = rhs.m_columns;

        m_ptr = rhs.m_ptr;

        rhs.m_ptr = nullptr;
        rhs.m_rows = 0;
        rhs.m_columns = 0;

        return *this;
    }


    Matrix<T> & operator+=(const Matrix<T> & rhs) {
        if (m_rows != rhs.m_rows || m_columns != rhs.m_columns) {
            throw std::logic_error("dimensions are not equal\n");
        }
        for (int i = 0; i < rhs.m_rows; ++i) {
            for (int j = 0; j < rhs.m_columns; ++j) {
                this->m_ptr[i * m_columns + j] += rhs.m_ptr[i * m_columns + j];
            }
        }
        return *this;
    }

    friend Matrix<T> operator+(Matrix<T> lhs, const Matrix<T> & rhs) {
        lhs += rhs;
        return lhs;
    }
    // m = m1 + m2;

    Matrix<T> & operator-=(const Matrix<T> & rhs) {
        if (m_rows != rhs.m_rows || m_columns != rhs.m_columns) {
            throw std::logic_error("dimensions are not equal\n");
        }
        for (int i = 0; i < rhs.m_rows; ++i) {
            for (int j = 0; j < rhs.m_columns; ++j) {
                this->m_ptr[i * m_columns + j] -= rhs.m_ptr[i * m_columns + j];
            }
        }
        return *this;
    }

    friend Matrix<T> operator-(Matrix<T> lhs, const Matrix<T> & rhs) {
        lhs -= rhs;
        return lhs;
    }

    Matrix<T> & operator*=(T rhs) {
        for (int i = 0; i < rhs.m_rows; ++i) {
            for (int j = 0; j < rhs.m_columns; ++j) {
                this->m_ptr[i * m_columns + j] *= rhs;
            }
        }
        return *this;
    }

    friend Matrix<T> operator*(Matrix<T> lhs, T rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Matrix<T> operator*(T lhs, Matrix<T> rhs) {
        rhs *= lhs;
        return rhs;
    }

    friend Matrix<T> operator*(Matrix<T> lhs, const Matrix<T> & rhs) {
        if (lhs.m_columns != rhs.m_rows) {
            throw std::logic_error("dimensions are not equal\n");
        }
        Matrix<T> tmp(lhs.m_rows, rhs.m_columns);
        for (int i = 0; i < lhs.m_rows; ++i) {
            for (int j = 0; j < rhs.m_columns; ++j) {
                for (int k = 0; k < lhs.m_columns; ++k) {
                    tmp(i, j) += lhs.m_ptr[i * lhs.m_columns + k] * rhs.m_ptr[k * rhs.m_columns + j];
                }
            }
        }
        return tmp;
    }


    double det() {

    }

    double trace() {
        double t = 0;

        for (int i = 0; i < std::min(m_rows, m_columns); ++i) {
            t += m_ptr[i * m_columns + i];
        }

        return t;
    }



    friend Matrix<T> pow(const Matrix<T> & m, unsigned p) {
        Matrix<T> tmp(m);

        for (int i = 0; i < p - 1; ++i) {
            tmp = tmp * m;
        }

        return tmp;
    }

};
