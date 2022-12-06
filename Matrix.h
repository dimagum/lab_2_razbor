#pragma once

template <class T = double>
class Matrix {
    T** m_ptr;
    unsigned m_rows;
    unsigned m_columns;
public:
    Matrix(unsigned r = 0, unsigned c = 1) {
        m_rows = r;
        m_columns = c;

        m_ptr = new T*[r];

        for (int i = 0; i < r; ++i) {
            m_ptr[i] = new T[c];
            for (int j = 0; j < c; ++j) {
                m_ptr[i][j] = 0;
            }
        }
    }
    Matrix(const Matrix<T> & other) {
        m_rows = other.m_rows;
        m_columns = other.m_columns;

        m_ptr = new T*[other.m_rows];

        for (int i = 0; i < other.m_rows; ++i) {
            m_ptr[i] = new T[other.m_columns];
            for (int j = 0; j < other.m_columns; ++j) {
                m_ptr[i][j] = other.m_ptr[i][j];
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
        m_ptr = new T*[m_rows];
        for (int i = 0; i < m_rows; ++i) {
            m_ptr[i] = new T[1];
            m_ptr[i][0] = *(lst.begin() + i);
        }
    }
    Matrix(const std::initializer_list<std::initializer_list<T>> & lst) {
        m_rows = lst.size();
        m_columns = lst.begin()->size();

        m_ptr = new T*[m_rows];
        for (int i = 0; i < m_rows; ++i) {
            m_ptr[i] = new T[m_columns];
            for (int j = 0; j < m_columns; ++j) {
                m_ptr[i][j] = *((lst.begin() + i)->begin() + j);
            }
        }

    }
    ~Matrix() {
        for (int i = 0; i < m_rows; ++i) {
            delete [] m_ptr[i];
        }
        delete [] m_ptr;
    }

    /*
    class Proxy {
        T* ptr;
    public:
        Proxy(T* other) {
            ptr = other;
        }
        T & operator[](unsigned idx) {
            return ptr[idx];
        }
    };

     Proxy operator[](unsigned idx) {
        return Proxy(m_ptr[idx]);
    }
    */

    T * operator[](unsigned idx) {
        return m_ptr[idx];
    }
    // m[0][0] = 1;
    // a = m[0][0];

    T * operator[](unsigned idx) const {
        return m_ptr[idx];
    }

    Matrix<T> & operator=(const Matrix<T> & rhs) {
        if (this == &rhs) {
            return *this;
        }

        for (int i = 0; i < m_rows; ++i) {
            delete [] m_ptr[i];
        }
        delete [] m_ptr;

        m_rows = rhs.m_rows;
        m_columns = rhs.m_columns;

        m_ptr = new T*[m_rows];
        for (int i = 0; i < m_rows; ++i) {
            m_ptr[i] = new T[m_columns];
            for (int j = 0; j < m_columns; ++j) {
                m_ptr[i][j] = rhs.m_ptr[i][j];
            }
        }

        return *this;
    }

    Matrix<T> & operator=(Matrix<T> && rhs)  noexcept {
        if (this == &rhs) {
            return *this;
        }

        for (int i = 0; i < m_rows; ++i) {
            delete [] m_ptr[i];
        }
        delete [] m_ptr;

        m_rows = rhs.m_rows;
        m_columns = rhs.m_columns;

        m_ptr = new T*[m_rows];
        for (int i = 0; i < m_rows; ++i) {
            m_ptr[i] = new T[m_columns];
            for (int j = 0; j < m_columns; ++j) {
                m_ptr[i][j] = rhs.m_ptr[i][j];
            }
        }

        for (int i = 0; i < rhs.m_rows; ++i) {
            delete [] rhs.m_ptr[i];
        }
        delete [] rhs.m_ptr;

        return *this;
    }

    T & operator()(int i, int j){
        return m_ptr[i][j];
    }

    const T & operator()(int i, int j) const {
        return m_ptr[i][j];
    }

    Matrix<T> & operator+=(const Matrix<T> & rhs) {
        for (int i = 0; i < rhs.m_rows; ++i) {
            for (int j = 0; j < rhs.m_columns; ++j) {
                this->m_ptr[i][j] += rhs.m_ptr[i][j];
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
        for (int i = 0; i < rhs.m_rows; ++i) {
            for (int j = 0; j < rhs.m_columns; ++j) {
                this->m_ptr[i][j] -= rhs.m_ptr[i][j];
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
                this->m_ptr[i][j] *= rhs;
            }
        }
        return *this;
    }

    friend Matrix<T> operator*(Matrix<T> lhs, T rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Matrix<T> operator*(Matrix<T> lhs, const Matrix<T> & rhs) {
        Matrix<T> tmp(lhs.m_rows, rhs.m_columns);
        for (int i = 0; i < lhs.m_rows; ++i) {
            for (int j = 0; j < rhs.m_columns; ++j) {
                for (int k = 0; k < lhs.m_columns; ++k) {
                    tmp[i][j] += lhs.m_ptr[i][k] * rhs.m_ptr[k][j];
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
            t += m_ptr[i][i];
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
