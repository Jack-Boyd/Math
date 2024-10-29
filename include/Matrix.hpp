#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>

class Matrix {
private:
    std::vector<std::vector<double>> values;
    size_t rows, cols;

    Matrix getMinor(const Matrix& mat, size_t row, size_t col, size_t size) const;
    double determinant(const Matrix& mat, size_t size) const;
    Matrix adjoint() const;
public:
    Matrix(size_t rows, size_t cols);
    Matrix(const std::vector<std::vector<double>>& values);

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;

    static Matrix identity(size_t size);

    void display() const;
    
    bool isSquare() const;
    Matrix inverse() const;
};

#endif // MATRIX_HPP