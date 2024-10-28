#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>

class Matrix {
private:
    std::vector<std::vector<double>> values;
    size_t rows, cols;
public:
    Matrix(size_t rows, size_t cols);
    Matrix(const std::vector<std::vector<double>>& values);

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    
    static Matrix identity(size_t size);
    void display() const;
    bool isSquare() const;
    size_t getRows() const;
    size_t getCols() const;
};

#endif // MATRIX_HPP