#include "Matrix.hpp"
#include <stdexcept>
#include <iomanip>

Matrix::Matrix(size_t rows, size_t cols) : rows(rows), cols(cols), values(rows, std::vector<double>(cols, 0.0)) {}
Matrix::Matrix(const std::vector<std::vector<double>>& values) : values(values), rows(values.size()), cols(values[0].size()) {}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must agree for addition.");
    }

    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result.values[i][j] = values[i][j] + other.values[i][j];
        }
    }
    return result;
}
Matrix Matrix::operator-(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must agree for addition.");
    }

    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result.values[i][j] = values[i][j] - other.values[i][j];
        }
    }
    return result;
}
Matrix Matrix::operator*(const Matrix& other) const {
    if (cols != other.rows) {
        throw std::invalid_argument("Matrix dimensions must agree for multiplication.");
    }

    Matrix result(rows, other.cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < other.cols; ++j) {
            for (size_t k = 0; k < cols; ++k) {
                result.values[i][j] += values[i][k] * other.values[k][j];
            }
        }
    }
    return result;
}

Matrix Matrix::identity(size_t size) {
    Matrix identity(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            identity.values[i][j] = (i == j) ? 1 : 0;
        }
    }
    return identity;
}
void Matrix::display() const {
    for (const auto& row : values) {
        for (double val : row) {
            std::cout << std::setw(4) << val << " ";
        }
        std::cout << std::endl;
    }
}
bool Matrix::isSquare() const { return rows == cols; }
size_t Matrix::getRows() const { return rows; }
size_t Matrix::getCols() const { return cols; }
