#pragma once

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

class Matrix {
private:
  std::vector<std::vector<double>> values;
  size_t rows, cols;

  static Matrix getMinor(const Matrix &matrix, size_t row, size_t col);

public:
  Matrix(size_t rows, size_t cols);
  Matrix(const std::vector<std::vector<double>> &values);

  Matrix(const Matrix &) = default;
  Matrix(Matrix &&) noexcept = default;
  Matrix &operator=(const Matrix &) = default;
  Matrix &operator=(Matrix &&) noexcept = default;
  ~Matrix() = default;

  Matrix operator+(const Matrix &other) const;
  Matrix operator-(const Matrix &other) const;
  Matrix operator*(const Matrix &other) const;

  Matrix operator*(double scalar) const;
  friend Matrix operator*(double scalar, const Matrix &matrix);

  bool operator==(const Matrix& other) const;
  bool operator!=(const Matrix& other) const;

  friend std::ostream& operator<<(std::ostream& os, const Matrix& m);

  static Matrix identity(size_t size);
  bool isSquare() const;
  size_t numRows() const { return rows; }
  size_t numCols() const { return cols; }
  
  double determinantRecursive() const;
  double determinantFast() const;

  Matrix adjoint() const;
  Matrix inverse() const;
  
  std::vector<double> solveGaussianElimination(const std::vector<double> &rhs) const;
};