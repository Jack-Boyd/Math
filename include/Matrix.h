#pragma once

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

class Matrix {
private:
  std::vector<std::vector<double>> values;
  size_t rows, cols;

  static Matrix getMinor(const Matrix &mat, size_t row, size_t col, size_t size);
  static double determinant(const Matrix &mat, size_t size);
  Matrix adjoint() const;

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
  friend Matrix operator*(double scalar, const Matrix &mat);

  bool operator==(const Matrix& other) const;
  bool operator!=(const Matrix& other) const;

  static Matrix identity(size_t size);

  bool isSquare() const;
  size_t numRows() const { return rows; }
  size_t numCols() const { return cols; }
  
  Matrix inverse() const;
  std::vector<double> solveGaussianElimination(const std::vector<double> &mat) const;

  friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
};