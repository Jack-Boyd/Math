#include "Matrix.h"
#include <stdexcept>
#include <iomanip>

Matrix::Matrix(size_t rows, size_t cols) : rows(rows), cols(cols), values(rows, std::vector<double>(cols, 0.0)) {}
Matrix::Matrix(const std::vector<std::vector<double>> &values) : values(values), rows(values.size()), cols(values[0].size()) {}

Matrix Matrix::operator+(const Matrix &other) const {
  if (rows != other.rows || cols != other.cols)
    throw std::invalid_argument("Matrix dimensions must agree for addition.");

  Matrix result(rows, cols);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      result.values[i][j] = values[i][j] + other.values[i][j];
    }
  }
  return result;
}

Matrix Matrix::operator-(const Matrix &other) const {
  if (rows != other.rows || cols != other.cols)
    throw std::invalid_argument("Matrix dimensions must agree for addition.");

  Matrix result(rows, cols);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      result.values[i][j] = values[i][j] - other.values[i][j];
    }
  }
  return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
  if (cols != other.rows)
    throw std::invalid_argument("Matrix dimensions must agree for multiplication.");

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
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      identity.values[i][j] = (i == j) ? 1 : 0;
    }
  }
  return identity;
}

void Matrix::display() const {
  for (const auto &row : values) {
    for (double val : row) {
      std::cout << std::setw(10) << val << " ";
    }
    std::cout << std::endl;
  }
}

bool Matrix::isSquare() const { return rows == cols; }
Matrix Matrix::getMinor(const Matrix &mat, size_t row, size_t col, size_t size) const {
  Matrix temp(size - 1, size - 1);
  int i = 0, j = 0;

  for (size_t r = 0; r < size; ++r) {
    for (size_t c = 0; c < size; ++c) {
      if (r != row && c != col) {
        temp.values[i][j++] = mat.values[r][c];
        if (j == size - 1) {
          j = 0;
          ++i;
        }
      }
    }
  }
  return temp;
}

double Matrix::determinant(const Matrix &matrix, size_t size) const {
  if (size == 1)
    return matrix.values[0][0];

  double det = 0;
  int sign = 1;

  for (size_t i = 0; i < size; ++i) {
    Matrix cf = getMinor(matrix, 0, i, size);
    det += sign * matrix.values[0][i] * determinant(cf, size - 1);
    sign = -sign;
  }
  return det;
}

Matrix Matrix::adjoint() const {
  if (!isSquare())
    throw std::invalid_argument("Matrix must be square.");

  Matrix adj(rows, cols);
  if (rows == 1) {
    adj.values[0][0] = 1;
    return adj;
  }

  int sign = 1;
  Matrix cf(rows - 1, cols - 1);

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      cf = getMinor(*this, i, j, rows);
      sign = ((i + j) % 2 == 0) ? 1 : -1;
      adj.values[j][i] = sign * determinant(cf, rows - 1);
    }
  }
  return adj;
}

Matrix Matrix::inverse() const {
  if (!isSquare())
    throw std::invalid_argument("Matrix must be square to find its inverse.");

  double det = determinant(*this, rows);
  if (det == 0)
    throw std::runtime_error("Singular matrix, can't find its inverse.");
  Matrix adj = adjoint();
  Matrix inv(rows, cols);

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      inv.values[i][j] = adj.values[i][j] / det;
    }
  }
  return inv;
}

std::vector<double> Matrix::solveGaussianElimination(const std::vector<double> &mat) const {
  if (!isSquare() && mat.size() != rows)
    throw std::invalid_argument("Matrix must be square and match the size of the vector.");

  Matrix augmented(*this);
  std::vector<double> result(mat);
  for (size_t i = 0; i < rows; ++i) {
    if (augmented.values[i][i] == 0)
      throw std::runtime_error("Singular matrix, cannot solve.");
    for (size_t j = i + 1; j < rows; ++j) {
      double ratio = augmented.values[j][i] / augmented.values[i][i];
      for (size_t k = i; k < cols; ++k) {
        augmented.values[j][k] -= augmented.values[i][k] * ratio;
      }
      result[j] -= result[i] * ratio;
    }
  }

  for (size_t i = rows - 1; i >= 0; --i) {
    for (size_t j = i + 1; j < rows; ++j) {
      result[i] -= augmented.values[i][j] * result[j];
    }
    result[i] /= augmented.values[i][i];
  }

  return result;
}