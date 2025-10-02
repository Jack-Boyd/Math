#include "Matrix.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <cmath>

Matrix::Matrix(size_t rows, size_t cols)
  : m_values(rows, std::vector<double>(cols, 0.0)), m_rows(rows), m_cols(cols) {}
Matrix::Matrix(const std::vector<std::vector<double>>& values)
  : m_values(values), m_rows(values.size()), m_cols(values[0].size()) {}

Matrix Matrix::operator+(const Matrix& other) const {
  if (m_rows != other.m_rows || m_cols != other.m_cols)
    throw std::invalid_argument("Matrix addition requires dimensions "
                                "to match: (" +
      std::to_string(m_rows) + "x" + std::to_string(m_cols) + ") vs (" +
      std::to_string(other.m_rows) + "x" + std::to_string(other.m_cols) +
      ")");

  Matrix result(m_rows, m_cols);
  for (size_t i = 0; i < m_rows; ++i)
    for (size_t j = 0; j < m_cols; ++j)
      result.m_values[i][j] = m_values[i][j] + other.m_values[i][j];
  return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
  if (m_rows != other.m_rows || m_cols != other.m_cols)
    throw std::invalid_argument("Matrix subtraction requires dimensions "
                                "to match: (" +
      std::to_string(m_rows) + "x" + std::to_string(m_cols) + ") vs (" +
      std::to_string(other.m_rows) + "x" + std::to_string(other.m_cols) +
      ")");

  Matrix result(m_rows, m_cols);
  for (size_t i = 0; i < m_rows; ++i)
    for (size_t j = 0; j < m_cols; ++j)
      result.m_values[i][j] = m_values[i][j] - other.m_values[i][j];
  return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
  if (m_rows != other.m_rows || m_cols != other.m_cols)
    throw std::invalid_argument("Matrix multiplication requires dimensions "
                                "to match: (" +
      std::to_string(m_rows) + "x" + std::to_string(m_cols) + ") vs (" +
      std::to_string(other.m_rows) + "x" + std::to_string(other.m_cols) +
      ")");

  Matrix result(m_rows, other.m_cols);
    for (size_t i = 0; i < m_rows; ++i)
      for (size_t j = 0; j < other.m_cols; ++j)
        for (size_t k = 0; k < m_cols; ++k)
          result.m_values[i][j] += m_values[i][k] * other.m_values[k][j];
  return result;
}

Matrix Matrix::operator*(double scalar) const {
  Matrix result(*this);
  for (auto& row : result.m_values)
    for (auto& val : row)
      val *= scalar;
  return result;
}

Matrix operator*(double scalar, const Matrix& mat) {
  return mat * scalar;
}

bool Matrix::operator==(const Matrix& other) const {
  return m_rows == other.m_rows && m_cols == other.m_cols && m_values == other.m_values;
}

bool Matrix::operator!=(const Matrix& other) const { return !(*this == other); }

Matrix Matrix::identity(size_t size) {
  Matrix id(size, size);
  for (size_t i = 0; i < size; ++i)
    id.m_values[i][i] = 1.0;
  return id;
}

bool Matrix::isSquare() const { return m_rows == m_cols; }

Matrix Matrix::getMinor(const Matrix& matrix, size_t row, size_t col) {
  if (!matrix.isSquare()) 
    throw std::invalid_argument("Minor requires a square matrix.");

  size_t size = matrix.numRows();
  Matrix temp(size - 1, size - 1);

  for (size_t r = 0, i = 0; r < size; ++r) {
    if (r == row) continue;
    for (size_t c = 0, j = 0; c < size; ++c) {
      if (c == col) continue;
      temp.m_values[i][j++] = matrix.m_values[r][c];
    }
    ++i;
  }
  return temp;
}

double Matrix::determinantRecursive() const {
  if (!isSquare())
    throw std::invalid_argument("Determinant requires square matrix");
  
  size_t n = m_rows;
  if (n == 1) return m_values[0][0];
  if (n == 2)
    return m_values[0][0] * m_values[1][1] -
           m_values[0][1] * m_values[1][0];

  double det = 0.0;
  int sign = 1;
  for (size_t j = 0; j < n; ++j) {
    Matrix minor = getMinor(*this, 0, j);
    det += sign * m_values[0][j] * minor.determinantRecursive();
    sign = -sign;
  }
  return det;
}

double Matrix::determinantFast() const {
  if (!isSquare())
    throw std::invalid_argument("Determinant requires square matrix");
  
  size_t n = m_rows;
  Matrix temp = *this;
  double det = 1.0;

  for (size_t i = 0; i < n; ++i) {
    size_t pivot = i;
    for (size_t r = i + 1; r < n; ++r) {
      if (std::fabs(temp.m_values[r][i]) > std::fabs(temp.m_values[pivot][i]))
        pivot = r;
    }

    if (std::fabs(temp.m_values[pivot][i]) < 1e-12)
      return 0.0;

    if (pivot != i) {
      std::swap(temp.m_values[i], temp.m_values[pivot]);
      det = -det;
    }

    det *= temp.m_values[i][i];

    for (size_t r = i + 1; r < n; ++r) {
      double factor = temp.m_values[r][i] / temp.m_values[i][i];
      for (size_t c = i; c < n; ++c) {
        temp.m_values[r][c] -= factor * temp.m_values[i][c];
      }
    }
  }

  return det;
}

Matrix Matrix::adjoint() const {
  if (!isSquare())
    throw std::invalid_argument("Adjoint requires a square matrix.");
  
  size_t n = m_rows;
  Matrix adj(n, n);

  if (n == 1) {
    adj.m_values[0][0] = 1;
    return adj;
  }

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      Matrix cf = getMinor(*this, i, j);
      int sign = ((i + j) % 2 == 0) ? 1 : -1;
      adj.m_values[j][i] = sign * cf.determinantRecursive();
    }
  }

  return adj;
}

Matrix Matrix::inverse() const {
  if (!isSquare())
    throw std::invalid_argument("Inverse requires a square matrix.");
  
  double det = determinantFast();
  if (std::fabs(det) < 1e-12) 
    throw std::runtime_error("Matrix is singular");
  
  Matrix adj = adjoint();
  return adj * (1.0 / det);
}

std::vector<double> Matrix::solveGaussianElimination(const std::vector<double>& rhs) const {
  if (!isSquare() || rhs.size() != m_rows)
    throw std::invalid_argument("System must be square and dimensions "
                                "consistent.");

  Matrix augmented(*this);
  std::vector<double> result(rhs);

  for (size_t i = 0; i < m_rows; ++i) {
    if (augmented.m_values[i][i] == 0)
      throw std::runtime_error("Zero pivot encountered.");

    for (size_t j = i + 1; j < m_rows; ++j) {
      double ratio = augmented.m_values[j][i] / augmented.m_values[i][i];
      for (size_t k = i; k < m_cols; ++k)
        augmented.m_values[j][k] -= augmented.m_values[i][k] * ratio;
      result[j] -= result[i] * ratio;
    }
  }

  for (std::ptrdiff_t i = m_rows - 1; i >= 0; --i) {
    for (size_t j = i + 1; j < m_rows; ++j)
      result[i] -= augmented.m_values[i][j] * result[j];
    result[i] /= augmented.m_values[i][i];
  }

  return result;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
  for (const auto& row : m.m_values) {
    for (double val : row)
      os << std::setw(10) << val << " ";
    os << '\n';
  }
  return os;
}