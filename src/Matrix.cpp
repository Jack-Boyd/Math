#include "Matrix.h"

Matrix::Matrix(size_t rows, size_t cols)
  : values(rows, std::vector<double>(cols, 0.0)), rows(rows), cols(cols) {}

Matrix::Matrix(const std::vector<std::vector<double>>& values)
  : values(values), rows(values.size()), cols(values[0].size()) {}

Matrix Matrix::operator+(const Matrix& other) const {
  if (rows != other.rows || cols != other.cols)
    throw std::invalid_argument("Matrix addition requires dimensions "
                                "to match: (" +
      std::to_string(rows) + "x" + std::to_string(cols) + ") vs (" +
      std::to_string(other.rows) + "x" + std::to_string(other.cols) +
      ")");

  Matrix result(rows, cols);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j)
      result.values[i][j] = values[i][j] + other.values[i][j];
  return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
  if (rows != other.rows || cols != other.cols)
    throw std::invalid_argument("Matrix subtraction requires dimensions "
                                "to match: (" +
      std::to_string(rows) + "x" + std::to_string(cols) + ") vs (" +
      std::to_string(other.rows) + "x" + std::to_string(other.cols) +
      ")");

  Matrix result(rows, cols);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j)
      result.values[i][j] = values[i][j] - other.values[i][j];
  return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
  if (rows != other.rows || cols != other.cols)
    throw std::invalid_argument("Matrix multiplication requires dimensions "
                                "to match: (" +
      std::to_string(rows) + "x" + std::to_string(cols) + ") vs (" +
      std::to_string(other.rows) + "x" + std::to_string(other.cols) +
      ")");

  Matrix result(rows, other.cols);
    for (size_t i = 0; i < rows; ++i)
      for (size_t j = 0; j < other.cols; ++j)
        for (size_t k = 0; k < cols; ++k)
          result.values[i][j] += values[i][k] * other.values[k][j];
  return result;
}

Matrix Matrix::operator*(double scalar) const {
  Matrix result(*this);
  for (auto& row : result.values)
    for (auto& val : row)
      val *= scalar;
  return result;
}

Matrix operator*(double scalar, const Matrix& mat) {
  return mat * scalar;
}

bool Matrix::operator==(const Matrix& other) const {
    return rows == other.rows && cols == other.cols && values == other.values;
}

bool Matrix::operator!=(const Matrix& other) const { return !(*this == other); }

Matrix Matrix::identity(size_t size) {
  Matrix id(size, size);
  for (size_t i = 0; i < size; ++i)
    id.values[i][i] = 1.0;
  return id;
}

bool Matrix::isSquare() const { return rows == cols; }

Matrix Matrix::getMinor(const Matrix& mat, size_t row, size_t col, size_t size) {
  Matrix temp(size - 1, size - 1);
  size_t i = 0, j = 0;
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

double Matrix::determinant(const Matrix& matrix, size_t size) {
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
    throw std::invalid_argument("Adjoint requires a square matrix.");

  Matrix adj(rows, cols);
  if (rows == 1) {
    adj.values[0][0] = 1;
    return adj;
  }

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      Matrix cf = getMinor(*this, i, j, rows);
      int sign = ((i + j) % 2 == 0) ? 1 : -1;
      adj.values[j][i] = sign * determinant(cf, rows - 1);
    }
  }
  return adj;
}

Matrix Matrix::inverse() const {
  if (!isSquare())
    throw std::invalid_argument("Inverse requires a square matrix.");

  double det = determinant(*this, rows);
  if (det == 0)
    throw std::runtime_error("Singular matrix, cannot be inverted.");

  return adjoint() * (1.0 / det);
}

std::vector<double> Matrix::solveGaussianElimination(const std::vector<double>& rhs) const {
  if (!isSquare() || rhs.size() != rows)
    throw std::invalid_argument("System must be square and dimensions "
                                "consistent.");

  Matrix augmented(*this);
  std::vector<double> result(rhs);

  for (size_t i = 0; i < rows; ++i) {
    if (augmented.values[i][i] == 0)
      throw std::runtime_error("Zero pivot encountered.");

    for (size_t j = i + 1; j < rows; ++j) {
      double ratio = augmented.values[j][i] / augmented.values[i][i];
      for (size_t k = i; k < cols; ++k)
        augmented.values[j][k] -= augmented.values[i][k] * ratio;
      result[j] -= result[i] * ratio;
    }
  }

  for (int i = static_cast<int>(rows) - 1; i >= 0; --i) {
    for (size_t j = i + 1; j < rows; ++j)
      result[i] -= augmented.values[i][j] * result[j];
    result[i] /= augmented.values[i][i];
  }

  return result;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
  for (const auto& row : m.values) {
    for (double val : row)
      os << std::setw(10) << val << " ";
    os << '\n';
  }
  return os;
}