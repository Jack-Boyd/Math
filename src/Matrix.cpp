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

Matrix Matrix::getMinor(const Matrix& matrix, size_t row, size_t col) {
  if (!matrix.isSquare()) 
    throw std::invalid_argument("Minor requires a square matrix.");

  size_t size = matrix.numRows();
  Matrix temp(size - 1, size - 1);

  for (size_t r = 0, i = 0; r < size; ++r) {
    if (r == row) continue;
    for (size_t c = 0, j = 0; c < size; ++c) {
      if (c == col) continue;
      temp.values[i][j++] = matrix.values[r][c];
    }
    ++i;
  }
  return temp;
}

double Matrix::determinantRecursive() const {
  if (!isSquare())
    throw std::invalid_argument("Determinant requires square matrix");
  
  size_t n = rows;
  if (n == 1) return values[0][0];
  if (n == 2)
    return values[0][0] * values[1][1] -
           values[0][1] * values[1][0];

  double det = 0.0;
  int sign = 1;
  for (size_t j = 0; j < n; ++j) {
    Matrix minor = getMinor(*this, 0, j);
    det += sign * values[0][j] * minor.determinantRecursive();
    sign = -sign;
  }
  return det;
}

double Matrix::determinantFast() const {
  if (!isSquare())
    throw std::invalid_argument("Determinant requires square matrix");
  
  size_t n = rows;
  Matrix temp = *this;
  double det = 1.0;

  for (size_t i = 0; i < n; ++i) {
    size_t pivot = i;
    for (size_t r = i + 1; r < n; ++r) {
      if (std::fabs(temp.values[r][i]) > std::fabs(temp.values[pivot][i]))
        pivot = r;
    }

    if (std::fabs(temp.values[pivot][i]) < 1e-12)
      return 0.0;

    if (pivot != i) {
      std::swap(temp.values[i], temp.values[pivot]);
      det = -det;
    }

    det *= temp.values[i][i];

    for (size_t r = i + 1; r < n; ++r) {
      double factor = temp.values[r][i] / temp.values[i][i];
      for (size_t c = i; c < n; ++c) {
        temp.values[r][c] -= factor * temp.values[i][c];
      }
    }
  }

  return det;
}

Matrix Matrix::adjoint() const {
  if (!isSquare())
    throw std::invalid_argument("Adjoint requires a square matrix.");
  
  size_t n = rows;
  Matrix adj(n, n);

  if (n == 1) {
    adj.values[0][0] = 1;
    return adj;
  }

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      Matrix cf = getMinor(*this, i, j);
      int sign = ((i + j) % 2 == 0) ? 1 : -1;
      adj.values[j][i] = sign * cf.determinantRecursive();
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
      os << std::setw(5) << val << " ";
    os << '\n';
  }
  return os;
}