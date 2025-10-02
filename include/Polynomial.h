#pragma once

#include <complex>
#include <vector>

class Polynomial {
private:
  std::vector<std::complex<double>> m_coefficients;
public:
  Polynomial();
  Polynomial(const std::vector<std::complex<double>>& coefficients);

  Polynomial(const Polynomial &) = default;
  Polynomial(Polynomial &&) noexcept = default;
  Polynomial &operator=(const Polynomial &) = default;
  Polynomial &operator=(Polynomial &&) noexcept = default;
  ~Polynomial() = default;

  Polynomial operator+(const Polynomial& rhs) const;
  Polynomial operator-(const Polynomial& rhs) const;
  Polynomial operator*(const Polynomial& rhs) const;
  std::complex<double> operator[](const size_t index) const { return m_coefficients[index]; }

  int getSize() const { return m_coefficients.size(); }
  int getDegree() const { return m_coefficients.size() - 1; }
  std::complex<double> getCoefficient(size_t index) const { return m_coefficients[index]; }
  void normalize();

  std::complex<double> evaluate(const std::complex<double>& x) const;
  Polynomial derivative() const;
  Polynomial integral(std::complex<double> constant = 0.0) const;

  std::complex<double> solveLinear();
  std::vector<std::complex<double>> solveQuadratic();
  std::vector<std::complex<double>> solveCubic();
  std::vector<std::complex<double>> solveQuartic();
  std::vector<std::complex<double>> solve();
  
  friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);
};