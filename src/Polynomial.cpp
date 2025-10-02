#include "Polynomial.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

Polynomial::Polynomial() : m_coefficients(1, {0.0, 0.0}) {}
Polynomial::Polynomial(const std::vector<std::complex<double>>& coefficients) 
  : m_coefficients(coefficients) {}

void Polynomial::normalize() {
  while(m_coefficients.size() > 1 && 
    m_coefficients.back() == std::complex<double>(0.0, 0.0)) {
    m_coefficients.pop_back();
  }
}

Polynomial Polynomial::operator+(const Polynomial& rhs) const {
  size_t maxSize = std::max(getSize(), rhs.getSize());
  std::vector<std::complex<double>> resultCoeffs(maxSize, {0.0, 0.0});

  for (size_t i = 0; i < getSize(); ++i) {
    resultCoeffs[i] += m_coefficients[i];
  }

  for (size_t i = 0; i < rhs.getSize(); ++i) {
    resultCoeffs[i] += rhs.m_coefficients[i];
  }

  Polynomial result(resultCoeffs);
  result.normalize();
  return result;
}

Polynomial Polynomial::operator-(const Polynomial& rhs) const {
  size_t maxSize = std::max(getSize(), rhs.getSize());
  std::vector<std::complex<double>> resultCoeffs(maxSize, {0.0, 0.0});

  for (size_t i = 0; i < getSize(); ++i) {
    resultCoeffs[i] += m_coefficients[i];
  }

  for (size_t i = 0; i < rhs.getSize(); ++i) {
    resultCoeffs[i] -= rhs.m_coefficients[i];
  }

  Polynomial result(resultCoeffs);
  result.normalize();
  return result;
}

Polynomial Polynomial::operator*(const Polynomial& rhs) const {
  if ((m_coefficients.size() == 1 && m_coefficients[0] == std::complex<double>(0.0, 0.0)) ||
    (rhs.m_coefficients.size() == 1 && rhs.m_coefficients[0] == std::complex<double>(0.0, 0.0))) {
    return Polynomial();
  }

  size_t resultSize = m_coefficients.size() + rhs.m_coefficients.size() - 1;
  std::vector<std::complex<double>> resultCoeffs(resultSize, {0.0, 0.0});

  for (size_t i = 0; i < getSize(); ++i) {
    for (size_t j = 0; j < rhs.getSize(); ++j) {
      resultCoeffs[i + j] += m_coefficients[i] * rhs.m_coefficients[j];
    }
  }

  Polynomial result(resultCoeffs);
  result.normalize();
  return result;
}

std::string toSuperscript(int num) {
  static const std::unordered_map<char, std::string> superscript = {
    {'0', "⁰"}, {'1', "¹"}, {'2', "²"}, {'3', "³"}, 
    {'4', "⁴"}, {'5', "⁵"}, {'6', "⁶"}, {'7', "⁷"}, 
    {'8', "⁸"}, {'9', "⁹"}, {'-', "⁻"}
  };

  std::string str = std::to_string(num);
  std::string result;
  result.reserve(str.size());

  for (char c : str) {
    auto it = superscript.find(c);
    if (it != superscript.end()) {
      result += it->second;
    }
    else {
      throw std::invalid_argument("Unsupported character for superscript");
    }
  }
  return result;
}

std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
  std::ostringstream ss;
  bool firstTerm = true;

  for (size_t i = p.getSize(); i -- > 0;) {
    auto coeff = p.m_coefficients[i];
    if (coeff == std::complex<double>(0.0, 0.0)) continue;

    double realCoeff = coeff.real();

    if (!firstTerm) {
      ss << (realCoeff < 0 ? " - " : " + ");
    }
    else {
      if (realCoeff < 0) ss << "-";
      firstTerm = false;
    }

    double absCoeff = std::abs(realCoeff);

    if (i == 0) {
      ss << absCoeff;
    }
    else {
      if (absCoeff != 1.0) ss << absCoeff;

      ss << "x";
      if (i > 1) ss << toSuperscript(static_cast<int>(i));
    }
  }

  if (firstTerm) ss << "0";
  os << ss.str();
  return os;
}


std::complex<double> Polynomial::evaluate(const std::complex<double>& x) const {
  std::complex<double> result(0.0, 0.0);
  for (size_t i = getSize(); i-- > 0;) {
    result = result * x + m_coefficients[i];
  }
  return result;
}

Polynomial Polynomial::derivative() const {
  std::vector<std::complex<double>> resultCoeffs(getSize() - 1, {0.0, 0.0});
  for (size_t i = 1; i < getSize(); ++i) {
    resultCoeffs[i - 1] = m_coefficients[i].real() * static_cast<double>(i);
  }

  return Polynomial(resultCoeffs);
}

Polynomial Polynomial::integral(std::complex<double> constant) const {
  std::vector<std::complex<double>> resultCoeffs(getSize() + 1, {0.0, 0.0});
  resultCoeffs[0] = constant;
  for (size_t i = 0; i < getSize(); ++i) {
    resultCoeffs[i + 1] = m_coefficients[i].real() / static_cast<double>(i + 1);
  }

  return Polynomial(resultCoeffs);
}

std::complex<double> Polynomial::solveLinear() {
  if (getSize() != 2) {
    throw std::invalid_argument("Polynomial is not linear");
  }

  std::complex<double> a0 = m_coefficients[0];
  std::complex<double> a1 = m_coefficients[1];

  if (a1 == std::complex<double>(0.0, 0.0)) {
    throw std::invalid_argument("Coefficient of x cannot be zero in a linear equation");
  }

  return -a0 / a1;
}

std::vector<std::complex<double>> Polynomial::solveQuadratic() {
  if (getSize() != 3) {
    throw std::invalid_argument("Polynomial is not a quadratic");
  }

  std::complex<double> a = m_coefficients[2];
  std::complex<double> b = m_coefficients[1];
  std::complex<double> c = m_coefficients[0];

  if (a == std::complex<double>(0.0, 0.0)) {
    throw std::invalid_argument("Coefficient 'a' cannot be 0");
  }

  std::complex<double> discriminant = b * b - std::complex<double>(4.0, 0.0) * a * c;
  std::complex<double> sqrtDiscriminant = std::sqrt(discriminant);

  std::complex<double> root1 = (-b + sqrtDiscriminant) / (std::complex<double>(2.0, 0.0) * a);
  std::complex<double> root2 = (-b - sqrtDiscriminant) / (std::complex<double>(2.0, 0.0) * a);

  return {root1, root2};
}

std::vector<std::complex<double>> Polynomial::solveCubic() {
  if (getSize() != 4) {
    throw std::invalid_argument("Polynomial is not a cubic");
  }
  std::complex<double> a = m_coefficients[3];
  std::complex<double> b = m_coefficients[2];
  std::complex<double> c = m_coefficients[1];
  std::complex<double> d = m_coefficients[0];

  if (a == std::complex<double>(0.0, 0.0)) {
    throw std::invalid_argument("Coefficient 'a' cannot be zero in a cubic");
  }

  std::complex<double> A = b / a;
  std::complex<double> B = c / a;
  std::complex<double> C = d / a;

  std::complex<double> p = B - A * A / std::complex<double>(3.0, 0.0);
  std::complex<double> q = 
    (2.0 * A * A * A / std::complex<double>(27.0, 0.0)) -
    (A * B / std::complex<double>(3.0, 0.0)) + C;

  std::complex<double> discriminant = 
    (q * q) / std::complex<double>(4.0, 0.0) + 
    (p * p * p) / std::complex<double>(27.0, 0.0);

  std::complex<double> sqrtDiscriminant = std::sqrt(discriminant);

  std::complex<double> u = std::pow(-q / std::complex<double>(2.0, 0.0) + sqrtDiscriminant, 1.0 / 3.0);
  std::complex<double> v = std::pow(-q / std::complex<double>(2.0, 0.0) - sqrtDiscriminant, 1.0 / 3.0);

  if (u == std::complex<double>(0.0, 0.0)) 
    u = std::pow(-q / std::complex<double>(2.0, 0.0) - sqrtDiscriminant, 1.0 / 3.0);

  std::complex<double> omega(-0.5, std::sqrt(3) / 2.0);

  std::vector<std::complex<double>> roots;
  roots.push_back(u + v - A / std::complex<double>(3.0, 0.0));
  roots.push_back(u * omega + v * std::conj(omega) - A / std::complex<double>(3.0, 0.0));
  roots.push_back(u * std::conj(omega) + v * omega - A / std::complex<double>(3.0, 0.0));

  return roots;
}

std::vector<std::complex<double>> Polynomial::solveQuartic() {
  if (getSize() != 5) {
    throw std::invalid_argument("Polynomial is not a quartic");
  }

  std::complex<double> a = m_coefficients[4];
  std::complex<double> b = m_coefficients[3];
  std::complex<double> c = m_coefficients[2];
  std::complex<double> d = m_coefficients[1];
  std::complex<double> e = m_coefficients[0];

  if (a == std::complex<double>(0.0, 0.0)) {
    throw std::invalid_argument("Coefficient 'a' cannot be zero in a quartic");
  }

  b /= a;
  c /= a;
  d /= a;
  e /= a;

  std::complex<double> alpha = -3.0*(b*b)/8.0 + c;
  std::complex<double> beta = b*b*b/8.0 - b*c/2.0 + d;
  std::complex<double> gamma = -3.0*(b*b*b*b)/256.0 + (c*b*b)/16.0 - (b*d)/4.0 + e;

  Polynomial resolvent({
    (beta*beta)/(-8.0),
    (gamma - alpha*alpha/4.0),
    (-alpha/2.0),
    1.0
  });

  auto cubicRoots = resolvent.solveCubic();
  std::complex<double> P = cubicRoots[0];

  std::complex<double> R = std::sqrt(0.25*b*b - c + P);
  std::complex<double> D, E;

  if (R == std::complex<double>(0.0, 0.0)) {
    D = std::sqrt(0.75*b*b - 2.0*c + 2.0*std::sqrt(P*P - 4.0*gamma));
    E = std::sqrt(0.75*b*b - 2.0*c - 2.0*std::sqrt(P*P - 4.0*gamma));
  }
  else {
    D = std::sqrt(0.75*b*b - R*R - 2.0*c + (0.25*beta)/R);
    E = std::sqrt(0.75*b*b - R*R - 2.0*c - (0.25*beta)/R);
  }

  std::complex<double> shift = -b/4.0;
  std::vector<std::complex<double>> roots;
  roots.push_back(shift + (R + D)/2.0);
  roots.push_back(shift + (R - D)/2.0);
  roots.push_back(shift + (-R + E)/2.0);
  roots.push_back(shift + (-R - E)/2.0);

  return roots;
}

std::vector<std::complex<double>> Polynomial::solve() {
  normalize();
  size_t deg = getDegree();

  switch(deg) {
    case 0:
      throw std::invalid_argument("Constant polynomial has no roots.");
    case 1:
      return { solveLinear() };
    case 2:
      return solveQuadratic();
    case 3:
      return solveCubic();
    case 4:
      return solveQuartic();
    default:
      throw std::runtime_error(
        "Polynomial degree too large"
      );
  }
}