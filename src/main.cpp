#include "Matrix.h"
#include "Polynomial.h"
#include <iostream>
#include <vector>

int main() {
  try {
    //MATRIX
    Matrix A({{1, 2, 3},
              {0, 1, 4},
              {5, 6, 0}});

    Matrix B({{2, 0, 1},
              {1, 3, 2},
              {0, 1, 1}});

    std::cout << "Matrix A:\n" << A << "\n";
    std::cout << "Matrix B:\n" << B << "\n";

    std::cout << "A + B:\n" << (A + B) << "\n";
    std::cout << "A - B:\n" << (A - B) << "\n";
    std::cout << "A * B:\n" << (A * B) << "\n";
    std::cout << "2 * A:\n" << (2 * A) << "\n";

    std::cout << "det(A) (recursive): " << A.determinantRecursive() << "\n";
    std::cout << "det(A) (fast): " << A.determinantFast() << "\n";
    std::cout << "Inverse of A:\n" << A.inverse() << "\n";
    std::cout << "Identity matrix (size 3):\n" << Matrix::identity(3) << "\n";

    std::vector<double> b = {1, 2, 3};
    std::vector<double> solution = A.solveGaussianElimination(b);

    std::cout << "Solution to A * x = b (where b = [1,2,3]):\n";
    for (double v : solution) {
      std::cout << v << " ";
    }
    std::cout << "\n";

    //POLYNOMIAL
    std::cout << "\n\n\n------------------------------------\n";
    std::cout << "POLYNOMIAL \n\n";

    Polynomial linear({2.0, -4.0});
    Polynomial quadratic({1.0, -3.0, 2.0});
    Polynomial cubic({-6.0 , 11.0, -6.0, 1.0});
    Polynomial quartic({1.0, 0.0, -5.0, 0.0, 4.0});
    std::cout << "A (Linear): " << linear << "\n";
    std::cout << "B (Quadratic): " << quadratic << "\n";
    std::cout << "C (Cubic): " << cubic << "\n";
    std::cout << "D (Quartic): " << quartic << "\n";

    std::cout << "A + B: " << (linear + quadratic) << "\n";
    std::cout << "B + C: " << (quadratic + cubic) << "\n";
    std::cout << "A - D: " << (linear - cubic) << "\n";
    std::cout << "B * D: " << (quadratic * cubic) << "\n\n";

    std::cout << "A(4): " << linear.evaluate(4).real() << "\n";
    std::cout << "B(4): " << quadratic.evaluate(4).real() << "\n";
    std::cout << "C(4): " << cubic.evaluate(4).real() << "\n";
    std::cout << "D(4): " << quartic.evaluate(4).real() << "\n\n";

    std::cout << "Derivatives \n";
    std::cout << "dB/dx: " << quadratic.derivative() << "\n";
    std::cout << "dC/dx: " << cubic.derivative() << "\n";
    std::cout << "dD/dx: " << quartic.derivative() << "\n\n";

    std::cout << "Integral A: " << linear.integral(4) << "\n";
    std::cout << "Integral B: " << quadratic.integral(4) << "\n";
    std::cout << "Integral C: " << cubic.integral(4) << "\n";
    std::cout << "Integral D: " << quartic.integral(4) << "\n";

    auto linearRoot = linear.solve();
    std::cout << "A Root" << "\n";
    for (auto root : linearRoot) {
      std::cout << root << " ";
    }
    std::cout << "\n";

    auto quadraticRoots = quadratic.solve();
    std::cout << "B Roots" << "\n";
    for (auto root : quadraticRoots) {
      std::cout << root << " ";
    }
    std::cout << "\n";

    auto cubicRoots = cubic.solve();
    std::cout << "C Roots" << "\n";
    for (auto root : cubicRoots) {
      std::cout << root << " ";
    }
    std::cout << "\n";

    auto quarticRoots = quartic.solve();
    std::cout << "D Roots" << "\n";
    for (auto root : quarticRoots) {
      std::cout << root << " ";
    }
    std::cout << "\n";
  }
  catch (const std::exception &ex) {
    std::cerr << "Error: " << ex.what() << "\n";
  }

  return 0;
}