#include "Matrix.h"
#include <iostream>
#include <vector>

int main() {
  try {
    Matrix A({{1, 2, 3},
              {0, 1, 4},
              {5, 6, 0}});

    Matrix B({{2, 0, 1},
              {1, 3, 2},
              {0, 1, 1}});

    std::cout << "Matrix A:\n"
              << A << "\n";
    std::cout << "Matrix B:\n"
              << B << "\n";

    std::cout << "A + B:\n"
              << (A + B) << "\n";

    std::cout << "A - B:\n"
              << (A - B) << "\n";

    std::cout << "A * B:\n"
              << (A * B) << "\n";

    std::cout << "2 * A:\n"
              << (2 * A) << "\n";

    std::cout << "det(A) (recursive): " << A.determinantRecursive() << "\n";
    std::cout << "det(A) (fast): " << A.determinantFast() << "\n";

    std::cout << "Inverse of A:\n"
              << A.inverse() << "\n";

    std::cout << "Identity matrix (size 3):\n"
              << Matrix::identity(3) << "\n";

    std::vector<double> b = {1, 2, 3};
    std::vector<double> solution = A.solveGaussianElimination(b);

    std::cout << "Solution to A * x = b (where b = [1,2,3]):\n";
    for (double v : solution) {
      std::cout << v << " ";
    }
    std::cout << "\n";
  }
  catch (const std::exception &ex) {
    std::cerr << "Error: " << ex.what() << "\n";
  }

  return 0;
}