#include "Matrix.hpp"

int main() {
  Matrix m({{9, 3, 4}, {4, 3, 4}, {1, 1, 1}});
  std::cout << "Matrix:" << std::endl;
  m.display();

  std::vector<double> vec({7, 8, 3});
  std::vector<double> sol = m.solveGaussianElimination(vec);
  for (double s : sol) {
    std::cout << s << std::endl;
  }
  return 0;
}
