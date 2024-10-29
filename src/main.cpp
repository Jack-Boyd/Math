#include "Matrix.hpp"

int main() {
    Matrix m({{9, 3, 4}, {4, 3, 4}, {1, 1, 1}});
    std::cout << "Matrix:" << std::endl;
    m.display();

    Matrix inv = m.inverse();
    std::cout << "Inverse:" << std::endl;
    inv.display();
    return 0;
}
