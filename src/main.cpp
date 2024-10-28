#include "Matrix.hpp"

int main() {
    Matrix m1({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    Matrix m2({{9, 8, 7}, {6, 5, 4}, {3, 2, 1}});

    std::cout << "Matrix 1:" << std::endl;
    m1.display();
    std::cout << "Matrix 2:" << std::endl;
    m2.display();
    Matrix m3 = m1 + m2;
    std::cout << "Matrix 1 + 2:" << std::endl;
    m3.display();
    Matrix m4 = m1 - m2;
    std::cout << "Matrix 1 - 2:" << std::endl;
    m4.display();

    Matrix m5 = m1 * m2;
    std::cout << "Matrix 1 * 2:" << std::endl;
    m5.display();
    return 0;
}
