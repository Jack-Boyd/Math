#include <iostream>
#include <vector>
#include <stdexcept>
#include <random>
#include <ctime>

using namespace std;

class Matrix {
private:
    vector<vector<int>> values;
    int rows, cols;
public:
    Matrix(int rows, int cols) : rows(rows), cols(cols) {
        values.resize(rows, vector<int>(cols, 0));
    }
    void populate() {
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> distr(1, 9);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                values[i][j] = distr(gen);
            }
        }
    }
    void print() const {
        for (const auto& row : values) {
            for (const auto& val : row) {
                cout << val << " ";
            }
            cout << endl;
        }
    }
    Matrix operator+(const Matrix& other) {
        if (rows != other.rows && cols != other.cols) {
            throw invalid_argument("Matrices incompatible with addition");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.values[i][j] = values[i][j] + other.values[i][j];
            }
        }
        return result;
    }
    Matrix operator-(const Matrix& other) {
        if (rows != other.rows && cols != other.cols) {
            throw invalid_argument("Matrices incompatible with subtraction");
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.values[i][j] = values[i][j] - other.values[i][j];
            }
        }
        return result;
    }
    Matrix operator*(const Matrix& other) {
        if (cols != other.rows) {
            throw invalid_argument("Matrices incompatible with multiplication");
        }
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                result.values[i][j] = 0;
                for (int k = 0; k < cols; k++) {
                    result.values[i][j] += (values[i][k] * other.values[k][j]);
                }
            }
        }
        return result;
    }
    static Matrix identity(int size) {
        Matrix identity(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                identity.values[i][j] = i == j ? 1 : 0;
            }
        }
        return identity;
    }

    int getRows() const { return rows; }
    int getCols() const { return cols; }
};

int main() {
    srand(static_cast<unsigned>(time(nullptr)));

    Matrix m(2, 2);
    Matrix n(2, 2);

    m.populate();
    n.populate();

    Matrix a = m + n;
    Matrix s = m - n;
    Matrix r = m * n;
    Matrix i = Matrix::identity(3);

    m.print();
    cout << "----" << endl;
    n.print();
    cout << "----" << endl;
    a.print();
    cout << "----" << endl;
    s.print();
    cout << "----" << endl;
    r.print();
    cout << "----" << endl;
    i.print();
    return 0;
}