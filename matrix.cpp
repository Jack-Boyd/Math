#include <iostream>
#include <vector>
#include <stdexcept>
#include <random>
#include <ctime>

using namespace std;

class Matrix {
private:
    vector<vector<double>> values;
    int rows, cols;
    static mt19937 gen;

public:
    Matrix(int rows, int cols) : rows(rows), cols(cols), values(rows, vector<double>(cols, 0)) {}
    void populate() {
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
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
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
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
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
    Matrix operator*(const Matrix& other) const {
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
                identity.values[i][j] = (i == j) ? 1 : 0;
            }
        }
        return identity;
    }
    ~Matrix() {}

    int getRows() const { return rows; }
    int getCols() const { return cols; }
};

mt19937 Matrix::gen(time(nullptr));

int main() {
    Matrix m(2, 2);
    Matrix n(2, 2);

    m.populate();
    n.populate();

    try {
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
    } catch (const invalid_argument& e) {
        cerr << e.what() << endl;
    }

    return 0;
}
