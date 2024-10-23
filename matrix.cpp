#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
    int** values;
    int rows;
    int cols;
} Matrix;

Matrix* createMatrix(int rows, int cols) {
    Matrix* matrix = (Matrix*)malloc(sizeof(Matrix));
    matrix->rows = rows;
    matrix->cols = cols;

    matrix->values = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++) {
        matrix->values[i] = (int*)malloc(cols * sizeof(int));
    }

    return matrix;
}
void populateMatrix(Matrix* matrix) {
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            matrix->values[i][j] = (rand() % 9) + 1;
        }
    }
}
void printMatrix(Matrix* matrix) {
    for (int r = 0; r < matrix->rows; r++) {
        for (int c = 0; c < matrix->cols; c++) {
            printf("%d ", matrix->values[r][c]);
        }
        printf("\n");
    }
}
void freeMatrix(Matrix* matrix) {
    for (int i = 0; i < matrix->rows; i++) {
        free(matrix->values[i]);
    }
    free(matrix->values);
    matrix->values = NULL;
}

Matrix* addMatrices(Matrix* m, Matrix* n) {
    if (m->rows != n->rows && m->cols != n->cols) {
        printf("Matrices incompatible with adding\n");
    }

    Matrix* res = createMatrix(m->rows, m->cols);

    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < n->cols; j++) {
            res->values[i][j] = m->values[i][j] + n->values[i][j];
        }
    }
    return res;
}
Matrix* subtractMatrices(Matrix* m, Matrix* n) {
    if (m->rows != n->rows && m->cols != n->cols) {
        printf("Matrices incompatible with subtracting\n");
    }
    
    Matrix* res = createMatrix(m->rows, m->cols);

    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < n->cols; j++) {
            res->values[i][j] = m->values[i][j] - n->values[i][j];
        }
    }
    return res;
}
Matrix* multiplyMatrices(Matrix* m, Matrix* n) {
    if (m->cols != n->rows) {
        printf("Matrices incompatible with multiplication\n");
        return NULL;
    }

    Matrix* res = createMatrix(m->rows, n->cols);

    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < n->cols; j++) {
            res->values[i][j] = 0;
            for (int k = 0; k < m->cols; k++) {
                res->values[i][j] += (m->values[i][k] * n->values[k][j]);
            }
        }
    }
    return res;
}
Matrix* identityMatrix(int size) {
    Matrix* identity = createMatrix(size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            identity->values[i][j] = i == j ? 1 : 0;
        }
    }
    return identity;
}

int main() {
    srand(time(NULL));

    Matrix* m = createMatrix(2, 2);
    Matrix* n = createMatrix(2, 2);

    populateMatrix(m);
    populateMatrix(n);

    Matrix* a = addMatrices(m, n);
    Matrix* s = subtractMatrices(m, n);
    Matrix* r = multiplyMatrices(m, n);
    Matrix* i = identityMatrix(3);

    printMatrix(m);
    printf("-----\n");
    printMatrix(n);
    printf("+\n");
    printMatrix(a);
    printf("-\n");
    printMatrix(s);
    printf("*\n");
    printMatrix(r);
    printf("-----\n");
    printMatrix(i);

    freeMatrix(m);
    freeMatrix(n);
    freeMatrix(a);
    freeMatrix(s);
    freeMatrix(r);
    freeMatrix(i);

    return 0;
}