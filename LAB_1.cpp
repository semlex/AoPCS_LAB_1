#include <iostream>
#include <omp.h>
// Задаем 6 потоков
#define NUM_THREADS 6

using namespace std;

// Последовательный алгоритм умножения матрицы на вектор
void mult_matrix_vector_serial(double** matrix, double* vector, double* result, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

// Распараллеливания циклов с директивой for
void mult_matrix_vertor_parralel_for(double** matrix, double* vector, double* result, int rows, int cols) {
#pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

// Распараллеливания циклов без директивы for (использование директивы sections) - “ручное” задание работ
void mult_matrix_vertor_parralel_manual(double** matrix, double* vector, double* result, int rows, int cols) {
// 6 секций потому что мы задали 6 потоков
#pragma omp parallel sections
    {
#pragma omp section
        {
            for (int i = 0; i < rows / 6; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = rows / 6; i < 2 * rows / 6; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = 2 * rows / 6; i < 3 * rows / 6; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = 3 * rows / 6; i < 4 * rows / 6; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = 4 * rows / 6; i < 5 * rows / 6; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = 5 * rows / 6; i < rows; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
    }
}

int main() {
    int rows = 100;
    int cols = 100;

    // Создаем матрицу с рандомными значениями
    double** matrix = new double* [rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = rand();
        }
    }

    // Создаем вектор с рандомными значениями
    double* vector = new double[cols];
    for (int i = 0; i < cols; i++) {
        vector[i] = rand();
    }

    double* serial_result = new double[rows];            // Результат для последовательного алгоритма
    double* parralel_for_result = new double[rows];      // Результат при распараллеливании циклов с директивой for
    double* parralel_manual_result = new double[rows];   // Распараллеливании циклов без директивы for

    mult_matrix_vector_serial(matrix, vector, serial_result, rows, cols);
    mult_matrix_vertor_parralel_for(matrix, vector, parralel_for_result, rows, cols);
    mult_matrix_vertor_parralel_manual(matrix, vector, parralel_manual_result, rows, cols);

    for (int i = 0; i < rows; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] vector;
    delete[] serial_result;
    delete[] parralel_for_result;
    delete[] parralel_manual_result;

    return 0;
}
