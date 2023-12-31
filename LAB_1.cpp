﻿// Алгоритм умножения матрицы на вектор

#include <iostream>
#include <ctime>
#include <omp.h>
#include <chrono>

// Задаем 5 потоков
#define NUM_THREADS 4

using namespace std;

// Последовательный алгоритм умножения матрицы на вектор
double* mult_matrix_vector_serial(double** matrix, double* vector, int rows, int cols) {
    double* result = new double[rows];

    for (int i = 0; i < rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

// Распараллеливание циклов с директивой for
double* mult_matrix_vertor_parallel_for(double** matrix, double* vector, int rows, int cols) {
    double* result = new double[rows];

    // Цикл выполняется несколькими потоками параллельно (число потоков задается в NUM_THEADS)
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i = 0; i < rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

// Распараллеливания циклов без директивы for (использование директивы sections) - “ручное” задание работ
double* mult_matrix_vertor_parallel_manual(double** matrix, double* vector, int rows, int cols) {
    double* result = new double[rows];

// 5 секций потому что мы задали 5 потоков (при мзменении числа потоков
// необходимо изменить число секций, а также переписать циклы внутри секций
// для равномерного распределения данных по потокам)
#pragma omp parallel sections
    {
#pragma omp section
        {
            // Цикл идет от 0 до rows / 5 потому что мы равномерно разделяем данные на 5 потоков
            for (int i = 0; i < rows / 5; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = rows / 5; i < 2 * rows / 5; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = 2 * rows / 5; i < 3 * rows / 5; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = 3 * rows / 5; i < 4 * rows / 5; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
#pragma omp section
        {
            for (int i = 4 * rows / 5; i < rows; i++) {
                result[i] = 0.0;
                for (int j = 0; j < cols; j++) {
                    result[i] += matrix[i][j] * vector[j];
                }
            }
        }
    }

    return result;
}

int main() {
    // Число строк
    int rows = 10000;
    // Число столбцов
    int cols = 10000;

    // Создаем матрицу с рандомными значениями
    double** matrix = new double* [rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = rand() / 1000;
        }
    }

    // Создаем вектор с рандомными значениями
    double* vector = new double[cols];
    for (int i = 0; i < cols; i++) {
        vector[i] = rand() / 1000;
    }

    // Переменные для замера времени
    chrono::steady_clock::time_point start, stop;
    std::chrono::duration<double> duration;

    // Время выполнения последовательного алгоритма
    double serial_time;
    // Время выполнения алгоритма с распараллеливанием циклов с директивой for
    double parallel_for_time;
    // Время выполнения алгоритма с распараллеливанием циклов без директивы for
    double parallel_manual_time;


    // Результат для последовательного алгоритма
    start = chrono::high_resolution_clock::now();
    double* serial_result = mult_matrix_vector_serial(matrix, vector, rows, cols);
    stop = chrono::high_resolution_clock::now();
    duration = stop - start;
    serial_time = duration.count();
    // Результат при распараллеливании циклов с директивой for
    start = chrono::high_resolution_clock::now();
    double* parallel_for_result = mult_matrix_vertor_parallel_for(matrix, vector, rows, cols);
    stop = chrono::high_resolution_clock::now();
    duration = stop - start;
    parallel_for_time = duration.count();
    // Распараллеливании циклов без директивы for
    start = chrono::high_resolution_clock::now();
    double* parallel_manual_result = mult_matrix_vertor_parallel_manual(matrix, vector, rows, cols);
    stop = chrono::high_resolution_clock::now();
    duration = stop - start;
    parallel_manual_time = duration.count();

    cout << "Running time of the sequential algorithm: " << serial_time << " seconds" << endl;
    cout << "Running time with directive for: " << parallel_for_time << " seconds" << endl;
    cout << "Running time without directive for: " << parallel_manual_time << " seconds" << endl;

    // Очищаем память
    for (int i = 0; i < rows; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] vector;
    delete[] serial_result;
    delete[] parallel_for_result;
    delete[] parallel_manual_result;

    return 0;
}
