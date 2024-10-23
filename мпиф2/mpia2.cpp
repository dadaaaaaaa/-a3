#include <iostream>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cstdlib>

// Константа для длины вектора
const int N = 10;

// Функция для последовательного вычисления скалярного произведения
double sequential_scalar_product(const std::vector<double>& v) {
    double sum = 0.0;
    for (double value : v) {
        sum += value * value;
    }
    return sum;
}

// Функция для параллельного вычисления скалярного произведения с использованием OpenMP
double parallel_scalar_product(const std::vector<double>& v) {
    double sum = 0.0;
#pragma omp parallel reduction(+:sum) num_threads(threads_num)
    {
        double local_sum = 0.0;
#pragma omp for
        for (int i = 0; i < v.size(); ++i) {
            local_sum += v[i] * v[i];
        }
        sum += local_sum;
    }
    return sum;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Генерация вектора только на процессе 0
    std::vector<double> vec;
    if (rank == 0) {
        vec.resize(N);
        for (int i = 0; i < N; ++i) {
            vec[i] = (double)(rank+1);
        }
    }

    // Определение количества элементов, которые будут обрабатываться каждым процессом
    int elements_per_process = N / size;
    int remainder = N % size; // Остаток от деления

    // Определение количества элементов для текущего процесса
    int local_size;
    if (rank < remainder) {
        local_size = elements_per_process + 1; // Процессам с рангом меньше remainder выделим на 1 элемент больше
    }
    else {
        local_size = elements_per_process; // Остальные процессы получают равное количество
    }
    std::vector<double> local_vec(local_size);

    // Распределение вектора по процессам
    MPI_Scatter(vec.data(), elements_per_process + (rank < remainder ? 1 : 0), MPI_DOUBLE,
        local_vec.data(), local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Последовательное вычисление на каждом процессе
    double sequential_result = sequential_scalar_product(local_vec);

    // Параллельное вычисление на каждом процессе
    double parallel_result = parallel_scalar_product(local_vec);

    // Сбор результатов на процессе 0
    double total_result = 0.0;
    MPI_Reduce(&parallel_result, &total_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Процесс 0 выводит результаты
    if (rank == 0) {
        std::cout << "Скалярное произведение (параллельно): " << total_result << std::endl;
        std::cout << "Скалярное произведение (последовательно, на процессе 0): "<< sequential_scalar_product(vec) << std::endl;
    }

    MPI_Finalize();
    return 0;
}