#include <iostream>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cstdlib>

// Константа для длины вектора
const int N = 100000;

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

    std::vector<double> vec;

    // Процесс 0 инициализирует вектор и распределяет количество элементов
    if (rank == 0) {
        vec.resize(N);
        for (int i = 0; i < N; ++i) {
            vec[i] = (double)(rank + 1); // Инициализация вектора значениями
        }
    }

    // Массив, в котором будем хранить количество элементов для каждого процесса
    std::vector<int> counts(size);
    std::vector<int> displs(size);

    // Неправильное распределение количества элементов (пример)
    for (int i = 0; i < size; ++i) {
        counts[i] = (i < N % size) ? N / size + 1 : N / size; // Ненормализованное распределение
        displs[i] = i == 0 ? 0 : displs[i - 1] + counts[i - 1]; // Смещение
    }

    // Локальный вектор для каждого процесса
    std::vector<double> local_vec(counts[rank]);

    // Распределение данных
    MPI_Scatterv(vec.data(), counts.data(), displs.data(), MPI_DOUBLE,
        local_vec.data(), counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Выполнение последовательного вычисления на каждом процессе
    double sequential_result = sequential_scalar_product(local_vec);

    // Выполнение параллельного вычисления на каждом процессе
    double parallel_result = parallel_scalar_product(local_vec);

    // Сбор результатов на процессе 0
    double total_result = 0.0;
    double total_result2 = 0.0;
    MPI_Reduce(&parallel_result, &total_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sequential_result, &total_result2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // Процесс 0 выводит результаты
    if (rank == 0) {
        std::cout << "Paralel: " << total_result << std::endl;
        std::cout << "Posled: " << total_result2 << std::endl;
    }

    MPI_Finalize();
    return 0;
}