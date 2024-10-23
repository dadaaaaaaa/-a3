#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define VECTOR_SIZE 10

double calculate_norm(double* vector) {
    double sum = 0.0;
#pragma omp parallel reduction(+:sum) num_threads(threads_num)
    {
#pragma omp for
        for (int i = 0; i < VECTOR_SIZE; i++) {
            sum += vector[i] * vector[i];
        }
    }
    return sqrt(sum);
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank != 0) {
        double* vector = (double*)malloc(VECTOR_SIZE * sizeof(double));
        for (int i = 0; i < VECTOR_SIZE; i++) {
            vector[i] = (double)(rank + 1);  // Пример заполнения
        }
        double norm = calculate_norm(vector);
        MPI_Send(&norm, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        free(vector);
    }
    else {
        for (int i = 1; i < size; i++) {
            double norm;
            MPI_Recv(&norm, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Proces %d: norm = %f\n", i, norm);
        }
    }

    MPI_Finalize();
    return 0;
}