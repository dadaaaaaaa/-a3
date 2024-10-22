#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define VECTOR_SIZE 100000

double calculate_norm(double* vector) {
    double sum = 0.0;

#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < VECTOR_SIZE; i++) {
        sum += vector[i] * vector[i];
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
        double norm;
        MPI_Request request;

        for (int i = 1; i < size; i++) {
            MPI_Irecv(&norm, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &request);
            int flag = 0;

            while (flag == 0) {
                // Проверяем завершение передачи сообщения
                MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
                // Здесь можно добавить какую-то дополнительную логику
            }
            printf("Процесс %d: Норма = %f\n", i, norm);
        }
    }

    MPI_Finalize();
    return 0;
}