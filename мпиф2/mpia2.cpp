#include <mpi.h>
#include <iostream>


int main(int argc, char* argv[]) {
    // Инициализация MPI
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Общее количество процессов

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Номер процесса

    // Вывод номера процесса и общего количества процессов
    std::cout << "Potok " << world_rank << " iz " << world_size << std::endl;

    // Завершение MPI
    MPI_Finalize();
    return 0;
}