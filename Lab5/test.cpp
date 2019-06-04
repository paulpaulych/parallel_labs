#include <cmath>
#include <cstdlib>
#include <iostream>
#include <pthread.h>
#include <mpi/mpi.h>
#include <queue>
#include <random>
#include <unistd.h>

int rank = 0, size = 0;

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
    MPI_Status status;
    int receive;
    MPI_Recv(&receive, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
