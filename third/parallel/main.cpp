//
// Created by Sergey Bogolepov on 5/11/15.
//
#include <mpi.h>
#include <iostream>
#include "poisson.h"

int main(int argc, char** argv) {

    int size;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Layer layer = init(rank, size);
    while (!iteration_step(layer));
    std::cout << "Rank" << layer.rank << " completed" << std::endl;
    MPI_Finalize();
}
