//
// Created by Sergey Bogolepov on 5/11/15.
//
#include <mpi.h>
#include "poisson.h"

int main(int argc, char** argv) {

    int size;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Input data initialization
    InputData::cells_x = 20;
    InputData::cells_y = 20;
    InputData::cells_z = 20;

    InputData::X = 2.0;
    InputData::Y = 2.0;
    InputData::Z = 2.0;

    InputData::a = 1;
    InputData::e = 0.00001;

    Layer layer = init(rank, size);

    while (!iteration_step(layer));

    MPI_Finalize();
}
