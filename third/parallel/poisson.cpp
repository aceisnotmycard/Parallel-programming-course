//
// Created by Sergey Bogolepov on 5/10/15.
//

#include <mpi.h>
#include <math.h>
#include <ostream>
#include <iostream>
#include "poisson.h"

int get_chunk_size(int given_rank, int rank);

double InputData::border(double x, double y, double z) {
    return x + y + z;
}

double InputData::ro(double x, double y, double z) {
    return -a*(x+y+z);
}

Layer init(int rank, int size) {
    int size_z = get_chunk_size(rank, size);
    Layer layer(InputData::cells_x, InputData::cells_y, size_z, rank, size);

    for (int i = 0; i < layer.size_x; i++) {
        for (int j = 0; j < layer.size_y; j++) {
            for (int k = 0; k < layer.size_z; k++) {
                if (i > 0 && i < layer.size_x && j > 0 && j < layer.size_y && k > 0 && k < layer.size_z) {
                    layer.at(i, j , k) = 0;
                    layer.prev->at(i, j, k) = 0;
                } else {
                    layer.at(i, j , k) = InputData::border(i * hx(), j * hy(), k * hz());
                    layer.prev->at(i, j, k) = InputData::border(i * hx(), j * hy(), k * hz());
                }
                if ((layer.is_lower() && k == layer.size_z) || (layer.is_upper() && k == 0)) {
                    layer.at(i, j , k) = InputData::border(i * hx(), j * hy(), k * hz());
                    layer.prev->at(i, j, k) = InputData::border(i * hx(), j * hy(), k * hz());
                }
            }
        }
    }
    return layer;
}

double hx() {
    return InputData::X / InputData::cells_x;
}

double hy() {
    return InputData::Y / InputData::cells_y;
}

double hz() {
    return InputData::Z / InputData::cells_z;
}

int get_chunk_size(int given_rank, int size) {
    int basic_chunk = InputData::cells_z / size;
    int rest = InputData::cells_z % size;
    return basic_chunk + (given_rank < rest ? 1 : 0);
}

//
bool iteration_step(Layer& layer) {
    bool is_over = true;

    double Fi, Fj, Fk;

    double powhx = pow(hx(), 2);
    double powhy = pow(hy(), 2);
    double powhz = pow(hz(), 2);
    double* received_from_above = new double[layer.size_x * layer.size_y];
    double* received_from_below = new double[layer.size_x * layer.size_y];
    double* send_to_below = layer.get_row(layer.size_z - 1);
    double* send_to_above = layer.get_row(0);

    MPI_Request reqs1, reqs2, reqr1, reqr2;

    if (!layer.is_lower()) {
        MPI_Isend(send_to_above, layer.size_x * layer.size_y, MPI_DOUBLE, (layer.rank + 1) % layer.size, TALK_ABOVE,
                  MPI_COMM_WORLD, &reqs2);

    }

    if (!layer.is_upper()) {
        MPI_Isend(send_to_below, layer.size_x * layer.size_y, MPI_DOUBLE, (layer.rank - 1) % layer.size, TALK_BELOW, MPI_COMM_WORLD, &reqs1);

    }
    if (!layer.is_lower()) {
        MPI_Irecv(received_from_below, layer.size_x * layer.size_y, MPI_DOUBLE, (layer.rank + 1) % layer.size, TALK_BELOW, MPI_COMM_WORLD,
                  &reqr2);

    }
    if (!layer.is_upper()) {
        MPI_Irecv(received_from_above, layer.size_x * layer.size_y, MPI_DOUBLE, (layer.rank - 1) % layer.size, TALK_ABOVE, MPI_COMM_WORLD,
                  &reqr1);

    }

    if (!layer.is_lower()) {
        MPI_Wait(&reqs2, MPI_STATUS_IGNORE);
        MPI_Wait(&reqr2, MPI_STATUS_IGNORE);
    }
    if (!layer.is_upper()) {
        MPI_Wait(&reqs1, MPI_STATUS_IGNORE);
        MPI_Wait(&reqr1, MPI_STATUS_IGNORE);
    }


    double c = 2/powhx + 2/powhy + 2/powhz + InputData::a;
    for (int i = 1; i < layer.size_x - 1; i++) {
        for (int j = 1; j < layer.size_y - 1; j++) {
            for (int k = 0; k < layer.size_z; k++) {
                Fi = (layer.prev->at(i+1, j, k) + layer.prev->at(i-1, j, k)) / powhx;
                Fj = (layer.prev->at(i, j+1, k) + layer.prev->at(i, j-1, k)) / powhy;
                if (!layer.is_upper() && k == 0) {
                    Fk = (layer.prev->at(i, j, k + 1) + received_from_above[layer.size_x * i + j]) / powhz;
                }
                if (!layer.is_lower() && k == layer.size_z - 1) {
                    Fk = (layer.prev->at(i, j, k + 1) + received_from_below[layer.size_x * i + j]) / powhz;
                }
                if (k != 0 && k != layer.size_z - 1) {
                    Fk = (layer.prev->at(i, j, k+1) + layer.prev->at(i, j, k-1)) / powhz;
                }
                // After all assignments
                layer.copy_data_to_prev();
                layer.at(i, j, k) = (Fi + Fj + Fk - InputData::ro(i*hx(), j*hy(), k*hz())) /  c;
                if (fabs(layer.at(i, j , k) - layer.prev->at(i, j, k)) > InputData::e) {
                    is_over = false;
                }
            }
        }
    }
    return is_over;
}
