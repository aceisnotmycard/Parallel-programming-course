//
// Created by Sergey Bogolepov on 5/10/15.
//

#include <math.h>
#include "poisson.h"

int get_chunk_size(int given_rank);
int calculate_shift(int given_rank);

double InputData::border(double x, double y, double z) {
    return x + y + z;
}

double InputData::ro(int x, int y, int z) {
    return -a*(x+y+z);
}

Layer init() {
    int size_z = get_chunk_size(rank);
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

int get_chunk_size(int given_rank) {
    int basic_chunk = InputData::cells_z / size;
    int rest = InputData::cells_z % size;
    return basic_chunk + (given_rank < rest ? 1 : 0);
}

int calculate_shift(int given_rank) {
    int result = 0;
    for (int i = 0; i < given_rank; i++) {
        result += get_chunk_size(i);
    }
    return result;
}

bool iteration_step(Layer& layer) {
    bool is_over = true;

    double Fi, Fj, Fk;

    double powhx = pow(hx(), 2);
    double powhy = pow(hy(), 2);
    double powhz = pow(hz(), 2);

    double c = 2/powhx + 2/powhy + 2/powhz + InputData::a;

    // TODO: Init with MPI
    double* received_from_above;
    double* received_from_below;

    for (int i = 0; i < layer.size_x; i++) {
        for (int j = 0; j < layer.size_y; j++) {
            for (int k = 0; k < layer.size_z; k++) {
                if (!layer.is_upper() && k == 0) {
                    Fi = (layer.prev->at(i+1, j, k) + layer.prev->at(i-1, j, k)) / powhx;
                    Fj = (layer.prev->at(i, j+1, k) + layer.prev->at(i-1, j-1, k)) / powhy;
                    Fk = (layer.prev->at(i, j, k + 1) + received_from_above[layer.size_x * i + j]) / powhz;
                }
            }
        }
    }

    return is_over;
}
