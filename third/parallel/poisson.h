//
// Created by Sergey Bogolepov on 5/10/15.
//

#ifndef PARALLEL_POISSON_H
#define PARALLEL_POISSON_H

#include "layer.h"

int rank, size;

namespace InputData {
    double X;
    double Y;
    double Z;

    int cells_x;
    int cells_y;
    int cells_z;

    double e;

    int a;

    double ro(int x, int y, int z);

    double border(double x, double y, double z);
};

double hz();
double hy();
double hx();

Layer init();

bool iteration_step();

#endif //PARALLEL_POISSON_H
