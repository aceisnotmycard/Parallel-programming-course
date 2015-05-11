//
// Created by Sergey Bogolepov on 5/10/15.
//

#ifndef PARALLEL_POISSON_H
#define PARALLEL_POISSON_H

#include "layer.h"

#define TALK_ABOVE 100
#define TALK_BELOW 101

namespace InputData {
    double X;
    double Y;
    double Z;

    int cells_x;
    int cells_y;
    int cells_z;

    double e;

    int a;

    double ro(double x, double y, double z);

    double border(double x, double y, double z);
};

double hz();
double hy();
double hx();

Layer init(int rank, int size);

bool iteration_step(Layer& layer);

#endif //PARALLEL_POISSON_H
