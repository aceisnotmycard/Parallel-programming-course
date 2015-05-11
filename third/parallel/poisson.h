//
// Created by Sergey Bogolepov on 5/10/15.
//

#ifndef PARALLEL_POISSON_H
#define PARALLEL_POISSON_H

#include "layer.h"

#define TALK_ABOVE 100
#define TALK_BELOW 101

namespace InputData {
    static double X = 2.0;
    static double Y = 2.0;
    static double Z = 2.0;

    static int cells_x = 20;
    static int cells_y = 20;
    static int cells_z = 20;

    static double e = 0.0001;

    static int a = -1;

    double ro(double x, double y, double z);

    double border(double x, double y, double z);
};

double hz();
double hy();
double hx();

Layer init(int rank, int size);

bool iteration_step(Layer& layer);

#endif //PARALLEL_POISSON_H
