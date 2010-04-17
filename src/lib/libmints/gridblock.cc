#include "gridblock.h"
#include "string.h"
#include <libciomr/libciomr.h>
#include <stdlib.h>

using namespace psi; 

GridBlock::GridBlock(int n) 
{
    max_points_ = n;
    x_ = init_array(n);
    y_ = init_array(n);
    z_ = init_array(n);
    w_ = init_array(n);
    true_points_ = 0;
}
GridBlock::~GridBlock()
{
    free(x_);
    free(y_);
    free(z_);
    free(w_);
}
void GridBlock::setGrid(double *x, double *y, double *z, double* w, int n)
{
    true_points_ = n;
    memcpy((void*) x_, (void*) x, n);
    memcpy((void*) y_, (void*) y, n);
    memcpy((void*) z_, (void*) z, n);
    memcpy((void*) w_, (void*) w, n);
}



