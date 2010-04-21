#include "gridblock.h"
#include "string.h"
#include <exception.h>
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
    if (n > max_points_) {
        throw SanityCheckError("GridBlock::setGrid: n > max_points_", __FILE__, __LINE__);
    }
    size_t length = n * sizeof(double);
    true_points_ = n;
    memcpy((void*) x_, (void*) x, length);
    memcpy((void*) y_, (void*) y, length);
    memcpy((void*) z_, (void*) z, length);
    memcpy((void*) w_, (void*) w, length);
}



