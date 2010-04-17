#include <libmints/vector3.h>
#include "integrator.h"
using namespace psi;
using namespace std;

namespace psi { namespace scf {

GridBlock::GridBlock(int max_points)
{
    max_points_ = max_points;
    w_ = init_array(max_points_);
    x_ = init_array(max_points_);
    y_ = init_array(max_points_);
    z_ = init_array(max_points_);
    true_points_ = 0;
}
GridBlock::~GridBlock()
{
    free(w_);
    free(x_);
    free(y_);
    free(z_);
}
void GridBlock::setGrid(double* x, double* y, double* z, double* w, int n)
{
    memcpy((void*) x_, (void*) x, n);   
    memcpy((void*) y_, (void*) y, n);   
    memcpy((void*) z_, (void*) z, n);   
    memcpy((void*) w_, (void*) w, n);
    true_points_ = n;   
}

}}

