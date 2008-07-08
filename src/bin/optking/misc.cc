
#include <libciomr/libciomr.h>

namespace psi { namespace optking {

// return unit matrix
double **unit_matrix(int dim) {
    double **u = block_matrix(dim,dim);
    for (int i=0; i<dim; ++i)
      u[i][i] = 1.0;
    return u;
}

// put inverse masses on diagonal
double **mass_matrix(int dim, double *masses) {
    double **u = block_matrix(dim,dim);
    for (int i=0; i<dim; i+=3) {
      if (masses[i] == 0.0)
        u[i][i] = 0.0;
      else
        u[i][i] = 1.0/masses[i];
    }
    return u;
}


}}
