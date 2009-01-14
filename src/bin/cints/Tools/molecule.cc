/*! \file molecule.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>
#include<libint/libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static void get_geometry(void);

void init_molecule()
{
  Molecule.label = chkpt_rd_label();
  Molecule.num_atoms = chkpt_rd_natom();
  Molecule.Rref = chkpt_rd_rref();
  get_geometry();

  return;
}


void cleanup_molecule()
{
  free(Molecule.centers);
  free(Molecule.Rref);
  
  return;
}

void get_geometry()
{
   int i;
   double *Z;  /* nuclear charges */
   double **g; /* cartesian geometry */
   
   Molecule.centers = (struct coordinates *)malloc(sizeof(struct coordinates)*Molecule.num_atoms);

   g = chkpt_rd_geom();
   Z = chkpt_rd_zvals();

   /*--- move it into the appropriate struct form ---*/
   for (i=0; i<Molecule.num_atoms; i++){
      Molecule.centers[i].x = g[i][0];
      Molecule.centers[i].y = g[i][1];
      Molecule.centers[i].z = g[i][2];
      Molecule.centers[i].Z_nuc = Z[i];
   }
   free_block(g);
   free(Z);

   return;
}


/*!---------------------------------------
  enuc computes nuclear repulsion energy
 ---------------------------------------*/

void compute_enuc()
{  
  int i, j;
  double Z1Z2, r2, oor;
  double E = 0.0;

  if(Molecule.num_atoms > 1)
    for(i=1; i<Molecule.num_atoms; i++)
      for(j=0; j<i; j++){
	r2 = 0.0;
	r2 += (Molecule.centers[i].x-Molecule.centers[j].x)*
	     (Molecule.centers[i].x-Molecule.centers[j].x);
	r2 += (Molecule.centers[i].y-Molecule.centers[j].y)*
	     (Molecule.centers[i].y-Molecule.centers[j].y);
	r2 += (Molecule.centers[i].z-Molecule.centers[j].z)*
	     (Molecule.centers[i].z-Molecule.centers[j].z);
	oor = 1.0/sqrt(r2);
        Z1Z2 = Molecule.centers[i].Z_nuc*Molecule.centers[j].Z_nuc;
#ifdef HAVE_FUNC_ISINF
        if (std::isnan(oor) || std::isinf(oor)) {
#elif HAVE_FUNC_FINITE
        if (std::isnan(oor) || !std::finite(oor)) {
#endif
          if (fabs(Z1Z2) != 0.0)
            throw std::domain_error("compute_enuc -- charges too close to each other");
        }
        else
	  E += Z1Z2 * oor;
      }
  Molecule.Enuc = E;

  return;
}

}}
