/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>

#include <cstring>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libint/libint.h>

#include "defines.h"
#define EXTERN
#include "global.h"
#include <Tools/prints.h>

namespace psi { namespace CINTS {
void file11()
{
  int i;
  char wfnstring[80];
  char tmpwfn[80];
  double **Geom, **GradRef, **GeomRef;
  double Etot;
  FILE *fp11;

  Etot = chkpt_rd_etot();

  /*--- Geometry in the canonical frame ---*/
  Geom = block_matrix(Molecule.num_atoms,3);
  /*--- Geometry and gradient in the reference frame ---*/
  GeomRef = block_matrix(Molecule.num_atoms,3);
  GradRef = block_matrix(Molecule.num_atoms,3);
  
  /*--- Rotate back to the reference frame ---*/
  for(i=0;i<Molecule.num_atoms;i++) {
      Geom[i][0] = Molecule.centers[i].x;
      Geom[i][1] = Molecule.centers[i].y;
      Geom[i][2] = Molecule.centers[i].z;
  }
  mmult(Geom,0,Molecule.Rref,0,GeomRef,0,Molecule.num_atoms,3,3,0);
  mmult(Grad,0,Molecule.Rref,0,GradRef,0,Molecule.num_atoms,3,3,0);

  if (UserOptions.empirical_dispersion) {
    sprintf(tmpwfn, "%s+D", UserOptions.wfn);
  }
  else 
    sprintf(tmpwfn, "%s", UserOptions.wfn);

  sprintf(wfnstring,"%s forces in the reference frame (a.u.)", tmpwfn);
  print_atomvec(wfnstring,GradRef);
  ffile(&fp11, "file11.dat", 1);

  fprintf(fp11,"%-59.59s %-10.10s%-8.8s\n",Molecule.label, tmpwfn,
    UserOptions.dertype);

  fprintf(fp11,"%5d",Molecule.num_atoms);
  if (strcmp(UserOptions.wfn,"SCF") == 0 && !UserOptions.empirical_dispersion)
    fprintf(fp11,"%20.10lf\n",MOInfo.Escf);
  else
    fprintf(fp11,"%20.10lf\n",Etot); 
  
  for(i=0;i<Molecule.num_atoms;i++)
    fprintf(fp11,"%20.10lf%20.10lf%20.10lf%20.10lf\n",
	    Molecule.centers[i].Z_nuc,GeomRef[i][0],GeomRef[i][1],GeomRef[i][2]);
  for(i=0;i<Molecule.num_atoms;i++)
    fprintf(fp11,"                    %20.10lf%20.10lf%20.10lf\n",
	    GradRef[i][0],GradRef[i][1],GradRef[i][2]);
  fclose(fp11);

  free_block(Geom);
  free_block(GeomRef);
  free_block(GradRef);
  
  return;
}
};};
