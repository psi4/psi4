/*! \file prints.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/

#include<cstdio>
#include<cstdlib>

#include<libint/libint.h>
#include<libqt/qt.h>
#include<libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi {
  namespace CINTS {

void print_intro()
{
  if (UserOptions.print_lvl >= PRINT_INTRO) {
     tstart(outfile);
     fprintf(outfile,"                  --------------------------------------------\n");
     fprintf(outfile,"                    CINTS: An integrals program written in C\n");
     fprintf(outfile,"                     Justin T. Fermann and Edward F. Valeev\n");
     fprintf(outfile,"                  --------------------------------------------\n\n");
  }

  return;
}


void print_scalars()
{
  int i;

  /*----------------------
     Print constants, etc.
    ----------------------*/
   if (UserOptions.print_lvl >= PRINT_OPTIONS) {
     fprintf(outfile,"\n  -OPTIONS:\n");
     fprintf(outfile,"    Print level                 = %d\n",UserOptions.print_lvl);
     if (UserOptions.restart)
       fprintf(outfile,"    Restart                     = yes\n");
     fprintf(outfile,"    Integral tolerance          = %1.0e\n",UserOptions.cutoff);
     fprintf(outfile,"    Max. memory to use          = %ld double words\n",UserOptions.max_memory);
     fprintf(outfile,"    Number of threads           = %d\n",UserOptions.num_threads);
     fprintf(outfile,"    LIBINT's real type length   = %d bit\n",sizeof(REALTYPE)*8);
     fprintf(outfile,"\n  -CALCULATION CONSTANTS:\n");
     fprintf(outfile,"    Label                       = %s\n", Molecule.label);
     fprintf(outfile,"    Number of atoms             = %d\n", Molecule.num_atoms);
     if (UserOptions.print_lvl >= PRINT_BASIS) {
       fprintf(outfile,"    Number of shells            = %d\n", BasisSet.num_shells);
       fprintf(outfile,"    Number of primitives        = %d\n", BasisSet.num_prims);
     }
     fprintf(outfile,"    Number of atomic orbitals   = %d\n", BasisSet.num_ao);
     fprintf(outfile,"    Number of symmetry orbitals = %d\n", Symmetry.num_so);
     fprintf(outfile,"    Maximum AM in the basis     = %d\n", BasisSet.max_am - 1);
     if (UserOptions.fine_structure_alpha != 1.0)
       fprintf(outfile,"    Fine-structure alpha scaling = %8.4lf\n", UserOptions.fine_structure_alpha);
     fprintf(outfile,"\n  -SYMMETRY INFORMATION;\n");
     fprintf(outfile,"    Computational point group        = %s\n", Symmetry.symlabel);
     fprintf(outfile,"    Number of irreps                 = %d\n", Symmetry.nirreps);
     if (UserOptions.print_lvl >= PRINT_BASIS) {
       fprintf(outfile,"    Number of symmetry unique shells = %d\n", Symmetry.num_unique_shells);
     }
   }
   if (UserOptions.print_lvl >= PRINT_GEOMETRY) {
     fprintf(outfile,"\n  -CARTESIAN COORDINATES (a.u.):\n");
     fprintf(outfile,"     Nuc. charge           X                  Y                   Z\n");
     fprintf(outfile,"    -------------   -----------------  -----------------  -----------------\n");

     for(i=0;i<Molecule.num_atoms;i++){
       fprintf(outfile,"      %8.3lf    ",Molecule.centers[i].Z_nuc);
       fprintf(outfile,"  %17.12lf  %17.12lf  %17.12lf",
	       Molecule.centers[i].x, Molecule.centers[i].y, Molecule.centers[i].z);
       fprintf(outfile,"\n");
     }
     fprintf(outfile,"\n");
   }
   fflush(outfile);

   return;
}


void print_basisset()
{
  int i, j;

  /*---------------------
     Print basis set info
    ---------------------*/
   if (UserOptions.print_lvl >= PRINT_BASIS) {
     fprintf(outfile,"\n  -BASIS SET INFORMATION:\n");
     fprintf(outfile,"    Prim#     Exponent   ");
     for(i=0;i<BasisSet.max_am;i++)
       fprintf(outfile,"  Norm. CCoeff. (L=%1d)",i);
     fprintf(outfile,"\n");
     fprintf(outfile,"    -----  --------------");
     for(i=0;i<BasisSet.max_am;i++)
       fprintf(outfile,"  -------------------");
     fprintf(outfile,"\n");
     for(i=0;i<BasisSet.num_prims;i++) {
       fprintf(outfile,"    %3d    %14.7lf",i+1,BasisSet.cgtos[i].exp);
       for(j=0;j<BasisSet.max_am;j++)
	 fprintf(outfile,"    %15.10lf  ",BasisSet.cgtos[i].ccoeff[j]);
       fprintf(outfile,"\n");
     }
     fprintf(outfile,"\n\n");
     fprintf(outfile,"    Shell#  Nuc#  L  SPRIM  SNUMG\n");
     fprintf(outfile,"    ------  ----  -  -----  -----\n");
     for(i=0;i<BasisSet.num_shells;i++)
       fprintf(outfile,"    %4d    %3d  %2d  %3d    %3d\n",i+1,
	       BasisSet.shells[i].center,BasisSet.shells[i].am-1,
	       BasisSet.shells[i].fprim,BasisSet.shells[i].n_prims);
     fprintf(outfile,"\n\n");
     fflush(outfile);
   }

   return;
}


void print_quote()
{
  /*--- It's about recurrence relations, well, sort of... ---*/
   if (UserOptions.print_lvl) {
     fprintf(outfile,"           ---------------------------------------------\n");
     fprintf(outfile,"                 'I will find a center in you.\n");
     fprintf(outfile,"                  I will chew it up and leave,\n");
     fprintf(outfile,"                   I will work to elevate you\n");
     fprintf(outfile,"                 just enough to bring you down.'\n");
     fprintf(outfile,"                                       'Sober', Tool\n");
     fprintf(outfile,"           ---------------------------------------------\n\n");
   }

   return;
}


void print_opdm()
{
  if (UserOptions.print_lvl >= PRINT_OPDM) {
    fprintf(outfile,"  -Total density matrix in AO basis :\n");
    print_mat(Dens,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"\n\n");
    fprintf(outfile,"  -Energy weighted density matrix in AO basis :\n");
    print_mat(Lagr,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"\n\n");
  }

  return;
}


/*!----------------------------------------------------
  May be used to print out gradients and other vector
  quantities associated with atoms
 ----------------------------------------------------*/
void print_atomvec(char *quantity, double **vecs)
{
  int i;

  fprintf(outfile,"\n  -%s:\n",quantity);
  fprintf(outfile,"     Atom            X                  Y                   Z\n");
  fprintf(outfile,"    ------   -----------------  -----------------  -----------------\n");

  for(i=0;i<Molecule.num_atoms;i++){
    fprintf(outfile,"    %4d   ",i+1);
    fprintf(outfile,"  %17.12lf  %17.12lf  %17.12lf",vecs[i][0], vecs[i][1], vecs[i][2]);
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  
  return;
}

/*!---------------------------------------------------------
  May be used to print out Hessian and other square matrix
  quantities associated with atoms

  N.B. for now uses print_mat()
 ---------------------------------------------------------*/
void print_atommat(char *quantity, double **mat)
{
  int i;
  
  fprintf(outfile,"\n  -%s:\n",quantity);

  /*  fprintf(outfile,"     Atom            X                  Y                   Z\n");
  fprintf(outfile,"    ------   -----------------  -----------------  -----------------\n");
  
  for(i=0;i<Molecule.num_atoms;i++){
    fprintf(outfile,"    %4d   ",i+1);
    fprintf(outfile,"  %17.12lf  %17.12lf  %17.12lf",vecs[i][0], vecs[i][1], vecs[i][2]);
    fprintf(outfile,"\n");
    }*/

  print_mat(mat,Molecule.num_atoms*3,Molecule.num_atoms*3,outfile);
  fprintf(outfile,"\n");
  
  return;
}

void print_moinfo_corr()
{
  int i;
  
  if (UserOptions.print_lvl >= PRINT_MOINFO_CORR) {
    fprintf(outfile,"\n  -MO information:\n");
    fprintf(outfile,"    Number of MOs    = %d\n",MOInfo.num_mo);
    fprintf(outfile,
	    "    Label\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
    fprintf(outfile,
	    "    -----\t------\t------\t------\t------\t------\n");
    for(i=0; i < Symmetry.nirreps; i++) {
	fprintf(outfile,
		"     %s\t   %d\t   %d\t    %d\t    %d\t    %d\n",
		Symmetry.irr_labels[i],MOInfo.frozen_docc[i],
		MOInfo.clsdpi[i],MOInfo.openpi[i],MOInfo.virtpi[i],
		MOInfo.frozen_uocc[i]);
    }
    fprintf(outfile,"\n");
  }
  fflush(outfile);

  return;
}
};};
