/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <libciomr/libciomr.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "input.h"
#include "global.h"
#include "defines.h"
#include <physconst.h>
#include <symmetry.h>

namespace {
  void alloc_irr_char();
  void permute_axes(double **, int, int);
}

namespace psi { namespace input {

void find_symmetry()
{
  int i, j;
  double tmp;
  double *v1;
  int **redun_atom_orbit;     /*Redundant (i.e. for all 8 operations) atomic orbits*/
  int c2x_flag = 1;
  int c2y_flag = 1;
  int c2z_flag = 1;
  int inv_flag = 1;
  int sig_xy_flag = 1;
  int sig_xz_flag = 1;
  int sig_yz_flag = 1;
  int naxes, nplanes;
  int *row;

  if (num_atoms == 1) { /* Atomic calculation */
/*    symmetry = strdup("D2h ");
    nirreps = 8;
    irr_labels = (char **) malloc(sizeof(char *)*nirreps);
    sym_oper = init_int_array(nirreps);
    alloc_irr_char();
    irr_labels[0] = strdup("Ag  ");
    irr_labels[1] = strdup("B1g ");
    irr_labels[2] = strdup("B2g ");
    irr_labels[3] = strdup("B3g ");
    irr_labels[4] = strdup("Au  ");
    irr_labels[5] = strdup("B1u ");
    irr_labels[6] = strdup("B2u ");
    irr_labels[7] = strdup("B3u ");
    sym_oper[1] = C2ZFLAG;
    sym_oper[2] = C2YFLAG;
    sym_oper[3] = C2XFLAG;
    sym_oper[4] = IFLAG;
    sym_oper[5] = SIGXYFLAG;
    sym_oper[6] = SIGXZFLAG;
    sym_oper[7] = SIGYZFLAG;
    atom_orbit = init_int_matrix(1,nirreps);
    atom_orbit[0][0] = atom_orbit[0][1] = atom_orbit[0][2] = atom_orbit[0][3] = atom_orbit[0][4] =
    atom_orbit[0][5] = atom_orbit[0][6] = atom_orbit[0][7] = 0;
    return; */
  }

  v1 = init_array(3);
  redun_atom_orbit = init_int_matrix(num_atoms,8);

  for(i=0;i<num_atoms;i++) {
    if (c2x_flag && !redun_atom_orbit[i][C2XFLAG]) {
      c2x(geometry[i],v1);
      c2x_flag = 0;
      for(j=0;j<num_atoms;j++)
        if ((int)nuclear_charges[i] == (int)nuclear_charges[j] && vect_equiv(geometry[j],v1)){
          redun_atom_orbit[i][C2XFLAG] = j;
	  redun_atom_orbit[j][C2XFLAG] = i;
	  c2x_flag = 1;
	  break;
	}
    }
    if (c2y_flag && !redun_atom_orbit[i][C2YFLAG]) {
      c2y(geometry[i],v1);
      c2y_flag = 0;
      for(j=0;j<num_atoms;j++)
        if ((int)nuclear_charges[i] == (int)nuclear_charges[j] && vect_equiv(geometry[j],v1)){
          redun_atom_orbit[i][C2YFLAG] = j;
	  redun_atom_orbit[j][C2YFLAG] = i;
	  c2y_flag = 1;
	  break;
	}
    }
    if (c2z_flag && !redun_atom_orbit[i][C2ZFLAG]) {
      c2z(geometry[i],v1);
      c2z_flag = 0;
      for(j=0;j<num_atoms;j++)
        if ((int)nuclear_charges[i] == (int)nuclear_charges[j] && vect_equiv(geometry[j],v1)){
          redun_atom_orbit[i][C2ZFLAG] = j;
	  redun_atom_orbit[j][C2ZFLAG] = i;
	  c2z_flag = 1;
	  break;
	}
    }
    if (inv_flag && !redun_atom_orbit[i][IFLAG]) {
      inversion(geometry[i],v1);
      inv_flag = 0;
      for(j=0;j<num_atoms;j++)
        if ((int)nuclear_charges[i] == (int)nuclear_charges[j] && vect_equiv(geometry[j],v1)){
          redun_atom_orbit[i][IFLAG] = j;
	  redun_atom_orbit[j][IFLAG] = i;
	  inv_flag = 1;
	  break;
	}
    }
    if (sig_xy_flag && !redun_atom_orbit[i][SIGXYFLAG]) {
      sig_xy(geometry[i],v1);
      sig_xy_flag = 0;
      for(j=0;j<num_atoms;j++)
        if ((int)nuclear_charges[i] == (int)nuclear_charges[j] && vect_equiv(geometry[j],v1)){
          redun_atom_orbit[i][SIGXYFLAG] = j;
	  redun_atom_orbit[j][SIGXYFLAG] = i;
	  sig_xy_flag = 1;
	  break;
	}
    }
    if (sig_xz_flag && !redun_atom_orbit[i][SIGXZFLAG]) {
      sig_xz(geometry[i],v1);
      sig_xz_flag = 0;
      for(j=0;j<num_atoms;j++)
        if ((int)nuclear_charges[i] == (int)nuclear_charges[j] && vect_equiv(geometry[j],v1)){
          redun_atom_orbit[i][SIGXZFLAG] = j;
	  redun_atom_orbit[j][SIGXZFLAG] = i;
	  sig_xz_flag = 1;
	  break;
	}
    }
    if (sig_yz_flag && !redun_atom_orbit[i][SIGYZFLAG]) {
      sig_yz(geometry[i],v1);
      sig_yz_flag = 0;
      for(j=0;j<num_atoms;j++)
        if ((int)nuclear_charges[i] == (int)nuclear_charges[j] && vect_equiv(geometry[j],v1)){
          redun_atom_orbit[i][SIGYZFLAG] = j;
	  redun_atom_orbit[j][SIGYZFLAG] = i;
	  sig_yz_flag = 1;
	  break;
	}
    }
  }

  free(v1);

  if (print_lvl >= DEBUGPRINT) {
    fprintf(outfile,"\n  -Symmetry flags:\n");
    fprintf(outfile,"    C2x = %d\n",c2x_flag);
    fprintf(outfile,"    C2y = %d\n",c2y_flag);
    fprintf(outfile,"    C2z = %d\n",c2z_flag);
    fprintf(outfile,"    I   = %d\n",inv_flag);
    fprintf(outfile,"    sxy = %d\n",sig_xy_flag);
    fprintf(outfile,"    sxz = %d\n",sig_xz_flag);
    fprintf(outfile,"    syz = %d\n",sig_yz_flag);
  }
  naxes = c2x_flag + c2y_flag + c2z_flag;
  nplanes = sig_xy_flag + sig_xz_flag + sig_yz_flag;
  nirreps = 1 + naxes + inv_flag + nplanes;


  /*-------------------
    Handling subgroups
   --------------------*/

  if (subgroup != NULL) {
    switch (nirreps) {
      case 8: /*Can be C2v, C2h, D2, etc.*/
	      if (!strcmp(subgroup,"D2")) {
		nirreps = 4; nplanes = 0; inv_flag = 0;
		sig_xy_flag = 0; sig_xz_flag = 0; sig_yz_flag = 0;
		break;
	      }
	      if (!strcmp(subgroup,"C2V")) {
		if (unique_axis == NULL) {
		  printf("  UNIQUE_AXIS must be specified for this symmetry.\n  The largest Abelian point group will be used.\n\n");
		  break;
		}
		nirreps = 4; naxes = 1; nplanes = 2; inv_flag = 0;
		if (!strcmp(unique_axis,"X")) {
		  c2y_flag = 0; c2z_flag = 0; sig_yz_flag = 0;
		}
		else if (!strcmp(unique_axis,"Y")) {
		  c2x_flag = 0; c2z_flag = 0; sig_xz_flag = 0;
		}
		else if (!strcmp(unique_axis,"Z")) {
		  c2x_flag = 0; c2y_flag = 0; sig_xy_flag = 0;
		}
		break;
	      }
	      if (!strcmp(subgroup,"C2H")) {
		if (unique_axis == NULL) {
		  printf("  UNIQUE_AXIS must be specified for this symmetry.\n  The largest Abelian point group will be used.\n\n");
		  break;
		}
		nirreps = 4; naxes = 1; nplanes = 1; inv_flag = 1;
		if (!strcmp(unique_axis,"X")) {
		  c2y_flag = 0; c2z_flag = 0; sig_xy_flag = 0; sig_xz_flag = 0;
		}
		else if (!strcmp(unique_axis,"Y")) {
		  c2x_flag = 0; c2z_flag = 0; sig_xy_flag = 0; sig_yz_flag = 0;
		}
		else if (!strcmp(unique_axis,"Z")) {
		  c2x_flag = 0; c2y_flag = 0; sig_xz_flag = 0; sig_yz_flag = 0;
		}
		break;
	      }
      case 4: /*Can be C2, Cs, Ci, etc.*/
	      if (inv_flag && !strcmp(subgroup,"CI")) {
		nirreps = 2; naxes = 0; nplanes = 0;
		c2x_flag = 0; c2y_flag = 0; c2z_flag = 0;
		sig_xy_flag = 0; sig_xz_flag = 0; sig_yz_flag = 0;
		break;
	      }
	      if ((naxes > 0) && !strcmp(subgroup,"C2")) {
		if (naxes == 3) {
		  if (unique_axis == NULL) {
		    printf("  UNIQUE_AXIS must be specified for this symmetry.\n  The largest Abelian point group will be used.\n\n");
		    break;
		  }
		  if (!strcmp(unique_axis,"X")) {
		    c2y_flag = 0; c2z_flag = 0;
		  }
		  else if (!strcmp(unique_axis,"Y")) {
		    c2x_flag = 0; c2z_flag = 0;
		  }
		  else if (!strcmp(unique_axis,"Z")) {
		    c2x_flag = 0; c2y_flag = 0;
		  }
		  naxes = 1;
		}
		nirreps = 2; nplanes = 0; inv_flag = 0;
		sig_xy_flag = 0; sig_xz_flag = 0; sig_yz_flag = 0;
		break;
	      }
	      if ((nplanes > 0) && !strcmp(subgroup,"CS")) {
		if (nplanes > 1) {
		  if (unique_axis == NULL) {
		    printf("  UNIQUE_AXIS must be specified for this symmetry.\n  The largest Abelian point group will be used.\n\n");
		    break;
		  }
		  if (!strcmp(unique_axis,"Z"))
		    if (sig_xy_flag) {
		      sig_xz_flag = 0; sig_yz_flag = 0;
		    }
		    else
		      punt("  UNIQUE_AXIS defined incorrectly.");
		  else if (!strcmp(unique_axis,"Y"))
		    if (sig_xz_flag) {
		      sig_xy_flag = 0; sig_yz_flag = 0;
		    }
		    else
		      punt("  UNIQUE_AXIS defined incorrectly.");
		  else if (!strcmp(unique_axis,"X"))
		    if (sig_yz_flag) {
		      sig_xy_flag = 0; sig_xz_flag = 0;
		    }
		    else
		      punt("  UNIQUE_AXIS defined incorrectly.");
		  nplanes = 1;
		}
		nirreps = 2; naxes = 0; inv_flag = 0;
		c2x_flag = 0; c2y_flag = 0; c2z_flag = 0;
		break;
	      }
      case 2: /*Can only be C1*/
  	      if (!strcmp(subgroup,"C1")) {
		nirreps = 1; naxes = 0; nplanes = 0;
		c2x_flag = 0; c2y_flag = 0; c2z_flag = 0; inv_flag = 0;
		sig_xy_flag = 0; sig_xz_flag = 0; sig_yz_flag = 0;
		break;
	      }
      case 1: /*What subgroup???*/
  	      if (strcmp(subgroup,"C1")) {
		printf("  Subgroup specification doesn't comply with the actual symmetry.\n");
		printf("  The largest Abelian point group will be used.\n\n");
	      }
	      break;
    }
    if (print_lvl >= DEBUGPRINT) {
      fprintf(outfile,"\n  -Flags within the subgroup:\n");
      fprintf(outfile,"    C2x = %d\n",c2x_flag);
      fprintf(outfile,"    C2y = %d\n",c2y_flag);
      fprintf(outfile,"    C2z = %d\n",c2z_flag);
      fprintf(outfile,"    I   = %d\n",inv_flag);
      fprintf(outfile,"    sxy = %d\n",sig_xy_flag);
      fprintf(outfile,"    sxz = %d\n",sig_xz_flag);
      fprintf(outfile,"    syz = %d\n",sig_yz_flag);
    }
  }

  
  /*-------------------------
    Allocate character table
   -------------------------*/

  alloc_irr_char();
  irr_labels = (char **) malloc(sizeof(char *)*nirreps);

  
  /*---------------------------------------------------------------------
    Figure out point group, reorient, if necessary, and get rid of zeros
    in the array of atomic orbits
   ---------------------------------------------------------------------*/
  
  atom_orbit = init_int_matrix(num_atoms,nirreps);
  for(i=0;i<num_atoms;i++)
    atom_orbit[i][EFLAG] = i;
  sym_oper = init_int_array(nirreps);
  sym_oper[0] = EFLAG;
  switch (nirreps) {
    case 1: /* There's one irreducible representation - it is C1 */
	    symmetry = strdup("C1  ");
	    irr_labels[0] = strdup("A   ");
	    break;
    case 2: /* Could be Ci, C2, or Cs */
	    if (inv_flag) {
	      symmetry = strdup("Ci  ");
	      irr_labels[0] = strdup("Ag  ");
	      irr_labels[1] = strdup("Au  ");
	      sym_oper[1] = IFLAG;
	      for(i=0;i<num_atoms;i++)
		atom_orbit[i][1] = redun_atom_orbit[i][IFLAG]; /* i */
            }
            else if (naxes) {
	      symmetry = strdup("C2  ");
	      irr_labels[0] = strdup("A   ");
	      irr_labels[1] = strdup("B   ");
	      sym_oper[1] = C2ZFLAG;
	      /* Making sure that the symmetry axis is along Z */
	      if (c2z_flag)
		for(i=0;i<num_atoms;i++)
		  atom_orbit[i][1] = redun_atom_orbit[i][C2ZFLAG]; /* C2 = C2z */
	      else
		if (c2x_flag) { /* Unique axis is along X */
		  permute_axes(geometry,XAXIS,ZAXIS);
		  for(i=0;i<num_atoms;i++) {
		    atom_orbit[i][1] = redun_atom_orbit[i][C2XFLAG]; /* C2 = C2x */
		  }
		}
		else { /* Unique axis is along Y */
		  permute_axes(geometry,YAXIS,ZAXIS);
		  for(i=0;i<num_atoms;i++) {
		    atom_orbit[i][1] = redun_atom_orbit[i][C2YFLAG]; /* C2 = C2y */
		  }
		}
	    }
            else if (nplanes) { /* Cs symmetry */
	      symmetry = strdup("Cs  ");
	      irr_labels[0] = strdup("Ap  ");
	      irr_labels[1] = strdup("App ");
	      sym_oper[1] = SIGXYFLAG;
	      if (sig_xy_flag)
		for(i=0;i<num_atoms;i++)
		  atom_orbit[i][1] = redun_atom_orbit[i][SIGXYFLAG]; /* sig_xy */
	      else
		if (sig_yz_flag) { /* Unique axis is along X */
		  permute_axes(geometry,XAXIS,ZAXIS);
		  for(i=0;i<num_atoms;i++) {
		    atom_orbit[i][1] = redun_atom_orbit[i][SIGYZFLAG]; /* sig_xy = sig_yz */
		  }
		}
		else { /* Unique axis is along Y */
		  permute_axes(geometry,YAXIS,ZAXIS);
		  for(i=0;i<num_atoms;i++) {
		    atom_orbit[i][1] = redun_atom_orbit[i][SIGXZFLAG]; /* sig_xy = sig_xz */
		  }
		}
	    }
	    else
	      punt("Unidentified symmetry.");
	    break;
	    
    case 4: /* It could be C2v, C2h, or D2 */
	    if (inv_flag) {
	      symmetry = strdup("C2h ");
	      /* Swap rows 1 and 2 in irr_char since ordering of irreps in C2h is different from that in C2v and D2*/
	      row = irr_char[1];
	      irr_char[1] = irr_char[2];
	      irr_char[2] = row;
	      irr_labels[0] = strdup("Ag  ");
	      irr_labels[1] = strdup("Bg  ");
	      irr_labels[2] = strdup("Au  ");
	      irr_labels[3] = strdup("Bu  ");
	      sym_oper[1] = C2ZFLAG;
	      sym_oper[2] = IFLAG;
	      sym_oper[3] = SIGXYFLAG;
	      if (c2z_flag)
		for(i=0;i<num_atoms;i++) {
		  atom_orbit[i][1] = redun_atom_orbit[i][C2ZFLAG]; /* C2 */
		  atom_orbit[i][2] = redun_atom_orbit[i][IFLAG]; /* i */
		  atom_orbit[i][3] = redun_atom_orbit[i][SIGXYFLAG]; /* sig_xy */
		}
	      else
		if (c2x_flag) { /* Unique axis is along X */
		  permute_axes(geometry,XAXIS,ZAXIS);
		  for(i=0;i<num_atoms;i++) {
		    atom_orbit[i][1] = redun_atom_orbit[i][C2XFLAG]; /* C2 = C2x */
		    atom_orbit[i][2] = redun_atom_orbit[i][IFLAG]; /* i */
		    atom_orbit[i][3] = redun_atom_orbit[i][SIGYZFLAG]; /* sig_xy = sig_yz */
		  }
		}
		else { /* Unique axis is along Y */
		  permute_axes(geometry,YAXIS,ZAXIS);
		  for(i=0;i<num_atoms;i++) {
		    atom_orbit[i][1] = redun_atom_orbit[i][C2YFLAG]; /* C2 = C2y */
		    atom_orbit[i][2] = redun_atom_orbit[i][IFLAG]; /* i */
		    atom_orbit[i][3] = redun_atom_orbit[i][SIGXZFLAG]; /* sig_xy = sig_xz */
		  }
		}
	    }
	    else if (nplanes == 2) {
	      symmetry = strdup("C2v ");
	      irr_labels[0] = strdup("A1  ");
	      irr_labels[1] = strdup("A2  ");
	      irr_labels[2] = strdup("B1  ");
	      irr_labels[3] = strdup("B2  ");
	      sym_oper[1] = C2ZFLAG;
	      sym_oper[2] = SIGXZFLAG;
	      sym_oper[3] = SIGYZFLAG;
	      if (c2z_flag)
		for(i=0;i<num_atoms;i++) {
		  atom_orbit[i][1] = redun_atom_orbit[i][C2ZFLAG]; /* C2 */
		  atom_orbit[i][2] = redun_atom_orbit[i][SIGXZFLAG]; /* sig_xz */
		  atom_orbit[i][3] = redun_atom_orbit[i][SIGYZFLAG]; /* sig_yz */
		}
	      else
		if (c2x_flag) { /* Unique axis is along X */
		  permute_axes(geometry,XAXIS,ZAXIS);
		  for(i=0;i<num_atoms;i++) {
		    atom_orbit[i][1] = redun_atom_orbit[i][C2XFLAG]; /* C2 = C2x */
		    atom_orbit[i][2] = redun_atom_orbit[i][SIGXYFLAG]; /* sig_xz = sig_yx */
		    atom_orbit[i][3] = redun_atom_orbit[i][SIGXZFLAG]; /* sig_yz = sig_zx */
		  }
		}
		else { /* Unique axis is along Y */
		  permute_axes(geometry,YAXIS,ZAXIS);
		  for(i=0;i<num_atoms;i++) {
		    atom_orbit[i][1] = redun_atom_orbit[i][C2YFLAG]; /* C2 = C2y */
		    atom_orbit[i][2] = redun_atom_orbit[i][SIGYZFLAG]; /* sig_xz = sig_zy */
		    atom_orbit[i][3] = redun_atom_orbit[i][SIGXYFLAG]; /* sig_yz = sig_xy */
		  }
		}
	    }
	    else if (naxes == 3) { /* D2 symmetry */
	      symmetry = strdup("D2  ");
	      irr_labels[0] = strdup("A   ");
	      irr_labels[1] = strdup("B1  ");
	      irr_labels[2] = strdup("B2  ");
	      irr_labels[3] = strdup("B3  ");
	      sym_oper[1] = C2ZFLAG;
	      sym_oper[2] = C2YFLAG;
	      sym_oper[3] = C2XFLAG;
	      for(i=0;i<num_atoms;i++) {
		atom_orbit[i][1] = redun_atom_orbit[i][C2ZFLAG]; /* C2z */
		atom_orbit[i][2] = redun_atom_orbit[i][C2YFLAG]; /* C2y */
		atom_orbit[i][3] = redun_atom_orbit[i][C2XFLAG]; /* C2x */
	      }
	    }
	    else
	      punt("Unidentified symmetry.");
	    break;
	    
    case 8: /* D2h symmetry */
            if (nplanes == 3 && naxes == 3 && inv_flag) {
	      symmetry = strdup("D2h ");
	      irr_labels[0] = strdup("Ag  ");
	      irr_labels[1] = strdup("B1g ");
	      irr_labels[2] = strdup("B2g ");
	      irr_labels[3] = strdup("B3g ");
	      irr_labels[4] = strdup("Au  ");
	      irr_labels[5] = strdup("B1u ");
	      irr_labels[6] = strdup("B2u ");
	      irr_labels[7] = strdup("B3u ");
	      sym_oper[1] = C2ZFLAG;
	      sym_oper[2] = C2YFLAG;
	      sym_oper[3] = C2XFLAG;
	      sym_oper[4] = IFLAG;
	      sym_oper[5] = SIGXYFLAG;
	      sym_oper[6] = SIGXZFLAG;
	      sym_oper[7] = SIGYZFLAG;
	      for(i=0;i<num_atoms;i++) {
		atom_orbit[i][1] = redun_atom_orbit[i][C2ZFLAG]; /* C2z */
		atom_orbit[i][2] = redun_atom_orbit[i][C2YFLAG]; /* C2y */
		atom_orbit[i][3] = redun_atom_orbit[i][C2XFLAG]; /* C2x */
		atom_orbit[i][4] = redun_atom_orbit[i][IFLAG]; /* i */
		atom_orbit[i][5] = redun_atom_orbit[i][SIGXYFLAG]; /* sig_xy */
		atom_orbit[i][6] = redun_atom_orbit[i][SIGXZFLAG]; /* sig_xz */
		atom_orbit[i][7] = redun_atom_orbit[i][SIGYZFLAG]; /* sig_yz */
	      }
	    }
	    else
	      punt("Unidentified symmetry.");
	    break;
	    
    default: /* Error */
	    punt("Unidentified symmetry.");
  } /* end of switch(nirreps) */

  free_int_matrix(redun_atom_orbit);
  
  return;
}
  
}} // namespace psi::input

namespace {

using namespace psi;
using namespace psi::input;

void alloc_irr_char()
{
  int i;

  /* Allocate global arrays */
  irr_char = init_int_matrix(nirreps,nirreps);

  for(i=0;i<nirreps;i++)
    irr_char[0][i] = 1;
  switch (nirreps) {
    case 1:
      break;
    case 2:
      irr_char[1][0] =  1; irr_char[1][1] = -1;
      break;
    case 4:
      irr_char[1][0] =  1; irr_char[1][1] =  1; irr_char[1][2] = -1; irr_char[1][3] = -1;
      irr_char[2][0] =  1; irr_char[2][1] = -1; irr_char[2][2] =  1; irr_char[2][3] = -1;
      irr_char[3][0] =  1; irr_char[3][1] = -1; irr_char[3][2] = -1; irr_char[3][3] =  1;
      break;
    case 8:
      irr_char[1][0] =  1; irr_char[1][1] =  1; irr_char[1][2] = -1; irr_char[1][3] = -1;
      irr_char[1][4] =  1; irr_char[1][5] =  1; irr_char[1][6] = -1; irr_char[1][7] = -1;
      irr_char[2][0] =  1; irr_char[2][1] = -1; irr_char[2][2] =  1; irr_char[2][3] = -1;
      irr_char[2][4] =  1; irr_char[2][5] = -1; irr_char[2][6] =  1; irr_char[2][7] = -1;
      irr_char[3][0] =  1; irr_char[3][1] = -1; irr_char[3][2] = -1; irr_char[3][3] =  1;
      irr_char[3][4] =  1; irr_char[3][5] = -1; irr_char[3][6] = -1; irr_char[3][7] =  1;
      irr_char[4][0] =  1; irr_char[4][1] =  1; irr_char[4][2] =  1; irr_char[4][3] =  1;
      irr_char[4][4] = -1; irr_char[4][5] = -1; irr_char[4][6] = -1; irr_char[4][7] = -1;
      irr_char[5][0] =  1; irr_char[5][1] =  1; irr_char[5][2] = -1; irr_char[5][3] = -1;
      irr_char[5][4] = -1; irr_char[5][5] = -1; irr_char[5][6] =  1; irr_char[5][7] =  1;
      irr_char[6][0] =  1; irr_char[6][1] = -1; irr_char[6][2] =  1; irr_char[6][3] = -1;
      irr_char[6][4] = -1; irr_char[6][5] =  1; irr_char[6][6] = -1; irr_char[6][7] =  1;
      irr_char[7][0] =  1; irr_char[7][1] = -1; irr_char[7][2] = -1; irr_char[7][3] =  1;
      irr_char[7][4] = -1; irr_char[7][5] =  1; irr_char[7][6] =  1; irr_char[7][7] = -1;
      break;
    default:
      punt("Number of irreps in find_symmetry is invalid.");
  }

  return;
}

/* This function permutes axes so that axis1 becomes axis2
   This is done by a rotation about C3 axis relating x, y, and z
   unit vectors */
void permute_axes(double **geom, int axis1, int axis2)
{
  using namespace psi;
  int atom,i,j;
  int axis3;
  double tmp;
  double **R;    /* The matrix representation of the forward rotation */

  if (axis1 == axis2)
      punt("Called permute_axes with axis1=axis2");

  for(i=0;i<3;i++)
    if ( (i != axis1) && (i != axis2)) {
      axis3 = i;
      break;
    }

  R = block_matrix(3,3);
  R[axis1][axis2] = 1.0;
  R[axis2][axis3] = 1.0;
  R[axis3][axis1] = 1.0;
  rotate_full_geom(R);
  free_block(R);

  return;
}

} // namespace

