/*! \file dcr.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<libciomr/libciomr.h>
#include<symmetry.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"small_fns.h"

namespace psi { namespace cints {

/*!------------------------------------------------------------------------------------------
  Generate sets of Double Coset Representatives for each combination of nuclear stabilizers
  (which are just subgroups of the point group). The subgroup-number correspondence is
  arbitrary, so I just hardwired it hoping that noone will want to go beyond D2h
 ------------------------------------------------------------------------------------------*/

void init_dcr()
{
  int atom;
  int u,v;
  /*--- Local equivalents for equivalent members of Symmetry ---*/
  int ***dcr;
  int **dcr_dim;
  int **dcr_deg;
  int **GnG;

  /*-----------------------------
    Check if need to do anything
   -----------------------------*/
  if (!strcmp(Symmetry.symlabel,"C1 "))
    return;
  /*-------------------------------------------
    C2, Ci, and Cs cases are almost identical:
    C1        =   0
    C2/Ci/Cs  =   1
   -------------------------------------------*/
  else if (!strcmp(Symmetry.symlabel,"C2 ") || !strcmp(Symmetry.symlabel,"Cs ") || !strcmp(Symmetry.symlabel,"Ci ")) {
    /*--- build DCR-related arrays ---*/
    dcr_dim = init_int_matrix(2,2);
    dcr_dim[0][0] = 2;
    dcr_dim[1][0] = dcr_dim[0][1] = dcr_dim[1][1] = 1;

    dcr = (int ***) malloc(sizeof(int **)*2);
    for(u=0;u<2;u++) {
      dcr[u] = (int **) malloc(sizeof(int *)*2);
      for(v=0;v<2;v++)
	dcr[u][v] = (int *) malloc(sizeof(int)*dcr_dim[u][v]);
    }
    dcr[0][0][0] = 0; dcr[0][0][1] = 1;
    dcr[1][0][0] = dcr[0][1][0] = dcr[1][1][0] = 0;

    dcr_deg = init_int_matrix(2,2);
    dcr_deg[0][0] = dcr_deg[1][0] = dcr_deg[0][1] = 1;
    dcr_deg[1][1] = 2;

    GnG = init_int_matrix(2,2);
    GnG[0][0] = GnG[0][1] = GnG[1][0] = 0;
    GnG[1][1] = 1;

    /*--- convert symmetry positions to stabilizers ---*/
    if (!strcmp(Symmetry.symlabel,"Cs "))
      for(atom=0;atom<Molecule.num_atoms;atom++)
	if (Symmetry.atom_positions[atom] == SIGXYCODE)
	  Symmetry.atom_positions[atom] = 1;
	else if (Symmetry.atom_positions[atom] == ECODE)
	  Symmetry.atom_positions[atom] = 0;
	else
	  throw std::domain_error("Unrecognized symmetry position");
    else if (!strcmp(Symmetry.symlabel,"C2 "))
      for(atom=0;atom<Molecule.num_atoms;atom++)
	if (Symmetry.atom_positions[atom] == C2ZCODE)
	  Symmetry.atom_positions[atom] = 1;
	else if (Symmetry.atom_positions[atom] == ECODE)
	  Symmetry.atom_positions[atom] = 0;
	else
	  throw std::domain_error("Unrecognized symmetry position");
    else /*--- Ci case ---*/
      for(atom=0;atom<Molecule.num_atoms;atom++)
	if (Symmetry.atom_positions[atom] == ICODE)
	  Symmetry.atom_positions[atom] = 1;
	else if (Symmetry.atom_positions[atom] == ECODE)
	  Symmetry.atom_positions[atom] = 0;
	else
	  throw std::domain_error("Unrecognized symmetry position");
  }
  /*------------------
    subgroups of C2v:
    C1         = 0
    Cs(sig_xz) = 1
    Cs(sig_yz) = 2
    C2v        = 3
   ------------------*/
  else if (!strcmp(Symmetry.symlabel,"C2v")) {
    /*--- build DCR-related arrays ---*/
    dcr_dim = init_int_matrix(4,4);
    dcr_dim[0][0] = 4;
    dcr_dim[1][0] = dcr_dim[0][1] = dcr_dim[1][1] = dcr_dim[2][0] = dcr_dim[0][2] = dcr_dim[2][2] = 2;
    dcr_dim[3][0] = dcr_dim[0][3] = dcr_dim[1][2] = dcr_dim[2][1] = dcr_dim[1][3] = dcr_dim[3][1] =
		    dcr_dim[2][3] = dcr_dim[3][2] = dcr_dim[3][3] = 1;

    dcr = (int ***) malloc(sizeof(int **)*4);
    for(u=0;u<4;u++) {
      dcr[u] = (int **) malloc(sizeof(int *)*4);
      for(v=0;v<4;v++)
	dcr[u][v] = (int *) malloc(sizeof(int)*dcr_dim[u][v]);
    }
    /* DCR(C1,C1)         = E,C2,sig_xz,sig_yz */
    dcr[0][0][0] = 0; dcr[0][0][1] = 1; dcr[0][0][2] = 2; dcr[0][0][3] = 3;
    /* DCR(C1,Cs(xz))     = E,C2 */
    dcr[0][1][0] = 0; dcr[0][1][1] = 1;
    dcr[1][0][0] = 0; dcr[1][0][1] = 1;
    /* DCR(C1,Cs(yz))     = E,C2 */
    dcr[0][2][0] = 0; dcr[0][2][1] = 1;
    dcr[2][0][0] = 0; dcr[2][0][1] = 1;
    /* DCR(C1,C2v)        = E */
    dcr[0][3][0] = 0;
    dcr[3][0][0] = 0;
    /* DCR(Cs(xz),Cs(xz)) = E,C2 */
    dcr[1][1][0] = 0; dcr[1][1][1] = 1;
    /* DCR(Cs(xz),Cs(yz)) = E */
    dcr[1][2][0] = 0;
    dcr[2][1][0] = 0;
    /* DCR(Cs(xz),C2v)    = E */
    dcr[1][3][0] = 0;
    dcr[3][1][0] = 0;
    /* DCR(Cs(yz),Cs(yz)) = E,C2 */
    dcr[2][2][0] = 0; dcr[2][2][1] = 1; 
    /* DCR(Cs(yz),C2v)    = E */
    dcr[2][3][0] = 0;
    dcr[3][2][0] = 0;
    /* DCR(C2v,C2v)       = E */
    dcr[3][3][0] = 0;

    dcr_deg = init_int_matrix(4,4);
    dcr_deg[0][0] = dcr_deg[1][0] = dcr_deg[0][1] = dcr_deg[2][0] = dcr_deg[0][2] = dcr_deg[3][0] =
		    dcr_deg[0][3] = dcr_deg[1][2] = dcr_deg[2][1] = 1;
    dcr_deg[1][1] = dcr_deg[2][2] = dcr_deg[1][3] = dcr_deg[3][1] = dcr_deg[2][3] = dcr_deg[3][2] = 2;
    dcr_deg[3][3] = 4;

    GnG = init_int_matrix(4,4);
    GnG[0][0] = GnG[0][1] = GnG[1][0] = GnG[0][2] = GnG[2][0] = GnG[0][3] = GnG[3][0] = GnG[1][2] = GnG[2][1] = 0;
    GnG[1][1] = GnG[1][3] = GnG[3][1] = 1;
    GnG[2][2] = GnG[2][3] = GnG[3][2] = 2;
    GnG[3][3] = 3;

    /*--- convert symmetry positions to stabilizers ---*/
    for(atom=0;atom<Molecule.num_atoms;atom++)
      if (Symmetry.atom_positions[atom] == C2ZCODE)
	Symmetry.atom_positions[atom] = 3;
      else if (Symmetry.atom_positions[atom] == SIGYZCODE)
	Symmetry.atom_positions[atom] = 2;
      else if (Symmetry.atom_positions[atom] == SIGXZCODE)
	Symmetry.atom_positions[atom] = 1;
      else if (Symmetry.atom_positions[atom] == ECODE)
	Symmetry.atom_positions[atom] = 0;
      else
	throw std::domain_error("Unrecognized symmetry position");
  }
  /*------------------
    subgroups of C2h:
    C1         = 0
    Cs         = 1
    C2         = 2
    C2h        = 3
   ------------------*/
  else if (!strcmp(Symmetry.symlabel,"C2h")) {
    /*--- build DCR-related arrays ---*/
    dcr_dim = init_int_matrix(4,4);
    dcr_dim[0][0] = 4;
    dcr_dim[1][0] = dcr_dim[0][1] = dcr_dim[1][1] = dcr_dim[2][0] = dcr_dim[0][2] = dcr_dim[2][2] = 2;
    dcr_dim[3][0] = dcr_dim[0][3] = dcr_dim[1][2] = dcr_dim[2][1] = dcr_dim[1][3] = dcr_dim[3][1] =
		    dcr_dim[2][3] = dcr_dim[3][2] = dcr_dim[3][3] = 1;

    dcr = (int ***) malloc(sizeof(int **)*4);
    for(u=0;u<4;u++) {
      dcr[u] = (int **) malloc(sizeof(int *)*4);
      for(v=0;v<4;v++)
	dcr[u][v] = (int *) malloc(sizeof(int)*dcr_dim[u][v]);
    }
    /* DCR(C1,C1)     = E,C2,i,sig */
    dcr[0][0][0] = 0; dcr[0][0][1] = 1; dcr[0][0][2] = 2; dcr[0][0][3] = 3;
    /* DCR(C1,Cs)     = E,C2 */
    dcr[0][1][0] = 0; dcr[0][1][1] = 1;
    dcr[1][0][0] = 0; dcr[1][0][1] = 1;
    /* DCR(C1,C2)     = E,i */
    dcr[0][2][0] = 0; dcr[0][2][1] = 2;
    dcr[2][0][0] = 0; dcr[2][0][1] = 2;
    /* DCR(C1,C2h)    = E */
    dcr[0][3][0] = 0;
    dcr[3][0][0] = 0;
    /* DCR(Cs,Cs)     = E,C2 */
    dcr[1][1][0] = 0; dcr[1][1][1] = 1;
    /* DCR(Cs,C2)     = E */
    dcr[1][2][0] = 0;
    dcr[2][1][0] = 0;
    /* DCR(Cs,C2h)    = E */
    dcr[1][3][0] = 0;
    dcr[3][1][0] = 0;
    /* DCR(C2,C2)     = E,i */
    dcr[2][2][0] = 0; dcr[2][2][1] = 2; 
    /* DCR(C2,C2h)    = E */
    dcr[2][3][0] = 0;
    dcr[3][2][0] = 0;
    /* DCR(C2h,C2h)   = E */
    dcr[3][3][0] = 0;

    dcr_deg = init_int_matrix(4,4);
    dcr_deg[0][0] = dcr_deg[1][0] = dcr_deg[0][1] = dcr_deg[2][0] = dcr_deg[0][2] = dcr_deg[3][0] =
		    dcr_deg[0][3] = dcr_deg[1][2] = dcr_deg[2][1] = 1;
    dcr_deg[1][1] = dcr_deg[2][2] = dcr_deg[1][3] = dcr_deg[3][1] = dcr_deg[2][3] = dcr_deg[3][2] = 2;
    dcr_deg[3][3] = 4;

    GnG = init_int_matrix(4,4);
    GnG[0][0] = GnG[0][1] = GnG[1][0] = GnG[0][2] = GnG[2][0] = GnG[0][3] = GnG[3][0] = GnG[1][2] = GnG[2][1] = 0;
    GnG[1][1] = GnG[1][3] = GnG[3][1] = 1;
    GnG[2][2] = GnG[2][3] = GnG[3][2] = 2;
    GnG[3][3] = 3;

    /*--- convert symmetry positions to stabilizers ---*/
    for(atom=0;atom<Molecule.num_atoms;atom++)
      if (Symmetry.atom_positions[atom] == ICODE)
	Symmetry.atom_positions[atom] = 3;
      else if (Symmetry.atom_positions[atom] == C2ZCODE)
	Symmetry.atom_positions[atom] = 2;
      else if (Symmetry.atom_positions[atom] == SIGXYCODE)
	Symmetry.atom_positions[atom] = 1;
      else if (Symmetry.atom_positions[atom] == ECODE)
	Symmetry.atom_positions[atom] = 0;
      else
	throw std::domain_error("Unrecognized symmetry position");
  }
  /*------------------
    subgroups of D2:
    C1         = 0
    C2z        = 1
    C2y        = 2
    C2x        = 3
    D2         = 4
   ------------------*/
  else if (!strcmp(Symmetry.symlabel,"D2 ")) {
    /*--- build DCR-related arrays ---*/
    dcr_dim = init_int_matrix(5,5);
    dcr_dim[0][0] = 4;
    dcr_dim[1][0] = dcr_dim[0][1] = dcr_dim[2][0] = dcr_dim[0][2] = dcr_dim[3][0] = dcr_dim[0][3] =
		    dcr_dim[1][1] = dcr_dim[2][2] = dcr_dim[3][3] = 2;
    dcr_dim[4][0] = dcr_dim[0][4] = dcr_dim[1][2] = dcr_dim[2][1] = dcr_dim[1][3] = dcr_dim[3][1] =
		    dcr_dim[1][4] = dcr_dim[4][1] = dcr_dim[2][3] = dcr_dim[3][2] = dcr_dim[2][4] =
		    dcr_dim[4][2] = dcr_dim[3][4] = dcr_dim[4][3] = dcr_dim[4][4] = 1;

    dcr = (int ***) malloc(sizeof(int **)*5);
    for(u=0;u<5;u++) {
      dcr[u] = (int **) malloc(sizeof(int *)*5);
      for(v=0;v<5;v++)
	dcr[u][v] = (int *) malloc(sizeof(int)*dcr_dim[u][v]);
    }
    /* DCR(C1,C1)     = E,C2z,C2y,C2x */
    dcr[0][0][0] = 0; dcr[0][0][1] = 1; dcr[0][0][2] = 2; dcr[0][0][3] = 3;
    /* DCR(C1,C2z)    = E,C2y */
    dcr[0][1][0] = 0; dcr[0][1][1] = 2;
    dcr[1][0][0] = 0; dcr[1][0][1] = 2;
    /* DCR(C1,C2y)    = E,C2z */
    dcr[0][2][0] = 0; dcr[0][2][1] = 1;
    dcr[2][0][0] = 0; dcr[2][0][1] = 1;
    /* DCR(C1,C2x)    = E,C2z */
    dcr[0][3][0] = 0; dcr[0][3][1] = 1;
    dcr[3][0][0] = 0; dcr[3][0][1] = 1;
    /* DCR(C1,D2)     = E */
    dcr[0][4][0] = 0;
    dcr[4][0][0] = 0;
    /* DCR(C2z,C2z)   = E,C2y */
    dcr[1][1][0] = 0; dcr[1][1][1] = 2;
    /* DCR(C2z,C2y)   = E */
    dcr[1][2][0] = 0;
    dcr[2][1][0] = 0;
    /* DCR(C2z,C2x)   = E */
    dcr[1][3][0] = 0;
    dcr[3][1][0] = 0;
    /* DCR(C2z,D2)    = E */
    dcr[1][4][0] = 0;
    dcr[4][1][0] = 0;
    /* DCR(C2y,C2y)   = E,C2z */
    dcr[2][2][0] = 0; dcr[2][2][1] = 1;
    /* DCR(C2y,C2x)   = E */
    dcr[2][3][0] = 0;
    dcr[3][2][0] = 0;
    /* DCR(C2y,D2)    = E */
    dcr[2][4][0] = 0;
    dcr[4][2][0] = 0;
    /* DCR(C2x,C2x)   = E,C2z */
    dcr[3][3][0] = 0; dcr[3][3][1] = 1;
    /* DCR(C2x,D2)    = E */
    dcr[3][4][0] = 0;
    dcr[4][3][0] = 0;
    /* DCR(D2,D2)     = E */
    dcr[4][4][0] = 0;

    dcr_deg = init_int_matrix(5,5);
    dcr_deg[0][0] = dcr_deg[0][1] = dcr_deg[1][0] = dcr_deg[0][2] = dcr_deg[2][0] = dcr_deg[0][3] = dcr_deg[3][0] = 
		    dcr_deg[0][4] = dcr_deg[4][0] = dcr_deg[1][2] = dcr_deg[2][1] = dcr_deg[1][3] = dcr_deg[3][1] =
		    dcr_deg[2][3] = dcr_deg[3][2] = 1;
    dcr_deg[1][1] = dcr_deg[1][4] = dcr_deg[4][1] = dcr_deg[2][2] = dcr_deg[2][4] = dcr_deg[4][2] = dcr_deg[3][3] =
		    dcr_deg[3][4] = dcr_deg[4][3] = 2;
    dcr_deg[4][4] = 4;

    GnG = init_int_matrix(5,5);
    GnG[0][0] = GnG[0][1] = GnG[1][0] = GnG[0][2] = GnG[2][0] = GnG[0][3] = GnG[3][0] = 
		GnG[0][4] = GnG[4][0] = GnG[1][2] = GnG[2][1] = GnG[1][3] = GnG[3][1] =
		GnG[2][3] = GnG[3][2] = 0;
    GnG[1][1] = GnG[1][4] = GnG[4][1] = 1;
    GnG[2][2] = GnG[2][4] = GnG[4][2] = 2;
    GnG[3][3] = GnG[3][4] = GnG[4][3] = 3;
    GnG[4][4] = 4;

    /*--- convert symmetry positions to stabilizers ---*/
    for(atom=0;atom<Molecule.num_atoms;atom++)
      if (Symmetry.atom_positions[atom] == ICODE)
	Symmetry.atom_positions[atom] = 4;
      else if (Symmetry.atom_positions[atom] == C2XCODE)
	Symmetry.atom_positions[atom] = 3;
      else if (Symmetry.atom_positions[atom] == C2YCODE)
	Symmetry.atom_positions[atom] = 2;
      else if (Symmetry.atom_positions[atom] == C2ZCODE)
	Symmetry.atom_positions[atom] = 1;
      else if (Symmetry.atom_positions[atom] == ECODE)
	Symmetry.atom_positions[atom] = 0;
      else
	throw std::domain_error("Unrecognized symmetry position");
  }
  /*------------------
    subgroups of D2h:
    C1         = 0
    Cs(xy)     = 1
    Cs(xz)     = 2
    Cs(yz)     = 3
    C2z        = 4
    C2y        = 5
    C2x        = 6
    D2h        = 7
   ------------------*/
  else if (!strcmp(Symmetry.symlabel,"D2h")) {
    /*--- build DCR-related arrays ---*/
    dcr_dim = init_int_matrix(8,8);
    dcr_dim[0][0] = 8;
    dcr_dim[1][0] = dcr_dim[0][1] = dcr_dim[2][0] = dcr_dim[0][2] = dcr_dim[3][0] = dcr_dim[0][3] =
		    dcr_dim[1][1] = dcr_dim[2][2] = dcr_dim[3][3] = 4;
    dcr_dim[4][0] = dcr_dim[0][4] = dcr_dim[5][0] = dcr_dim[0][5] = dcr_dim[6][0] = dcr_dim[0][6] = 
		    dcr_dim[1][2] = dcr_dim[2][1] = dcr_dim[1][3] = dcr_dim[3][1] =
		    dcr_dim[1][5] = dcr_dim[5][1] = dcr_dim[1][6] = dcr_dim[6][1] =
		    dcr_dim[2][3] = dcr_dim[3][2] = dcr_dim[2][4] = dcr_dim[4][2] =
		    dcr_dim[2][6] = dcr_dim[6][2] = dcr_dim[3][4] = dcr_dim[4][3] =
		    dcr_dim[3][5] = dcr_dim[5][3] = dcr_dim[4][4] = dcr_dim[5][5] = dcr_dim[6][6] = 2;
    dcr_dim[0][7] = dcr_dim[7][0] = dcr_dim[1][4] = dcr_dim[4][1] = dcr_dim[2][5] = dcr_dim[5][2] =
		    dcr_dim[3][6] = dcr_dim[6][3] = dcr_dim[1][7] = dcr_dim[7][1] = dcr_dim[2][7] =
		    dcr_dim[7][2] = dcr_dim[3][7] = dcr_dim[7][3] = dcr_dim[4][5] = dcr_dim[5][4] =
		    dcr_dim[4][6] = dcr_dim[6][4] = dcr_dim[5][6] = dcr_dim[6][5] =
		    dcr_dim[4][7] = dcr_dim[7][4] = dcr_dim[5][7] = dcr_dim[7][5] =
		    dcr_dim[6][7] = dcr_dim[7][6] = dcr_dim[7][7] = 1;

    dcr = (int ***) malloc(sizeof(int **)*8);
    for(u=0;u<8;u++) {
      dcr[u] = (int **) malloc(sizeof(int *)*8);
      for(v=0;v<8;v++)
	dcr[u][v] = (int *) malloc(sizeof(int)*dcr_dim[u][v]);
    }
    /* DCR(C1,C1)          = E,C2z,C2y,C2x,i,sig_xy,sig_xz,sig_yz */
    dcr[0][0][0] = 0; dcr[0][0][1] = 1; dcr[0][0][2] = 2; dcr[0][0][3] = 3;
    dcr[0][0][4] = 4; dcr[0][0][5] = 5; dcr[0][0][6] = 6; dcr[0][0][7] = 7;
    /* DCR(C1,Cs(xy))      = E,C2z,C2y,C2x */
    dcr[0][1][0] = 0; dcr[0][1][1] = 1; dcr[0][1][2] = 2; dcr[0][1][3] = 3; 
    dcr[1][0][0] = 0; dcr[1][0][1] = 1; dcr[1][0][2] = 2; dcr[1][0][3] = 3;
    /* DCR(C1,Cs(xz))      = E,C2z,C2y,C2x */
    dcr[0][2][0] = 0; dcr[0][2][1] = 1; dcr[0][2][2] = 2; dcr[0][2][3] = 3; 
    dcr[2][0][0] = 0; dcr[2][0][1] = 1; dcr[2][0][2] = 2; dcr[2][0][3] = 3;
    /* DCR(C1,Cs(yz))      = E,C2z,C2y,C2x */
    dcr[0][3][0] = 0; dcr[0][3][1] = 1; dcr[0][3][2] = 2; dcr[0][3][3] = 3; 
    dcr[3][0][0] = 0; dcr[3][0][1] = 1; dcr[3][0][2] = 2; dcr[3][0][3] = 3;
    /* DCR(C1,C2v(z))      = E,i */
    dcr[0][4][0] = 0; dcr[0][4][1] = 4;
    dcr[4][0][0] = 0; dcr[4][0][1] = 4;
    /* DCR(C1,C2v(y))      = E,i */
    dcr[0][5][0] = 0; dcr[0][5][1] = 4;
    dcr[5][0][0] = 0; dcr[5][0][1] = 4;
    /* DCR(C1,C2v(x))      = E,i */
    dcr[0][6][0] = 0; dcr[0][6][1] = 4;
    dcr[6][0][0] = 0; dcr[6][0][1] = 4;
    /* DCR(C1,D2h)         = E */
    dcr[0][7][0] = 0;
    dcr[7][0][0] = 0;
    /* DCR(Cs(xy),Cs(xy))  = E,C2z,C2y,C2x */
    dcr[1][1][0] = 0; dcr[1][1][1] = 1; dcr[1][1][2] = 2; dcr[1][1][3] = 3;
    /* DCR(Cs(xy),Cs(xz))  = E,i */
    dcr[1][2][0] = 0; dcr[1][2][1] = 4;
    dcr[2][1][0] = 0; dcr[2][1][1] = 4;
    /* DCR(Cs(xy),Cs(yz))  = E,i */
    dcr[1][3][0] = 0; dcr[1][3][1] = 4;
    dcr[3][1][0] = 0; dcr[3][1][1] = 4;
    /* DCR(Cs(xy),C2v(z))  = E */
    dcr[1][4][0] = 0;
    dcr[4][1][0] = 0;
    /* DCR(Cs(xy),C2v(y))  = E,i */
    dcr[1][5][0] = 0; dcr[1][5][1] = 4;
    dcr[5][1][0] = 0; dcr[5][1][1] = 4;
    /* DCR(Cs(xy),C2v(x))  = E,i */
    dcr[1][6][0] = 0; dcr[1][6][1] = 4;
    dcr[6][1][0] = 0; dcr[6][1][1] = 4;
    /* DCR(Cs(xy),D2h)     = E */
    dcr[1][7][0] = 0;
    dcr[7][1][0] = 0;
    /* DCR(Cs(xz),Cs(xz))  = E,C2z,C2y,C2x */
    dcr[2][2][0] = 0; dcr[2][2][1] = 1; dcr[2][2][2] = 2; dcr[2][2][3] = 3;
    /* DCR(Cs(xz),Cs(yz))  = E,i */
    dcr[2][3][0] = 0; dcr[2][3][1] = 4;
    dcr[3][2][0] = 0; dcr[3][2][1] = 4;
    /* DCR(Cs(xz),C2v(z))  = E,i */
    dcr[2][4][0] = 0; dcr[2][4][1] = 4;
    dcr[4][2][0] = 0; dcr[4][2][1] = 4;
    /* DCR(Cs(xz),C2v(y))  = E */
    dcr[2][5][0] = 0;
    dcr[5][2][0] = 0;
    /* DCR(Cs(xz),C2v(x))  = E,i */
    dcr[2][6][0] = 0; dcr[2][6][1] = 4;
    dcr[6][2][0] = 0; dcr[6][2][1] = 4;
    /* DCR(Cs(xz),D2h)     = E */
    dcr[2][7][0] = 0;
    dcr[7][2][0] = 0;
    /* DCR(Cs(yz),Cs(yz))  = E,C2z,C2y,C2x */
    dcr[3][3][0] = 0; dcr[3][3][1] = 1; dcr[3][3][2] = 2; dcr[3][3][3] = 3;
    /* DCR(Cs(yz),C2v(z))  = E,i */
    dcr[3][4][0] = 0; dcr[3][4][1] = 4;
    dcr[4][3][0] = 0; dcr[4][3][1] = 4;
    /* DCR(Cs(yz),C2v(y))  = E,i */
    dcr[3][5][0] = 0; dcr[3][5][1] = 4;
    dcr[5][3][0] = 0; dcr[5][3][1] = 4;
    /* DCR(Cs(yz),C2v(x))  = E */
    dcr[3][6][0] = 0;
    dcr[6][3][0] = 0;
    /* DCR(Cs(yz),D2h)     = E */
    dcr[3][7][0] = 0;
    dcr[7][3][0] = 0;
    /* DCR(C2v(z),C2v(z))  = E,i */
    dcr[4][4][0] = 0; dcr[4][4][1] = 4;
    /* DCR(C2v(z),C2v(y))  = E */
    dcr[4][5][0] = 0;
    dcr[5][4][0] = 0;
    /* DCR(C2v(z),C2v(x))  = E */
    dcr[4][6][0] = 0;
    dcr[6][4][0] = 0;
    /* DCR(C2v(z),D2h)     = E */
    dcr[4][7][0] = 0;
    dcr[7][4][0] = 0;
    /* DCR(C2v(y),C2v(y))  = E,i */
    dcr[5][5][0] = 0; dcr[5][5][1] = 4;
    /* DCR(C2v(y),C2v(x))  = E */
    dcr[5][6][0] = 0;
    dcr[6][5][0] = 0;
    /* DCR(C2v(y),D2h)     = E */
    dcr[5][7][0] = 0;
    dcr[7][5][0] = 0;
    /* DCR(C2v(x),C2v(x))  = E,i */
    dcr[6][6][0] = 0; dcr[6][6][1] = 4;
    /* DCR(C2v(x),D2h)     = E */
    dcr[6][7][0] = 0;
    dcr[7][6][0] = 0;
    /* DCR(D2h,D2h)        = E */
    dcr[7][7][0] = 0;
    
    dcr_deg = init_int_matrix(8,8);
    dcr_deg[0][0] = dcr_deg[0][1] = dcr_deg[1][0] = dcr_deg[0][2] = dcr_deg[2][0] = dcr_deg[0][3] = dcr_deg[3][0] = 
		    dcr_deg[0][4] = dcr_deg[4][0] = dcr_deg[0][5] = dcr_deg[5][0] = dcr_deg[0][6] = dcr_deg[6][0] =
		    dcr_deg[0][7] = dcr_deg[7][0] = dcr_deg[1][2] = dcr_deg[2][1] = dcr_deg[1][3] = dcr_deg[3][1] =
		    dcr_deg[2][3] = dcr_deg[3][2] = dcr_deg[1][4] = dcr_deg[4][1] = dcr_deg[2][5] = dcr_deg[5][2] =
		    dcr_deg[3][6] = dcr_deg[6][3] = 1;
    dcr_deg[1][1] = dcr_deg[1][5] = dcr_deg[5][1] = dcr_deg[1][6] = dcr_deg[6][1] = dcr_deg[1][7] = dcr_deg[7][1] =
		    dcr_deg[5][6] = dcr_deg[6][5] = dcr_deg[2][2] = dcr_deg[2][4] = dcr_deg[4][2] =
		    dcr_deg[2][6] = dcr_deg[6][2] = dcr_deg[2][7] = dcr_deg[7][2] = dcr_deg[4][6] = dcr_deg[6][4] =
		    dcr_deg[3][3] = dcr_deg[3][4] = dcr_deg[4][3] = dcr_deg[3][5] = dcr_deg[5][3] =
		    dcr_deg[3][7] = dcr_deg[7][3] = dcr_deg[4][5] = dcr_deg[5][4] = 2;
    dcr_deg[4][4] = dcr_deg[4][7] = dcr_deg[7][4] = dcr_deg[5][5] = dcr_deg[5][7] = dcr_deg[7][5] = dcr_deg[6][6] =
		    dcr_deg[6][7] = dcr_deg[7][6] = 4;
    dcr_deg[7][7] = 8;

    GnG = init_int_matrix(8,8);
    GnG[0][0] = GnG[0][1] = GnG[1][0] = GnG[0][2] = GnG[2][0] = GnG[0][3] = GnG[3][0] = 
		GnG[0][4] = GnG[4][0] = GnG[0][5] = GnG[5][0] = GnG[0][6] = GnG[6][0] =
		GnG[0][7] = GnG[7][0] = GnG[1][2] = GnG[2][1] = GnG[1][3] = GnG[3][1] =
		GnG[2][3] = GnG[3][2] = GnG[1][4] = GnG[4][1] = GnG[2][5] = GnG[5][2] =
		GnG[3][6] = GnG[6][3] = 0;
    GnG[1][1] = GnG[1][5] = GnG[5][1] = GnG[1][6] = GnG[6][1] = GnG[1][7] = GnG[7][1] =
		GnG[5][6] = GnG[6][5] = 1;
    GnG[2][2] = GnG[2][4] = GnG[4][2] = GnG[2][6] = GnG[6][2] = GnG[2][7] = GnG[7][2] =
		GnG[4][6] = GnG[6][4] = 2;
    GnG[3][3] = GnG[3][4] = GnG[4][3] = GnG[3][5] = GnG[5][3] = GnG[3][7] = GnG[7][3] =
		GnG[4][5] = GnG[5][4] = 3;
    GnG[4][4] = GnG[4][7] = GnG[7][4] = 4;
    GnG[5][5] = GnG[5][7] = GnG[7][5] = 5;
    GnG[6][6] = GnG[6][7] = GnG[7][6] = 6;
    GnG[7][7] = 7;
    

    /*--- convert symmetry positions to stabilizers ---*/
    for(atom=0;atom<Molecule.num_atoms;atom++)
      if (Symmetry.atom_positions[atom] == ICODE)
	Symmetry.atom_positions[atom] = 7;
      else if (Symmetry.atom_positions[atom] == C2XCODE)
	Symmetry.atom_positions[atom] = 6;
      else if (Symmetry.atom_positions[atom] == C2YCODE)
	Symmetry.atom_positions[atom] = 5;
      else if (Symmetry.atom_positions[atom] == C2ZCODE)
	Symmetry.atom_positions[atom] = 4;
      else if (Symmetry.atom_positions[atom] == SIGYZCODE)
	Symmetry.atom_positions[atom] = 3;
      else if (Symmetry.atom_positions[atom] == SIGXZCODE)
	Symmetry.atom_positions[atom] = 2;
      else if (Symmetry.atom_positions[atom] == SIGXYCODE)
	Symmetry.atom_positions[atom] = 1;
      else if (Symmetry.atom_positions[atom] == ECODE)
	Symmetry.atom_positions[atom] = 0;
      else
	throw std::domain_error("Unrecognized symmetry position");
  }

  Symmetry.dcr = dcr;
  Symmetry.dcr_dim = dcr_dim;
  Symmetry.dcr_deg = dcr_deg;
  Symmetry.GnG = GnG;

  return;
}


void cleanup_dcr(void)
{
  /*--- It's too messy to clean it up ---*/
  return;
}
}}
