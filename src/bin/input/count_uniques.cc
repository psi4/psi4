/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <cmath>
#include <symmetry.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

/*-----------------------------------------------------------------------------------------------------------------
  This function determines number of unique atoms (num_uniques), fills in the reference array (u2a),
  computes degeneracy of each unique atom (unique_degen, number of symmetry equivalent atoms generated from the unique + 1),
  atom's classes and positions.
 -----------------------------------------------------------------------------------------------------------------*/

void count_uniques()
{
  int i,j,k,l,m,un,atom;
  int last_class;
  int last_unique_class;
  int sym_oper_code;
  int incr;
  int last_gener_atom;
  int class_cnt;
    /*Maps global ordering into local ordering so that the class numbering is neat */
  int map_global2local[] = {IFLAG, C2ZFLAG, C2YFLAG, C2XFLAG, SIGXYFLAG, SIGXZFLAG, SIGYZFLAG, EFLAG};

   /*Initialize global variables*/
  num_uniques = 0;
  ap_flag = 0;
  u2a = init_int_array(num_atoms);
  unique_degen = init_int_array(num_atoms);
  atom_position = init_int_array(num_atoms);
  atom_class = init_int_array(num_atoms);
  for(i=0;i<num_atoms;i++)
    atom_class[i] = -1;

  if (!strcmp(symmetry,"C1  ")) {
    num_uniques = num_atoms;
    ap_flag = ECODE;
    for(i=0;i<num_uniques;i++) { /* all positions are general in C1 */
      u2a[i] = i;
      unique_degen[i] = 1;
      atom_class[i] = 1;
      atom_position[i] = ECODE;
    }
  }
  else if (!strcmp(symmetry,"Ci  ")) {
    for(i=0;i<num_atoms;i++)
      if ( (geometry[i][2] > ZERO) || 
	   (fabs(geometry[i][2]) < ZERO && geometry[i][0] > ZERO) ||
	   (fabs(geometry[i][2]) < ZERO && fabs(geometry[i][0]) < ZERO && geometry[i][1] > ZERO) ) { /* General position */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | ECODE;
	atom_position[i] = ECODE;
      }
      else if (inv_related(geometry[i],geometry[i])){ /* Atom in the inversion center */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 1;
	ap_flag = ap_flag | ICODE;
	atom_position[i] = ICODE;
      }
  }
  else if (!strcmp(symmetry,"Cs  ")) {
    for(i=0;i<num_atoms;i++)
      if (geometry[i][2] > ZERO) { /* General position */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | ECODE;
	atom_position[i] = ECODE;
      }
      else if (fabs(geometry[i][2]) < ZERO) { /* Atom in the plane */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 1;
	ap_flag = ap_flag | SIGXYCODE;
	atom_position[i] = SIGXYCODE;
      }
  }
  else if (!strcmp(symmetry,"C2  ")) {
    for(i=0;i<num_atoms;i++)
      if ( (geometry[i][0] > ZERO) ||
	   (fabs(geometry[i][0]) < ZERO && geometry[i][1] > ZERO) ) { /* General position */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | ECODE;
	atom_position[i] = ECODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && fabs(geometry[i][1]) < ZERO) { /* Atom on the axis */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 1;
	ap_flag = ap_flag | C2ZCODE;
	atom_position[i] = C2ZCODE;
      }
  }
  else if (!strcmp(symmetry,"C2h ")) {
    for(i=0;i<num_atoms;i++)
      if ( (geometry[i][0] > ZERO && geometry[i][2] > ZERO) || /* General position */
	   (fabs(geometry[i][0]) < ZERO && geometry[i][1] > ZERO && geometry[i][2] > ZERO) ) {
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 4;
	ap_flag = ap_flag | ECODE;
	atom_position[i] = ECODE;
      }
      else if ( (geometry[i][0] > ZERO && fabs(geometry[i][2]) < ZERO) ||
		(fabs(geometry[i][0]) < ZERO && geometry[i][1] > ZERO && fabs(geometry[i][2]) < ZERO) ) { /* Atom in the plane */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | SIGXYCODE;
	atom_position[i] = SIGXYCODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && fabs(geometry[i][1]) < ZERO && geometry[i][2] > ZERO) { /* Atom on the axis */
        u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | C2ZCODE;
	atom_position[i] = C2ZCODE;
      }
      else if (inv_related(geometry[i],geometry[i])) { /* Atom in the inversion center */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 1;
	ap_flag = ap_flag | ICODE;
	atom_position[i] = ICODE;
      }
  }
  else if (!strcmp(symmetry,"C2v ")) {
    for(i=0;i<num_atoms;i++)
      if (geometry[i][0] > ZERO && geometry[i][1] > ZERO) { /* General position */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 4;
	ap_flag = ap_flag | ECODE;
	atom_position[i] = ECODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && geometry[i][1] > ZERO) { /* Atom in xz-plane */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | SIGXZCODE;
	atom_position[i] = SIGXZCODE;
      }
      else if (geometry[i][0] > ZERO && fabs(geometry[i][1]) < ZERO) { /* Atom in yz-plane */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | SIGYZCODE;
	atom_position[i] = SIGYZCODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && fabs(geometry[i][1]) < ZERO) { /* Atom on the axis */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 1;
	ap_flag = ap_flag | C2ZCODE;
	atom_position[i] = C2ZCODE;
      }
  }
  else if (!strcmp(symmetry,"D2  ")) {
    for(i=0;i<num_atoms;i++)
      if ( (geometry[i][0] > ZERO && geometry[i][2] > ZERO) || /* General position */
           (fabs(geometry[i][0]) < ZERO && geometry[i][1] > ZERO && geometry[i][2] > ZERO) ||
           (geometry[i][0] > ZERO && geometry[i][1] > ZERO && fabs(geometry[i][2]) < ZERO) ) {
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 4;
	ap_flag = ap_flag | ECODE;
	atom_position[i] = ECODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && fabs(geometry[i][1]) < ZERO && geometry[i][2] > ZERO) { /* Atom on the Z-axis */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | C2ZCODE;
	atom_position[i] = C2ZCODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && fabs(geometry[i][2]) < ZERO && geometry[i][1] > ZERO) { /* Atom on the Y-axis */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | C2YCODE;
	atom_position[i] = C2YCODE;
      }
      else if (fabs(geometry[i][1]) < ZERO && fabs(geometry[i][2]) < ZERO && geometry[i][0] > ZERO) { /* Atom on the X-axis */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | C2XCODE;
	atom_position[i] = C2XCODE;
      }
      else if (inv_related(geometry[i],geometry[i])) { /* Atom in the origin */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 1;
	ap_flag = ap_flag | ICODE;
	atom_position[i] = ICODE;
      }
  }
  else if (!strcmp(symmetry,"D2h ")) {
    for(i=0;i<num_atoms;i++)
      if (geometry[i][0] > ZERO && geometry[i][1] > ZERO && geometry[i][2] > ZERO) { /* General position */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 8;
	ap_flag = ap_flag | ECODE;
	atom_position[i] = ECODE;
      }
      else if (geometry[i][0] > ZERO && geometry[i][1] > ZERO && fabs(geometry[i][2]) < ZERO) { /* In xy-plane */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 4;
	ap_flag = ap_flag | SIGXYCODE;
	atom_position[i] = SIGXYCODE;
      }
      else if (geometry[i][0] > ZERO && fabs(geometry[i][1]) < ZERO && geometry[i][2] > ZERO) { /* In xz-plane */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 4;
	ap_flag = ap_flag | SIGXZCODE;
	atom_position[i] = SIGXZCODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && geometry[i][1] > ZERO && geometry[i][2] > ZERO) { /* In yz-plane */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 4;
	ap_flag = ap_flag | SIGYZCODE;
	atom_position[i] = SIGYZCODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && fabs(geometry[i][1]) < ZERO && geometry[i][2] > ZERO) { /* On Z-axis */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | C2ZCODE;
	atom_position[i] = C2ZCODE;
      }
      else if (geometry[i][0] > ZERO && fabs(geometry[i][1]) < ZERO && fabs(geometry[i][2]) < ZERO) { /* On X-axis */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | C2XCODE;
	atom_position[i] = C2XCODE;
      }
      else if (fabs(geometry[i][0]) < ZERO && geometry[i][1] > ZERO && fabs(geometry[i][2]) < ZERO) { /* On Y-axis */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 2;
	ap_flag = ap_flag | C2YCODE;
	atom_position[i] = C2YCODE;
      }
      else if (inv_related(geometry[i],geometry[i])) { /* Atom in the inversion center */
	u2a[num_uniques] = i;
	unique_degen[num_uniques++] = 1;
	ap_flag = ap_flag | ICODE;
	atom_position[i] = ICODE;
      }
  }
  else
    punt("Unrecognized symmetry in count_uniques().");


  /*Computing atom classes atom_class and reduced unique atom's orbits red_unique_orbit:
    each row contains the unique atom itself plus all non_unique_symmetry_equiv_atoms
    generated from a given unique atom*/

  num_unique_classes = 0;
  num_unique_classes += (ap_flag & ECODE) ? 1 : 0;
  num_unique_classes += (ap_flag & ICODE) ? 1 : 0;
  num_unique_classes += (ap_flag & C2XCODE) ? 1 : 0;
  num_unique_classes += (ap_flag & C2YCODE) ? 1 : 0;
  num_unique_classes += (ap_flag & C2ZCODE) ? 1 : 0;
  num_unique_classes += (ap_flag & SIGXYCODE) ? 1 : 0;
  num_unique_classes += (ap_flag & SIGXZCODE) ? 1 : 0;
  num_unique_classes += (ap_flag & SIGYZCODE) ? 1 : 0;
  uc2c = init_int_array(num_unique_classes);
  unique_class_degen = init_int_array(num_unique_classes);
  red_unique_orbit = init_int_matrix(num_uniques,nirreps);
  for(i=0;i<num_uniques;i++)
    red_unique_orbit[i][0] = u2a[i];
  last_class = 0;
  last_unique_class = 0;
  for(i=0;i<8;i++) {
    sym_oper_code = 1 << map_global2local[i];
    if (ap_flag & sym_oper_code) {
      uc2c[last_unique_class] = last_class;
      for(j=0;j<num_uniques;j++) {
	un = u2a[j];
	if (atom_position[un] == sym_oper_code) {
	  atom_class[un] = last_class;
	  incr = unique_degen[j];
	  unique_class_degen[last_unique_class] = unique_degen[j];
	  last_gener_atom = unique_degen[j] - 1;
	  class_cnt = 1;
	  for(k=0;k<last_gener_atom;k++) {
	    m = atom_orbit[un][k+1];
	    if (atom_class[m] == -1) {
	      red_unique_orbit[j][class_cnt] = m;
	      atom_class[m] = last_class + class_cnt;
	      atom_position[m] = sym_oper_code;
	      class_cnt++;
	    }
	    else
	      last_gener_atom++;
	  }
	}
      }
      last_class += incr;
      last_unique_class++;
    }
  }
  /* Compute total number of classes and a symmetry orbit for each class*/
  num_classes = last_class;
  class_orbit = init_int_matrix(num_classes,nirreps);
  for(i=0;i<num_uniques;i++)
    for(k=0;k<unique_degen[i];k++) {
      atom = red_unique_orbit[i][k];
      m = atom_class[atom];
      for(l=0;l<nirreps;l++)
	class_orbit[m][l] = atom_class[atom_orbit[atom][l]];
    }

}

}} // namespace psi::input
