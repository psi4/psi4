#ifndef _psi3_libchkpt_config_h_
#define _psi3_libchkpt_config_h_

#ifndef MAX_ELEMNAME
#define MAX_ELEMNAME 13
#endif

#ifndef CHKPT_PREFIX_LEN
#define CHKPT_PREFIX_LEN 32
#endif

/*--- Z-matrix entry type ---*/
struct z_entry {
  int bond_atom;            /* first reference atom (bond) */
  int angle_atom;           /* second reference atom (angle) */
  int tors_atom;            /* third reference atom (torsion) */
  int bond_opt;             /* flags indicating to optimize values */
  int angle_opt; 
  int tors_opt; 
  double bond_val;          /* coordinate values */
  double angle_val; 
  double tors_val; 
  char bond_label[20];      /* variable labels, if any */
  char angle_label[20]; 
  char tors_label[20]; 
};

/*--- Types of reference determinants ---*/
typedef enum {ref_rhf = 0, ref_uhf = 1, ref_rohf = 2, ref_tcscf = 3,
	      ref_rks = 4, ref_uks = 5} reftype;

#endif
