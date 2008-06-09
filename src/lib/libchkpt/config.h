#ifndef _psi3_libchkpt_config_h_
#define _psi3_libchkpt_config_h_

namespace psi {

#ifndef MAX_ELEMNAME
#define MAX_ELEMNAME 13
#endif

#ifndef CHKPT_PREFIX_LEN
#define CHKPT_PREFIX_LEN 32
#endif

#define CHKPT_ZMAT_LABEL_LEN 20
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
  char bond_label[CHKPT_ZMAT_LABEL_LEN];      /* variable labels, if any */
  char angle_label[CHKPT_ZMAT_LABEL_LEN]; 
  char tors_label[CHKPT_ZMAT_LABEL_LEN]; 
};

/*--- Types of reference determinants ---*/
typedef enum {ref_rhf = 0, ref_uhf = 1, ref_rohf = 2, ref_tcscf = 3,
	      ref_rks = 4, ref_uks = 5} reftype;

}

#endif
