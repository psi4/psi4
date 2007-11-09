#ifndef _psi_src_lib_libciomr_structs_h_
#define _psi_src_lib_libciomr_structs_h_

#ifdef __cplusplus
extern "C" {
#endif

struct integrals {
     double *smat;
     double *hmat;
     double *twomat;
     double *pmat;
     double **sahalf;
     double **fock_mat;
     double **cmat;
     double **fock_trans;
     double *fock_trns_pac;
     double *fock_eig_vals;
     double **ctrans;
     double *gmat;
     double *fock_pac;
     double *pmat_old;
     } ;

#ifdef __cplusplus
}
#endif

#endif /* header guard */
