/*! \file
    \ingroup OPTKING
    \brief opt.h : function declarations
*/

#ifndef _psi3_bin_optking_opt_h_
#define _psi3_bin_optking_opt_h_

#include <psi4-dec.h>

namespace psi { namespace optking {

void open_PSIF(void);
void close_PSIF(void);
void print_mat2(double **matrix, int rows, int cols, FILE *of);
void print_mat5(double **matrix, int rows, int cols, FILE *of);
int div_int(int big, int little);
double **symm_matrix_invert(double **A, int dim, int print_det, int redundant);
double energy_chkpt(void);
void dgeev_optking(int L, double **G, double *lambda, double **alpha);
void print_evects(double **evects, double *evals, int nrow, int ncol, FILE *out);
double **mass_mat(double *masses);
double **unit_matrix(int dim);
void symmetrize_geom(double *x);

double nuclear_repulsion(double *fatomic_num, double *geom);
double get_disp_nuclear_repulsion(const cartesians & carts, int x1, int x2, double disp1,
    double disp2, double disp_size);

// character-table related functions
int **get_char_table(char *ptgrp);    /* returns the character table */
int get_nirreps(char *ptgrp);         /* "     " number of irreps */
const char **get_irrep_labels(char *ptgrp); /* "     " number of symmetry operations */
const char **get_symm_ops(char *ptgrp);/* "     " symm operation labels */
int *get_ops_coeffs(char *ptgrp);      /* "     " coefficients of the symmetry operations */
int get_num_ops(char *ptgrp);          /* "     " number of operations */
int get_num_classes(char *ptgrp);      /* "     " number of classes of operations */
int *get_ops_in_class(char *ptgrp, int nirreps);
int get_irrep_xyz(double **cartrep, int xyz);

int check_coordinates(int natom, double *coord, double *masses, double *Zvals,
    int *ndisp_small, double ***disp_small);
double **irrep(const simples_class &simples, double **di_coord);

double **compute_B(const simples_class &simples, const salc_set &salcs);
double *compute_q(const simples_class &simples, const salc_set &symm);
double **compute_G(double **B, int num_intcos, const cartesians &carts);
void get_optinfo(void);
void get_syminfo(const simples_class &simples);

int disp_user(const cartesians &carts, simples_class & simples, const salc_set &all_salcs);

int make_disp_irrep(const cartesians &carts, simples_class &simples, const salc_set &all_salcs);
int make_disp_nosymm(const cartesians &carts, simples_class &simples, const salc_set &all_salcs);

int disp_fc_grad_selected(const cartesians &carts, simples_class &simples, const salc_set &symm);
void fc_grad_selected(const cartesians &carts, simples_class &simples, const salc_set &symm);


void freq_grad_irrep(const cartesians &carts, simples_class &simples, const salc_set &all_salcs);

void freq_grad_nosymm(const cartesians &carts, simples_class &simples, const salc_set &all_salcs);

void grad_energy(cartesians &carts, simples_class &simples, const salc_set &all_salcs);
void grad_save(const cartesians &carts);
void energy_save(void);

int opt_step(cartesians &carts, simples_class &simples, const salc_set &symm_salcs);

int opt_step_cart(cartesians &carts);

int *read_constraints(const simples_class &simples);
void opt_report(FILE *of);
int disp_freq_grad_cart(const cartesians &carts);
void freq_grad_cart(const cartesians &carts);
int disp_freq_energy_cart(const cartesians &carts);

void freq_energy_cart(void);

int test_B(const cartesians &carts, simples_class &simples, const salc_set &symm);

void H_update_cart(double **H, cartesians & carts);
void sort_evals_all(int nsalc_all, double *evals_all, int *evals_all_irrep);
void empirical_H(const simples_class & simples, const salc_set &symm, const cartesians &carts);

// functions involving z-matrices
void zmat_to_intco(void);
void compute_zmat(const cartesians &carts, int *unique_zvars);
void print_zmat(FILE *outfile, int *unique_zvars);


double **compute_H(simples_class & simples, const salc_set &symm,
    double **P, const cartesians &carts);
double **compute_H_cart(cartesians & carts, double **P);

void delocalize(const simples_class &simples, const cartesians &carts);
void rm_rotations(const simples_class & simples, const cartesians &carts,
    int &num_nonzero, double **evects);


bool new_geom(const cartesians &carts, simples_class &simples, const salc_set &all_salcs,
    double *dq, int print_flag, int restart_geom_file, char *disp_label, int disp_num,
    int last_disp, double *return_geom);

void step_limit(const simples_class &, const salc_set &, double *dq);
void check_zero_angles(const simples_class &, const salc_set &, double *dq);

void fconst_init(const cartesians &carts, const simples_class &simples, const salc_set &symm);
void fconst_init_cart(cartesians &carts);


}} /* namespace psi::optking */

#endif
