/*
** Declarations for functions found in libciomr.a
**
** C. David Sherrill and T. Daniel Crawford
**
**
*/

#ifndef _psi_src_lib_libciomr_libciomr_h_
#define _psi_src_lib_libciomr_libciomr_h_

#include <cstdio>

namespace psi {

int psi_start(FILE** infile, FILE** outfile, char** psi_file_prefix, int argc, char *argv[], int overwrite_output);
int psi_stop(FILE* infile, FILE* outfile, char* psi_file_prefix);
char* psi_ifname();
char* psi_ofname();
char* psi_fprefix();

void ffile(FILE **fptr, const char *suffix, int code);
void ffile_noexit(FILE **fptr, char *suffix, int code);
void ffileb(FILE **fptr, char *suffix, int code);
void ffileb_noexit(FILE **fptr, char *suffix, int code);

void add_arr(double *a, double *b, double *c, int n);
void add_mat(double **a,double **b,double **c,int n,int m);
void balance(double **a, int n);

/* Functions under dot.cc */
void dot_arr(double *a, double *b, int size, double *value) ;
double dot_mat(double **a,double **b,int n);

void eigout(double **a,double *b,double *c,int m,int n,FILE *out);
void eigsort(double *d,double **v,int n);
void eivout(double **a, double *b, int m, int n, FILE *out) ;
void mosort(double *d, double **v, int *sym, int nso, int nmo);

void flin(double **a,double *b,int in,int im,double *det);
void free_matrix(double **array, unsigned long int size) ;
double * init_array(unsigned long int size) ;
double ** init_matrix(unsigned long int rows, unsigned long int cols) ;

void lubksb(double **a,int n,int *indx,double *b);
void ludcmp(double **a,int n,int *indx,double *d);

/* Functions under mat_to_arr.c */
void mat_to_arr(double **a,double *b, int m, int n);
void arr_to_mat(double **a,double *b,int m,int n);


void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
           int nr, int nl, int nc, int add) ;
void mxmb(double **a,int ia,int ja,double **b,int ib,int jb,double **c,
          int ic,int jc, int nrow, int nlnk, int ncol);
void print_array(double *a, int m, FILE *out) ;
void print_mat(double **a, int rows, int cols, FILE *out) ;

void rsp(int nm, int n, int nv, double *array, double *evals, int matz,
         double **evecs, double toler) ;
void sq_rsp(int nm, int n, double **array, double *evals, int matz,
            double **evecs, double toler) ;
void sq_to_tri(double **bmat,double *amat,int size);

/* Functions under tri_to_block.c */
void tri_to_block(double *a,double **b,int num_ir,int *num_so,int *ioff);
void block_to_tri(double *a,double **b,int num_ir,int *num_so,int *ioff);

void tri_to_sq(double *amat,double **bmat,int size);

/* Functions under tstart.c */
void tstart() ;
void tstop() ;

/* Functions in zero.c */
void zero_arr(double *a, int size) ;
void zero_mat(double **a, int rows, int cols) ;

/* Functions in int_array.c */
int * init_int_array(int size) ;
void zero_int_array(int *a, int size);
int **init_int_matrix(int rows, int cols);
void free_int_matrix(int **array);
void zero_int_matrix(int **array, int rows, int cols);
void print_int_mat(int **a, int m, int n, FILE *out);

/* Functions in long_int_array.c */
long int * init_long_int_array(int size) ;
void zero_long_int_array(long int *a, int size);
long int **init_long_int_matrix(int rows, int cols);
void free_long_int_matrix(long int **array);
void zero_long_int_matrix(long int **array, int rows, int cols);
void print_long_int_mat(long int **a, int m, int n, FILE *out);

/* Functions in block_matrix.c */
double ** block_matrix(unsigned long int n, unsigned long int m, bool mlock = false);
void free_block(double **array);

/* Functions in fndcor */
void fndcor(long int *maxcrb, FILE *infile, FILE *outfile);

}

#endif /* header guard */
