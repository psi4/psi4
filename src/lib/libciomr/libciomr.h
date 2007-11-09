/* 
** Declarations for functions found in libciomr.a
**
** C. David Sherrill and T. Daniel Crawford
**
**
*/

#ifndef _psi_src_lib_libciomr_libciomr_h_
#define _psi_src_lib_libciomr_libciomr_h_

#include <stdio.h>
#include <libciomr/iomrparam.h>

#ifdef __cplusplus
extern "C" {
#endif

int psi_start(FILE** infile, FILE** outfile, char** psi_file_prefix, int argc, char *argv[], int overwrite_output);
int psi_stop(FILE* infile, FILE* outfile, char* psi_file_prefix);
char* psi_ifname();
char* psi_ofname();
char* psi_fprefix();

void ffile(FILE **fptr, char *suffix, int code);
void ffile_noexit(FILE **fptr, char *suffix, int code);
void ffileb(FILE **fptr, char *suffix, int code);
void ffileb_noexit(FILE **fptr, char *suffix, int code);

void add_arr(double *a, double *b, double *c, int n);
void add_mat(double **a,double **b,double **c,int n,int m);
void balance(double **a, int n);

/* Functions under block_alloc.c */
double *** block_mat_alloc(int n_so_typs,int num_ir,int *num_so);
void block_mat_dealloc(double ***array,int num_ir,int *num_so);
double ** block_arr_alloc(int n_so_typs,int num_ir,int *num_so);
void block_arr_dealloc(double **array,int n_so_typs);

/* Functions under dot.c */
void dot_arr(double *a, double *b, int size, double *value) ;
void dot_mat(double **a,double **b,int n,double *value);

void eigout(double **a,double *b,double *c,int m,int n,FILE *out);
void eigsort(double *d,double **v,int n);
void eivout(double **a, double *b, int m, int n, FILE *out) ;
void mosort(double *d, double **v, int *sym, int nso, int nmo);

/* Functions under errors.c */
void no_path_given(char *name);
void malloc_check(char *caller,char *data);
void fopen_check(char *caller,char *path,char *data);
void fread_error(char *caller);
void fwrite_error(char *caller);

void flin(double **a,double *b,int in,int im,double *det);
void free_matrix(double **array, unsigned long int size) ;
#ifdef DEC
int get_file_info(char *token,char *format,char *val);
int get_param(char *token,char *format,char *val);
#else
int get_file_info(char *token,char *format,void *val);
int get_param(char *token,char *format,void *val);
#endif
int i2sec(PSI_FPTR n);
double * init_array(unsigned long int size) ;
double ** init_matrix(unsigned long int rows, unsigned long int cols) ;
void init_ptrs();

/* Functions under int_pac.c */
void int_pac(int *i,int *ib,int *j,int *jb,unsigned int *ipak,
             unsigned int *jpak);
void int_unpac(int *i,int *ib,int *j,int *jb,unsigned int *ipak,
               unsigned int *jpak);

int io_getline(FILE *input,char *line);
int io_locate(FILE *input, char *loc_token);

/* Functions under ioopen.c */
void ioinit_();
void ioopen_(int *unit);
void ioclos_(int *unit, int *status);
void iowrr_(int *unit,char *buffer,PSI_FPTR *first, int *length);
void iordr_(int *unit,char *buffer,PSI_FPTR *first, int *length);
void ioabort();

void lubksb(double **a,int n,int *indx,double *b);
void ludcmp(double **a,int n,int *indx,double *d);

/* Functions under mat_to_arr.c */
void mat_to_arr(double **a,double *b, int m, int n);
void arr_to_mat(double **a,double *b,int m,int n);


void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
           int nr, int nl, int nc, int add) ;
void mxmb(double **a,int ia,int ja,double **b,int ib,int jb,double **c,
          int ic,int jc, int nrow, int nlnk, int ncol);
int oldstyleinput();
void print_array(double *a, int m, FILE *out) ;
void print_mat(double **a, int rows, int cols, FILE *out) ;
void rclose(int unit, int status) ;
void rfile(int unit) ;
void rgetsa(int unit,int *iadr);
void rread(int itape, char *array, int nlen, int irec) ;
void rsetsa(int unit,int address);
void rsp(int nm, int n, int nv, double *array, double *evals, int matz,
         double **evecs, double toler) ;
void rwrit(int itape, char *array, int nlen, int irec) ;
PSI_FPTR sec2i(int n);
void sq_rsp(int nm, int n, double **array, double *evals, int matz, 
            double **evecs, double toler) ;
void sq_to_tri(double **bmat,double *amat,int size);
void sread(int itape, char *array, int nlen) ;
void srew(int itape) ;
void swrit(int itape,char *array,int nlen);

/* Functions under tri_to_block.c */
void tri_to_block(double *a,double **b,int num_ir,int *num_so,int *ioff);
void block_to_tri(double *a,double **b,int num_ir,int *num_so,int *ioff);

void tri_to_sq(double *amat,double **bmat,int size);

/* Functions under tstart.c */
void tstart(FILE *outfile) ;
void tstop(FILE *outfile) ;

void wreadw(int tape, char *buffer, int size, PSI_FPTR fword, 
            PSI_FPTR *nxtwrd);
void wwritw(int unit, char *buffer, int nwords, PSI_FPTR fword, 
            PSI_FPTR *nxtwrd);

/* Functions in zero.c */
void zero_arr(double *a, int size) ;
void zero_mat(double **a, int rows, int cols) ;

PSI_FPTR get_file_ptr(int unit);

/* Functions in flen.c */
PSI_FPTR flen(int unit);

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
double ** block_matrix(unsigned long int n, unsigned long int m);
void free_block(double **array);

/* Functions in fndcor */
void fndcor(long int *maxcrb, FILE *infile, FILE *outfile);

#ifdef __cplusplus
}
#endif

#endif /* header guard */
