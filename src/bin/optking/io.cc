/*! \file
    \ingroup OPTKING
    \brief io/library interface
*/

#define EXTERN
#include "globals.h"
#undef EXTERN

#include <libpsio/psio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

namespace psi { namespace optking {

// open a text file for 
// code == 0 : writing
// code == 1 : writing/appending
// code == 2 : reading
void opt_ffile(FILE **fptr, const char *suffix, int code) {
  ffile(fptr, suffix, code);
}

// same as above but do not abort if requested file is not present
void opt_ffile_noexit(FILE **fptr, const char *suffix, int code) {
  ffile_noexit(fptr, suffix, code);
}

// determine input and output file pointers
int opt_psi_start(FILE** infile, FILE** outfile, char** psi_file_prefix,
int argc, char *argv[], int overwrite_output) {
    int i;
    i = psi_start(infile, outfile, psi_file_prefix, argc, argv, overwrite_output);
    return i;
}

int opt_psi_stop(FILE* infile, FILE* outfile, char* psi_file_prefix) {
  ip_done();
  free(psi_file_prefix);
  fclose(outfile);
  fclose(infile);
  return(PSI_RETURN_SUCCESS);
}

void exit_io(void) {
  fprintf(outfile,"\n******** OPTKING execution completed ********\n\n");
  psio_done();
  opt_psi_stop(infile,outfile,psi_file_prefix);
}

// matrix multiplication
void opt_mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
  int nr, int nl, int nc, int add) {
    mmult(AF, ta, BF, tb, CF, tc, nr, nl, nc, add);
}

/*
** sq_rsp(): diagomalize a symmetric square matrix ('array').
**
** \param nm     = rows of matrix
** \param n      = columns of matrix
** \param nv     = number of elements in lower triangle (n*(n+1)/2)
** \param array  = matrix to diagonalize
** \param e_vals = array to hold eigenvalues
** \param matz   = 0 (no eigenvectors, eigenvals in ascending order)
**               = 1 (eigenvectors and eigenvalues in ascending order)
**               = 2 (no eigenvectors, eigenvalues in descending order)
**               = 3 (eigenvectors and eigenvalues in descending order)
** \param e_vecs = matrix of eigenvectors (one column for each eigvector)
** \param toler  = tolerance for eigenvalues?  Often 1.0E-14.
*/
void opt_sq_rsp(int nm, int n, double **array, double *e_vals, int matz,
  double **e_vecs, double toler) {
    sq_rsp(nm, n, array, e_vals, matz, e_vecs, toler);
}



}} /* namespace psi::optking */

