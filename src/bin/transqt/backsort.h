/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here 
*/
#ifndef _psi3_bin_transqt_backsort_h_
#define _psi3_bin_transqt_backsort_h_

namespace psi { namespace transqt {

void backsort_prep(int uhf);
void backsort(int first_tmp_file, double tolerance, int uhf);
void backsort_write(int i, int j, double **A, int kfirst, int klast,
		    int lfirst, int llast, int printflag, FILE *outfile,
		       struct iwlbuf *twopdm_out, int uhf);

}} // end namespace psi::transqt
#endif // header guard
