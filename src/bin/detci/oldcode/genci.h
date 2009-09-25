/*
** GENCI.H
** 
** Contains the function prototypes for the C routines in the GENCI
** module.
**
*/

#ifndef _psi_src_bin_detci_genci_h
#define _psi_src_bin_detci_genci_h

namespace psi { namespace detci {

int schmidt_addoc(double *buffer4, double *buffer5, int buf_size, 
                  int extra_buf, int num_buf, PSI_FPTR d_index, 
                  int N, int L, int b_file, int d_file);
void v_normalize(double *A, PSI_FPTR index, int buf_size, 
                 int extra_buf, int num_buf, int d_file);
double *v_schmidt(double *buffer4, double *buffer5, int buf_size, 
                  int extra_buf, int num_buf, int N, int L, int b_file);
void det2strings(BIGINT det, int *alp_code, int *alp_idx,
                 int *bet_code, int *bet_idx);
BIGINT strings2det(int alp_code, int alp_idx, int bet_code, int bet_idx);
void unit_guess(int alp_code, int alp_idx, int bet_code, int bet_idx,
                int switch_buf3, double *buffer, int buf_size,
                int num_buf, int extra_buf, PSI_FPTR b_file,
                PSI_FPTR b_writ, int M, int N);
void max_element(double *buffer, int num_elements, double *max, int *max_num);
void min_element(double *buffer, int num_elements, double *min, int *min_num);
void read_c(int switch_buf3, double *buffer, int buf_size, int num_buf,
            int extra_buf, int b_file, PSI_FPTR b_writ,
            int c_file, PSI_FPTR c_index);

}} // namespace psi::detci

#endif // header guard

