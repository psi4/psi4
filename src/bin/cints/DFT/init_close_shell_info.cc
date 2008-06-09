/*! \file init_close_shell_info.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
namespace psi { namespace CINTS {
struct close_shell_info_s init_close_shell_info(void){
    
    struct close_shell_info_s close;
    close.shells_close_to_chunk = (int *)malloc(sizeof(int)*BasisSet.num_shells);
    close.aos_close_to_chunk = (int *)malloc(sizeof(int)*BasisSet.num_ao);
    close.close_shells_per_am = (int *)malloc(sizeof(int)*BasisSet.max_am);
    /*close.close_COCC = (double **)malloc(sizeof(double *)*BasisSet.num_ao);*/
    if(UserOptions.reftype == uhf){
	close.close_COCC_a = block_matrix(BasisSet.num_ao,MOInfo.alpha_occ);
	close.close_COCC_b = block_matrix(BasisSet.num_ao,MOInfo.alpha_occ);
    }
    else{
	close.close_COCC = block_matrix(BasisSet.num_ao,MOInfo.ndocc);
    }
    
    return close;
}

void free_close_shell_info(struct close_shell_info_s close_shell_info){
    
    free(close_shell_info.shells_close_to_chunk);
    free(close_shell_info.close_shells_per_am);
    if(UserOptions.reftype == uhf){
	free(close_shell_info.close_COCC_a);
	free(close_shell_info.close_COCC_b);
    }
    else{
	free(close_shell_info.close_COCC);
    }
}

void print_close_shell_info(struct close_shell_info_s close){
    
    int i;

    fprintf(outfile,"\nClose Shell Info Data Structure");
    
    fprintf(outfile,"\nNumber of close AO's = %d",close.num_close_aos);
    fprintf(outfile,"\nNumber of close Shells = %d",close.num_close_shells);
    for(i=0;i<BasisSet.num_shells;i++)
	fprintf(outfile,"\nshell_close_to_chunk[%d] = %d",i,close.shells_close_to_chunk[i]);
    for(i=0;i<BasisSet.num_ao;i++)
	fprintf(outfile,"\naos_close_to_chunk[%d] = %d",i,close.aos_close_to_chunk[i]);
    for(i=0;i<BasisSet.max_am;i++)
	fprintf(outfile,"\nclose_shells_per_am[%d] = %d",i,close.close_shells_per_am[i]);
    if(UserOptions.reftype == uhf){
	fprintf(outfile,"\nClose Alpha Occupied Eigenvector");
	print_mat(close.close_COCC_a,close.num_close_aos,MOInfo.alpha_occ,outfile);
	fprintf(outfile,"\n\n");
	fprintf(outfile,"\nClose Beta Occupied Eigenvector");
	print_mat(close.close_COCC_b,close.num_close_aos,MOInfo.beta_occ,outfile);
	fprintf(outfile,"\n\n");
    }
    else{
	fprintf(outfile,"\nClose Occupied Eigenvector");
	print_mat(close.close_COCC,close.num_close_aos,MOInfo.ndocc,outfile);
	fprintf(outfile,"\n\n");
    }
    fflush(outfile);
}
};};
