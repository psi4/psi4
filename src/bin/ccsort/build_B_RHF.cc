/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

/* build_B_RHF_singlet(): Builds the RHF orbital Hessian for computing
** spatial- and spin-symmetry instabilities.  In spin orbitals the Hessian
** matrix is:
**
** B(ai,bj) = delta_ij f_ab - delta_ab f_ij + <ij||ab> - <ja||ib>
**
** RHF references and singlet eigenstates:
**  B(AI,BJ) = delta_IJ f_AB - delta_AB f_IJ + <IJ|BA> - <IA|JB>
**
** This routine will build each spin component of the entire matrix 
** for later in-core diagonalization.
**
** TDC, March 2003
**
** Adapted for use in local-CC response calculations by TDC, April
** 2004.  Note that this version of the MO Hessian is correct for
** non-canonical Hartree-Fock orbitals.
*/

void build_B_RHF(double omega)
{
  int h, nirreps;
  int a, b, i, j, ai, bj, A, B, I, J, Asym, Bsym, Isym, Jsym;
  dpdbuf4 C, D, Bmat;
  dpdfile2 fIJ, fij, fAB, fab;
	
  nirreps = moinfo.nirreps;
	
  psio_open(PSIF_MO_HESS, 0);
	
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&D, PSIF_MO_HESS, sprq, 11, 11, "B(AI,BJ)");
  dpd_buf4_close(&D);
	
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_sort_axpy(&C, PSIF_MO_HESS, qpsr, 11, 11, "B(AI,BJ)", -1);
  dpd_buf4_close(&C);
	
  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
  dpd_file2_mat_init(&fab);
  dpd_file2_mat_rd(&fab);
	
  dpd_buf4_init(&Bmat, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "B(AI,BJ)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Bmat, h);
    dpd_buf4_mat_irrep_rd(&Bmat, h);
    for(ai=0; ai < Bmat.params->rowtot[h]; ai++) {
      a = Bmat.params->roworb[h][ai][0];
      i = Bmat.params->roworb[h][ai][1];
      A = fAB.params->rowidx[a];
      I = fIJ.params->rowidx[i];
      Asym = fAB.params->psym[a];
      Isym = fIJ.params->psym[i];
      for(bj=0; bj < Bmat.params->coltot[h]; bj++) {
	b = Bmat.params->colorb[h][bj][0];
	j = Bmat.params->colorb[h][bj][1];
	B = fAB.params->colidx[b];
	J = fIJ.params->colidx[j];
	Bsym = fAB.params->qsym[b];
	Jsym = fIJ.params->qsym[j];
	if((A==B) && (Isym==Jsym)) Bmat.matrix[h][ai][bj] -= fIJ.matrix[Isym][I][J];
	if((I==J) && (Asym==Bsym)) Bmat.matrix[h][ai][bj] += fAB.matrix[Asym][A][B];
	if(ai==bj) Bmat.matrix[h][ai][bj] -= omega;
      }
    }
    dpd_buf4_mat_irrep_wrt(&Bmat, h);
    dpd_buf4_mat_irrep_close(&Bmat, h);
  }
  dpd_buf4_close(&Bmat);
	
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fab);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fAB);
  dpd_file2_mat_close(&fij);
  dpd_file2_close(&fij);
  dpd_file2_mat_close(&fIJ);
  dpd_file2_close(&fIJ);
	
  psio_close(PSIF_MO_HESS, 1);
}

}} // namespace psi::ccsort
