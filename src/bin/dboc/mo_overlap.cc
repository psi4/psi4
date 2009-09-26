/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include "defines.h"
#include "float.h"
#include "params.h"
#include "linalg.h"
#include "moinfo.h"
#include "mo_overlap.h"
#include "hfwfn.h"
#include <libbasis/basisset.h>
#include <libbasis/overlap.h>
#include <libbasis/rotation.h>
#include <psi4-dec.h>

namespace psi { namespace DBOC {

extern void done(const char * message);
extern BasisSet* BasisSets[MAX_NUM_DISP];
extern HFWavefunction* HFVectors[MAX_NUM_DISP];
extern Params_t Params;

FLOAT **eval_S_alpha(DisplacementIndex LDisp, DisplacementIndex RDisp)
{
  
  HFWavefunction* HFWfn_R = HFVectors[RDisp];
  HFWavefunction* HFWfn_L = HFVectors[LDisp];
  int num_ao = HFWfn_R->num_ao();
  double** hf_evec_l = HFWfn_L->alpha_evec();
  double** aotoso_l = HFWfn_L->aotoso();
  double** rref_l = HFWfn_L->rref();
  double** hf_evec_r = HFWfn_R->alpha_evec();
  double** aotoso_r = HFWfn_R->aotoso();
  double** rref_r = HFWfn_R->rref();

#if USE_MOINFO
  extern MOInfo_t MOInfo;
  int num_mo = MOInfo.num_mo;
  int num_so = MOInfo.num_so;
#else
  int num_mo = HFWfn_R->num_mo();
  int num_so = HFWfn_R->num_so();
#endif

  //
  // Convert matrices of doubles into matrices of FLOAT's
  //
  FLOAT** hf_evec_r_FLOAT = convert_matrix(hf_evec_r, num_so, num_mo, 0);
  FLOAT** hf_evec_l_FLOAT_transp = convert_matrix(hf_evec_l, num_so, num_mo, 1);

  // Compute plus/minus overlap
  OverlapEngine overlap(BasisSets[LDisp],BasisSets[RDisp]);
  double** Slr_AO_AO = overlap.compute_full_matrix();

  // Rotate bases to the original coordinate systems (prior to reorientation into principal axis system)
  RotationOp Rop_r(BasisSets[RDisp]);
  RotationOp Rop_l(BasisSets[LDisp]);
  double** basisRref_r = Rop_r.full_rotation_mat(rref_r);
  double** basisRref_l = Rop_l.full_rotation_mat(rref_l);

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  -Rotation matrix for AO basis (disp = %d)\n", RDisp);
    psi::DBOC::print_mat(basisRref_r, num_ao, num_ao, outfile);
    
    fprintf(outfile, "  -Rotation matrix for AO basis (disp = %d)\n", LDisp);
    psi::DBOC::print_mat(basisRref_l, num_ao, num_ao, outfile);

    double** tmp_l1 = block_matrix(num_ao, num_mo);
    double** tmp_l2 = block_matrix(num_ao, num_mo);
    double** tmp_r1 = block_matrix(num_ao, num_mo);
    double** tmp_r2 = block_matrix(num_ao, num_mo);
    mmult(aotoso_r,1,hf_evec_r,0,tmp_r1,0,num_ao,num_so,num_mo,0);
    mmult(basisRref_r,0,tmp_r1,0,tmp_r2,0,num_ao,num_ao,num_mo,0);
    mmult(aotoso_l,1,hf_evec_l,0,tmp_l1,0,num_ao,num_so,num_mo,0);
    mmult(basisRref_l,0,tmp_l1,0,tmp_l2,0,num_ao,num_ao,num_mo,0);
    
    fprintf(outfile, "  -Original alpha eigenvector (disp = %d)\n", RDisp);
    psi::DBOC::print_mat(tmp_r1, num_ao, num_mo, outfile);
    fprintf(outfile, "  -Rotated alpha eigenvector (disp = %d)\n", RDisp);
    psi::DBOC::print_mat(tmp_r2, num_ao, num_mo, outfile);
    
    fprintf(outfile, "  -Original alpha eigenvector (disp = %d)\n", LDisp);
    psi::DBOC::print_mat(tmp_l1, num_ao, num_mo, outfile);
    fprintf(outfile, "  -Rotated alpha eigenvector (disp = %d)\n", LDisp);
    psi::DBOC::print_mat(tmp_l2, num_ao, num_mo, outfile);
  }

  double** tmpmat = block_matrix(num_ao,num_ao);
  mmult(Slr_AO_AO,0,basisRref_r,0,tmpmat,0,num_ao,num_ao,num_ao,0);
  mmult(basisRref_l,1,tmpmat,0,Slr_AO_AO,0,num_ao,num_ao,num_ao,0);
  free_block(tmpmat);
  free_block(basisRref_l);
  free_block(basisRref_r);

  double** Slr_AO_SO = block_matrix(num_ao,num_so);
  mmult(Slr_AO_AO, 0, aotoso_r, 1, Slr_AO_SO, 0, num_ao, num_ao, num_so, 0);
  free_block(Slr_AO_AO);
  double** Slr_SO_SO = block_matrix(num_so,num_so);
  mmult(aotoso_l, 0, Slr_AO_SO, 0, Slr_SO_SO, 0, num_so, num_ao, num_so, 0);
  free_block(Slr_AO_SO);
  FLOAT** Slr_FLOAT = convert_matrix(Slr_SO_SO, num_so, num_so, 0);
  free_block(Slr_SO_SO);

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, " (%d/%d) overlap matrix (SO basis)\n", LDisp, RDisp);
    psi::DBOC::print_mat(Slr_FLOAT, num_so, num_so, outfile);
  }

  FLOAT** tmpmat1 = create_matrix(num_mo,num_so);
  if (matrix_mult(hf_evec_l_FLOAT_transp, num_mo, num_so, Slr_FLOAT, num_so, num_so, tmpmat1))
    done("matrix_mult failed. Report the problem to the author.");
  FLOAT** S = create_matrix(num_mo, num_mo);
  if (matrix_mult(tmpmat1, num_mo, num_so, hf_evec_r_FLOAT, num_so, num_mo, S))
    done("matrix_mult failed. Report the problem to the author.");

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  (%d/%d) alpha overlap matrix (MO basis)\n", LDisp, RDisp);
    psi::DBOC::print_mat(S, num_mo, num_mo, outfile);
  }

  delete_matrix(tmpmat1);
  delete_matrix(Slr_FLOAT);
  delete_matrix(hf_evec_l_FLOAT_transp);
  delete_matrix(hf_evec_r_FLOAT);

  return S;
}


FLOAT **eval_S_beta(DisplacementIndex LDisp, DisplacementIndex RDisp)
{
  HFWavefunction* HFWfn_R = HFVectors[RDisp];
  HFWavefunction* HFWfn_L = HFVectors[LDisp];
  int num_ao = HFWfn_R->num_ao();
  double** hf_evec_l = HFWfn_L->beta_evec();
  double** aotoso_l = HFWfn_L->aotoso();
  double** rref_l = HFWfn_L->rref();
  double** hf_evec_r = HFWfn_R->beta_evec();
  double** aotoso_r = HFWfn_R->aotoso();
  double** rref_r = HFWfn_R->rref();

#if USE_MOINFO
  extern MOInfo_t MOInfo;
  int num_mo = MOInfo.num_mo;
  int num_so = MOInfo.num_so;
#else
  int num_mo = HFWfn_R->num_mo();
  int num_so = HFWfn_R->num_so();
#endif

//
  // Convert matrices of doubles into matrices of FLOAT's
  //
  FLOAT** hf_evec_r_FLOAT = convert_matrix(hf_evec_r, num_so, num_mo, 0);
  FLOAT** hf_evec_l_FLOAT_transp = convert_matrix(hf_evec_l, num_so, num_mo, 1);

  // Compute plus/minus overlap
  OverlapEngine overlap(BasisSets[LDisp],BasisSets[RDisp]);
  double** Slr_AO_AO = overlap.compute_full_matrix();

  // Rotate bases to the original coordinate systems (prior to reorientation into principal axis system)
  RotationOp Rop_r(BasisSets[RDisp]);
  RotationOp Rop_l(BasisSets[LDisp]);
  double** basisRref_r = Rop_r.full_rotation_mat(rref_r);
  double** basisRref_l = Rop_l.full_rotation_mat(rref_l);

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  -Rotation matrix for AO basis (disp = %d)\n", RDisp);
    psi::DBOC::print_mat(basisRref_r, num_ao, num_ao, outfile);
    
    fprintf(outfile, "  -Rotation matrix for AO basis (disp = %d)\n", LDisp);
    psi::DBOC::print_mat(basisRref_l, num_ao, num_ao, outfile);

    double** tmp_l1 = block_matrix(num_ao, num_mo);
    double** tmp_l2 = block_matrix(num_ao, num_mo);
    double** tmp_r1 = block_matrix(num_ao, num_mo);
    double** tmp_r2 = block_matrix(num_ao, num_mo);
    mmult(aotoso_r,1,hf_evec_r,0,tmp_r1,0,num_ao,num_so,num_mo,0);
    mmult(basisRref_r,0,tmp_r1,0,tmp_r2,0,num_ao,num_ao,num_mo,0);
    mmult(aotoso_l,1,hf_evec_l,0,tmp_l1,0,num_ao,num_so,num_mo,0);
    mmult(basisRref_l,0,tmp_l1,0,tmp_l2,0,num_ao,num_ao,num_mo,0);
    
    fprintf(outfile, "  -Original beta eigenvector (disp = %d)\n", RDisp);
    psi::DBOC::print_mat(tmp_r1, num_ao, num_mo, outfile);
    fprintf(outfile, "  -Rotated beta eigenvector (disp = %d)\n", RDisp);
    psi::DBOC::print_mat(tmp_r2, num_ao, num_mo, outfile);
    
    fprintf(outfile, "  -Original beta eigenvector (disp = %d)\n", LDisp);
    psi::DBOC::print_mat(tmp_l1, num_ao, num_mo, outfile);
    fprintf(outfile, "  -Rotated beta eigenvector (disp = %d)\n", LDisp);
    psi::DBOC::print_mat(tmp_l2, num_ao, num_mo, outfile);
  }

  double** tmpmat = block_matrix(num_ao,num_ao);
  mmult(Slr_AO_AO,0,basisRref_r,0,tmpmat,0,num_ao,num_ao,num_ao,0);
  mmult(basisRref_l,1,tmpmat,0,Slr_AO_AO,0,num_ao,num_ao,num_ao,0);
  free_block(tmpmat);
  free_block(basisRref_l);
  free_block(basisRref_r);

  double** Slr_AO_SO = block_matrix(num_ao,num_so);
  mmult(Slr_AO_AO, 0, aotoso_r, 1, Slr_AO_SO, 0, num_ao, num_ao, num_so, 0);
  free_block(Slr_AO_AO);
  double** Slr_SO_SO = block_matrix(num_so,num_so);
  mmult(aotoso_l, 0, Slr_AO_SO, 0, Slr_SO_SO, 0, num_so, num_ao, num_so, 0);
  free_block(Slr_AO_SO);
  FLOAT** Slr_FLOAT = convert_matrix(Slr_SO_SO, num_so, num_so, 0);
  free_block(Slr_SO_SO);

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, " (%d/%d) overlap matrix (SO basis)\n", LDisp, RDisp);
    psi::DBOC::print_mat(Slr_FLOAT, num_so, num_so, outfile);
  }

  FLOAT** tmpmat1 = create_matrix(num_mo,num_so);
  if (matrix_mult(hf_evec_l_FLOAT_transp, num_mo, num_so, Slr_FLOAT, num_so, num_so, tmpmat1))
    done("matrix_mult failed. Report the problem to the author.");
  FLOAT** S = create_matrix(num_mo, num_mo);
  if (matrix_mult(tmpmat1, num_mo, num_so, hf_evec_r_FLOAT, num_so, num_mo, S))
    done("matrix_mult failed. Report the problem to the author.");

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  (%d/%d) beta overlap matrix (MO basis)\n", LDisp, RDisp);
    psi::DBOC::print_mat(S, num_mo, num_mo, outfile);
  }

  delete_matrix(tmpmat1);
  delete_matrix(Slr_FLOAT);
  delete_matrix(hf_evec_l_FLOAT_transp);
  delete_matrix(hf_evec_r_FLOAT);

  return S;
}

}} // namespace psi::DBOC
