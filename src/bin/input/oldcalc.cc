/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>

#include "input.h"
#include "defines.h"
#define EXTERN
#include "global.h"

#include "float.h"
#include "linalg.h"

namespace psi { namespace input {

void init_oldcalc()
{
  int s, ncalcs;

  /*--- initialize the old checkpoint file ---*/
  chkpt_init(PSIO_OPEN_OLD);

  /*--- check if the symmetry is the same ---*/
  Oldcalc.symmetry = chkpt_rd_sym_label();
  if (strncmp(Oldcalc.symmetry,symmetry,3)) {
      fprintf(outfile,"  Old chkpt file uses different point group from the current one\n");
      fprintf(outfile,"  Will continue without the MO projection\n");
      fprintf(stderr,"  Old chkpt file uses different point group from the current one\n");
      fprintf(stderr,"  Will continue without the MO projection\n");
      chkpt_mos = 0;
      chkpt_close();
      return;
  }
  Oldcalc.nirreps = chkpt_rd_nirreps();

  /*--- read the old geometry ---*/
  Oldcalc.natom = chkpt_rd_natom();
  Oldcalc.geometry = chkpt_rd_geom();
  Oldcalc.Rref = chkpt_rd_rref();

  /*--- read the old basis ---*/
  Oldcalc.num_ao = chkpt_rd_nao();
  Oldcalc.num_so = chkpt_rd_nso();
  Oldcalc.sopi = chkpt_rd_sopi();
  Oldcalc.num_shells = chkpt_rd_nshell();
  Oldcalc.num_prims = chkpt_rd_nprim();
  Oldcalc.exponents = chkpt_rd_exps();
  Oldcalc.contr_coeff = chkpt_rd_contr_full();
  Oldcalc.shell_nucleus = chkpt_rd_snuc();
  Oldcalc.first_prim_shell = chkpt_rd_sprim();
  Oldcalc.shell_ang_mom = chkpt_rd_stype();
  Oldcalc.nprim_in_shell = chkpt_rd_snumg();
  Oldcalc.first_ao_shell = chkpt_rd_sloc();
  Oldcalc.usotao = chkpt_rd_usotao();

  /*--- Convert to C-numbering ---*/
  Oldcalc.max_angmom = 0;
  for(s=0;s<Oldcalc.num_shells;s++) {
      Oldcalc.shell_nucleus[s]--;
      Oldcalc.shell_ang_mom[s]--;
      Oldcalc.first_ao_shell[s]--;
      Oldcalc.first_prim_shell[s]--;
      if (Oldcalc.shell_ang_mom[s] > Oldcalc.max_angmom)
	  Oldcalc.max_angmom = Oldcalc.shell_ang_mom[s];
  }

  /*--- read the old SCF eigenvector, RHF case so far ---*/
  Oldcalc.ref = (reftype)chkpt_rd_ref();
  if (Oldcalc.ref == ref_rhf   ||
      Oldcalc.ref == ref_rohf  ||
      Oldcalc.ref == ref_tcscf ||
      Oldcalc.ref == ref_rks)
      Oldcalc.spinrestr_ref = 1;
  else if (Oldcalc.ref == ref_uhf   ||
	   Oldcalc.ref == ref_uks)
      Oldcalc.spinrestr_ref = 0;
  else {
      fprintf(outfile,"  Input is not aware of reference type %d\n",Oldcalc.ref);
      fprintf(outfile,"  Will continue without the MO projection\n");
      chkpt_mos = 0;
      chkpt_close();
      return;
  }

  Oldcalc.iopen = chkpt_rd_iopen();
  Oldcalc.num_mo = chkpt_rd_nmo();
  Oldcalc.orbspi = chkpt_rd_orbspi();
  Oldcalc.clsdpi = chkpt_rd_clsdpi();
  Oldcalc.openpi = chkpt_rd_openpi();
  Oldcalc.escf = chkpt_rd_escf();
  if (Oldcalc.spinrestr_ref) { /* Spin-restricted case */
      Oldcalc.scf_evect_so = chkpt_rd_scf();
  }
  else { /* Spin-unrestricted case */
      Oldcalc.scf_evect_so_alpha = chkpt_rd_alpha_scf();
      Oldcalc.scf_evect_so_beta = chkpt_rd_beta_scf();
  }

  /* read the old localized MOs, if present */
  Oldcalc.local = chkpt_rd_local_scf();

  /*--- Done with the checkpoint file ---*/
  chkpt_close();

  if (print_lvl >= DEBUGPRINT) {
    if (Oldcalc.spinrestr_ref) {
      fprintf(outfile,"  -Old eigenvector (in SO basis):\n");
      print_mat(Oldcalc.scf_evect_so,Oldcalc.num_so,Oldcalc.num_mo,outfile);
    }
    else {
      fprintf(outfile,"  -Old alpha eigenvector (in SO basis):\n");
      print_mat(Oldcalc.scf_evect_so_alpha,Oldcalc.num_so,Oldcalc.num_mo,outfile);
      fprintf(outfile,"  -Old beta eigenvector (in SO basis):\n");
      print_mat(Oldcalc.scf_evect_so_beta,Oldcalc.num_so,Oldcalc.num_mo,outfile);
    }
  }
  
  return;
}


void cleanup_oldcalc()
{
  free_block(Oldcalc.geometry);
  free_block(Oldcalc.Rref);
  free(Oldcalc.exponents);
  free_block(Oldcalc.contr_coeff);
  free(Oldcalc.shell_nucleus);
  free(Oldcalc.first_prim_shell);
  free(Oldcalc.shell_ang_mom);
  free(Oldcalc.nprim_in_shell);
  free(Oldcalc.first_ao_shell);
  free_block(Oldcalc.usotao);

  free(Oldcalc.orbspi);
  free(Oldcalc.clsdpi);
  free(Oldcalc.openpi);
  if (Oldcalc.spinrestr_ref)
      free_block(Oldcalc.scf_evect_so);
  else {
      free_block(Oldcalc.scf_evect_so_alpha);
      free_block(Oldcalc.scf_evect_so_beta);
  }

  if(Oldcalc.local != NULL) free_block(Oldcalc.local);

  return;
}


void store_oldcalc()
{
  FLOAT **S12;

  if (max_angmom > Oldcalc.max_angmom)
    init_gto(max_angmom);
  else
    init_gto(Oldcalc.max_angmom);

  S12 = overlap_new_old_float();
  
  /*--- write things out ---*/
  psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_NEW);
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Num. atoms", (char*) &Oldcalc.natom,
		   sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Geometry", (char*) Oldcalc.geometry[0],
		   Oldcalc.natom*3*sizeof(double));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Transmat to reference frame", (char*) Oldcalc.Rref[0],
		   9*sizeof(double));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Num. irreps", (char*) &Oldcalc.nirreps,
		   sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:SO's per irrep", (char*) Oldcalc.sopi,
		   Oldcalc.nirreps*sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Reference", (char*) &Oldcalc.ref,
		   sizeof(reftype));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Num. SO", (char*) &Oldcalc.num_so,
		   sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Num. MO", (char*) &Oldcalc.num_mo,
		   sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:MO's per irrep", (char*) Oldcalc.orbspi,
		   Oldcalc.nirreps*sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Closed shells per irrep", (char*) Oldcalc.clsdpi,
		   Oldcalc.nirreps*sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:Open shells per irrep", (char*) Oldcalc.openpi,
		   Oldcalc.nirreps*sizeof(int));
  if (Oldcalc.spinrestr_ref) { /* Spin-restricted case */
    psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs alpha", (char*) Oldcalc.scf_evect_so[0],
		     Oldcalc.num_so*Oldcalc.num_mo*sizeof(double));
  }
  else { /* Spin-unrestricted case */
    psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs alpha", (char*) Oldcalc.scf_evect_so_alpha[0],
		     Oldcalc.num_so*Oldcalc.num_mo*sizeof(double));
    psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs beta", (char*) Oldcalc.scf_evect_so_beta[0],
		     Oldcalc.num_so*Oldcalc.num_mo*sizeof(double));
  }
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:AO->SO", (char *) Oldcalc.usotao[0],
		   Oldcalc.num_ao*Oldcalc.num_so*sizeof(double));
  psio_write_entry(PSIF_OLD_CHKPT, ":PrevCalc:New-Old basis overlap", (char*) S12[0],
		   num_so*Oldcalc.num_so*sizeof(FLOAT));
  psio_close(PSIF_OLD_CHKPT, 1);

  delete_matrix(S12);
}

}} // namespace psi::input
