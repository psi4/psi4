/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

extern int *correlate(char *ptgroup, int irrep, int *nirreps_ref, int *nirreps);

void read_chkpt_geom()
{
  int i, j, errcod, disp_irrep, *correlation, *clsdpi_ref, *clsdpi, mm;
  int atom, nirreps_ref, nirreps, h, *openpi, *openpi_ref, *states_per_irrep;
  int nso, nmo, h_ref, cnt, *cnt_orb_irr, *orbspi, *orbspi_ref, *orbs_off_ref, *so_off;
  int ref, *sopi, *sopi_ref, *orbs_off;
  double Z = 0.0, **scf, **scf_col, **scf_ref, escf_ref;
  double tmp = 0.0;
  char *atom_label, *ptgrp_ref, *save_prefix;

  /*** read geometry with Zvals from chkpt file */
  chkpt_init(PSIO_OPEN_OLD);
  num_atoms = chkpt_rd_natom();
  num_allatoms = chkpt_rd_nallatom();
  if (num_atoms == 0)
    punt("GEOMETRY in the checkpoint file is empty!");
  else if (num_atoms > MAXATOM)
    punt("There are more atoms than allowed!");
  full_geom = chkpt_rd_fgeom();
  geometry = (double **) malloc(num_atoms*sizeof(double *));
  atom_dummy = chkpt_rd_atom_dummy();
  atom = 0;
  for(i=0;i<num_allatoms;i++)
    if (!atom_dummy[i]) {
      geometry[atom] = full_geom[i];
      ++atom;
    }
  nuclear_charges = chkpt_rd_zvals();
  full_element = chkpt_rd_felement();
  /* Grab subgroup and get rid of a possible blank
  subgroup = chkpt_rd_sym_label();
  if (subgroup[2] == ' ') subgroup[2] = '\0';
  */

  /* these are set by optking */
  disp_irrep = chkpt_rd_disp_irrep();
  save_prefix = chkpt_rd_prefix();

  /*** read symmetry info and MOs for undisplaced geometry from
       root section of checkpoint file ***/
  chkpt_reset_prefix();
  chkpt_commit_prefix();

  ptgrp_ref = chkpt_rd_sym_label();
  clsdpi_ref = chkpt_rd_clsdpi(); /*closed MOs per irrep*/
  openpi_ref = chkpt_rd_openpi(); /*open MOs per irrep*/

  /* Lookup irrep correlation table */
  correlation = correlate(ptgrp_ref, disp_irrep, &nirreps_ref, &nirreps);

  if (print_lvl > 2) {
    fprintf(outfile,"Reference point group is %s.\n", ptgrp_ref);
    free(ptgrp_ref);
    fprintf(outfile,"Irrep of this displacement is %d.\n", disp_irrep);
    fprintf(outfile,"Irrep correlation:");
    for (i=0; i<nirreps_ref; ++i)
      fprintf(outfile," %d",correlation[i]);
    fprintf(outfile,"\n");
  }

  /* build orbital information for current point group */
  clsdpi = init_int_array(nirreps);
  openpi = init_int_array(nirreps);
  for (h=0; h < nirreps_ref; ++h) {
    clsdpi[ correlation[h] ] += clsdpi_ref[h];
    openpi[ correlation[h] ] += openpi_ref[h];
  }

  chkpt_set_prefix(save_prefix);
  chkpt_commit_prefix();
  free(save_prefix);

  /* write orbital information to chkpt file */
  chkpt_wt_nirreps(nirreps);
  chkpt_wt_clsdpi(clsdpi);
  chkpt_wt_openpi(openpi);

  if (cc_wfn(wfn)) {
    if (cc_excited(wfn)) {
      states_per_irrep = init_int_array(nirreps);
      for (h=0; h < nirreps_ref; ++h) {
        errcod = ip_data("STATES_PER_IRREP","%d",&(i),1,h);
        states_per_irrep[ correlation[h] ] += i;
      }
      chkpt_wt_statespi(states_per_irrep);
      if (print_lvl > 2) {
        fprintf(outfile,"states_per_irrep");
        for (h=0; h < nirreps; ++h)
          fprintf(outfile, " %d",states_per_irrep[h]);
        fprintf(outfile,"\n");
      }
      free(states_per_irrep);
    }
  }

  if (print_lvl > 2) {
    fprintf(outfile,"clsdpi");
    for (h=0; h < nirreps; ++h)
      fprintf(outfile, " %d",clsdpi[h]);
    fprintf(outfile,"\n");
    fprintf(outfile,"openpi");
    for (h=0; h < nirreps; ++h)
      fprintf(outfile, " %d",openpi[h]);
    fprintf(outfile,"\n");
    fprintf(outfile,"orbspi");
    for (h=0; h < nirreps; ++h)
      fprintf(outfile, " %d",orbspi[h]);
    fprintf(outfile,"\n");
  }

  chkpt_close();

  free(clsdpi); free(clsdpi_ref);
  free(openpi); free(openpi_ref);

  element = (char **) malloc(sizeof(char *)*num_atoms);
  elemsymb_charges = init_array(num_atoms);
  for(i=0;i<num_atoms;i++) {
    element[i] = elem_name[(int)nuclear_charges[i]];
    elemsymb_charges[i] = nuclear_charges[i];
  }

  return;
}

}} // namespace psi::input
