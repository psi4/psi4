/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

using namespace psi;

/*-----------------------------------------------------------------------------------------------------------------
  This function writes information out to checkpoint file
 -----------------------------------------------------------------------------------------------------------------*/

void write_to_chkpt(double repulsion, const char *basiskey)
{

  char *calc_label;
  int i,j,k,l;
  int atom,ua,shell,irr,symop,shell_first,shell_last,us;
  int eq_atom;
  int *arr_int;
  double *arr_double;
  double **mat_double;
  int *ict, **ict_tmp;
  int *scf_pointers;
  double *cspd;     /*Array of contraction coefficients in chkpt file format*/
  char *atom_label, **tmp_atom_label;
  int **shell_transm;
  int **ioff_irr;
  int max_num_prims;
  int max_atom_degen;
  int max_angmom_unique;
  psio_address chkptr;

  /* Check the max_angmom. If it's >= MAXANGMOM from chkpt_params.h - die */
  if (max_angmom >= MAXANGMOM)
    punt("Angular momentum is too high to be handled by your version of Psi (chkpt file)");


  /*----------------------------------
    Write out the label then 80 zeros
    ----------------------------------*/
  calc_label = init_char_array(80);
  strncpy(calc_label,label,80);

  /*----------------------------------
    Write out basic info to chkpt
   ----------------------------------*/
  if (basiskey[0] == '\0') {
    chkpt_init(keep_chkpt ? PSIO_OPEN_OLD : PSIO_OPEN_NEW);
  }
  else {
    chkpt_init(PSIO_OPEN_OLD);
  }

  chkpt_wt_label(calc_label);
  chkpt_wt_num_unique_atom(num_uniques);
  chkpt_wt_num_unique_shell(num_unique_shells, basiskey);
  chkpt_wt_rottype(rotor);
  chkpt_wt_max_am(max_angmom, basiskey);
  chkpt_wt_puream(puream, basiskey);
  chkpt_wt_nso(num_so, basiskey);
  chkpt_wt_nao(num_ao, basiskey);
  chkpt_wt_nshell(num_shells, basiskey);
  chkpt_wt_nirreps(nirreps);
  chkpt_wt_nprim(num_prims, basiskey);
  chkpt_wt_natom(num_atoms);
  chkpt_wt_nallatom(num_allatoms);
  chkpt_wt_nfzc(nfzc);
  chkpt_wt_nfzv(nfzv);
  chkpt_wt_nfragment(nfragments);
  if (num_atoms > 1) {
    chkpt_wt_rotconst(rotconst);
  }
  chkpt_wt_rot_symm_num(rot_symm_num);
  if (nfragments > 1) {
    chkpt_wt_natom_per_fragment(frag_num_atoms);
    chkpt_wt_nallatom_per_fragment(frag_num_allatoms);
    chkpt_wt_nref_per_fragment(nref_per_fragment);
    chkpt_wt_fragment_coeff(ref_pts_lc);
  }
  free(calc_label);


  /*-----------------------------------
    Start writing data to the file
    -----------------------------------*/

  /* Nuclear charges */
  chkpt_wt_zvals(nuclear_charges);

  /* Transformation table for atoms - just atom_orbit transposed */
  ict = init_int_array(num_atoms);
  ict_tmp = init_int_matrix(nirreps, num_atoms);
  for(i=0;i<nirreps;i++) {
    for(j=0;j<num_atoms;j++) {
      ict[j] = atom_orbit[j][i]+1;
      ict_tmp[i][j] = ict[j];
    }
  }
  chkpt_wt_ict(ict_tmp);
  free_int_matrix(ict_tmp);
  free(ict);

  /* Exponents of primitive gaussians */
  chkpt_wt_exps(exponents, basiskey);

  /* Contraction coefficients */
  /*------This piece of code is for segmented contractions ONLY------*/
  cspd = init_array(num_prims*MAXANGMOM);
  for(j=0;j<num_shells;j++)
    for(k=0;k<nprim_in_shell[j];k++)
      /*---
	Pitzer normalization of Psi 2 is NOT used - cc's for d-functions used to be
	multiplied by sqrt(3), f - by sqrt(15), g - sqrt(105), etc
	---*/
      cspd[shell_ang_mom[j]*num_prims+first_prim_shell[j]+k] = contr_coeff[first_prim_shell[j]+k];
  chkpt_wt_contr(cspd, basiskey);
  free(cspd);

  /* Pointer to primitives for a shell */
  arr_int = init_int_array(num_shells);
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_prim_shell[i]+1;
  chkpt_wt_sprim(arr_int, basiskey);

  /* Atom on which nth shell is centered */
  for(i=0;i<num_shells;i++)
    arr_int[i] = shell_nucleus[i]+1;
  chkpt_wt_snuc(arr_int, basiskey);

  /* Angular momentum of a shell */
  for(i=0;i<num_shells;i++)
    arr_int[i] = shell_ang_mom[i]+1;
  chkpt_wt_stype(arr_int, basiskey);

  /* Number of contracted functions (primitives) in a shell */
  chkpt_wt_snumg(nprim_in_shell, basiskey);

  /* Pointer to the first AO in shell */
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_ao_shell[i]+1;
  chkpt_wt_sloc(arr_int, basiskey);
  free(arr_int);

  /* Labels of irreps */
  chkpt_wt_irr_labs(irr_labels);

  /* Transformation matrices for coordinates (or p-functions) */
  mat_double = block_matrix(nirreps,9);
  for(symop=0;symop<nirreps;symop++) {
    mat_double[symop][0] = ao_type_transmat[1][symop][0];
    mat_double[symop][4] = ao_type_transmat[1][symop][1];
    mat_double[symop][8] = ao_type_transmat[1][symop][2];
  }
  chkpt_wt_cartrep(mat_double);
  free(mat_double);

  /* Transformation matrix for shells */
  shell_transm = init_int_matrix(num_shells,nirreps);;
  for(atom=0;atom<num_atoms;atom++)
    for(symop=0;symop<nirreps;symop++) {
      eq_atom = atom_orbit[atom][symop];
      for(i=0;i<nshells_per_atom[atom];i++)
	shell_transm[first_shell_on_atom[atom]+i][symop] = first_shell_on_atom[eq_atom] + i + 1;
    }
  chkpt_wt_shell_transm(shell_transm, basiskey);
  free_int_matrix(shell_transm);

  /* Labels of atoms including dummy atoms */
  tmp_atom_label = (char **) malloc(num_allatoms*sizeof(char *));
  for(atom=0; atom<num_allatoms; atom++) {
    tmp_atom_label[atom] = init_char_array(MAX_ELEMNAME);

    if(strlen(full_element[atom]) > MAX_ELEMNAME)
      punt("Element name exceeds limit, MAX_ELEMNAME.\n");

    strncpy(tmp_atom_label[atom],full_element[atom],strlen(full_element[atom])+1);
  }
  chkpt_wt_felement(tmp_atom_label);
  for(atom=0; atom<num_allatoms; atom++) {
    free(tmp_atom_label[atom]);
  }
  free(tmp_atom_label);

  /* Orbitals per irrep */
  chkpt_wt_sopi(num_so_per_irrep, basiskey);

  /* Symmetry label */
  chkpt_wt_sym_label(symmetry);

  /* Symmetry positions of atoms - for more info see count_uniques.c */
  chkpt_wt_atom_position(atom_position);

  /* Unique shell number to full shell number mapping array */
  arr_int = init_int_array(num_unique_shells);
  us = 0;
  for(ua=0;ua<num_uniques;ua++) {
    atom = u2a[ua];
    shell = first_shell_on_atom[atom];
    for(i=0;i<nshells_per_atom[atom];i++,shell++,us++)
      arr_int[us] = shell;
  }
  chkpt_wt_us2s(arr_int, basiskey);
  free(arr_int);

  /* SO to AO transformation matrix */
  chkpt_wt_usotao(usotao, basiskey);

  /* SO to basis functions transformation matrix */
  if (puream) {
    chkpt_wt_usotbf(usotbf, basiskey);
  }

  /* Pointers to first basis functions from shells */
  arr_int = init_int_array(num_shells);
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_basisfn_shell[i]+1;
  chkpt_wt_sloc_new(arr_int, basiskey);
  free(arr_int);

  /* Unique atom number to full atom number mapping array */
  chkpt_wt_ua2a(u2a);

  /* Mapping between canonical Cotton ordering of symmetry operations
     in the point group to the symmetry.h-defined ordering */
  chkpt_wt_symoper(sym_oper);

  /* write z_mat if it exists, see global.h for info about z_entry structure */
  if(!cartOn) {
    chkpt_wt_zmat(z_geom);
  }

  /* Number of shells in each angmom block */
  chkpt_wt_shells_per_am(shells_per_am, basiskey);

  /* Mapping array from the am-blocked to the canonical (in the order of
     appearance) ordering of shells */
  chkpt_wt_am2canon_shell_order(am2canon_shell_order, basiskey);

  /* Matrix representation of rotation back to the reference frame */
  chkpt_wt_rref(Rref);

  /* write full_geom, cartesian geometry with dummy atoms included */
  chkpt_wt_fgeom(full_geom);

  /* write array of flags that indicate whether atoms in full_geom are dummy or not */
  chkpt_wt_atom_dummy(atom_dummy);

  /* nuclear repulsion energy */
  chkpt_wt_enuc(repulsion);

  /* Cartesian displacement SALCs */
  chkpt_wt_cdsalc2cd(const_cast<const double**>(cdsalc2cd));

  /* cartdisp SALCs per irrep */
  chkpt_wt_cdsalcpi(cdsalc_pi);

  chkpt_close();
  return;
}

}} // namespace psi::input
