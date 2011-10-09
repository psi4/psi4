/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

#include <libmints/wavefunction.h>

namespace psi { namespace ccsort {

/* get_moinfo(): Routine to obtain basic orbital information from
** chkpt and compute the associated lookup arrays.
**
** Several loookup arrays are written to CC_INFO at the end of this
** routine for use in other CC programs.  I do this for two reasons:
** (1) I don't want to have to use redundant code for re-computing the
** arrays and (2) I don't want to allow certain types of user input
** (e.g. orbital occupations) to change between CC programs.
**
** T. Daniel Crawford, October 1996
** Modified by TDC March, 1999 */

void get_moinfo(void)
{
  int i, j, h, count, ocount, vcount;
  int active_count, all_count;
  int active_a_count, active_b_count, all_a_count, all_b_count;
  int fc_offset, offset, cl_offset, op_offset, vr_offset, fv_offset;
  int nfzc, nuocc, nopen, nclsd;
  double escf;
  double ***evects, ***scf_vector;
  double ***evects_A, ***scf_vector_A;
  double ***evects_B, ***scf_vector_B;
  double ***C;
  int *pitz2qt, *qt2pitz, J, qt_j, pitz_j, pitz_J, *pitz_offset;
  int *pitz2qt_A, *qt2pitz_A, *pitz2qt_B, *qt2pitz_B;
  psio_address next;

  chkpt_init(PSIO_OPEN_OLD);
  moinfo.nirreps = chkpt_rd_nirreps();
  moinfo.nmo = chkpt_rd_nmo();
  moinfo.nso = chkpt_rd_nso();
  moinfo.nao = chkpt_rd_nao();
  moinfo.iopen = chkpt_rd_iopen();
  moinfo.labels = chkpt_rd_irr_labs();
  moinfo.enuc = chkpt_rd_enuc();
  escf = chkpt_rd_escf();
  moinfo.sopi = chkpt_rd_sopi();
  moinfo.orbspi = chkpt_rd_orbspi();
  moinfo.clsdpi = chkpt_rd_clsdpi();
  moinfo.openpi = chkpt_rd_openpi();
  // TODO: Uncomment the following line.
  //moinfo.usotao = chkpt_rd_usotao();
  if(params.ref == 0) moinfo.scf = chkpt_rd_scf();

  /* Dump the reference wave function ID to CC_INFO */
  psio_write_entry(CC_INFO, "Reference Wavefunction",
           (char *) &(params.ref), sizeof(int));

  moinfo.frdocc = Process::environment.reference_wavefunction()->frzcpi();
  moinfo.fruocc = Process::environment.reference_wavefunction()->frzvpi();
  chkpt_close();

  moinfo.nfzc = moinfo.nfzv = 0;
  for(i=0; i < moinfo.nirreps; i++) {
    moinfo.nfzc += moinfo.frdocc[i];
    moinfo.nfzv += moinfo.fruocc[i];
  }

  if(moinfo.nfzc) {
    chkpt_init(PSIO_OPEN_OLD);
    moinfo.efzc = chkpt_rd_efzc();
    chkpt_close();
  }
  else moinfo.efzc = 0.0;

  /* Calculation consistency check */
  if(moinfo.nfzc && (fabs(moinfo.efzc) < 1e-7)) {
    fprintf(outfile, "\tCCSORT Error: Orbitals are frozen in input,\n");
    fprintf(outfile, "\tbut frozen core energy is small!\n");
    fprintf(outfile, "\tCalculation will be aborted...\n");
    fflush(outfile);
    exit(PSI_RETURN_FAILURE);
  }
  else if(!moinfo.nfzc && fabs(moinfo.efzc)) {
    fprintf(outfile, "\tCCSORT Warning: No orbitals are frozen,\n");
    fprintf(outfile, "\tbut the frozen-core energy in chkpt is non-zero.\n");
    fprintf(outfile, "\tCalculation will continue with zero efzc...\n");
    fflush(outfile);
    moinfo.efzc = 0.0;
  }

  /* Dump the frozen orbital arrays to CC_INFO */
  psio_write_entry(CC_INFO, "Frozen Core Orbs Per Irrep",
           (char *) moinfo.frdocc, sizeof(int)*moinfo.nirreps);
  psio_write_entry(CC_INFO, "Frozen Virt Orbs Per Irrep",
           (char *) moinfo.fruocc, sizeof(int)*moinfo.nirreps);

  /* Get the Pitzer->QT reordering array and compute its inverse */
  if(params.ref == 2) { /* UHF */
    moinfo.pitz2qt_A = init_int_array(moinfo.nmo);
    moinfo.qt2pitz_A = init_int_array(moinfo.nmo);
    moinfo.pitz2qt_B = init_int_array(moinfo.nmo);
    moinfo.qt2pitz_B = init_int_array(moinfo.nmo);
    reorder_qt_uhf(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
           moinfo.pitz2qt_A, moinfo.pitz2qt_B, moinfo.orbspi, moinfo.nirreps);

    for(i=0; i < moinfo.nmo; i++) {
      moinfo.qt2pitz_A[moinfo.pitz2qt_A[i]] = i;
      moinfo.qt2pitz_B[moinfo.pitz2qt_B[i]] = i;
    }
  }
  else { /* RHF and ROHF */
    moinfo.pitz2qt = init_int_array(moinfo.nmo);
    moinfo.qt2pitz = init_int_array(moinfo.nmo);
    reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc,
           moinfo.fruocc, moinfo.pitz2qt,  moinfo.orbspi, moinfo.nirreps);
    for(i=0; i < moinfo.nmo; i++) moinfo.qt2pitz[moinfo.pitz2qt[i]] = i;
  }

  /* Compute spatial-orbital reordering arrays */
  moinfo.pitzer2qt = init_int_array(moinfo.nmo);
  moinfo.qt2pitzer = init_int_array(moinfo.nmo);
  reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
         moinfo.pitzer2qt, moinfo.orbspi, moinfo.nirreps);
  for(i=0; i < moinfo.nmo; i++) {
    j = moinfo.pitzer2qt[i];
    moinfo.qt2pitzer[j] = i;
  }

  /* Adjust clsdpi array for frozen orbitals */
  for(i=0; i < moinfo.nirreps; i++)
    moinfo.clsdpi[i] -= moinfo.frdocc[i];

  moinfo.uoccpi = init_int_array(moinfo.nirreps);
  for(i=0; i < moinfo.nirreps; i++)
    moinfo.uoccpi[i] = moinfo.orbspi[i] - moinfo.clsdpi[i] -
      moinfo.openpi[i] - moinfo.fruocc[i] -
      moinfo.frdocc[i];

  nclsd = nopen = nuocc = 0;
  for(i=0; i < moinfo.nirreps; i++) {
    nclsd += moinfo.clsdpi[i];
    nopen += moinfo.openpi[i];
    nuocc += moinfo.uoccpi[i];
  }
  nfzc = moinfo.nfzc;

  /* Number of non-frozen orbitals */
  moinfo.nactive = nclsd + nopen + nuocc;
  psio_write_entry(CC_INFO, "No. of Active Orbitals", (char *) &(moinfo.nactive),
                   sizeof(int));

  /* Number of occupied and virtual orbirals per irrep including
     open-shells */
  if(params.ref == 2) {
    moinfo.aoccpi = init_int_array(moinfo.nirreps);
    moinfo.boccpi = init_int_array(moinfo.nirreps);
    moinfo.avirtpi = init_int_array(moinfo.nirreps);
    moinfo.bvirtpi = init_int_array(moinfo.nirreps);
    moinfo.all_aoccpi = init_int_array(moinfo.nirreps);
    moinfo.all_boccpi = init_int_array(moinfo.nirreps);
    moinfo.all_avirtpi = init_int_array(moinfo.nirreps);
    moinfo.all_bvirtpi = init_int_array(moinfo.nirreps);
    for(h=0; h < moinfo.nirreps; h++) {
      moinfo.aoccpi[h] = moinfo.clsdpi[h] + moinfo.openpi[h];
      moinfo.boccpi[h] = moinfo.clsdpi[h];
      moinfo.avirtpi[h] = moinfo.uoccpi[h];
      moinfo.bvirtpi[h] = moinfo.uoccpi[h] + moinfo.openpi[h];

      moinfo.all_aoccpi[h] = moinfo.frdocc[h] + moinfo.clsdpi[h] + moinfo.openpi[h];
      moinfo.all_boccpi[h] = moinfo.frdocc[h] + moinfo.clsdpi[h];
      moinfo.all_avirtpi[h] = moinfo.fruocc[h] + moinfo.uoccpi[h];
      moinfo.all_bvirtpi[h] = moinfo.fruocc[h] + moinfo.uoccpi[h] + moinfo.openpi[h];
    }
    /* Dump occpi and virtpi arrays to CC_INFO */
    psio_write_entry(CC_INFO, "Active Alpha Occ Orbs Per Irrep",
             (char *) moinfo.aoccpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "Active Beta Occ Orbs Per Irrep",
             (char *) moinfo.boccpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "Active Alpha Virt Orbs Per Irrep",
             (char *) moinfo.avirtpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "Active Beta Virt Orbs Per Irrep",
             (char *) moinfo.bvirtpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Alpha Occ Orbs Per Irrep",
             (char *) moinfo.all_aoccpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Beta Occ Orbs Per Irrep",
             (char *) moinfo.all_boccpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Alpha Virt Orbs Per Irrep",
             (char *) moinfo.all_avirtpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Beta Virt Orbs Per Irrep",
             (char *) moinfo.all_bvirtpi, sizeof(int)*moinfo.nirreps);

  }
  else {
    moinfo.occpi = init_int_array(moinfo.nirreps);
    moinfo.virtpi = init_int_array(moinfo.nirreps);
    moinfo.all_occpi = init_int_array(moinfo.nirreps);
    moinfo.all_virtpi = init_int_array(moinfo.nirreps);
    for(h=0; h < moinfo.nirreps; h++) {
      moinfo.occpi[h] = moinfo.clsdpi[h] + moinfo.openpi[h];
      moinfo.virtpi[h] = moinfo.uoccpi[h] + moinfo.openpi[h];
      moinfo.all_occpi[h] = moinfo.frdocc[h] + moinfo.clsdpi[h] +
    moinfo.openpi[h];
      moinfo.all_virtpi[h] = moinfo.fruocc[h] + moinfo.uoccpi[h] +
    moinfo.openpi[h];
    }
    /* Dump occpi and virtpi arrays to CC_INFO */
    psio_write_entry(CC_INFO, "Active Occ Orbs Per Irrep",
             (char *) moinfo.occpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "Active Virt Orbs Per Irrep",
             (char *) moinfo.virtpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Occ Orbs Per Irrep",
             (char *) moinfo.all_occpi, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Virt Orbs Per Irrep",
             (char *) moinfo.all_virtpi, sizeof(int)*moinfo.nirreps);

  }

  /* Build the boolean arrays for the orbital classification routines.
     The occ and vir arrays identify active occupied or virtual
     orbitals (including open shells).  The socc array identifies only
     the singly occupied orbitals.  The all_occ and all_vir arrays
     identify any occupied or virtual orbitals (include open shells
     and frozen orbitals).  The all_socc array identifies only singly
     occupied orbitals within the full list of orbitals.  The frozen
     array identifies only frozen orbitals.  The argument for all
     arrays much be a QT-ordered index. Most of these are used in the
     integral sorting routines found in CCSORT and CCDENSITY. */

  if(params.ref == 2) {
    moinfo.aocc = init_int_array(moinfo.nactive);
    moinfo.bocc = init_int_array(moinfo.nactive);
    moinfo.avir = init_int_array(moinfo.nactive);
    moinfo.bvir = init_int_array(moinfo.nactive);
    moinfo.all_aocc = init_int_array(moinfo.nmo);
    moinfo.all_bocc = init_int_array(moinfo.nmo);
    moinfo.all_avir = init_int_array(moinfo.nmo);
    moinfo.all_bvir = init_int_array(moinfo.nmo);
    moinfo.frozen = init_int_array(moinfo.nmo);

    active_a_count = active_b_count = all_a_count = all_b_count = 0;
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.frdocc[i]; j++, all_a_count++, all_b_count++) {
    moinfo.all_aocc[all_a_count] = 1;
    moinfo.all_bocc[all_b_count] = 1;
    moinfo.frozen[all_a_count] = 1;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.clsdpi[i] + moinfo.openpi[i]; j++) {
    moinfo.aocc[active_a_count] = 1;
    moinfo.all_aocc[all_a_count] = 1;
    active_a_count++;
    all_a_count++;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.clsdpi[i]; j++) {
    moinfo.bocc[active_b_count] = 1;
    moinfo.all_bocc[all_b_count] = 1;
    active_b_count++;
    all_b_count++;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.openpi[i] + moinfo.uoccpi[i]; j++) {
    moinfo.bvir[active_b_count] = 1;
    moinfo.all_bvir[all_b_count] = 1;
    active_b_count++;
    all_b_count++;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.uoccpi[i]; j++) {
    moinfo.avir[active_a_count] = 1;
    moinfo.all_avir[all_a_count] = 1;
    active_a_count++;
    all_a_count++;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.fruocc[i]; j++, all_a_count++,all_b_count++) {
    moinfo.all_avir[all_a_count] = 1;
    moinfo.all_bvir[all_b_count] = 1;
    moinfo.frozen[all_a_count] = 1;
      }
    }

    /* Dump the Boolean arrays to CC_INFO */
    psio_write_entry(CC_INFO, "Active Alpha Occ Orbital Boolean",
             (char *) moinfo.aocc, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Beta Occ Orbital Boolean",
             (char *) moinfo.bocc, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Alpha Virt Orbital Boolean",
             (char *) moinfo.avir, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Beta Virt Orbital Boolean",
             (char *) moinfo.bvir, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "All Alpha Occ Orbital Boolean",
             (char *) moinfo.all_aocc, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "All Beta Occ Orbital Boolean",
             (char *) moinfo.all_bocc, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "All Alpha Virt Orbital Boolean",
             (char *) moinfo.all_avir, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "All Beta Virt Orbital Boolean",
             (char *) moinfo.all_bvir, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "Frozen Orbital Boolean",
             (char *) moinfo.frozen, sizeof(int)*moinfo.nmo);
  }
  else {
    moinfo.occ = init_int_array(moinfo.nactive);
    moinfo.vir = init_int_array(moinfo.nactive);
    moinfo.socc = init_int_array(moinfo.nactive);
    moinfo.all_occ = init_int_array(moinfo.nmo);
    moinfo.all_vir = init_int_array(moinfo.nmo);
    moinfo.all_socc = init_int_array(moinfo.nmo);
    moinfo.frozen = init_int_array(moinfo.nmo);
    moinfo.orbsym = init_int_array(moinfo.nmo);

    active_count = all_count = 0;
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.frdocc[i]; j++, all_count++) {
    moinfo.all_occ[all_count] = 1;
    moinfo.frozen[all_count] = 1;
    moinfo.orbsym[all_count] = i;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.clsdpi[i]; j++, active_count++, all_count++) {
    moinfo.occ[active_count] = 1;
    moinfo.all_occ[all_count] = 1;
    moinfo.orbsym[all_count] = i;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.openpi[i]; j++, active_count++, all_count++) {
    moinfo.occ[active_count] = 1;
    moinfo.vir[active_count] = 1;
    moinfo.socc[active_count] = 1;
    moinfo.all_occ[all_count] = 1;
    moinfo.all_vir[all_count] = 1;
    moinfo.all_socc[all_count] = 1;
    moinfo.orbsym[all_count] = i;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.uoccpi[i]; j++, active_count++, all_count++) {
    moinfo.vir[active_count] = 1;
    moinfo.all_vir[all_count] = 1;
    moinfo.orbsym[all_count] = i;
      }
    }
    for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.fruocc[i]; j++, all_count++) {
    moinfo.all_vir[all_count] = 1;
    moinfo.frozen[all_count] = 1;
    moinfo.orbsym[all_count] = i;
      }
    }

    /* Dump the Boolean arrays to CC_INFO */
    psio_write_entry(CC_INFO, "Active Occ Orbital Boolean",
             (char *) moinfo.occ, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Virt Orbital Boolean",
             (char *) moinfo.vir, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Socc Orbital Boolean",
             (char *) moinfo.socc, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "All Occ Orbital Boolean",
             (char *) moinfo.all_occ, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "All Virt Orbital Boolean",
             (char *) moinfo.all_vir, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "All Socc Orbital Boolean",
             (char *) moinfo.all_socc, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "Frozen Orbital Boolean",
             (char *) moinfo.frozen, sizeof(int)*moinfo.nmo);
  }

  if(params.ref == 2) {

    /**** UHF references ****/

    /*

    Build relative index arrays which map a QT Standard orbital index
    into the CC occupied or virtual relative index.  Since the QT
    Standard ordering is different for alpha and beta UHF indices, the
    CC occupied and virtual orderings are also necessarily different
    for the different spins.  The cc_aocc, cc_bocc, cc_avir, and
    cc_bvir arrays do not include frozen orbitals and are of length
    moinfo.nactive.  The cc_allaocc, cc_allbocc, cc_allavir,
    cc_allbvir do include frozen orbitals and are of length
    moinfo.nmo.  In these latter arrays, the frozen orbitals are given
    FIRST, even in the cc_allavir and cc_allbvir orderings.  (This is
    just like the RHF/ROHF case given below.)  These arrays are needed
    by classification routines used in sorting two-electron integrals
    (in CCSORT and CCDENSITY) where integrals are divided up by their
    occupied and virtual index patterns (see classify.c and
    distribute.c).


    For example, for 3B1 CH2 in a DZ basis with one frozen core
    orbital and three frozen virtual orbitals:

    QT Alpha:    0   1   2   3   4   5   6   7   8   9  10  11  12  13
    Space:      fc   d   s   d   s   u   u   u   u   u   u  fv  fv  fv
    Symm:       a1  a1  a1  b2  b1  a1  a1  a1  b1  b2  b2  a1  a1  b2
    ------------------------------------------------------------------
    cc_aocc:     X   0   1   2   3  -1  -1  -1  -1  -1  -1   X   X   X
    cc_avir:     X  -1  -1  -1  -1   0   1   2   3   4   5   X   X   X
    cc_allaocc:  0   1   2   3   4  -1  -1  -1  -1  -1  -1  -1  -1  -1
    cc_allavir: -1  -1  -1  -1  -1   2   3   4   5   7   8   0   1   6

    QT Beta:     0   1   2   3   4   5   6   7   8   9  10  11  12  13
    Space:      fc   d   d   s   u   u   u   s   u   u   u  fv  fv  fv
    Symm:       a1  a1  b2  a1  a1  a1  a1  b1  b1  b2  b2  a1  a1  b2
    ------------------------------------------------------------------
    cc_bocc:     X   0   1  -1  -1  -1  -1  -1  -1  -1  -1   X   X   X
    cc_bvir:     X  -1  -1   0   1   2   3   4   5   6   7   X   X   X
    cc_allbocc:  0   1   2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
    cc_allbvir: -1  -1  -1   2   3   4   5   6   7   9  10   0   1   8

    A "-1" is used to mark QT indices irrelevant to the given CC
    ordering array, and an "X" indicates that the element doesn't
    appear in the CC array at all (i.e., cc_occ and cc_vir are of
    length nactive and cc_allocc and cc_allvir are of length nmo).

    Note that I also compute the CC->QT versions of the lookup arrays
    (qt_aocc, qt_bocc, and qt_avir, and qt_bvir) simultaneously.
    Hence, if one has an array in CC ordering without frozen orbitals,
    and needs to convert the array to CC ordering WITH frozen orbitals
    (i.e. add zeroes in all the right places), then one could loop
    over the CC active orbitals, compute the QT index for the current
    orbital, and finally compute the index for the full CC list.

    */

    /* First the active-orbital CC ordering and symmetry arrays */
    moinfo.cc_aocc = init_int_array(moinfo.nactive);
    moinfo.cc_bocc = init_int_array(moinfo.nactive);
    moinfo.cc_avir = init_int_array(moinfo.nactive);
    moinfo.cc_bvir = init_int_array(moinfo.nactive);
    moinfo.qt_aocc = init_int_array(moinfo.nactive);
    moinfo.qt_bocc = init_int_array(moinfo.nactive);
    moinfo.qt_avir = init_int_array(moinfo.nactive);
    moinfo.qt_bvir = init_int_array(moinfo.nactive);
    for(i=0; i < moinfo.nactive; i++) {
      moinfo.cc_aocc[i] = moinfo.cc_avir[i] = -1;
      moinfo.cc_bocc[i] = moinfo.cc_bvir[i] = -1;
      moinfo.qt_aocc[i] = moinfo.qt_avir[i] = -1;
      moinfo.qt_bocc[i] = moinfo.qt_bvir[i] = -1;
    }

    moinfo.aocc_sym = init_int_array(moinfo.nactive);
    moinfo.bocc_sym = init_int_array(moinfo.nactive);
    moinfo.avir_sym = init_int_array(moinfo.nactive);
    moinfo.bvir_sym = init_int_array(moinfo.nactive);

    count = 0;
    offset = 0;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) offset += moinfo.clsdpi[h-1] + moinfo.openpi[h-1];
      for(i=0; i < moinfo.clsdpi[h] + moinfo.openpi[h]; i++,count++) {
    moinfo.cc_aocc[offset+i] = count;
    moinfo.qt_aocc[count] = nfzc + offset + i;
    moinfo.aocc_sym[count] = h;
      }
    }

    count = 0;
    offset = 0;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) offset += moinfo.clsdpi[h-1];
      for(i=0; i < moinfo.clsdpi[h]; i++,count++) {
    moinfo.cc_bocc[offset+i] = count;
    moinfo.qt_bocc[count] = nfzc + offset + i;
    moinfo.bocc_sym[count] = h;
      }
    }

    count = 0;
    offset = nclsd + nopen;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) offset += moinfo.uoccpi[h-1];
      for(i=0; i < moinfo.uoccpi[h]; i++,count++) {
    moinfo.cc_avir[offset+i] = count;
    moinfo.qt_avir[count] = nfzc + offset + i;
    moinfo.avir_sym[count] = h;
      }
    }

    count = 0;
    offset = nclsd;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) offset += moinfo.uoccpi[h-1] + moinfo.openpi[h-1];
      for(i=0; i < moinfo.uoccpi[h] + moinfo.openpi[h]; i++,count++) {
    moinfo.cc_bvir[offset+i] = count;
    moinfo.qt_bvir[count] = nfzc + offset + i;
    moinfo.bvir_sym[count] = h;
      }
    }

    /* Dump the active-orbital CC and QT ordering and symmetry arrays to CC_INFO */
    psio_write_entry(CC_INFO, "Active Alpha Occ Orb Symmetry",
             (char *) moinfo.aocc_sym, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Beta Occ Orb Symmetry",
             (char *) moinfo.bocc_sym, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Alpha Virt Orb Symmetry",
             (char *) moinfo.avir_sym, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Beta Virt Orb Symmetry",
             (char *) moinfo.bvir_sym, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "QT->CC Alpha Active Occ Order",
             (char *) moinfo.cc_aocc, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "QT->CC Beta Active Occ Order",
             (char *) moinfo.cc_bocc, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "QT->CC Alpha Active Virt Order",
             (char *) moinfo.cc_avir, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "QT->CC Beta Active Virt Order",
             (char *) moinfo.cc_bvir, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "CC->QT Alpha Active Occ Order",
             (char *) moinfo.qt_aocc, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "CC->QT Beta Active Occ Order",
             (char *) moinfo.qt_bocc, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "CC->QT Alpha Active Virt Order",
             (char *) moinfo.qt_avir, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "CC->QT Beta Active Virt Order",
             (char *) moinfo.qt_bvir, sizeof(int)*moinfo.nactive);

    /* Now the all-orbital CC ordering and symmetry arrays */
    moinfo.cc_allaocc = init_int_array(moinfo.nmo);
    moinfo.cc_allbocc = init_int_array(moinfo.nmo);
    moinfo.cc_allavir = init_int_array(moinfo.nmo);
    moinfo.cc_allbvir = init_int_array(moinfo.nmo);
    moinfo.qt_allaocc = init_int_array(moinfo.nmo);
    moinfo.qt_allbocc = init_int_array(moinfo.nmo);
    moinfo.qt_allavir = init_int_array(moinfo.nmo);
    moinfo.qt_allbvir = init_int_array(moinfo.nmo);
    for(i=0; i < moinfo.nmo; i++) {
      moinfo.cc_allaocc[i] = moinfo.cc_allavir[i] = -1;
      moinfo.cc_allbocc[i] = moinfo.cc_allbvir[i] = -1;
      moinfo.qt_allaocc[i] = moinfo.qt_allavir[i] = -1;
      moinfo.qt_allbocc[i] = moinfo.qt_allbvir[i] = -1;
    }

    moinfo.allaocc_sym = init_int_array(moinfo.nmo);
    moinfo.allbocc_sym = init_int_array(moinfo.nmo);
    moinfo.allavir_sym = init_int_array(moinfo.nmo);
    moinfo.allbvir_sym = init_int_array(moinfo.nmo);

    count = 0;
    fc_offset = 0;
    offset = nfzc;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) fc_offset += moinfo.frdocc[h-1];
      for(i=0; i < moinfo.frdocc[h]; i++,count++) {
    moinfo.cc_allaocc[fc_offset+i] = count;
    moinfo.qt_allaocc[count] = fc_offset+i;
    moinfo.allaocc_sym[count] = h;
      }
      if(h) offset += moinfo.clsdpi[h-1] + moinfo.openpi[h-1];
      for(i=0; i < moinfo.clsdpi[h] + moinfo.openpi[h]; i++,count++) {
    moinfo.cc_allaocc[offset+i] = count;
    moinfo.qt_allaocc[count] = offset+i;
    moinfo.allaocc_sym[count] = h;
      }
    }

    count = 0;
    fc_offset = 0;
    offset = nfzc;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) fc_offset += moinfo.frdocc[h-1];
      for(i=0; i < moinfo.frdocc[h]; i++,count++) {
    moinfo.cc_allbocc[fc_offset+i] = count;
    moinfo.qt_allbocc[count] = fc_offset+i;
    moinfo.allbocc_sym[count] = h;
      }
      if(h) offset += moinfo.clsdpi[h-1];
      for(i=0; i < moinfo.clsdpi[h]; i++,count++) {
    moinfo.cc_allbocc[offset+i] = count;
    moinfo.qt_allbocc[count] = offset+i;
    moinfo.allbocc_sym[count] = h;
      }
    }

    count = 0;
    fv_offset = nfzc + nclsd + nopen + nuocc;
    offset = nfzc + nclsd + nopen;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) fv_offset += moinfo.fruocc[h-1];
      for(i=0; i < moinfo.fruocc[h]; i++,count++) {
    moinfo.cc_allavir[fv_offset+i] = count;
    moinfo.qt_allavir[count] = fv_offset+i;
    moinfo.allavir_sym[count] = h;
      }
      if(h) offset += moinfo.uoccpi[h-1];
      for(i=0; i < moinfo.uoccpi[h]; i++,count++) {
    moinfo.cc_allavir[offset+i] = count;
    moinfo.qt_allavir[count] = offset+i;
    moinfo.allavir_sym[count] = h;
      }
    }

    count = 0;
    fv_offset = nfzc + nclsd + nopen + nuocc;
    offset = nfzc + nclsd;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) fv_offset += moinfo.fruocc[h-1];
      for(i=0; i < moinfo.fruocc[h]; i++,count++) {
    moinfo.cc_allbvir[fv_offset+i] = count;
    moinfo.qt_allbvir[count] = fv_offset+i;
    moinfo.allbvir_sym[count] = h;
      }
      if(h) offset += moinfo.uoccpi[h-1] + moinfo.openpi[h-1];
      for(i=0; i < moinfo.uoccpi[h] + moinfo.openpi[h]; i++,count++) {
    moinfo.cc_allbvir[offset+i] = count;
    moinfo.qt_allbvir[count] = offset+i;
    moinfo.allbvir_sym[count] = h;
      }
    }

    /* Dump the all-orbital CC ordering and symmetry arrays to CC_INFO */
    psio_write_entry(CC_INFO, "All Alpha Occ Orb Symmetry",
             (char *) moinfo.allaocc_sym, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "All Beta Occ Orb Symmetry",
             (char *) moinfo.allbocc_sym, sizeof(int)*moinfo.nmo);

    psio_write_entry(CC_INFO, "All Alpha Virt Orb Symmetry",
             (char *) moinfo.allavir_sym, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "All Beta Virt Orb Symmetry",
             (char *) moinfo.allbvir_sym, sizeof(int)*moinfo.nmo);

    psio_write_entry(CC_INFO, "QT->CC All Alpha Occ Order",
             (char *) moinfo.cc_allaocc, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "QT->CC All Beta Occ Order",
             (char *) moinfo.cc_allbocc, sizeof(int)*moinfo.nmo);

    psio_write_entry(CC_INFO, "QT->CC All Alpha Virt Order",
             (char *) moinfo.cc_allavir, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "QT->CC All Beta Virt Order",
             (char *) moinfo.cc_allbvir, sizeof(int)*moinfo.nmo);

    psio_write_entry(CC_INFO, "CC->QT All Alpha Occ Order",
             (char *) moinfo.qt_allaocc, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "CC->QT All Beta Occ Order",
             (char *) moinfo.qt_allbocc, sizeof(int)*moinfo.nmo);

    psio_write_entry(CC_INFO, "CC->QT All Alpha Virt Order",
             (char *) moinfo.qt_allavir, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "CC->QT All Beta Virt Order",
             (char *) moinfo.qt_allbvir, sizeof(int)*moinfo.nmo);

  }  /*** UHF references ***/
  else {

    /**** RHF and ROHF references ****/

    /*

    Build relative index arrays which map a QT Standard orbital index
    into the CC occupied or virtual relative index.  In the CC
    occupied ordering (given by "cc_allocc"), the orbitals are
    organized with all frozen core orbitals first, followed by doubly
    occupied orbitals, followed by singly occupied orbitals.  In the
    CC virtual ordering (given by "cc_allvir"), the orbitals are
    organized with all frozen virtual orbitals first, followed by
    active unoccupied orbitals, followed by the singly occupied
    orbitals.  The "cc_occ" and "cc_vir" arrays are defined similarly,
    but without frozen orbitals.  These arrays are needed by
    classification routines used in sorting two-electron integrals
    (in CCSORT and CCDENSITY) where integrals are divided up by their
    occupied and virtual index patterns (see classify.c and
    distribute.c).

    For example, for 3B1 CH2 in a DZ basis with one frozen core
    orbital and three frozen virtual orbitals:

    QT:         0   1   2   3   4   5   6   7   8   9  10  11  12  13
    Space:     fc   d   d   s   s   u   u   u   u   u   u  fv  fv  fv
    Symm:      a1  a1  b2  a1  b1  a1  a1  a1  b1  b2  b2  a1  a1  b2
    -----------------------------------------------------------------
    cc_occ:     X   0   3   1   2  -1  -1  -1  -1  -1  -1   X   X   X
    cc_vir:     X  -1  -1   3   5   0   1   2   4   6   7   X   X   X
    cc_allocc:  0   1   4   2   3  -1  -1  -1  -1  -1  -1  -1  -1  -1
    cc_allvir: -1  -1  -1   5   7   2   3   4   6   9  10   0   1   8

    A "-1" is used to mark QT indices irrelevant to the given CC
    ordering array, and an "X" indicates that the element doesn't
    appear in the CC array at all (i.e., cc_occ and cc_vir are of
    length nactive and cc_allocc and cc_allvir are of length nmo).

    Note that I also compute the CC->QT versions of the cc_occ and
    cc_vir arrays (qt_occ and qt_vir, respectively) simultaneously.
    Hence, if one has an array in CC ordering without frozen orbitals,
    and needs to convert the array to CC ordering WITH frozen orbitals
    (i.e. add zeroes in all the right places), then one could loop
    over the CC active orbitals, compute the QT index for the current
    orbital (using qt_occ or qt_vir), and finally compute the index
    for the full CC list (using cc_allocc or cc_allvir).

    */

    /* First the active-orbital CC ordering and symmetry arrays */
    moinfo.cc_occ = init_int_array(moinfo.nactive);
    moinfo.cc_vir = init_int_array(moinfo.nactive);
    moinfo.qt_occ = init_int_array(moinfo.nactive);
    moinfo.qt_vir = init_int_array(moinfo.nactive);
    for(i=0; i < moinfo.nactive; i++) {
      moinfo.cc_occ[i] = moinfo.cc_vir[i] = -1;
      moinfo.qt_occ[i] = moinfo.qt_vir[i] = -1;
    }

    moinfo.occ_sym = init_int_array(moinfo.nactive);
    moinfo.vir_sym = init_int_array(moinfo.nactive);

    count = 0;
    cl_offset = 0;
    op_offset = nclsd;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) cl_offset += moinfo.clsdpi[h-1];
      for(i=0; i < moinfo.clsdpi[h]; i++,count++) {
    moinfo.cc_occ[cl_offset+i] = count;
    moinfo.qt_occ[count] = nfzc + cl_offset+i;
    moinfo.occ_sym[count] = h;
      }
      if(h) op_offset += moinfo.openpi[h-1];
      for(i=0; i < moinfo.openpi[h]; i++,count++) {
    moinfo.cc_occ[op_offset+i] = count;
    moinfo.qt_occ[count] = nfzc + op_offset+i;
    moinfo.occ_sym[count] = h;
      }
    }

    count = 0;
    vr_offset = nclsd + nopen;
    op_offset = nclsd;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) vr_offset += moinfo.uoccpi[h-1];
      for(i=0; i < moinfo.uoccpi[h]; i++,count++) {
    moinfo.cc_vir[vr_offset+i] = count;
    moinfo.qt_vir[count] = nfzc + vr_offset+i;
    moinfo.vir_sym[count] = h;
      }
      if(h) op_offset += moinfo.openpi[h-1];
      for(i=0; i < moinfo.openpi[h]; i++,count++) {
    moinfo.cc_vir[op_offset+i] = count;
    moinfo.qt_vir[count] = nfzc + op_offset+i;
    moinfo.vir_sym[count] = h;
      }
    }

    /* Dump the active-orbital CC and QT ordering and symmetry arrays to CC_INFO */
    psio_write_entry(CC_INFO, "Active Occ Orb Symmetry",
             (char *) moinfo.occ_sym, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "Active Virt Orb Symmetry",
             (char *) moinfo.vir_sym, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "QT->CC Active Occ Order",
             (char *) moinfo.cc_occ, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "QT->CC Active Virt Order",
             (char *) moinfo.cc_vir, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "CC->QT Active Occ Order",
             (char *) moinfo.qt_occ, sizeof(int)*moinfo.nactive);
    psio_write_entry(CC_INFO, "CC->QT Active Virt Order",
             (char *) moinfo.qt_vir, sizeof(int)*moinfo.nactive);

    /* Now the all-orbital CC ordering and symmetry arrays */
    moinfo.cc_allocc = init_int_array(moinfo.nmo);
    moinfo.cc_allvir = init_int_array(moinfo.nmo);
    moinfo.qt_allocc = init_int_array(moinfo.nmo);
    moinfo.qt_allvir = init_int_array(moinfo.nmo);
    for(i=0; i < moinfo.nmo; i++) {
      moinfo.cc_allocc[i] = moinfo.cc_allvir[i] = -1;
      moinfo.qt_allocc[i] = moinfo.qt_allvir[i] = -1;
    }

    moinfo.allocc_sym = init_int_array(moinfo.nmo);
    moinfo.allvir_sym = init_int_array(moinfo.nmo);

    count = 0;
    fc_offset = 0;
    cl_offset = nfzc;
    op_offset = nfzc + nclsd;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) fc_offset += moinfo.frdocc[h-1];
      for(i=0; i < moinfo.frdocc[h]; i++,count++) {
    moinfo.cc_allocc[fc_offset+i] = count;
    moinfo.qt_allocc[count] = fc_offset+i;
    moinfo.allocc_sym[count] = h;
      }
      if(h) cl_offset += moinfo.clsdpi[h-1];
      for(i=0; i < moinfo.clsdpi[h]; i++,count++) {
    moinfo.cc_allocc[cl_offset+i] = count;
    moinfo.qt_allocc[count] = cl_offset+i;
    moinfo.allocc_sym[count] = h;
      }
      if(h) op_offset += moinfo.openpi[h-1];
      for(i=0; i < moinfo.openpi[h]; i++,count++) {
    moinfo.cc_allocc[op_offset+i] = count;
    moinfo.qt_allocc[count] = op_offset+i;
    moinfo.allocc_sym[count] = h;
      }
    }

    count = 0;
    fv_offset = nfzc + nclsd + nopen + nuocc;
    vr_offset = nfzc + nclsd + nopen;
    op_offset = nfzc + nclsd;
    for(h=0; h < moinfo.nirreps; h++) {
      if(h) fv_offset += moinfo.fruocc[h-1];
      for(i=0; i < moinfo.fruocc[h]; i++,count++) {
    moinfo.cc_allvir[fv_offset+i] = count;
    moinfo.qt_allvir[count] = fv_offset+i;
    moinfo.allvir_sym[count] = h;
      }
      if(h) vr_offset += moinfo.uoccpi[h-1];
      for(i=0; i < moinfo.uoccpi[h]; i++,count++) {
    moinfo.cc_allvir[vr_offset+i] = count;
    moinfo.qt_allvir[count] = vr_offset+i;
    moinfo.allvir_sym[count] = h;
      }
      if(h) op_offset += moinfo.openpi[h-1];
      for(i=0; i < moinfo.openpi[h]; i++,count++) {
    moinfo.cc_allvir[op_offset+i] = count;
    moinfo.qt_allvir[count] = op_offset+i;
    moinfo.allvir_sym[count] = h;
      }
    }

    /* Dump the all-orbital CC ordering and symmetry arrays to CC_INFO */
    psio_write_entry(CC_INFO, "All Occ Orb Symmetry",
             (char *) moinfo.allocc_sym, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "All Virt Orb Symmetry",
             (char *) moinfo.allvir_sym, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "QT->CC All Occ Order",
             (char *) moinfo.cc_allocc, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "QT->CC All Virt Order",
             (char *) moinfo.cc_allvir, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "CC->QT All Occ Order",
             (char *) moinfo.qt_allocc, sizeof(int)*moinfo.nmo);
    psio_write_entry(CC_INFO, "CC->QT All Virt Order",
             (char *) moinfo.qt_allvir, sizeof(int)*moinfo.nmo);

  } /*** RHF and ROHF references ***/

  /* Calculate occupied and virtual orbital offsets within each irrep */

  if(params.ref == 2) {   /*** UHF references ***/
    moinfo.aocc_off = init_int_array(moinfo.nirreps);
    moinfo.bocc_off = init_int_array(moinfo.nirreps);
    moinfo.avir_off = init_int_array(moinfo.nirreps);
    moinfo.bvir_off = init_int_array(moinfo.nirreps);

    ocount = moinfo.aoccpi[0]; vcount = moinfo.avirtpi[0];
    for(h=1; h < moinfo.nirreps; h++) {
      moinfo.aocc_off[h] = ocount;
      ocount += moinfo.aoccpi[h];
      moinfo.avir_off[h] = vcount;
      vcount += moinfo.avirtpi[h];
    }

    ocount = moinfo.boccpi[0]; vcount = moinfo.bvirtpi[0];
    for(h=1; h < moinfo.nirreps; h++) {
      moinfo.bocc_off[h] = ocount;
      ocount += moinfo.boccpi[h];
      moinfo.bvir_off[h] = vcount;
      vcount += moinfo.bvirtpi[h];
    }

    moinfo.all_aocc_off = init_int_array(moinfo.nirreps);
    moinfo.all_bocc_off = init_int_array(moinfo.nirreps);
    moinfo.all_avir_off = init_int_array(moinfo.nirreps);
    moinfo.all_bvir_off = init_int_array(moinfo.nirreps);

    ocount = moinfo.all_aoccpi[0]; vcount = moinfo.all_avirtpi[0];
    for(h=1; h < moinfo.nirreps; h++) {
      moinfo.all_aocc_off[h] = ocount;
      ocount += moinfo.all_aoccpi[h];
      moinfo.all_avir_off[h] = vcount;
      vcount += moinfo.all_avirtpi[h];
    }

    ocount = moinfo.all_boccpi[0]; vcount = moinfo.all_bvirtpi[0];
    for(h=1; h < moinfo.nirreps; h++) {
      moinfo.all_bocc_off[h] = ocount;
      ocount += moinfo.all_boccpi[h];
      moinfo.all_bvir_off[h] = vcount;
      vcount += moinfo.all_bvirtpi[h];
    }

    /* Dump occ_off and vir_off arrays to CC_INFO */
    psio_write_entry(CC_INFO, "Active Alpha Occ Orb Offsets",
             (char *) moinfo.aocc_off, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "Active Beta Occ Orb Offsets",
             (char *) moinfo.bocc_off, sizeof(int)*moinfo.nirreps);

    psio_write_entry(CC_INFO, "Active Alpha Virt Orb Offsets",
             (char *) moinfo.avir_off,  sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "Active Beta Virt Orb Offsets",
             (char *) moinfo.bvir_off,  sizeof(int)*moinfo.nirreps);

    psio_write_entry(CC_INFO, "All Alpha Occ Orb Offsets",
             (char *) moinfo.all_aocc_off, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Beta Occ Orb Offsets",
             (char *) moinfo.all_bocc_off, sizeof(int)*moinfo.nirreps);

    psio_write_entry(CC_INFO, "All Alpha Virt Orb Offsets",
             (char *) moinfo.all_avir_off, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Beta Virt Orb Offsets",
             (char *) moinfo.all_bvir_off, sizeof(int)*moinfo.nirreps);

  }
  else {   /*** RHF/ROHF references ***/
    moinfo.occ_off = init_int_array(moinfo.nirreps);
    moinfo.vir_off = init_int_array(moinfo.nirreps);
    ocount = moinfo.occpi[0]; vcount = moinfo.virtpi[0];
    for(h=1; h < moinfo.nirreps; h++) {
      moinfo.occ_off[h] = ocount;
      ocount += moinfo.occpi[h];
      moinfo.vir_off[h] = vcount;
      vcount += moinfo.virtpi[h];
    }
    moinfo.all_occ_off = init_int_array(moinfo.nirreps);
    moinfo.all_vir_off = init_int_array(moinfo.nirreps);
    ocount = moinfo.all_occpi[0]; vcount = moinfo.all_virtpi[0];
    for(h=1; h < moinfo.nirreps; h++) {
      moinfo.all_occ_off[h] = ocount;
      ocount += moinfo.all_occpi[h];
      moinfo.all_vir_off[h] = vcount;
      vcount += moinfo.all_virtpi[h];
    }

    /* Dump occ_off and vir_off arrays to CC_INFO */
    psio_write_entry(CC_INFO, "Active Occ Orb Offsets",
             (char *) moinfo.occ_off, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "Active Virt Orb Offsets",
             (char *) moinfo.vir_off,  sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Occ Orb Offsets",
             (char *) moinfo.all_occ_off, sizeof(int)*moinfo.nirreps);
    psio_write_entry(CC_INFO, "All Virt Orb Offsets",
             (char *) moinfo.all_vir_off, sizeof(int)*moinfo.nirreps);

  }

  fprintf(outfile,"\n\tChkpt Parameters:\n");
  fprintf(outfile,"\t--------------------\n");
  fprintf(outfile,"\tNumber of irreps     = %d\n",moinfo.nirreps);
  fprintf(outfile,"\tNumber of MOs        = %d\n",moinfo.nmo);
  fprintf(outfile,"\tNumber of active MOs = %d\n\n",moinfo.nactive);
  fprintf(outfile,
      "\tLabel\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
  fprintf(outfile,
      "\t-----\t-----\t------\t------\t------\t------\t------\n");
  for(i=0; i < moinfo.nirreps; i++) {
    fprintf(outfile,
        "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
        moinfo.labels[i],moinfo.orbspi[i],moinfo.frdocc[i],
        moinfo.clsdpi[i],moinfo.openpi[i],moinfo.uoccpi[i],
        moinfo.fruocc[i]);
  }
  fprintf(outfile,"\n\tNuclear Rep. energy (chkpt) =  %20.14f\n", moinfo.enuc);
  fprintf(outfile,  "\tSCF energy          (chkpt) =  %20.14f\n", escf);

  /* Lastly, build the active virtual orbital SCF eigenvector array for
     the AO-basis code (see CCENERGY) */

  if(params.ref == 2) { /*** UHF references ***/

    pitz_offset = init_int_array(moinfo.nirreps);
    pitz_offset[0] = 0;
    for(h=1; h < moinfo.nirreps; h++) {
      pitz_offset[h] = pitz_offset[h-1] + moinfo.orbspi[h-1];
    }

    chkpt_init(PSIO_OPEN_OLD);

    evects_A = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    scf_vector_A = (double ***) malloc(moinfo.nirreps * sizeof(double **));

    next = PSIO_ZERO;
    for(h=0; h < moinfo.nirreps; h++) {
      evects_A[h] = chkpt_rd_alpha_scf_irrep(h); /* Pitzer-ordered MO's */

      if(moinfo.avirtpi[h]) {
    scf_vector_A[h] = block_matrix(moinfo.sopi[h], moinfo.avirtpi[h]);

    for(i=0; i < moinfo.sopi[h]; i++) {
      for(j=0; j < moinfo.avirtpi[h]; j++) { /* CC-ordered relative index */
        J = moinfo.avir_off[h] + j;          /* CC-ordered absolute index */
        qt_j = moinfo.qt_avir[J];            /* QT-ordered absolute index */
        pitz_J = moinfo.qt2pitz_A[qt_j];            /* Pitzer-ordered absolute index */
            pitz_j = pitz_J - pitz_offset[h];    /* Pitzer-ordered relative index */
            scf_vector_A[h][i][j] = evects_A[h][i][pitz_j];
      }
    }

    psio_write(CC_INFO, "UHF Active Alpha Virtual Orbs", (char *) scf_vector_A[h][0],
           moinfo.sopi[h]*moinfo.avirtpi[h]*sizeof(double), next, &next);

        free_block(evects_A[h]);
        free_block(scf_vector_A[h]);

      }

    }

    free(evects_A);  free(scf_vector_A);

    evects_B = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    scf_vector_B = (double ***) malloc(moinfo.nirreps * sizeof(double **));

    next = PSIO_ZERO;
    for(h=0; h < moinfo.nirreps; h++) {
      evects_B[h] = chkpt_rd_beta_scf_irrep(h); /* Pitzer-ordered MO's */

      if(moinfo.bvirtpi[h]) {
    scf_vector_B[h] = block_matrix(moinfo.sopi[h], moinfo.bvirtpi[h]);

    for(i=0; i < moinfo.sopi[h]; i++) {
      for(j=0; j < moinfo.bvirtpi[h]; j++) { /* CC-ordered relative index */
        J = moinfo.bvir_off[h] + j;          /* CC-ordered absolute index */
        qt_j = moinfo.qt_bvir[J];            /* QT-ordered absolute index */
        pitz_J = moinfo.qt2pitz_B[qt_j];            /* Pitzer-ordered absolute index */
            pitz_j = pitz_J - pitz_offset[h];    /* Pitzer-ordered relative index */
            scf_vector_B[h][i][j] = evects_B[h][i][pitz_j];
      }
    }

    psio_write(CC_INFO, "UHF Active Beta Virtual Orbs", (char *) scf_vector_B[h][0],
           moinfo.sopi[h]*moinfo.bvirtpi[h]*sizeof(double), next, &next);

    free_block(evects_B[h]);
    free_block(scf_vector_B[h]);

      }

    }

    free(evects_B);  free(scf_vector_B);

    free(pitz_offset);
    chkpt_close();
  }
  else { /*** RHF/ROHF references ***/

    pitz_offset = init_int_array(moinfo.nirreps);
    pitz_offset[0] = 0;
    for(h=1; h < moinfo.nirreps; h++) {
      pitz_offset[h] = pitz_offset[h-1] + moinfo.orbspi[h-1];
    }

    chkpt_init(PSIO_OPEN_OLD);
    evects = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    scf_vector = (double ***) malloc(moinfo.nirreps * sizeof(double **));

    next = PSIO_ZERO;
    for(h=0; h < moinfo.nirreps; h++) {
      evects[h] = chkpt_rd_scf_irrep(h); /* Pitzer-order MO's */

      /* Are there active virtuals in this irrep? */
      if(moinfo.virtpi[h]) {

    scf_vector[h] = block_matrix(moinfo.sopi[h],moinfo.virtpi[h]);

    for(i=0; i < moinfo.sopi[h]; i++) {
      for(j=0; j < moinfo.virtpi[h]; j++) { /* CC-ordered relative index */
        J = moinfo.vir_off[h] + j;          /* CC-ordered absolute index */
        qt_j = moinfo.qt_vir[J];            /* QT-ordered absolute index */
        pitz_J = moinfo.qt2pitz[qt_j];             /* Pitzer-ordered absolute index */
        pitz_j = pitz_J - pitz_offset[h];   /* Pitzer-order relative index */
        scf_vector[h][i][j] = evects[h][i][pitz_j];
      }
    }

    psio_write(CC_INFO, "RHF/ROHF Active Virtual Orbitals",
           (char *) scf_vector[h][0],
           moinfo.sopi[h]*moinfo.virtpi[h]*sizeof(double),
           next, &next);

    free_block(scf_vector[h]);

      }
      free_block(evects[h]);
    }

    for(h=0; h < moinfo.nirreps; h++) {
      evects[h] = chkpt_rd_scf_irrep(h); /* Pitzer-order MO's */
      next = PSIO_ZERO;
      if(moinfo.occpi[h]) {
    scf_vector[h] = block_matrix(moinfo.sopi[h],moinfo.occpi[h]);

    for(i=0; i < moinfo.sopi[h]; i++) {
      for(j=0; j < moinfo.occpi[h]; j++) { /* CC-ordered relative index */
        J = moinfo.occ_off[h] + j;          /* CC-ordered absolute index */
        qt_j = moinfo.qt_occ[J];            /* QT-ordered absolute index */
        pitz_J = moinfo.qt2pitz[qt_j];             /* Pitzer-ordered absolute index */
        pitz_j = pitz_J - pitz_offset[h];   /* Pitzer-order relative index */
        scf_vector[h][i][j] = evects[h][i][pitz_j];
      }
    }

    psio_write(CC_INFO, "RHF/ROHF Active Occupied Orbitals",
           (char *) scf_vector[h][0],
           moinfo.sopi[h]*moinfo.occpi[h]*sizeof(double),
           next, &next);

    free_block(scf_vector[h]);

      }

      free_block(evects[h]);

    /*
      fprintf(outfile, "\n\tOriginal SCF Eigenvectors:\n");
      print_mat(evects[h], moinfo.sopi[h], moinfo.orbspi[h], outfile);

      fprintf(outfile, "\n\tRe-ordered Virtual SCF Eigenvectors:\n");
      print_mat(scf_vector[h], moinfo.sopi[h], moinfo.virtpi[h], outfile);
    */
    }

    C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    next = PSIO_ZERO;
    for(h=0; h < moinfo.nirreps; h++) {
      if(moinfo.sopi[h] && moinfo.virtpi[h]) {
        C[h] = block_matrix(moinfo.sopi[h],moinfo.virtpi[h]);
        psio_read(CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *) C[h][0],
                  moinfo.sopi[h]*moinfo.virtpi[h]*sizeof(double), next, &next);
      }
    }
    moinfo.C = C;

    free(evects);  free(scf_vector);
    free(pitz_offset);
    chkpt_close();
  }
}

}} // namespace psi::ccsort
