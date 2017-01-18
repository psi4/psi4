/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_ccenergy_moinfo_h
#define _psi_src_bin_ccenergy_moinfo_h

namespace psi { namespace ccenergy {

struct MOInfo {
  int nirreps;           /* no. of irreducible representations */
  int nmo;               /* no. of molecular orbitals */
  int nso;               /* no. of symmetry orbitals */
  int nao;               /* no. of atomic orbitals */
  int *sopi;             /* no. of SOs per irrep (only used in AO-based algorithm) */
  int *sosym;            /* SO symmetry (Pitzer) */
  int *orbspi;           /* no. of MOs per irrep */
  int *clsdpi;           /* no. of closed-shells per irrep excl. frdocc */
  int *openpi;           /* no. of open-shells per irrep */
  int *uoccpi;           /* no. of unoccupied orbitals per irr. ex. fruocc */
  int *frdocc;           /* no. of frozen core orbitals per irrep */
  int *fruocc;           /* no. of frozen unoccupied orbitals per irrep */
  int nvirt;             /* total no. of virtual orbitals */
  char **labels;         /* irrep labels */
  int *occpi;            /* no. of occupied orbs. (incl. open) per irrep */
  int *aoccpi;           /* no. of alpha occupied orbs. (incl. open) per irrep */
  int *boccpi;           /* no. of beta occupied orbs. (incl. open) per irrep */
  int *virtpi;           /* no. of virtual orbs. (incl. open) per irrep */
  int *avirtpi;          /* no. of alpha virtual orbs. (incl. open) per irrep */
  int *bvirtpi;          /* no. of beta virtual orbs. (incl. open) per irrep */
  int *occ_sym;          /* relative occupied index symmetry */
  int *aocc_sym;         /* relative alpha occupied index symmetry */
  int *bocc_sym;         /* relative beta occupied index symmetry */
  int *vir_sym;          /* relative virtual index symmetry */
  int *avir_sym;         /* relative alpha virtual index symmetry */
  int *bvir_sym;         /* relative beta virtual index symmetry */

  int *occ_off;       /* occupied orbital offsets within each irrep */
  int *aocc_off;      /* alpha occupied orbital offsets within each irrep */
  int *bocc_off;      /* beta occupied orbital offsets within each irrep */
  int *vir_off;       /* virtual orbital offsets within each irrep */
  int *avir_off;      /* alpha virtual orbital offsets within each irrep */
  int *bvir_off;      /* beta virtual orbital offsets within each irrep */
  int *cc_occ;        /* QT->CC active occupied reordering array */
  int *cc_aocc;       /* QT->CC alpha active occupied reordering array */
  int *cc_bocc;       /* QT->CC beta active occupied reordering array */
  int *cc_vir;        /* QT->CC active virtiual reordering array */
  int *cc_avir;       /* QT->CC alpha active virtiual reordering array */
  int *cc_bvir;       /* QT->CC beta active virtiual reordering array */
  int *qt_occ;        /* CC->QT active occupied reordering array */
  int *qt_aocc;       /* CC->QT alpha active occupied reordering array */
  int *qt_bocc;       /* CC->QT beta active occupied reordering array */
  int *qt_vir;        /* CC->QT active virtiual reordering array */
  int *qt_avir;       /* CC->QT alpha active virtiual reordering array */
  int *qt_bvir;       /* CC->QT beta active virtiual reordering array */

  int *pitzer2qt;     /* Pitzer -> QT translation array */
  int *qt2pitzer;     /* QT -> Pitzer translation array */

  int *pitzer2qt_a;   /* Pitzer -> QT translation array for alpha orbitals */
  int *qt2pitzer_a;   /* QT -> Pitzer translation array for alpha orbitals */
  int *pitzer2qt_b;   /* Pitzer -> QT translation array for beta orbitals */
  int *qt2pitzer_b;   /* QT -> Pitzer translation array for beta orbitals */

  int iter;              /* Current CCSD iteration */
  double conv;           /* Current convergence level */
  double enuc;           /* Nuclear repulsion energy */
  double emp2;           /* MP2 energy */
  double emp2_ss;        /* Same-spin MP2 correlation energy*/
  double emp2_os;        /* Opposite-spin MP2 correlation energy*/
  double escf;           /* SCF energy (from wfn) */
  double eref;           /* Reference energy (file100) */
  double ecc;            /* Current coupled cluster correlation energy */
  double ecc_ss;	 /* Same-spin coupled cluster correlation energy*/
  double ecc_os;         /* Opposite-spin coupled cluster energy*/
  double t1diag;         /* Standard open- or closed-shell T1 diagnostic */
  double d1diag;         /* Janssen and Nielsen's D1 Diagnostic */
  double new_d1diag;     /* Lee's modified D1 Diagnostic */
  double d2diag;         /* Nielsen and Janssen's D2 Diagnostic */
  double ***Cv;           /* Virtual orbital transformation matrix (for AO-basis B terms) */
  double ***Cav;          /* UHF alpha virtual orbital transformation matrix (for AO-basis B terms) */
  double ***Cbv;          /* UHF beta virtual orbital transformation matrix (for AO-basis B terms) */
  double ***Co;           /* Occupied orbital transformation matrix (for AO-basis B terms) */
  double ***Cao;          /* UHF alpha occupied orbital transformation matrix (for AO-basis B terms) */
  double ***Cbo;          /* UHF beta occupied orbital transformation matrix (for AO-basis B terms) */
};

}} // namespace psi::ccenergy

#endif //  _psi_src_bin_ccenergy_moinfo_h