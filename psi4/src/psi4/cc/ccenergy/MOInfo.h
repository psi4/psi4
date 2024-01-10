/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
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

#include <string>
#include <vector>
#include "psi4/libmints/dimension.h"

namespace psi {
namespace ccenergy {

struct MOInfo {
    int nirreps;                     /* no. of irreducible representations */
    int nmo;                         /* no. of molecular orbitals */
    int nso;                         /* no. of symmetry orbitals */
    int nao;                         /* no. of atomic orbitals */
    int *sopi;                       /* no. of SOs per irrep (only used in AO-based algorithm) */
    int *sosym;                      /* SO symmetry (Pitzer) */
    int *orbspi;                     /* no. of MOs per irrep */
    Dimension clsdpi;                /* no. of closed-shells per irrep excl. frdocc */
    Dimension openpi;                /* no. of open-shells per irrep */
    int *uoccpi;                     /* no. of unoccupied orbitals per irr. ex. fruocc */
    int *frdocc;                     /* no. of frozen core orbitals per irrep */
    int *fruocc;                     /* no. of frozen unoccupied orbitals per irrep */
    int nvirt;                       /* total no. of virtual orbitals */
    std::vector<std::string> labels; /* irrep labels */
    int *occpi;                      /* no. of occupied orbs. (incl. open) per irrep */
    int *aoccpi;                     /* no. of alpha occupied orbs. (incl. open) per irrep */
    int *boccpi;                     /* no. of beta occupied orbs. (incl. open) per irrep */
    int *virtpi;                     /* no. of virtual orbs. (incl. open) per irrep */
    int *avirtpi;                    /* no. of alpha virtual orbs. (incl. open) per irrep */
    int *bvirtpi;                    /* no. of beta virtual orbs. (incl. open) per irrep */
    int *occ_sym;                    /* relative occupied index symmetry */
    int *aocc_sym;                   /* relative alpha occupied index symmetry */
    int *bocc_sym;                   /* relative beta occupied index symmetry */
    int *vir_sym;                    /* relative virtual index symmetry */
    int *avir_sym;                   /* relative alpha virtual index symmetry */
    int *bvir_sym;                   /* relative beta virtual index symmetry */

    int *occ_off;  /* occupied orbital offsets within each irrep */
    int *aocc_off; /* alpha occupied orbital offsets within each irrep */
    int *bocc_off; /* beta occupied orbital offsets within each irrep */
    int *vir_off;  /* virtual orbital offsets within each irrep */
    int *avir_off; /* alpha virtual orbital offsets within each irrep */
    int *bvir_off; /* beta virtual orbital offsets within each irrep */

    int iter;          /* Current CCSD iteration */
    double conv;       /* Current convergence level */
    double enuc;       /* Nuclear repulsion energy */
    double emp2;       /* MP2 energy */
    double emp2_ss;    /* Same-spin MP2 correlation energy*/
    double emp2_os;    /* Opposite-spin MP2 correlation energy*/
    double emp2_s;     /* Singles MP2 correlation energy*/
    double escf;       /* SCF energy (from wfn) */
    double eref;       /* Reference energy (file100) */
    double ecc;        /* Current coupled cluster correlation energy */
    double ecc_ss;     /* Same-spin coupled cluster correlation energy*/
    double ecc_os;     /* Opposite-spin coupled cluster energy*/
    double ecc_s;      /* Singles coupled cluster energy */
    double t1diag;     /* Standard open- or closed-shell T1 diagnostic */
    double d1diag;     /* Janssen and Nielsen's D1 Diagnostic */
    double new_d1diag; /* Lee's modified D1 Diagnostic */
    double d2diag;     /* Nielsen and Janssen's D2 Diagnostic */
    double ***Cv;      /* Virtual orbital transformation matrix (for AO-basis B terms) */
    double ***Cav;     /* UHF alpha virtual orbital transformation matrix (for AO-basis B terms) */
    double ***Cbv;     /* UHF beta virtual orbital transformation matrix (for AO-basis B terms) */
    double ***Co;      /* Occupied orbital transformation matrix (for AO-basis B terms) */
    double ***Cao;     /* UHF alpha occupied orbital transformation matrix (for AO-basis B terms) */
    double ***Cbo;     /* UHF beta occupied orbital transformation matrix (for AO-basis B terms) */
};

}  // namespace ccenergy
}  // namespace psi

#endif  //  _psi_src_bin_ccenergy_moinfo_h
