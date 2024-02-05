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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/

#ifndef _psi_src_bin_cceom_moinfo_h
#define _psi_src_bin_cceom_moinfo_h

#include <string>
#include <vector>

namespace psi {
namespace cceom {

struct MOInfo {
    int nirreps;                       /* no. of irreducible representations */
    int nmo;                           /* no. of molecular orbitals */
    int nso;                           /* no. of symmetry orbitals */
    int iopen;                         /* 0=closed shell; >0=open shell */
    int *sopi;                         /* no. of SOs per irrep */
    int *sosym;                        /* orbital symmetry (Pitzer/SO) */
    int *orbspi;                       /* no. of MOs per irrep */
    int *clsdpi;                       /* no. of closed-shells per irrep excl. frdocc */
    int *openpi;                       /* no. of open-shells per irrep */
    int *uoccpi;                       /* no. of unoccupied orbitals per irr. ex. fruocc */
    int *frdocc;                       /* no. of frozen core orbitals per irrep */
    int *fruocc;                       /* no. of frozen unoccupied orbitals per irrep */
    int nvirt;                         /* total no. of (active) virtual orbitals */
    std::vector<std::string> irr_labs; /* irrep labels */
    char **irr_labs_lowercase;         /* irrep labels */
    int *occpi;                        /* no. of occupied orbs. (incl. open) per irrep */
    int *aoccpi;                       /* no. of alpha occupied orbs. (incl. open) per irrep */
    int *boccpi;                       /* no. of beta occupied orbs. (incl. open) per irrep */
    int *virtpi;                       /* no. of virtual orbs. (incl. open) per irrep */
    int *avirtpi;                      /* no. of alpha virtual orbs. (incl. open) per irrep */
    int *bvirtpi;                      /* no. of beta virtual orbs. (incl. open) per irrep */
    int *occ_sym;                      /* relative occupied index symmetry */
    int *aocc_sym;                     /* relative alpha occupied index symmetry */
    int *bocc_sym;                     /* relative beta occupied index symmetry */
    int *vir_sym;                      /* relative virtual index symmetry */
    int *avir_sym;                     /* relative alpha virtual index symmetry */
    int *bvir_sym;                     /* relative beta virtual index symmetry */
    int iter;                          /* Current CCSD iteration */
    int sym;                           /* symmetry of converged CCSD state */
    int *occ_off;                      /* occupied orbital offsets within each irrep */
    int *aocc_off;                     /* alpha occupied orbital offsets within each irrep */
    int *bocc_off;                     /* beta occupied orbital offsets within each irrep */
    int *vir_off;                      /* virtual orbital offsets within each irrep */
    int *avir_off;                     /* alpha virtual orbital offsets within each irrep */
    int *bvir_off;                     /* beta virtual orbital offsets within each irrep */
    double conv;                       /* Current convergence level */
    double enuc;                       /* Nuclear repulsion energy */
    double escf;                       /* SCF energy (from wfn) */
    double eref;                       /* Reference energy (file100) */
    double ecc;                        /* Current coupled cluster energy */
    double t1diag;                     /* Standard open- or closed-shell T1 diagnostic */
    double d1diag;                     /* Janssen and Nielsen's D1 Diagnostic */
    double ***C;                       /* Virtual orbital transformation matrix (for AO-basis B terms) */
    double ***Ca;                      /* UHF alpha virtual orbital transformation matrix (for AO-basis B terms) */
    double ***Cb;                      /* UHF beta virtual orbital transformation matrix (for AO-basis B terms) */
};

}  // namespace cceom
}  // namespace psi

#endif // _psi_src_bin_cceom_moinfo_h
