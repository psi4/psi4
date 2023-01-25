/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#ifndef CCDENSITY_MOINFO_H
#define CCDENSITY_MOINFO_H

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/

#include <string>
#include <vector>
#include "psi4/libmints/matrix.h"

namespace psi {
namespace ccdensity {

struct MOInfo {
    int nirreps;                     /* no. of irreducible representations */
    int nmo;                         /* no. of molecular orbitals */
    int nso;                         /* no. of symmetry orbitals */
    int nactive;                     /* no. of active orbitals */
    Dimension orbspi;                /* no. of MOs per irrep */
    Dimension clsdpi;                /* no. of closed-shells per irrep excl. frdocc */
    Dimension openpi;                /* no. of open-shells per irrep */
    Dimension uoccpi;                /* no. of unoccupied orbitals per irrep excl. fruocc */
    Dimension frdocc;                /* no. of frozen core orbitals per irrep */
    Dimension fruocc;                /* no. of frozen unoccupied orbitals per irrep */
    std::vector<std::string> labels; /* irrep labels */
    int nfzc;                        /* total no. of frozen core orbitals */
    int nfzv;                        /* total no. of frozen virtual orbitals */
    int nclsd;                       /* total no. of closd shells excl. frdocc */
    int nopen;                       /* total no. of open shells  */
    int nuocc;                       /* total no. of unoccupied shells excl. fruocc */
    std::vector<int> occ_sym;        /* active occupied index symmetry */
    std::vector<int> aocc_sym;       /* alpha active occupied index symmetry */
    std::vector<int> bocc_sym;       /* beta active occupied index symmetry */
    std::vector<int> vir_sym;        /* active virtual index symmetry */
    std::vector<int> avir_sym;       /* alpha active virtual index symmetry */
    std::vector<int> bvir_sym;       /* beta active virtual index symmetry */
    int sym;                         /* symmetry of converged CCSD state */
    Dimension occpi;                 /* no. of active occ. orbs. (incl. open) per irrep */
    Dimension aoccpi;                /* no. of alpha active occ. orbs. (incl. open) per irrep */
    Dimension boccpi;                /* no. of beta active occ. orbs. (incl. open) per irrep */
    Dimension virtpi;                /* no. of active virt. orbs. (incl. open) per irrep */
    Dimension avirtpi;               /* no. of alpha active virt. orbs. (incl. open) per irrep */
    Dimension bvirtpi;               /* no. of beta active virt. orbs. (incl. open) per irrep */
    Dimension occ_off;               /* occupied orbital offsets within each irrep */
    Dimension aocc_off;              /* alpha occupied orbital offsets within each irrep */
    Dimension bocc_off;              /* beta occupied orbital offsets within each irrep */
    Dimension vir_off;               /* virtual orbital offsets within each irrep */
    Dimension avir_off;              /* alpha virtual orbital offsets within each irrep */
    Dimension bvir_off;              /* beta virtual orbital offsets within each irrep */
    std::vector<int> cc_occ;         /* QT->CC active occupied reordering array */
    std::vector<int> cc_aocc;        /* QT->CC alpha active occupied reordering array */
    std::vector<int> cc_bocc;        /* QT->CC beta active occupied reordering array */
    std::vector<int> cc_vir;         /* QT->CC active virtiual reordering array */
    std::vector<int> cc_avir;        /* QT->CC alpha active virtiual reordering array */
    std::vector<int> cc_bvir;        /* QT->CC beta active virtiual reordering array */
    std::vector<int> qt_occ;         /* CC->QT active occupied reordering array */
    std::vector<int> qt_aocc;        /* CC->QT alpha active occupied reordering array */
    std::vector<int> qt_bocc;        /* CC->QT beta active occupied reordering array */
    std::vector<int> qt_vir;         /* CC->QT active virtiual reordering array */
    std::vector<int> qt_avir;        /* CC->QT alpha active virtiual reordering array */
    std::vector<int> qt_bvir;        /* CC->QT beta active virtiual reordering array */
    double enuc;                     /* Nuclear repulsion energy */
    double escf;                     /* SCF energy from wfn */
    double eref;                     /* Reference energy */
    double ecc;                      /* CC energy (CC2, CCSD, or CC3) from ccenergy */
    double et;                       /* (T) energy from cctriples */
    Matrix opdm;                     /* Onepdm in the full (fzc+clsd+socc+uocc) space */
    Matrix opdm_a;                   /* Alpha Onepdm in the full (fzc+clsd+socc+uocc) space */
    Matrix opdm_b;                   /* Beta Onepdm in the full (fzc+clsd+socc+uocc) space */
    Matrix ltd_mat;                  /* <0|O|n> Left transition density */
    Matrix ltd_a_mat;                /* <0|O|n> Left transition alpha density */
    Matrix ltd_b_mat;                /* <0|O|n> Left transition beta density */
    Matrix rtd_mat;                  /* <n|O|0> Right transition density */
    Matrix rtd_a_mat;                /* <n|O|0> Right transition alpha density */
    Matrix rtd_b_mat;                /* <n|O|0> Right transition beta density */
    std::vector<int> pitzer2qt;      /* Pitzer to QT re-ordering array */
    std::vector<int> qt2pitzer;      /* QT to Pitzer re-ordering array */
    SharedMatrix Ca;                 /* SCF orbitals (standard ordering) */
    std::vector<SharedMatrix> L;
    std::vector<SharedMatrix> nabla;
    std::vector<SharedMatrix> dip;
};

}  // namespace ccdensity
}  // namespace psi

#endif
