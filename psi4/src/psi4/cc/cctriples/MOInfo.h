/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#ifndef CCTRIPLES_MOINFO_H
#define CCTRIPLES_MOINFO_H

/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here
*/

#include <string>
#include <vector>
#include "psi4/libmints/dimension.h"

namespace psi {
namespace cctriples {

struct MOInfo {
    int nirreps;                     /* no. of irreducible representations */
    int nmo;                         /* no. of molecular orbitals */
    Dimension orbspi;                /* no. of MOs per irrep */
    Dimension clsdpi;                /* no. of closed-shells per irrep excl. frdocc */
    Dimension openpi;                /* no. of open-shells per irrep */
    Dimension uoccpi;                /* no. of unoccupied orbitals per irrep excl. fruocc */
    Dimension frdocc;                /* no. of frozen core orbitals per irrep */
    Dimension fruocc;                /* no. of frozen unoccupied orbitals per irrep */
    std::vector<std::string> labels; /* irrep labels */
    Dimension occpi;                 /* no. of occupied orbs. (incl. open) per irrep */
    Dimension aoccpi;                /* no. of alpha occupied orbs. (incl. open) per irrep */
    Dimension boccpi;                /* no. of beta occupied orbs. (incl. open) per irrep */
    Dimension virtpi;                /* no. of virtual orbs. (incl. open) per irrep */
    Dimension avirtpi;               /* no. of alpha virtual orbs. (incl. open) per irrep */
    Dimension bvirtpi;               /* no. of beta virtual orbs. (incl. open) per irrep */
    int *occ_sym;                    /* relative occupied index symmetry */
    int *aocc_sym;                   /* relative alpha occupied index symmetry */
    int *bocc_sym;                   /* relative beta occupied index symmetry */
    int *vir_sym;                    /* relative virtual index symmetry */
    int *avir_sym;                   /* relative alpha virtual index symmetry */
    int *bvir_sym;                   /* relative beta virtual index symmetry */
    int *occ_off;                    /* occupied orbital offsets within each irrep */
    int *aocc_off;                   /* alpha occupied orbital offsets within each irrep */
    int *bocc_off;                   /* beta occupied orbital offsets within each irrep */
    int *vir_off;                    /* virtual orbital offsets within each irrep */
    int *avir_off;                   /* alpha virtual orbital offsets within each irrep */
    int *bvir_off;                   /* beta virtual orbital offsets within each irrep */
    int iter;                        /* Current lambda iteration */
    double conv;                     /* Current convergence level */
    double enuc;                     /* Nuclear repulsion energy */
    double escf;                     /* SCF energy from wfn */
    double eref;                     /* Reference energy (file100) */
    double ecc;                      /* CC energy from ccenergy */
};

}  // namespace cctriples
}  // namespace psi

#endif
