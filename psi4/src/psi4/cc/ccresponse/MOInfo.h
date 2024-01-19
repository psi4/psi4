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

#ifndef CCRESPONSE_MOINFO_H
#define CCRESPONSE_MOINFO_H

/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here
*/

#include <string>
#include <vector>
#include "psi4/libmints/matrix.h"

namespace psi {
namespace ccresponse {

struct MOInfo {
    int nirreps;                     /* no. of irreducible representations */
    int nmo;                         /* no. of molecular orbitals */
    int nso;                         /* no. of symmetry orbitals */
    int nao;                         /* no. of atomic orbitals */
    int noei;                        /* no. of elements in SOxSO lower triangle */
    int ntri;                        /* no. of elements in MOxMO lower triangle */
    int noei_ao;                     /* no. of elements in AOxAO lower triangle */
    int nactive;                     /* no. of active MO's */
    int nfzc;                        /* no. of frozen core orbitals */
    int *sopi;                       /* no. of SOs per irrep */
    int *orbspi;                     /* no. of MOs per irrep */
    int *clsdpi;                     /* no. of closed-shells per irrep  */
    int *openpi;                     /* no. of open-shells per irrep */
    int *uoccpi;                     /* no. of unoccupied orbitals per irrep  */
    int *frdocc;                     /* no. of frozen core orbitals per irrep */
    int *fruocc;                     /* no. of frozen unoccupied orbitals per irrep */
    int nvirt;                       /* total no. of (active) virtual orbitals */
    int *actpi;                      /* no. of active orbitals per irrep */
    Dimension act_pi;                /* Dimension form of actpi */
    std::vector<std::string> labels; /* irrep labels */
    int *occpi;                      /* no. of occupied orbs. (incl. open) per irrep */
    Dimension act_occpi;             /* Dimension form of occpi */
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
    int *occ_off;                    /* occupied orbital offsets within each irrep */
    int *aocc_off;                   /* occupied alpha orbital offsets within each irrep */
    int *bocc_off;                   /* occupied beta orbital offsets within each irrep */
    int *vir_off;                    /* virtual orbital offsets within each irrep */
    int *avir_off;                   /* virtual alpha orbital offsets within each irrep */
    int *bvir_off;                   /* virtual beta orbital offsets within each irrep */
    std::shared_ptr<Matrix> Ca;
    int *mu_irreps;                  /* irreps of x,y,z dipole components */
    int *l_irreps;                   /* irreps of x,y,z angular momentum components */
    int natom;                       /* number of atoms */
    double *zvals;                   /* atomic zvals */
    double ***C;                     /* Virtual orbital transformation matrix */
};

}  // namespace ccresponse
}  // namespace psi
#endif
