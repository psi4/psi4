/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "psi4/libmints/dimension.h"
#include "psi4/libmints/matrix.h"

namespace psi {

class Wavefunction;

namespace cc {

enum class Reference;

// FIXME Duplicated from cctransort/pitzer2qt.cc
std::vector<int> pitzer2qt(std::vector<Dimension>& spaces);

struct CCMOInfo final {
    CCMOInfo() {}
    CCMOInfo(std::shared_ptr<Wavefunction> wfn, Reference ref);

    /*! no. of irreducible representations */
    int nirreps;
    /*! no. of molecular orbitals */
    int nmo;
    /*! no. of symmetry orbitals */
    int nso;
    /*! no. of atomic orbitals */
    int nao;
    /*! total no. of virtual orbitals */
    int nvirt;
    /*! no. of active orbitals */
    int nactive;

    /*! irrep labels
     * FIXME this is currently only used for printing
     * it can probably be trashed.
     */
    std::vector<std::string> labels;

    /*! Nuclear repulsion energy */
    double enuc;
    /*! SCF energy */
    double escf;

    /*! @{ Orbital spaces dimensions */
    /*! no. of MOs per irrep */
    Dimension orbspi;
    /*! no. of SOs per irrep (only used in AO-based algorithm) */
    Dimension sopi;
    /*! no. of closed-shells per irrep excl. frdocc */
    Dimension clsdpi;
    /*! no. of open-shells per irrep */
    Dimension openpi;
    /*! no. of unoccupied orbitals per irr. ex. fruocc */
    Dimension uoccpi;
    /*! no. of frozen core orbitals per irrep */
    Dimension frdocc;
    /*! no. of frozen unoccupied orbitals per irrep */
    Dimension fruocc;
    /*! no. of occupied orbs. (incl. open) per irrep */
    Dimension occpi;
    /*! no. of alpha occupied orbs. (incl. open) per irrep */
    Dimension aoccpi;
    /*! no. of beta occupied orbs. (incl. open) per irrep */
    Dimension boccpi;
    /*! no. of virtual orbs. (incl. open) per irrep */
    Dimension virtpi;
    /*! no. of alpha virtual orbs. (incl. open) per irrep */
    Dimension avirtpi;
    /*! no. of beta virtual orbs. (incl. open) per irrep */
    Dimension bvirtpi;
    /*! SO symmetry (Pitzer) */
    std::vector<int> sosym;
    /*! @}*/

    /*! @{ Relative index symmetry */
    /*! relative occupied index symmetry */
    std::vector<int> occ_sym;
    /*! relative alpha occupied index symmetry */
    std::vector<int> aocc_sym;
    /*! relative beta occupied index symmetry */
    std::vector<int> bocc_sym;
    /*! relative virtual index symmetry */
    std::vector<int> vir_sym;
    /*! relative alpha virtual index symmetry */
    std::vector<int> avir_sym;
    /*! relative beta virtual index symmetry */
    std::vector<int> bvir_sym;
    /*! @}*/

    /*! @{ Offsets within each irrep */
    /*! occupied orbital offsets within each irrep */
    std::vector<int> occ_off;
    /*! alpha occupied orbital offsets within each irrep */
    std::vector<int> aocc_off;
    /*! beta occupied orbital offsets within each irrep */
    std::vector<int> bocc_off;
    /*! virtual orbital offsets within each irrep */
    std::vector<int> vir_off;
    /*! alpha virtual orbital offsets within each irrep */
    std::vector<int> avir_off;
    /*! beta virtual orbital offsets within each irrep */
    std::vector<int> bvir_off;
    /*! @}*/

    /*! @{ Index reordeding arrays */
    /*! QT->CC active occupied reordering array */
    std::vector<int> cc_occ;
    /*! QT->CC alpha active occupied reordering array */
    std::vector<int> cc_aocc;
    /*! QT->CC beta active occupied reordering array */
    std::vector<int> cc_bocc;
    /*! QT->CC active virtual reordering array */
    std::vector<int> cc_vir;
    /*! QT->CC alpha active virtual reordering array */
    std::vector<int> cc_avir;
    /*! QT->CC beta active virtual reordering array */
    std::vector<int> cc_bvir;
    /*! CC->QT active occupied reordering array */
    std::vector<int> qt_occ;
    /*! CC->QT alpha active occupied reordering array */
    std::vector<int> qt_aocc;
    /*! CC->QT beta active occupied reordering array */
    std::vector<int> qt_bocc;
    /*! CC->QT active virtual reordering array */
    std::vector<int> qt_vir;
    /*! CC->QT alpha active virtual reordering array */
    std::vector<int> qt_avir;
    /*! CC->QT beta active virtual reordering array */
    std::vector<int> qt_bvir;
    /*! Pitzer -> QT translation array */
    std::vector<int> pitzer2qt;
    /*! QT -> Pitzer translation array */
    std::vector<int> qt2pitzer;
    /*! Pitzer -> QT translation array for alpha orbitals */
    std::vector<int> pitzer2qt_a;
    /*! QT -> Pitzer translation array for alpha orbitals */
    std::vector<int> qt2pitzer_a;
    /*! Pitzer -> QT translation array for beta orbitals */
    std::vector<int> pitzer2qt_b;
    /*! QT -> Pitzer translation array for beta orbitals */
    std::vector<int> qt2pitzer_b;
    /*! @}*/

    /*! @{ MO coefficient matrices */
    /* Virtual orbital transformation matrix (for AO-basis B terms) */
    std::shared_ptr<Matrix> Cv;
    /* UHF alpha virtual orbital transformation matrix (for AO-basis B terms) */
    std::shared_ptr<Matrix> Cav;
    /* UHF beta virtual orbital transformation matrix (for AO-basis B terms) */
    std::shared_ptr<Matrix> Cbv;
    /* Occupied orbital transformation matrix (for AO-basis B terms) */
    std::shared_ptr<Matrix> Co;
    /* UHF alpha occupied orbital transformation matrix (for AO-basis B terms) */
    std::shared_ptr<Matrix> Cao;
    /* UHF beta occupied orbital transformation matrix (for AO-basis B terms) */
    std::shared_ptr<Matrix> Cbo;
    /*! @}*/
};

void print_ccmoinfo(const CCMOInfo& info);

}  // namespace cc
}  // namespace psi
