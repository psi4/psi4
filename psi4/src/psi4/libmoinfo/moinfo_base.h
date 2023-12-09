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

#ifndef _psi_src_lib_libmoinfo_moinfo_base_h_
#define _psi_src_lib_libmoinfo_moinfo_base_h_

/*! \file    moinfo_base.h
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding MOs
*/

#include <string>
#include "psi4/libpsi4util/libpsi4util.h"

namespace psi {
using intvec = std::vector<int>;
using boolvec = std::vector<bool>;

class Options;
class Wavefunction;

/// @brief This class stores all the basic info regarding MOs
class MOInfoBase {
   public:
    MOInfoBase(Wavefunction& ref_wfn_, Options& options_, bool silent_ = false);

    /// @brief Get the nuclear energy stored in an MOInfoBase object (or an object derived from MOInfoBase).
    /// @return The nuclear repulsion energy
    double get_nuclear_energy() const { return (nuclear_energy); }

    /// @brief Get one of the irrep labels that are stored in an MOInfoBase object (or an object derived from MOInfoBase). Not bounds-checked.
    /// @param i : Index of the irrep label
    /// @return The selected irrep label
    std::string get_irr_lab(size_t i) const { return (irr_labs[i]); }

    /// @brief Get the # of irreps an MOInfoBase object (or an object derived from MOInfoBase) has been constructed with.
    /// @return The # of irreps
    int get_nirreps() const { return (nirreps); }

    /// @brief Get the PSI nso value (# of symmetry-adapted atomic orbitals) that is stored in an MOInfoBase object (or an object derived from MOInfoBase).
    /// @return The # of symmetry-adapted atomic orbitals
    int get_nso() const { return (nso); }

    /// @brief Get a copy of the array holding the numbers of SOs per irrep, from an MOInfoBase object (or an object derived from MOInfoBase).
    /// @return A copy of the array holding the numbers of SOs per irrep
    intvec get_sopi() const { return (sopi); }

    /// @brief Get a copy of the array holding the numbers of doubly occupied orbitals (DOCC) per irrep, from an MOInfoBase object (or an object derived from MOInfoBase).
    /// @return A copy of the array holding the numbers of doubly occupied orbitals (DOCC) per irrep
    intvec get_docc() const { return (docc); }

    /// @brief Get a copy of the array holding the numbers of active orbitals per irrep, from an MOInfoBase object (or an object derived from MOInfoBase).
    /// @return A copy of the array holding the numbers of active orbitals per irrep
    intvec get_actv() const { return (actv); }

    bool get_guess_occupation() const { return (guess_occupation); }

    /// @brief Get the total # of active orbitals that is stored in an MOInfoBase object (or an object derived from MOInfoBase).
    /// @return Total # of active orbitals across all irreps
    int get_nactv() const { return (nactv); }

    /// @brief Get the # of alpha electrons (including frozen) that is stored in an MOInfoBase object (or an object derived from MOInfoBase).
    /// @return The # of alpha electrons (including frozen)
    int get_nael() const { return (nael); }

    /// @brief Get the # of beta electrons (including frozen) that is stored in an MOInfoBase object (or an object derived from MOInfoBase).
    /// @return The # of beta electrons (including frozen)
    int get_nbel() const { return (nbel); }

    double** get_scf_mos() const { return (scf); }

   protected:
    void read_data();
    void compute_number_of_electrons();
    void read_mo_space(int nirreps_ref, int& n, intvec& mo, std::string labels);
    void print_mo_space(int nmo, intvec& mo, std::string labels);
    intvec convert_int_array_to_vector(int n, const int* array);

    Wavefunction& ref_wfn;
    Options& options;
    const int nirreps; //The # of irreps this object has been constructed with
    int wfn_sym;
    const int charge; //The charge this object has been constructed with
    int multiplicity;

    int nso;  // PSI nso (# of symmetry-adapted atomic orbitals)
    int nmo;  // Psi nmo (# of molecular orbitals, including frozen core and frozen virtual)
    int ndocc;
    int nactv; //Total # of active orbitals across all irreps
    int nael; // The # of alpha electrons (including frozen)
    int nbel; // The # of beta electrons (including frozen)
    int nactive_ael;
    int nactive_bel;

    intvec sopi; //Array holding the numbers of SOs per irrep
    intvec docc; //Array holding the numbers of doubly occupied orbitals (DOCC) per irrep
    intvec actv; //Array holding the numbers of active orbitals per irrep
    bool guess_occupation;
    bool silent;

    double nuclear_energy; //The nuclear repulsion energy

    double** scf;         // MO coefficients
    double*** scf_irrep;  // MO coefficients

    std::vector<std::string> irr_labs;
};

}  // namespace psi

#endif  // _psi_src_lib_libmoinfo_moinfo_base_h_
