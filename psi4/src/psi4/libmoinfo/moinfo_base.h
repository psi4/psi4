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

#ifndef _psi_src_lib_libmoinfo_moinfo_base_h_
#define _psi_src_lib_libmoinfo_moinfo_base_h_

/*! \file    moinfo_base.h
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding MOs
*/

#include <string>
#include <vector>

namespace psi {

class Options;
class Wavefunction;
using intvec = std::vector<int>;
using boolvec = std::vector<bool>;

/// @brief Stores all the basic info regarding MOs, at least for the parts of Psi4 that are using it.
/// Mostly used for deriving other classes.
class MOInfoBase {
   public:
    MOInfoBase(Wavefunction& ref_wfn_, Options& options_);

    /// @brief Get the nuclear energy stored in an MOInfoBase object (or derived object).
    /// @return The nuclear repulsion energy
    double get_nuc_E() const { return nuclear_energy; }

    /// @brief Get a const ref. to one of the irrep labels that are stored in an MOInfoBase object (or derived object).
    /// Not bounds-checked!
    /// @param i : Index of the irrep label
    /// @return Const& of the selected irrep label
    const std::string& get_irr_lab(size_t i) const { return irr_labs[i]; }

    /// @brief Get the nirreps value (# of irreducible representations) that is stored in an MOInfoBase object (or
    /// derived object).
    /// @return The # of irreps
    int get_nirreps() const { return (nirreps); }

    /// @brief Get the nso value (# of symmetry-adapted atomic orbitals) that is stored in an MOInfoBase object (or
    /// derived object).
    /// @return The # of symmetry-adapted atomic orbitals
    int get_nso() const { return (nso); }

    /// @brief Get a const ref. to the array holding the numbers of SOs per irrep, from an MOInfoBase object (or derived
    /// object)
    /// @return Const& of the SOs per irrep array
    const intvec& get_sopi() const { return sopi; }

    /// @brief Get a const ref. to the array holding the numbers of doubly occupied orbitals (DOCC) per irrep, from an
    /// MOInfoBase object (or derived object).
    /// @return Const& of the array holding the numbers of doubly occupied orbitals (DOCC) per irrep
    const intvec& get_docc() const { return docc; }

    /// @brief Get a const ref. to the array holding the numbers of active orbitals per irrep, from an MOInfoBase object
    /// (or derived object).
    /// @return Const& of the array holding the numbers of active orbitals per irrep
    const intvec& get_actv() const { return actv; }
    bool get_guess_occupation() const { return (guess_occupation); }

    /// @brief Get the total # of active orbitals that is stored in an MOInfoBase object (or derived object).
    /// @return Total # of active orbitals across all irreps
    int get_nactv() const { return (nactv); }

    /// @brief Get the # of alpha electrons (including frozen) that is stored in an MOInfoBase object (or derived
    /// object).
    /// @return The # of alpha electrons (including frozen)
    int get_nael() const { return (nael); }

    /// @brief Get the # of beta electrons (including frozen) that is stored in an MOInfoBase object (or derived
    /// object).
    /// @return The # of beta electrons (including frozen)
    int get_nbel() const { return (nbel); }

   protected:
    /// @brief Get a const ref. to an element of the array holding the numbers of SOs per irrep, from an MOInfoBase
    /// object (or derived object). Not bounds-checked!
    /// @return Const& of the selected number of SOs per irrep
    const int& get_sopi(size_t i) const { return sopi[i]; }

    void read_data();
    void compute_number_of_electrons();
    void read_mo_space(const int nirreps_ref, int& n, intvec& mo, const std::string& labels);
    void print_mo_space(int nmo, const intvec& mo, const std::string& labels);

    Wavefunction& ref_wfn;
    Options& options;
    int nirreps;
    int wfn_sym;
    int multiplicity;

    int ndocc;
    int nactv;  // Total # of active orbitals across all irreps
    int nael;   // The # of alpha electrons (including frozen)
    int nbel;   // The # of beta electrons (including frozen)
    int nactive_ael;
    int nactive_bel;

    intvec docc;  // Array holding the numbers of doubly occupied orbitals (DOCC) per irrep
    intvec actv;  // Array holding the numbers of active orbitals per irrep
    bool guess_occupation;

    double*** scf_irrep;  // MO coefficients

   private:
    double nuclear_energy;  // The nuclear repulsion energy
    const int charge;       // The molecular charge this object has been constructed with

    std::vector<std::string> irr_labs;  // Array of the irrep labels
    int nso;                            // The # of symmetry-adapted atomic orbitals
    intvec sopi;                        // Array holding the numbers of SOs per irrep
};

}  // namespace psi

#endif  // _psi_src_lib_libmoinfo_moinfo_base_h_
