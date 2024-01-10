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

#ifndef CDSALCLIST_H
#define CDSALCLIST_H

#include <cstdio>
#include <vector>
#include "psi4/pragma.h"
#include <memory>
#include "psi4/libmints/typedefs.h"
#include "psi4/libpsi4util/exception.h"
#include "typedefs.h"

namespace psi {

class MatrixFactory;

class PSI_API CdSalc {
   public:
    class PSI_API Component {
       public:
        double coef;
        int atom;
        int xyz;

        Component(double coef_, int atom_, int xyz_) : coef(coef_), atom(atom_), xyz(xyz_) {}
    };

   private:
    /// All the components needed for this transformation
    std::vector<Component> components_;

    /// Symmetry of this SALC compatible for ^ (A1 = 0)
    char irrep_;

   public:
    CdSalc(char irrep) : irrep_(irrep) {}

    void add(double coef, int atom, int xyz) { components_.push_back(Component(coef, atom, xyz)); }

    size_t ncomponent() const { return components_.size(); }

    // made const function to access coef,atom,xyz through CdSalc[i]
    const Component& component(int com) const { return components_[com]; }
    const std::vector<Component>& get_components() const { return components_; }

    char irrep() const { return irrep_; }

    void print() const;
};

class PSI_API CdSalcWRTAtom {
   public:
    class PSI_API Component {
       public:
        double coef;
        int irrep;
        int salc;

        Component(double coef_, int irrep_, int salc_) : coef(coef_), irrep(irrep_), salc(salc_) {}
    };

   private:
    /// We split the components up for simplicity
    std::vector<Component> x_;
    std::vector<Component> y_;
    std::vector<Component> z_;

   public:
    CdSalcWRTAtom() {}

    void add(int xyz, double coef, int irrep, int salc) {
        if (xyz == 0)
            x_.push_back(Component(coef, irrep, salc));
        else if (xyz == 1)
            y_.push_back(Component(coef, irrep, salc));
        else if (xyz == 2)
            z_.push_back(Component(coef, irrep, salc));
    }

    int nx() const { return x_.size(); }
    int ny() const { return y_.size(); }
    int nz() const { return z_.size(); }

    const Component& x(int n) const { return x_[n]; }
    const Component& y(int n) const { return y_[n]; }
    const Component& z(int n) const { return z_[n]; }

    void print() const;
};

class PSI_API CdSalcList {
    SharedMolecule molecule_;

    char needed_irreps_;
    bool project_out_translations_;
    bool project_out_rotations_;

    int ncd_;
    int cdsalcpi_[8];
    int nirrep_;

    /// Vector of all requested SALCs
    std::vector<CdSalc> salcs_;
    std::vector<CdSalcWRTAtom> atom_salcs_;

   public:
    CdSalcList() { throw PSIEXCEPTION("CdSalcList(): U R STOOOPID!!!"); }

    /*! Constructor for generating Cartesian displacement symmetry adapted
     *  linear combinations.
     *
     *  For a gradient calculation, you only need totally symmetric SALCs
     *  (needed_irreps=1).
     *
     *  The project_* ability is not coded, yet.
     *
     *  \param mol Molecule to form SALCs for
     *  \param needed_irreps bit representation for desired irreps (default all)
     *  \param project_out_translations Project out translational SALCs
     *  \param project_out_rotations Project out rotational SALCs
     */
    CdSalcList(SharedMolecule mol, int needed_irreps = 0xFF, bool project_out_translations = true,
               bool project_out_rotations = true);
    ~CdSalcList();

    /*! Returns the number of SALCs. It may not be 3n-5 or 3n-6. The value
     *  returned depends on needed_irreps and the project_out* settings.
     */
    size_t ncd() const { return salcs_.size(); }

    std::vector<SharedMatrix> create_matrices(const std::string& basename, const MatrixFactory& factory) const;
    std::string salc_name(int index) const;

    char needed_irreps() const { return needed_irreps_; }
    int nirrep() const { return nirrep_; }
    int cdsalcpi(int h) const { return cdsalcpi_[h]; }
    bool project_out_translations() const { return project_out_translations_; }
    bool project_out_rotations() const { return project_out_rotations_; }

    const CdSalc& operator[](int i) const { return salcs_[i]; }

    const CdSalcWRTAtom& atom_salc(int i) const { return atom_salcs_[i]; }

    const std::vector<CdSalc>& get_salcs() const { return salcs_; }

    SharedMatrix matrix() const;
    SharedMatrix matrix_irrep(int h) const;  // return only salcs of a given irrep
    // SharedMatrix matrix_projected_out() const;

    void print() const;
};

}  // namespace psi

#endif  // CDSALCLIST_H
