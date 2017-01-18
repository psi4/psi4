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

#include "psi4/psi4-dec.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/cdsalclist.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libqt/qt.h"

#include <algorithm>

;

namespace {

static char direction(int xyz)
{
    switch (xyz) {
    case 0:
        return 'x';
    case 1:
        return 'y';
    case 2:
        return 'z';
    default:
        return '?';
    }
}

}

namespace psi {

void CdSalc::print() const
{
    outfile->Printf( "\tirrep = %d, ncomponent = %ld\n", irrep_, ncomponent());
    for (size_t i=0; i<ncomponent(); ++i) {
        outfile->Printf( "\t\t%d: atom %d, direction %c, coef %lf\n",
                i,
                components_[i].atom,
                direction(components_[i].xyz),
                components_[i].coef);
    }
}

void CdSalcWRTAtom::print() const
{
    outfile->Printf( "\tx component, size = %ld\n", x_.size());
    for (size_t i=0; i<x_.size(); ++i) {
        outfile->Printf( "\t\t%d: salc %d, irrep %d, coef %lf\n",
                i,
                x_[i].salc,
                x_[i].irrep,
                x_[i].coef);
    }

    outfile->Printf( "\ty component, size = %ld\n", y_.size());
    for (size_t i=0; i<y_.size(); ++i) {
        outfile->Printf( "\t\t%d: salc %d, irrep %d, coef %lf\n",
                i,
                y_[i].salc,
                y_[i].irrep,
                y_[i].coef);
    }

    outfile->Printf( "\tz component, size = %ld\n", z_.size());
    for (size_t i=0; i<z_.size(); ++i) {
        outfile->Printf( "\t\t%d: salc %d, irrep %d, coef %lf\n",
                i,
                z_[i].salc,
                z_[i].irrep,
                z_[i].coef);
    }
}

CdSalcList::CdSalcList(std::shared_ptr<Molecule> mol,
                       std::shared_ptr<MatrixFactory> fact,
                       int needed_irreps,
                       bool project_out_translations,
                       bool project_out_rotations)
    : molecule_(mol), factory_(fact), needed_irreps_(needed_irreps),
      project_out_translations_(project_out_translations),
      project_out_rotations_(project_out_rotations)
{
    // Ensure point group has been set.
    if (!molecule_->point_group()) {
        throw PSIEXCEPTION("CdSalcList::CdSalcList: Molecule point group has not been set.");
    }

    int natom = molecule_->natom();
    ncd_ = 3 * natom;

    // Immediately create the rotation and translation vectors to be projected out.
    Matrix constraints("COM & Rotational Constraints", 6, ncd_);

    SharedMatrix pI(molecule_->inertia_tensor());
    Vector ev(3);
    Matrix X(3, 3);
    pI->diagonalize(X, ev);

//    outfile->Printf( "pI[0][1] = %20.14lf\n", pI->get(0, 1));
//    molecule_->inertia_tensor()->print();
//    X.eivprint(ev);

    // Pull out data to local variables to reduce memory lookup
    double X00 = X(0, 0), X01 = X(0, 1), X02 = X(0, 2);
    double X10 = X(1, 0), X11 = X(1, 1), X12 = X(1, 2);
    double X20 = X(2, 0), X21 = X(2, 1), X22 = X(2, 2);

//    Matrix constraints("COM & Rotational Constraints", 6, 3*natom);
    double tval0, tval1, tval2;
    for (int i=0; i < natom; ++i) {
        // Local lookups
        double atomx = molecule_->x(i);
        double atomy = molecule_->y(i);
        double atomz = molecule_->z(i);
        double smass = sqrt(molecule_->mass(i));

        // COM constraints
        if (project_out_translations_) {
            constraints(0, 3*i+0) = smass;
            constraints(1, 3*i+1) = smass;
            constraints(2, 3*i+2) = smass;
        }

        // Rotational constraints
        if (project_out_rotations_) {
            tval0 = (atomx * X00) + (atomy * X10) + (atomz * X20);
            tval1 = (atomx * X01) + (atomy * X11) + (atomz * X21);
            tval2 = (atomx * X02) + (atomy * X12) + (atomz * X22);

            constraints(3, 3*i+0) = (tval1 * X02 - tval2 * X01) * smass;
            constraints(3, 3*i+1) = (tval1 * X12 - tval2 * X11) * smass;
            constraints(3, 3*i+2) = (tval1 * X22 - tval2 * X21) * smass;

            constraints(4, 3*i+0) = (tval2 * X00 - tval0 * X02) * smass;
            constraints(4, 3*i+1) = (tval2 * X10 - tval0 * X12) * smass;
            constraints(4, 3*i+2) = (tval2 * X20 - tval0 * X22) * smass;

            constraints(5, 3*i+0) = (tval0 * X01 - tval1 * X00) * smass;
            constraints(5, 3*i+1) = (tval0 * X11 - tval1 * X10) * smass;
            constraints(5, 3*i+2) = (tval0 * X21 - tval1 * X20) * smass;
        }
    }

    //constraints.print();

    // Remove NULL constraint (if present) and normalize the rest of them
    for (int i=0; i<6; ++i) {
        double normval = C_DDOT(ncd_, constraints[0][i], 1, constraints[0][i], 1);
        if (normval > 1.0E-10)
            constraints.scale_row(0, i, 1.0 / sqrt(normval));
        else
            constraints.scale_row(0, i, 0.0);
    }

    //constraints.print();

    Matrix constraints_ortho("Orthogonalized COM & Rotational constraints", 6, 3*natom);
    // Ensure rotations and translations are exactly orthogonal
    int count = 0;
    for (int i=0; i<6; ++i)
        count += constraints_ortho.schmidt_add_row(0, i, constraints[0][i]);

    //constraints_ortho.print();

    double *salc = new double[ncd_];

    // Obtain handy reference to point group.
    PointGroup& pg = *molecule_->point_group().get();
    CharacterTable char_table = pg.char_table();
    nirrep_ = char_table.nirrep();

    // I don't know up front how many I have per irrep
    // I'll resize afterwards
    Dimension dim(nirrep_);
    for (int i=0; i<nirrep_; ++i)
        dim[i] = ncd_;

    Matrix salcs("Requested SALCs", dim, dim);
    int *salcirrep = new int[ncd_];

    // We know how many atom_salcs_ we have.
    for (int i=0; i<natom; ++i)
        atom_salcs_.push_back(CdSalcWRTAtom());

    // Obtain atom mapping of atom * symm op to atom
    int **atom_map = compute_atom_map(molecule_);
    memset(cdsalcpi_, 0, sizeof(int)*8);

    int nsalc = 0;
    for (int uatom=0; uatom < molecule_->nunique(); ++uatom) {
        int atom = molecule_->unique(uatom);

        // Project each displacement
        for (int xyz=0; xyz<3; ++xyz) {

            // on each irrep
            for (int irrep=0; irrep<nirrep_; ++irrep) {
                IrreducibleRepresentation gamma = char_table.gamma(irrep);
                memset(salc, 0, sizeof(double)*ncd_);

                // This is the order of the atom stabilizer
                // ...how many times the symmetry operation keeps the atom the same
                int stab_order = 0;

                // Apply the projection
                for (int G=0; G<nirrep_; ++G) {
                    SymmetryOperation so = char_table.symm_operation(G);
                    int Gatom = atom_map[atom][G];
                    if (Gatom == atom)
                        ++stab_order;

                    // compute position in the salc
                    int Gcd = 3*Gatom + xyz;

                    // so(xyz, xyz) tells us how x, y, or z transforms in this
                    // symmetry operation, then we multiply by the character of the
                    // operation in the irrep
                    double coeff = so(xyz, xyz) * gamma.character(G);

                    // Add this contribution to the salc.
                    salc[Gcd] += coeff;
                }

                if (stab_order == 0)
                    throw PSIEXCEPTION("CdSalcList::CdSalcList: Stabilizer order is 0 this is not possible.");

                int nonzero=0;
                for (int cd=0; cd<ncd_; ++cd) {
                    // Normalize the salc
                    salc[cd] /= sqrt((double)nirrep_*stab_order);

                    // Count number of nonzeros
                    if (fabs(salc[cd]) > 1e-10)
                        ++nonzero;
                }

                // We're only interested in doing the following if there are nonzeros
                // AND the irrep that we're on is what the user wants.
                //if (nonzero && (1 << irrep) & needed_irreps) {
                if (nonzero && (1 << irrep) & needed_irreps) {
                    // Store the salc so we can project out constraints below
                    salcs.copy_to_row(irrep, cdsalcpi_[irrep], salc);
                    salcirrep[nsalc] = irrep;
                    ++cdsalcpi_[irrep];
                    nsalc++;
                }
            }
        }
    }

    // Raw - non-projected cartesian displacements
    //salcs.print("outfile", "Raw, Nonprojected Cartesian Displacements");

    // Project out any constraints
    salcs.project_out(constraints_ortho);
    salcs.set_name("Resulting SALCs after projections");
    //salcs.print();

    // Walk through the new salcs and populate our sparse vectors.
    int new_cdsalcpi[8];
    memset(new_cdsalcpi, 0, sizeof(int)*8);
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<cdsalcpi_[h]; ++i) {
            bool added = false;
            CdSalc new_salc(h);
            for (int cd=0; cd < ncd_; ++cd) {
                if (fabs(salcs(h, i, cd)) > 1.0e-10) {
                    added = true;
                    new_salc.add(salcs(h, i, cd), cd/3, cd % 3);
                    atom_salcs_[cd/3].add(cd % 3, salcs(h, i, cd), h, salcs_.size());
                }
            }
            if (added) {
                salcs_.push_back(new_salc);
                new_cdsalcpi[h]++;
            }
        }
    }
    ncd_ = salcs_.size();
    memcpy(cdsalcpi_, new_cdsalcpi, sizeof(int)*8);

    // Free memory.
    delete[] salcirrep;
    delete[] salc;
    delete_atom_map(atom_map, molecule_);
}

CdSalcList::~CdSalcList()
{
}

std::vector<SharedMatrix > CdSalcList::create_matrices(const std::string &basename)
{
    std::vector<SharedMatrix > matrices;
    std::string name;

    for (size_t i=0; i<salcs_.size(); ++i) {
        name = basename + " " + name_of_component(i);
        matrices.push_back(factory_->create_shared_matrix(name, salcs_[i].irrep()));
    }

    return matrices;
}

std::string CdSalcList::name_of_component(int component)
{
    std::string name;
    CdSalc& salc = salcs_[component];

    for (size_t i=0; i<salc.ncomponent(); ++i) {
        const CdSalc::Component& com = salc.component(i);

        name += com.coef > 0.0 ? "+" : "-";
        name += to_string(fabs(com.coef)) + " ";
        name += molecule_->label(com.atom);
        if (com.xyz == 0)
            name += "-x";
        else if (com.xyz == 1)
            name += "-y";
        else if (com.xyz == 2)
            name += "-z";
        name += " ";
    }

    return name;
}

SharedMatrix CdSalcList::matrix()
{
    SharedMatrix temp(new Matrix("Cartesian/SALC transformation", ncd(), 3*molecule_->natom()));

    for (size_t i=0; i<ncd(); ++i) {
        int nc = salcs_[i].ncomponent();
        for (int c=0; c<nc; ++c) {
            int a       = salcs_[i].component(c).atom;
            int xyz     = salcs_[i].component(c).xyz;
            double coef = salcs_[i].component(c).coef;
            temp->set(i, 3*a+xyz, coef);
        }
    }

    return temp;
}

SharedMatrix CdSalcList::matrix_irrep(int h)
{
    // cdsalcpi_ does not get updated after projected out translation and rotations
    // why?  if it ever is, I can use it.
    //SharedMatrix temp(new Matrix("Cartesian/SALC transformation", cdsalcpi_[h], 3*molecule_->natom()));

    int cnt = 0;
    for (size_t i=0; i<ncd(); ++i)
        if (salcs_[i].irrep() == h) ++cnt;

    SharedMatrix temp(new Matrix("Cartesian/SALC transformation", cnt, 3*molecule_->natom()));

    cnt = 0;
    for (size_t i=0; i<ncd(); ++i) {
        if (salcs_[i].irrep() == h) {
            int nc = salcs_[i].ncomponent();
            for (int c=0; c<nc; ++c) {
                int a       = salcs_[i].component(c).atom;
                int xyz     = salcs_[i].component(c).xyz;
                double coef = salcs_[i].component(c).coef;
                temp->set(cnt, 3*a+xyz, coef);
            }
            ++cnt;
        }
    }

    return temp;
}

void CdSalcList::print() const
{
    outfile->Printf( "  Cartesian Displacement SALCs\n  By SALC:\n");
    outfile->Printf( "  Number of SALCs: %ld, nirreps: %d\n"
            "  Project out translations: %s\n"
            "  Project out rotations: %s\n",
            salcs_.size(), needed_irreps_,
            project_out_translations_ ? "True" : "False",
            project_out_rotations_ ? "True" : "False");

    for (size_t i=0; i<salcs_.size(); ++i)
        salcs_[i].print();

    outfile->Printf( "\n  By Atomic Center:\n");
    outfile->Printf( "  Number of atomic centers: %ld\n", atom_salcs_.size());
    for (size_t i=0; i<atom_salcs_.size(); ++i) {
        outfile->Printf( "   Atomic Center %d:\n", i);
        atom_salcs_[i].print();
    }
    outfile->Printf( "\n");
}

// Generate and return those translations or rotations that
// have been projected out
/*
SharedMatrix CdSalcList::matrix_projected_out() const
{
    int natom = molecule_->natom();

    Matrix constraints("COM & Rotational Constraints", 6, 3*natom);

    SharedMatrix pI(molecule_->inertia_tensor());
    Vector ev(3);
    Matrix X(3, 3);
    pI->diagonalize(X, ev);

molecule_->inertia_tensor()->print();

    // Pull out data to local variables to reduce memory lookup
    double X00 = X(0, 0), X01 = X(0, 1), X02 = X(0, 2);
    double X10 = X(1, 0), X11 = X(1, 1), X12 = X(1, 2);
    double X20 = X(2, 0), X21 = X(2, 1), X22 = X(2, 2);

    double tval0, tval1, tval2;
    for (int i=0; i < natom; ++i) {
        // Local lookups
        double atomx = molecule_->x(i);
        double atomy = molecule_->y(i);
        double atomz = molecule_->z(i);
        double smass = sqrt(molecule_->mass(i));

        // COM constraints
        if (project_out_translations_) {
            constraints(0, 3*i+0) = smass;
            constraints(1, 3*i+1) = smass;
            constraints(2, 3*i+2) = smass;
        }

        // Rotational constraints
        if (project_out_rotations_) {
            tval0 = (atomx * X00) + (atomy * X10) + (atomz * X20);
            tval1 = (atomx * X01) + (atomy * X11) + (atomz * X21);
            tval2 = (atomx * X02) + (atomy * X12) + (atomz * X22);

            constraints(3, 3*i+0) = (tval1 * X02 - tval2 * X01) * smass;
            constraints(3, 3*i+1) = (tval1 * X12 - tval2 * X11) * smass;
            constraints(3, 3*i+2) = (tval1 * X22 - tval2 * X21) * smass;

            constraints(4, 3*i+0) = (tval2 * X00 - tval0 * X02) * smass;
            constraints(4, 3*i+1) = (tval2 * X10 - tval0 * X12) * smass;
            constraints(4, 3*i+2) = (tval2 * X20 - tval0 * X22) * smass;

            constraints(5, 3*i+0) = (tval0 * X01 - tval1 * X00) * smass;
            constraints(5, 3*i+1) = (tval0 * X11 - tval1 * X10) * smass;
            constraints(5, 3*i+2) = (tval0 * X21 - tval1 * X20) * smass;
        }
    }

    // Remove NULL constraint (if present) and normalize the rest of them
    std::vector<int> non_zero;
    for (int i=0; i<6; ++i) {
        double normval = C_DDOT(3*natom, constraints[0][i], 1, constraints[0][i], 1);
        if (normval > 1.0E-10) {
            constraints.scale_row(0, i, 1.0 / sqrt(normval));
            non_zero.push_back(i);
        }
        else
            constraints.scale_row(0, i, 0.0);
    }

    // modifying the above code slightly to only allocate the necessary number of rows
    SharedMatrix constraints_ortho(new Matrix("Orthogonalized COM & Rotational constraints",
        non_zero.size(), 3*natom));

    for (int i=0; i<non_zero.size(); ++i)
        constraints_ortho->schmidt_add(0, i, constraints[0][non_zero[i]]);

    //constraints_ortho->print();

    return constraints_ortho;
}
*/


}
