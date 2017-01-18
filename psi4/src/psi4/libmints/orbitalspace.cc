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
#include "psi4/libqt/qt.h"
#include "psi4/libmints/view.h"
#include "psi4/libmints/orbitalspace.h"
#include "psi4/libmints/orthog.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/sobasis.h"

#include <tuple>

namespace psi {

OrbitalSpace::OrbitalSpace(const std::string &id,
                           const std::string &name,
                           const SharedMatrix &full_C,
                           const std::shared_ptr<Vector> &evals,
                           const std::shared_ptr<BasisSet> &basis,
                           const std::shared_ptr<IntegralFactory> &ints)
        : id_(id),
          name_(name),
          C_(full_C),
          evals_(evals),
          basis_(basis),
          ints_(ints),
          dim_(full_C->colspi())
{
}

OrbitalSpace::OrbitalSpace(const std::string &id,
                           const std::string &name,
                           const SharedMatrix &full_C,
                           const std::shared_ptr<BasisSet> &basis,
                           const std::shared_ptr<IntegralFactory> &ints)
        : id_(id),
          name_(name),
          C_(full_C),
          basis_(basis),
          ints_(ints),
          dim_(full_C->colspi())
{
}

OrbitalSpace::OrbitalSpace(const std::string &id,
                           const std::string &name,
                           const std::shared_ptr<Wavefunction> &wave)
        : id_(id),
          name_(name),
          C_(wave->Ca()),
          evals_(wave->epsilon_a()),
          basis_(wave->basisset()),
          ints_(wave->integral()),
          dim_(wave->Ca()->colspi())
{
}

int OrbitalSpace::nirrep() const
{
    return C_->nirrep();
}

const std::string &OrbitalSpace::id() const
{
    return id_;
}

const std::string &OrbitalSpace::name() const
{
    return name_;
}

const SharedMatrix &OrbitalSpace::C() const
{
    return C_;
}

const std::shared_ptr<Vector> &OrbitalSpace::evals() const
{
    return evals_;
}

const std::shared_ptr<BasisSet> &OrbitalSpace::basisset() const
{
    return basis_;
}

const std::shared_ptr<IntegralFactory> &OrbitalSpace::integral() const
{
    return ints_;
}

const Dimension &OrbitalSpace::dim() const
{
    return dim_;
}

OrbitalSpace OrbitalSpace::transform(const OrbitalSpace &A, const std::shared_ptr<BasisSet> &B)
{
    SharedMatrix SBA = overlap(B, A.basisset());
    SBA->set_name("Sba");
    SharedMatrix SBB = overlap(B, B);
    SBB->set_name("SBB");

    // Follows Werner's method from Mol. Phys. 102, 21-22, 2311
    // just like HF::dualBasisProjection

    // 1. Invert SBB
    SBB->invert();
    SBB->set_name("SBB^-1");

    // 2. Form T
    SharedMatrix I = Matrix::create("I = SAB SBB SBA", SBA->colspi(), SBA->colspi());
    I->transform(SBB, SBA);

    SharedMatrix T = Matrix::create("T", A.dim(), A.dim());
    T->transform(I, A.C());
    I.reset(); // release memory

    // 3. Form T^{-1/2}
    T->power(-0.5);

    // 4. Cb = [Sbb]^-1 Sba] Ca T^{-1/2}
    // 4a. Ca T^{-1/2}
    SharedMatrix CaT = Matrix::create("Ca*T^{-1/2}", A.C()->rowspi(), A.C()->colspi());
    CaT->gemm(false, false, 1.0, A.C(), T, 0.0);

    // 4b. Sba * 4a
    SharedMatrix SbaCaT = Matrix::create("SbaCaT", SBB->rowspi(), A.C()->colspi());
    SbaCaT->gemm(false, false, 1.0, SBA, CaT, 0.0);

    // 4c. [Sbb]^-1 * 4b
    SharedMatrix Cb = Matrix::create("Cb", SBB->rowspi(), A.C()->colspi());
    Cb->gemm(false, false, 1.0, SBB, SbaCaT, 0.0);

    std::shared_ptr<IntegralFactory> i(new IntegralFactory(B, B, B, B));

    return OrbitalSpace("p",
                        "Ca transformed into Cb",
                        Cb,
                        A.evals(),
                        B,
                        i);
}

SharedMatrix OrbitalSpace::overlap(const OrbitalSpace &space1, const OrbitalSpace &space2)
{
    IntegralFactory mix_ints(space1.basisset(), space2.basisset());

    PetiteList p1(space1.basisset(), space1.integral());
    PetiteList p2(space2.basisset(), space2.integral());

    SharedMatrix Smat(new Matrix("Overlap between space1 and space2",
                                 p1.SO_basisdim(), p2.SO_basisdim()));

    OneBodySOInt *S = mix_ints.so_overlap();
    S->compute(Smat);
    delete S;

    return Smat;
}

SharedMatrix OrbitalSpace::overlap(const std::shared_ptr<BasisSet> &basis1,
                                   const std::shared_ptr<BasisSet> &basis2)
{
    IntegralFactory mix_ints(basis1, basis2);
    SOBasisSet sobasis1(basis1, &mix_ints);
    SOBasisSet sobasis2(basis2, &mix_ints);

    SharedMatrix Smat(new Matrix("Overlap between space1 and space2",
                                 sobasis1.dimension(), sobasis2.dimension()));

    OneBodySOInt *S = mix_ints.so_overlap();
    S->compute(Smat);
    delete S;

    return Smat;
}

void OrbitalSpace::print() const
{
    outfile->Printf("    Orbital space %s (%s)\n", name_.c_str(), id_.c_str());
    outfile->Printf("        Basis: %s\n", basis_->name().c_str());
    basis_->print_summary();
    outfile->Printf("        Dimensions: ");
    dim_.print();
//    outfile->Printf( "        Transformation matrix:\n");
//    C_->print();
}

namespace { // anonymous
OrbitalSpace orthogonalize(const std::string &id, const std::string &name,
                           const std::shared_ptr<BasisSet> &bs,
                           double lindep_tol)
{
    outfile->Printf("    Orthogonalizing basis for space %s.\n", name.c_str());

    SharedMatrix overlap = OrbitalSpace::overlap(bs, bs);
    Dimension SODIM = overlap->rowspi();
    SharedMatrix evecs = Matrix::create("evecs", SODIM, SODIM);
    SharedMatrix sqrtm = Matrix::create("evecs", SODIM, SODIM);
    SharedVector evals = Vector::create("evals", SODIM);

    int nlindep = 0;
    overlap->diagonalize(evecs, evals);
    for (int h = 0; h < SODIM.n(); h++) {
        for (int i = 0; i < SODIM[h]; i++) {
            if (fabs(evals->get(h, i)) > lindep_tol) {
                sqrtm->set(h, i, i, 1.0 / sqrt(evals->get(h, i)));
            }
            else {
                sqrtm->set(h, i, i, 0.0);
                nlindep++;
            }
        }
    }

    sqrtm->back_transform(evecs);

    outfile->Printf("    %d linear dependencies will be \'removed\'.\n", nlindep);

    std::shared_ptr<IntegralFactory> localfactory(new IntegralFactory(bs));
    return OrbitalSpace(id, name, sqrtm, bs, localfactory);
}

OrbitalSpace orthogonal_compliment(const OrbitalSpace &space1, const OrbitalSpace &space2,
                                   const std::string &id, const std::string &name, const double &lindep_tol)
{
    outfile->Printf("    Projecting out '%s' from '%s' to obtain space '%s'\n",
                    space1.name().c_str(), space2.name().c_str(), name.c_str());

    // If space1 is empty, return a copy of the original space.
    if (space1.dim().sum() == 0) {
        outfile->Printf("    '%s' space is empty. Nothing to project out.\n", space1.name().c_str());
        return OrbitalSpace(id, name, space2.C(), space2.evals(), space2.basisset(), space2.integral());
    }

    // O12 = O12
    SharedMatrix O12 = OrbitalSpace::overlap(space1, space2);
    SharedMatrix C12 = Matrix::create("C12", space1.C()->colspi(), space2.C()->colspi());

    // C12 = C1t * S12 * C2
    // TODO: Try using the svd_a function of Matrix, however it calls dgesdd not dgesvd.
    C12->transform(space1.C(), O12, space2.C());
    C12->print();

#if SVD
    std::tuple<SharedMatrix, SharedMatrix, SharedMatrix> svd_temps = C12->svd_a_temps();

#else
    // We're interested in the right side vectors (V) of an SVD solution.
    // We don't need a full SVD solution just part of it. Do it by hand:
    SharedMatrix D11 = Matrix::create("D11", C12->colspi(), C12->colspi());
    D11->gemm(true, false, 1.0, C12, C12, 0.0);
//        D11->print();

    SharedMatrix V11 = Matrix::create("V11", D11->rowspi(), D11->colspi());
    SharedVector E1 = Vector::create("E1", D11->colspi());
    D11->diagonalize(V11, E1);
//        V11->eivprint(E1);

    // Count the number of eigenvalues < lindep_tol
    Dimension zeros(space1.nirrep());
    for (int h = 0; h < space1.nirrep(); ++h) {
        for (int i = 0; i < E1->dimpi()[h]; ++i) {
            if (E1->get(h, i) < lindep_tol)
                zeros[h]++;
        }
    }

    outfile->Printf("        Orbital space before projecting out: ");
    space2.dim().print();
    outfile->Printf("        Orbital space after projecting out:  ");
    zeros.print();
    outfile->Printf("\n");

    // Pull out the nullspace vectors
    View nullspace(V11, V11->rowspi(), zeros);
    SharedMatrix V = nullspace();
//        V->print();

    // Half-back transform to space2
    SharedMatrix newC = Matrix::create("Transformation matrix", space2.C()->rowspi(), zeros);
    newC->gemm(false, false, 1.0, space2.C(), V, 0.0);
//        newC->print();

    return OrbitalSpace(id, name, newC, space2.basisset(), space2.integral());
#endif
}
} // namespace anonymous

OrbitalSpace OrbitalSpace::build_cabs_space(const OrbitalSpace &orb_space, const OrbitalSpace &ri_space,
                                            double lindep_tol)
{
    return orthogonal_compliment(orb_space, ri_space, "p''", "CABS", lindep_tol);
}

OrbitalSpace OrbitalSpace::build_ri_space(const std::shared_ptr<Molecule> &molecule, const std::string &obs_key,
                                          const std::string &aux_key, double lindep_tol)
{
    // Construct a combined basis set.
    Options &options = Process::environment.options;
    std::vector <std::string> keys, targets, roles, others;
    keys.push_back(obs_key);
    keys.push_back(aux_key);
    targets.push_back(options.get_str(obs_key));
    targets.push_back(options.get_str(aux_key));
    roles.push_back(obs_key);
    roles.push_back("F12");
    others.push_back(options.get_str(obs_key));
    others.push_back(options.get_str(obs_key));
    throw PSIEXCEPTION("build_ri_space has not been updated to the new python based basis set construction scheme.");
    //std::shared_ptr<BasisSet> combined = BasisSet::pyconstruct_combined(molecule, keys, targets, roles, others);
    std::shared_ptr<BasisSet> combined = BasisSet::zero_ao_basis_set();

    // orthogonalize the basis set projecting out linear dependencies.
    return orthogonalize("p'", "RIBS", combined, lindep_tol);
}

} // namespace psi
