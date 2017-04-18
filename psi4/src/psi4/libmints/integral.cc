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
#include "psi4/libmints/integral.h"
#include "psi4/libmints/shellrotation.h"
#include "psi4/libmints/cartesianiter.h"
#include "psi4/libmints/rel_potential.h"
#include "psi4/libmints/electricfield.h"
#include "psi4/libmints/tracelessquadrupole.h"
#include "psi4/libmints/efpmultipolepotential.h"
#include "psi4/libmints/eri.h"
#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/quadrupole.h"
#include "psi4/libmints/angularmomentum.h"
#include "psi4/libmints/nabla.h"
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/pseudospectral.h"
#include "psi4/libmints/kinetic.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libmints/overlap.h"
#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/potentialint.h"
#include "psi4/libmints/ecpint.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/erd_eri.h"

#ifdef USING_simint
#include "psi4/libmints/siminteri.h"
#endif

#include <libint/libint.h>

;
using namespace psi;

IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2,
                                 std::shared_ptr<BasisSet> bs3,
                                 std::shared_ptr<BasisSet> bs4,
                                 std::shared_ptr<BasisSet> bsecp)
{
    set_basis(bs1, bs2, bs3, bs4, bsecp);
}
IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2,
                                 std::shared_ptr<BasisSet> bs3,
                                 std::shared_ptr<BasisSet> bs4):
    bs_ecp_(nullptr)
{

    set_basis(bs1, bs2, bs3, bs4, bs_ecp_);
}


IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bsecp)
{
    set_basis(bs1, bs1, bs1, bs1, bsecp);
}

IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1):
    bs_ecp_(nullptr)
{
    set_basis(bs1, bs1, bs1, bs1, bs_ecp_);
}

IntegralFactory::~IntegralFactory()
{

}

std::shared_ptr<BasisSet> IntegralFactory::basis1() const
{
    return bs1_;
}

std::shared_ptr<BasisSet> IntegralFactory::basis2() const
{
    return bs2_;
}

std::shared_ptr<BasisSet> IntegralFactory::basis3() const
{
    return bs3_;
}

std::shared_ptr<BasisSet> IntegralFactory::basis4() const
{
    return bs4_;
}

std::shared_ptr<BasisSet> IntegralFactory::basisECP() const
{
	return bs_ecp_;
}

bool IntegralFactory::hasECP() const
{
	return bs_ecp_ != nullptr; 
}

void IntegralFactory::set_basis(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4, std::shared_ptr<BasisSet> bsecp)
{
    bs1_ = bs1;
    bs2_ = bs2;
    bs3_ = bs3;
    bs4_ = bs4;
	bs_ecp_ = bsecp;

    // Use the max am from libint
    init_spherical_harmonics(LIBINT_MAX_AM+1);
}

OneBodyAOInt* IntegralFactory::ao_overlap(int deriv)
{
    return new OverlapInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodySOInt* IntegralFactory::so_overlap(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_overlap(deriv));
    return new OneBodySOInt(ao_int, this);
}

ThreeCenterOverlapInt* IntegralFactory::overlap_3c()
{
    return new ThreeCenterOverlapInt(spherical_transforms_, bs1_, bs2_, bs3_);
}

OneBodyAOInt* IntegralFactory::ao_kinetic(int deriv)
{
    return new KineticInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodySOInt* IntegralFactory::so_kinetic(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_kinetic(deriv));
    return new OneBodySOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_potential(int deriv)
{
    return new PotentialInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodySOInt* IntegralFactory::so_potential(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_potential(deriv));
    return new PotentialSOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_ecp(int deriv)
{
	return new ECPInt(spherical_transforms_, bs1_, bs2_, bs_ecp_, deriv);
}

OneBodySOInt* IntegralFactory::so_ecp(int deriv)
{
	std::shared_ptr<OneBodyAOInt> ao_int(ao_ecp(deriv));
	return new ECPSOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_rel_potential(int deriv)
{
    return new RelPotentialInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodySOInt* IntegralFactory::so_rel_potential(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_rel_potential(deriv));
    return new RelPotentialSOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_pseudospectral(int deriv)
{
    return new PseudospectralInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodySOInt* IntegralFactory::so_pseudospectral(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_pseudospectral(deriv));
    return new OneBodySOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::electrostatic()
{
    return new ElectrostaticInt(spherical_transforms_, bs1_, bs2_, 0);
}

OneBodyAOInt* IntegralFactory::pcm_potentialint()
{
    return new PCMPotentialInt(spherical_transforms_, bs1_, bs2_, 0);
}

OneBodyAOInt* IntegralFactory::ao_dipole(int deriv)
{
    return new DipoleInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodySOInt* IntegralFactory::so_dipole(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_dipole(deriv));
    return new OneBodySOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_nabla(int deriv)
{
    return new NablaInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodySOInt* IntegralFactory::so_nabla(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_nabla(deriv));
    return new OneBodySOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_angular_momentum(int deriv)
{
    return new AngularMomentumInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodySOInt* IntegralFactory::so_angular_momentum(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_angular_momentum(deriv));
    return new OneBodySOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_quadrupole()
{
    return new QuadrupoleInt(spherical_transforms_, bs1_, bs2_);
}

OneBodySOInt* IntegralFactory::so_quadrupole()
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_quadrupole());
    return new OneBodySOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_multipoles(int order)
{
    return new MultipoleInt(spherical_transforms_, bs1_, bs2_, order);
}

OneBodyAOInt* IntegralFactory::ao_efp_multipole_potential(int order)
{
    return new EFPMultipolePotentialInt(spherical_transforms_, bs1_, bs2_, order);
}

OneBodySOInt* IntegralFactory::so_efp_multipole_potential(int order)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_efp_multipole_potential(order));
    return new OneBodySOInt(ao_int, this);
}

OneBodySOInt* IntegralFactory::so_multipoles(int order)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_multipoles(order));
    return new OneBodySOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::ao_traceless_quadrupole()
{
    return new TracelessQuadrupoleInt(spherical_transforms_, bs1_, bs2_);
}

OneBodySOInt* IntegralFactory::so_traceless_quadrupole()
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_traceless_quadrupole());
    return new OneBodySOInt(ao_int, this);
}

OneBodyAOInt* IntegralFactory::electric_field()
{
    return new ElectricFieldInt(spherical_transforms_, bs1_, bs2_);
}

TwoBodyAOInt* IntegralFactory::erd_eri(int deriv, bool use_shell_pairs)
{
#ifdef USING_simint
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "SIMINT")
        return new SimintERI(this, deriv, use_shell_pairs);
#elif defined USING_erd
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "ERD")
        return new ERDERI(this, deriv, use_shell_pairs);
#endif
    return eri(deriv, use_shell_pairs);
}

TwoBodyAOInt* IntegralFactory::eri(int deriv, bool use_shell_pairs)
{
#ifdef USING_simint
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "SIMINT")
        return new SimintERI(this, deriv, use_shell_pairs);
#elif defined USING_erd
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "ERD")
        return new ERDERI(this, deriv, use_shell_pairs);
#endif
    return new ERI(this, deriv, use_shell_pairs);
}

TwoBodyAOInt* IntegralFactory::erf_eri(double omega, int deriv, bool use_shell_pairs)
{
    return new ErfERI(omega, this, deriv, use_shell_pairs);
}

TwoBodyAOInt* IntegralFactory::erf_complement_eri(double omega, int deriv, bool use_shell_pairs)
{
    return new ErfComplementERI(omega, this, deriv, use_shell_pairs);
}

TwoBodyAOInt* IntegralFactory::f12(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return new F12(cf, this, deriv, use_shell_pairs);
}

TwoBodyAOInt* IntegralFactory::f12_scaled(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return new F12Scaled(cf, this, deriv, use_shell_pairs);
}

TwoBodyAOInt* IntegralFactory::f12_squared(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return new F12Squared(cf, this, deriv, use_shell_pairs);
}

TwoBodyAOInt* IntegralFactory::f12g12(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return new F12G12(cf, this, deriv, use_shell_pairs);
}

TwoBodyAOInt* IntegralFactory::f12_double_commutator(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return new F12DoubleCommutator(cf, this, deriv, use_shell_pairs);
}

void IntegralFactory::init_spherical_harmonics(int max_am)
{
    spherical_transforms_.clear();
    ispherical_transforms_.clear();

    for (int i=0; i<=max_am; ++i) {
        spherical_transforms_.push_back(SphericalTransform(i));
        ispherical_transforms_.push_back(ISphericalTransform(i));
    }
}

AOShellCombinationsIterator IntegralFactory::shells_iterator()
{
    return AOShellCombinationsIterator(bs1_, bs2_, bs3_, bs4_);
}

AOShellCombinationsIterator* IntegralFactory::shells_iterator_ptr()
{
    return new AOShellCombinationsIterator(bs1_, bs2_, bs3_, bs4_);
}

AOIntegralsIterator IntegralFactory::integrals_iterator(int p, int q, int r, int s)
{
    return AOIntegralsIterator(bs1_->shell(p), bs2_->shell(q), bs3_->shell(r), bs4_->shell(s));
}

CartesianIter* IntegralFactory::cartesian_iter(int l) const
{
    return new CartesianIter(l);
}

RedundantCartesianIter* IntegralFactory::redundant_cartesian_iter(int l) const
{
    return new RedundantCartesianIter(l);
}

RedundantCartesianSubIter* IntegralFactory::redundant_cartesian_sub_iter(int l) const
{
    return new RedundantCartesianSubIter(l);
}

ShellRotation IntegralFactory::shell_rotation(int am, SymmetryOperation &so, int pure) const
{
    ShellRotation r(am, so, this, pure);
    return r;
}

SphericalTransformIter* IntegralFactory::spherical_transform_iter(int am, int inv, int subl) const
{
    if (subl != -1)
        throw NOT_IMPLEMENTED_EXCEPTION();

    if (inv) {
        return new SphericalTransformIter(ispherical_transforms_[am]);
    }
    return new SphericalTransformIter(spherical_transforms_[am]);
}
