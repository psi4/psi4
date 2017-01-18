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

#ifndef _psi_src_lib_libmints_eri_h
#define _psi_src_lib_libmints_eri_h

#include <libint/libint.h>
#include <libint/libderiv.h>
#include "psi4/libmints/twobody.h"

namespace psi {

class BasisSet;
class GaussianShell;
class TwoBodyAOInt;
class IntegralFactory;
class SphericalTransform;
class Fjt;
class AOShellCombinationsIterator;
class CorrelationFactor;

/**
  * \ingroup MINTS
  * Structure to hold precomputed shell pair information
  */
typedef struct ShellPair_typ {
    //! Shells for this information.
    int i, j;
    //! Matrix over primitives with x, y, z coordinate of average Gaussian
    double ***  P;
    //! Distance between shell i and shell j centers
    double AB[3];
    //! Distance between P and shell i center
    double ***  PA;
    //! Distance between P and shell j center
    double ***  PB;
    //! Array of alphas for both centers
    double *  ai, *  aj;
    //! Array of the gammas (ai + aj)
    double **  gamma;
    //! Contraction coefficients
    double *  ci, *  cj;
    //! Overlap between primitives on i and j
    double **  overlap;
} ShellPair;

/*! \ingroup MINTS
 *  \class ERI
 *  \brief Capable of computing two-electron repulsion integrals.
 */
class TwoElectronInt : public TwoBodyAOInt
{
protected:
    //! Libint object.
    Libint_t libint_;
    //! Libderiv object
    Libderiv_t libderiv_;

    //! Maximum cartesian class size.
    int max_cart_;

    //! Computes the fundamental
    Fjt *fjt_;

    //! Computes the ERIs between four shells.
    size_t compute_quartet(int, int, int, int);

    //! Computes the ERI derivatives between four shells.
    size_t compute_quartet_deriv1(int, int, int, int);

    //! Computes the ERI second derivative between four shells.
    size_t compute_quartet_deriv2(int, int, int, int);

    //! Form shell pair information. Must be smart enough to handle arbitrary basis sets
    void init_shell_pairs12();
    void init_shell_pairs34();

    //! Free shell pair information
    void free_shell_pairs12();
    void free_shell_pairs34();

    //! Should we use shell pair information?
    bool use_shell_pairs_;

    //! Stack memory pointer, used in init_shell_pairs, freed in destructor
    double *stack12_, *stack34_;

    //! Shell pair information
    ShellPair **pairs12_, **pairs34_;

    //! Evaluates how much memory (in doubles) is needed to store shell pair data
    size_t memory_to_store_shell_pairs(const std::shared_ptr<BasisSet>&, const std::shared_ptr<BasisSet>&);

    //! Original shell index requested
    int osh1_, osh2_, osh3_, osh4_;

    //! Were the indices permuted?
    bool p13p24_, p12_, p34_;


public:
    //! Constructor. Use an IntegralFactory to create this object.
    TwoElectronInt(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);

    virtual ~TwoElectronInt();

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator&);

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual size_t compute_shell(int, int, int, int);

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    virtual size_t compute_shell_deriv1(int, int, int, int);

    /// Compute ERI second derivatives between 4 sheels. Result is stored in buffer.
    virtual size_t compute_shell_deriv2(int, int, int, int);
};

class ERI : public TwoElectronInt
{
public:
    ERI(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~ERI();
};

class F12 : public TwoElectronInt
{
public:
    F12(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~F12();
};

class F12Scaled : public TwoElectronInt
{
public:
    F12Scaled(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~F12Scaled();
};

class F12Squared : public TwoElectronInt
{
public:
    F12Squared(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~F12Squared();
};

class F12G12 : public TwoElectronInt
{
public:
    F12G12(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~F12G12();
};

class F12DoubleCommutator : public TwoElectronInt
{
public:
    F12DoubleCommutator(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~F12DoubleCommutator();
};

class ErfERI : public TwoElectronInt
{
public:
    ErfERI(double omega, const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~ErfERI();

    void setOmega(double omega);
};

class ErfComplementERI : public TwoElectronInt
{
public:
    ErfComplementERI(double omega, const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~ErfComplementERI();

    void setOmega(double omega);
};

}

#endif
