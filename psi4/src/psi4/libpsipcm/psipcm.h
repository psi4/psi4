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

#ifndef PCM_H
#define PCM_H
#ifdef USING_PCMSolver

#include <memory>
#include <tuple>
#include <vector>

#include <PCMSolver/pcmsolver.h>

namespace psi {
class BasisSet;
class Matrix;
class Molecule;
class Options;
class PCMPotentialInt;
class Vector;

using SharedMatrix = std::shared_ptr<Matrix>;
using SharedVector = std::shared_ptr<Vector>;

class PCM {
   public:
    enum class CalcType : int { Total, NucAndEle, EleOnly };
    PCM() = default;
    PCM(int print_level, std::shared_ptr<BasisSet> basisset);
    ~PCM();
    double compute_E(SharedMatrix &D, CalcType type = CalcType::Total);
    SharedMatrix compute_V();
    SharedMatrix compute_V_electronic();  // This is needed by the CC code (and maybe the LR-SCF code)

   protected:
    /// The number of tesserae in PCMSolver.
    int ntess_;
    /// The number of irreducible tesserae in PCMSolver.
    int ntessirr_;
    /// A matrix to hold the charges and {x,y,z} coordinates of the tesserae
    SharedMatrix tess_Zxyz_;
    /// A scratch array to hold the electronic potential values at the tesserae
    double *tess_pot_e_;
    /// A scratch array to hold the nuclear potential values at the tesserae (unchanging)
    double *tess_pot_n_;
    /// A scratch array to hold the total potential values at the tesserae
    double *tess_pot_;
    /// A scratch array to hold the electronic charges at the tesserae
    double *tess_charges_e_;
    /// A scratch array to hold the nuclear charges at the tesserae
    double *tess_charges_n_;
    /// A scratch array to hold the charges at the tesserae
    double *tess_charges_;
    /// Calculate energy using total charges and potentials
    double compute_E_total(SharedMatrix &D);
    /// Calculate energy separating between charges and potentials
    double compute_E_separate(SharedMatrix &D);
    /// Calculate electronic polarization energy (U_ee) only (for CC step)
    double compute_E_electronic(SharedMatrix &D);

    /// Current basis set (for puream and nao/nso info)
    std::shared_ptr<BasisSet> basisset_;

    /// The AO->SO transformation matrix, which is used for transforming
    /// matrices between pure and Cartesian representations.
    SharedMatrix my_aotoso_;

    /// Factory for the electrostatic integrals
    PCMPotentialInt *potential_int_;

    /// Handle to stuff provided by PCMSolver
    std::shared_ptr<pcmsolver_context_t> context_;

    /// print level
    int pcm_print_;
};

namespace detail {
std::tuple<std::vector<double>, std::vector<double>> collect_atoms(std::shared_ptr<Molecule>);

PCMInput pcmsolver_input();

void host_writer(const char *);

}  // detail
}  // psi
#endif
#endif
