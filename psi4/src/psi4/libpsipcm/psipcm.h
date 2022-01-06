/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "psi4/libmints/dimension.h"
#include "psi4/libmints/typedefs.h"

#include <PCMSolver/pcmsolver.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace psi {
class BasisSet;
class Options;
class PCMPotentialInt;

class PCM final {
   public:
    enum class CalcType : int { Total, NucAndEle, EleOnly };
    PCM() = default;
    PCM(const std::string &pcmsolver_parsed_fname, int print_level, std::shared_ptr<BasisSet> basisset);
    PCM(const PCM *);
    ~PCM() {}
    /*! \brief Compute polarization energy and Fock matrix contribution
     *  \param[in] D density matrix
     *  \param[in] type how to treat MEP and ASC
     */
    std::pair<double, SharedMatrix> compute_PCM_terms(const SharedMatrix &D, CalcType type = CalcType::Total) const;
    SharedMatrix compute_V(const SharedMatrix &D);

   private:
    /// The number of tesserae in PCMSolver.
    int ntess_;
    /// The number of irreducible tesserae in PCMSolver.
    int ntessirr_;
    Dimension tesspi_;
    /// Charges and {x,y,z} coordinates of the cavity points
    SharedMatrix tess_Zxyz_;
    /// Nucler MEP at cavity points
    SharedVector MEP_n_;
    /// Computes electronic MEP at cavity points
    SharedVector compute_electronic_MEP(const SharedMatrix &D) const;
    /// Calculate energy using total charges and potentials
    double compute_E_total(const SharedVector &MEP_e) const;
    /// Calculate energy separating between charges and potentials
    double compute_E_separate(const SharedVector &MEP_e) const;
    /// Calculate electronic polarization energy (U_ee) only (for CC step)
    double compute_E_electronic(const SharedVector &MEP_e) const;
    /*! \brief Compute PCM potential
     *  \param[in] ASC the apparent surface charge to contract with
     *  charge-attraction integrals
     */
    SharedMatrix compute_Vpcm(const SharedVector &ASC) const;

    /// Current basis set (for puream and nao/nso info)
    std::shared_ptr<BasisSet> basisset_;

    /// The AO->SO transformation matrix, which is used for transforming
    /// matrices between pure and Cartesian representations.
    SharedMatrix my_aotoso_;

    /// Factory for the electrostatic integrals
    PCMPotentialInt *potential_int_;

    /// Handle to stuff provided by PCMSolver
    std::shared_ptr<pcmsolver_context_t> context_;

    /// Filename for the PCM input file as parsed by PCMSolver
    std::string pcmsolver_parsed_fname_;

    /// print level
    int pcm_print_;
};

namespace detail {
std::pair<std::vector<double>, std::vector<double>> collect_atoms(std::shared_ptr<Molecule>);

PCMInput pcmsolver_input();

void host_writer(const char *);

std::shared_ptr<pcmsolver_context_t> init_PCMSolver(const std::string &pcmsolver_parsed_fname,
                                                    const std::shared_ptr<Molecule> &molecule);
}  // namespace detail
}  // namespace psi
#endif
#endif
