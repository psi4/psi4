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

#ifndef FISAPT_H
#define FISAPT_H

#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/wavefunction.h"
#include <map>
#include <tuple>

namespace psi {

class JK;

namespace fisapt {

class FISAPT {

protected:

    // sSAPT0 exchange-scaling
    double sSAPT0_scale_;
    /// Global options object
    Options& options_;
    /// Memory in doubles
    size_t doubles_;
    /// Reference wavefunction
    std::shared_ptr<Wavefunction> reference_;

    /// Orbital Basis Set (full molecule)
    std::shared_ptr<BasisSet> primary_;
    std::shared_ptr<BasisSet> df_basis_scf_;

    /// Global JK object
    std::shared_ptr<JK> jk_;

    /// Map of scalars
    std::map<std::string, double> scalars_;
    /// Map of vectors
    std::map<std::string, std::shared_ptr<Vector> > vectors_;
    /// Map of matrices
    std::map<std::string, std::shared_ptr<Matrix> > matrices_;

    // => FISAPT 0th-Order Wavefunction <= //

    /// Common initialization (bases, orbitals, eigenvalues, etc)
    void common_init();
    /// Print header, bases, sizes, etc
    void print_header();
    /// Localize the active occupied orbitals via IBO2
    void localize();
    /// Partition the nuclei and electrons
    void partition();
    /// Build the overlap integrals S
    void overlap();
    /// Build the kinetic integrals T
    void kinetic();
    /// Build the nuclear potentials V and interaction energies
    void nuclear();
    /// Build the J/K potentials for C, D, and E
    void coulomb();
    /// Solve the relaxed SCF equations for A0 and B0
    void scf();
    /// Freeze the core orbitals
    void freeze_core();
    /// Produce unified matrices for A', B', and C'
    void unify();
    /// Plot some analysis files
    void plot();

    // => F-SAPT0 <= //

    /// Localize
    void flocalize();
    /// Electrostatics
    void felst();
    /// Exchange
    void fexch();
    /// Induction
    void find();
    /// Dispersion
    void fdisp();
    /// Output
    void fdrop();

    // => SAPT0 <= //

    /// Delta HF
    void dHF();
    /// Electrostatics
    void elst();
    /// Exchange
    void exch();
    /// Induction
    void ind();
    /// Dispersion
    void disp();
    /// Print SAPT results
    void print_trailer();

    // Build the ExchInd20 potential in the monomer A ov space
    std::shared_ptr<Matrix> build_exch_ind_pot(std::map<std::string, std::shared_ptr<Matrix> >& vars);
    // Build the Ind20 potential in the monomer A ov space
    std::shared_ptr<Matrix> build_ind_pot(std::map<std::string, std::shared_ptr<Matrix> >& vars);

    /// Helper to drop a matrix to filepath/A->name().dat
    void drop(std::shared_ptr<Matrix> A, const std::string& filepath);
    /// Helper to drop a vector to filepath/A->name().dat
    void drop(std::shared_ptr<Vector> A, const std::string& filepath);
    /// Helper to extract columns from a matrix
    static std::shared_ptr<Matrix> extract_columns(
        const std::vector<int>& cols,
        std::shared_ptr<Matrix> A);

public:
    /// Initialize an FISAPT object with an SCF reference
    FISAPT(std::shared_ptr<Wavefunction> scf, Options& options);
    virtual ~FISAPT();

    /// Gogo!
    void compute_energy();

};

class FISAPTSCF {

protected:


    /// Global options object
    Options& options_;

    /// Global JK object
    std::shared_ptr<JK> jk_;

    /// Map of scalars
    std::map<std::string, double> scalars_;
    /// Map of vectors
    std::map<std::string, std::shared_ptr<Vector> > vectors_;
    /// Map of matrices
    std::map<std::string, std::shared_ptr<Matrix> > matrices_;

    /// Print orbitals
    void print_orbitals(
        const std::string& header,
        int start,
        std::shared_ptr<Vector> eps
        );


public:

    FISAPTSCF(
        std::shared_ptr<JK> jk,    // JK object
        double enuc,                 // Nuclear repulsion energy
        std::shared_ptr<Matrix> S, // Overlap integrals
        std::shared_ptr<Matrix> X, // Restricted orthogonalization matrix [nbf x nmo]
        std::shared_ptr<Matrix> T, // Kinetic integrals
        std::shared_ptr<Matrix> V, // Potential integrals
        std::shared_ptr<Matrix> W, // External embedding potential
        std::shared_ptr<Matrix> C, // Guess for occupied orbitals [nbf x nocc]
        Options& options
        );
    virtual ~FISAPTSCF();

    void compute_energy();

    std::map<std::string, double>& scalars()                      { return scalars_; }
    std::map<std::string, std::shared_ptr<Vector> >& vectors()  { return vectors_; }
    std::map<std::string, std::shared_ptr<Matrix> >& matrices() { return matrices_; }

};

class CPHF_FISAPT {

friend class FISAPT;

protected:

    // => Global Data <= //


    // Convergence tolerance
    double delta_;
    // Maximum allowed iterations
    int maxiter_;
    // JK Object
    std::shared_ptr<JK> jk_;

    // => Monomer A Problem <= //

    // Perturbation applied to A
    std::shared_ptr<Matrix> w_A_;
    // Response of A
    std::shared_ptr<Matrix> x_A_;
    // Active occ orbital coefficients of A
    std::shared_ptr<Matrix> Cocc_A_;
    // Active vir orbital coefficients of A
    std::shared_ptr<Matrix> Cvir_A_;
    // Active occ orbital eigenvalues of A
    std::shared_ptr<Vector> eps_occ_A_;
    // Active vir orbital eigenvalues of A
    std::shared_ptr<Vector> eps_vir_A_;

    // => Monomer B Problem <= //

    // Perturbation applied to B
    std::shared_ptr<Matrix> w_B_;
    // Response of B
    std::shared_ptr<Matrix> x_B_;
    // Active occ orbital coefficients of B
    std::shared_ptr<Matrix> Cocc_B_;
    // Active vir orbital coefficients of B
    std::shared_ptr<Matrix> Cvir_B_;
    // Active occ orbital eigenvalues of B
    std::shared_ptr<Vector> eps_occ_B_;
    // Active vir orbital eigenvalues of B
    std::shared_ptr<Vector> eps_vir_B_;

    // Form the s = Ab product for the provided vectors b (may or may not need more iterations)
    std::map<std::string, std::shared_ptr<Matrix> > product(std::map<std::string, std::shared_ptr<Matrix> > b);
    // Apply the denominator from r into z
    void preconditioner(std::shared_ptr<Matrix> r,
                        std::shared_ptr<Matrix> z,
                        std::shared_ptr<Vector> o,
                        std::shared_ptr<Vector> v);

public:
    CPHF_FISAPT();
    virtual ~CPHF_FISAPT();

    void compute_cphf();
};

} // Namespace fisapt

} // Namespace psi

#endif
