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

#ifndef LIBFOCK_DFT_H
#define LIBFOCK_DFT_H
#include "psi4/libmints/typedefs.h"
#include <vector>
#include <map>

namespace psi {
class BasisSet;
class Options;
class DFTGrid;
class PointFunctions;
class SuperFunctional;

// => BASE CLASS <= //

/**
 * Class VBase
 *
 * Class to compute KS-V matrices and
 * K-matrix-vector products
 **/

class VBase {

protected:
    /// Debug flag
    int debug_;
    /// Print flag
    int print_;
    /// Options object, used to build grid
    Options& options_;
    /// Basis set used in the integration
    std::shared_ptr<BasisSet> primary_;
    /// Desired superfunctional kernal
    std::shared_ptr<SuperFunctional> functional_;
    /// Point function computer (densities, gammas, basis values)
    std::shared_ptr<PointFunctions> properties_;
    /// Integration grid, built by KSPotential
    std::shared_ptr<DFTGrid> grid_;
    /// Quadrature values obtained during integration
    std::map<std::string, double> quad_values_;

    /// AO2USO matrix (if not C1)
    SharedMatrix AO2USO_;

    /// Vector of V matrices (built by form_D)
    std::vector<SharedMatrix> V_;
    /// Vector of C1 V matrices (built by USO2AO)
    std::vector<SharedMatrix> V_AO_;

    /// Vector of occupied C matrices (used for D and KE density)
    std::vector<SharedMatrix> C_;
    /// Vector of D matrices (built by form_D)
    std::vector<SharedMatrix> D_;
    /// Vector of C1 C matrices (built by USO2AO)
    std::vector<SharedMatrix> C_AO_;
    /// Vector of C1 D matrices (built by USO2AO)
    std::vector<SharedMatrix> D_AO_;

    /// Vector of Caocc matrices (TDDFT)
    std::vector<SharedMatrix> Caocc_;
    /// Vector of Cavir matrices (TDDFT)
    std::vector<SharedMatrix> Cavir_;
    /// Vector of Perturbation matrices (TDDFT, ia)
    std::vector<SharedMatrix> P_;
    /// Vector of Perturbation matrices (TDDFT, SO)
    std::vector<SharedMatrix> P_SO_;
    /// Vector of Perturbation matrices (TDDFT, AO)
    std::vector<SharedMatrix> P_AO_;

    virtual void compute_D();
    virtual void USO2AO();
    virtual void AO2USO();

    /// Actually build V_AO
    virtual void compute_V() = 0;
    /// Set things up
    void common_init();
public:
     VBase(std::shared_ptr<SuperFunctional> functional,
           std::shared_ptr<BasisSet> primary, Options& options);
     virtual ~VBase();

    static std::shared_ptr<VBase> build_V(std::shared_ptr<BasisSet> primary,
                                          std::shared_ptr<SuperFunctional> functional,
                                          Options& options,
                                          const std::string& type = "RV");

    std::shared_ptr<BasisSet> basis() const { return primary_; }
    std::shared_ptr<SuperFunctional> functional() const { return functional_; }
    std::shared_ptr<PointFunctions> properties() const { return properties_; }
    std::shared_ptr<DFTGrid> grid() const { return grid_; }
    std::map<std::string, double>& quadrature_values() { return quad_values_; }

    /// Grab this, clear, and push Cocc matrices (with symmetry) to change GS density
    std::vector<SharedMatrix>& C() { return C_; }
    std::vector<SharedMatrix>& Caocc() { return Caocc_; }
    std::vector<SharedMatrix>& Cavir() { return Cavir_; }
    std::vector<SharedMatrix>& P() { return P_; }
    const std::vector<SharedMatrix>& V() const { return V_; }
    const std::vector<SharedMatrix>& D() const { return D_; }

    /// Throws by default
    virtual SharedMatrix compute_gradient();

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

    virtual void initialize();
    virtual void compute();
    virtual void finalize();

    virtual void print_header() const;
};

// => APPLIED CLASSES <= //

class RV : public VBase {

protected:

    // Actually build V_AO
    virtual void compute_V();

public:
    RV(std::shared_ptr<SuperFunctional> functional,
        std::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~RV();

    virtual void initialize();
    virtual void finalize();

    virtual SharedMatrix compute_gradient();

    virtual void print_header() const;
};

class UV : public VBase {

protected:

    // Actually build V_AO
    virtual void compute_V();

public:
    UV(std::shared_ptr<SuperFunctional> functional,
        std::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~UV();

    virtual void initialize();
    virtual void finalize();

    virtual SharedMatrix compute_gradient();

    virtual void print_header() const;
};

class RK : public RV {

protected:

    // Actually build V_AO
    virtual void compute_V();

public:
    RK(std::shared_ptr<SuperFunctional> functional,
        std::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~RK();

    virtual void print_header() const;
};

class UK : public UV {

protected:

    // Actually build V_AO
    virtual void compute_V();

public:
    UK(std::shared_ptr<SuperFunctional> functional,
        std::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~UK();

    virtual void print_header() const;
};

}
#endif
