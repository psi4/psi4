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
class BlockOPoints;

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
    /// Number of threads
    int num_threads_;
    /// Number of basis functions;
    int nbf_;
    /// Rho threshold for the second derivative;
    double v2_rho_cutoff_;
    /// VV10 interior kernel threshold
    double vv10_rho_cutoff_;
    /// Options object, used to build grid
    Options& options_;
    /// Basis set used in the integration
    std::shared_ptr<BasisSet> primary_;
    /// Desired superfunctional kernal
    std::shared_ptr<SuperFunctional> functional_;
    /// Desired superfunctional kernal
    std::vector<std::shared_ptr<SuperFunctional>> functional_workers_;
    /// Point function computer (densities, gammas, basis values)
    std::vector<std::shared_ptr<PointFunctions>> point_workers_;
    /// Integration grid, built by KSPotential
    std::shared_ptr<DFTGrid> grid_;
    /// Quadrature values obtained during integration
    std::map<std::string, double> quad_values_;

    /// AO2USO matrix (if not C1)
    SharedMatrix AO2USO_;
    SharedMatrix USO2AO_;

    /// Vector of C1 D matrices (built by USO2AO)
    std::vector<SharedMatrix> D_AO_;

    // GRAC data
    bool grac_initialized_;

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
    std::vector<std::shared_ptr<PointFunctions>> properties() const { return point_workers_; }
    std::shared_ptr<DFTGrid> grid() const { return grid_; }
    std::shared_ptr<BlockOPoints> get_block(int block);
    size_t nblocks();
    std::map<std::string, double>& quadrature_values() { return quad_values_; }

    // Set the D matrix, get it back if needed
    void set_D(std::vector<SharedMatrix> Dvec);
    const std::vector<SharedMatrix>& Dao() const { return D_AO_; }

    // Set the site of the grac shift
    void set_grac_shift(double value);

    /// Throws by default
    virtual void compute_V(std::vector<SharedMatrix> ret);
    virtual void compute_Vx(std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret);
    virtual SharedMatrix compute_gradient();
    virtual SharedMatrix compute_hessian();

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

    virtual void initialize();
    virtual void finalize();

    virtual void print_header() const;
};

// => APPLIED CLASSES <= //

class RV : public VBase {

protected:

public:
    RV(std::shared_ptr<SuperFunctional> functional,
        std::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~RV();

    virtual void initialize();
    virtual void finalize();

    virtual void compute_V(std::vector<SharedMatrix> ret);
    virtual void compute_Vx(std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret);
    virtual SharedMatrix compute_gradient();
    virtual SharedMatrix compute_hessian();

    virtual void print_header() const;
};

class UV : public VBase {

protected:

public:
    UV(std::shared_ptr<SuperFunctional> functional,
        std::shared_ptr<BasisSet> primary,
        Options& options);
    virtual ~UV();

    virtual void initialize();
    virtual void finalize();

    virtual void compute_V(std::vector<SharedMatrix> ret);
    virtual void compute_Vx(std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret);
    virtual SharedMatrix compute_gradient();

    virtual void print_header() const;
};


}
#endif
