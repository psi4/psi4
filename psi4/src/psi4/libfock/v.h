/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include "psi4/pragma.h"
#include "integrator_manager.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <string>

namespace psi {
class BasisSet;
class Options;
class DFTGrid;
class PointFunctions;
class SuperFunctional;
class BlockOPoints;
#ifdef USING_BrianQC
class BrianQCBase;
#endif
#ifdef USING_gauxc
class GauXCBase;
#endif

// => BASE CLASS <= //

/**
 * Class VBase
 *
 * Class to compute KS-V matrices and
 * K-matrix-vector products
 **/

class PSI_API VBase: public IntegratorManager {
   protected:
    /// Debug flag
    int debug_;
    /// Print flag
    int print_;
    /// Number of threads
    int num_threads_;
    /// Rho threshold for the second derivative;
    double v2_rho_cutoff_;
    /// VV10 interior kernel threshold
    double vv10_rho_cutoff_;
    /// Desired superfunctional kernel
    std::shared_ptr<SuperFunctional> functional_;
    /// Desired superfunctional kernel
    std::vector<std::shared_ptr<SuperFunctional>> functional_workers_;
    /// Point function computer (densities, gammas, basis values)
    std::vector<std::shared_ptr<PointFunctions>> point_workers_;
    /// Integration grid, built by KSPotential
    std::shared_ptr<DFTGrid> grid_;
#ifdef USING_BrianQC
    /// Integrator object for BrianQC based integration
    std::shared_ptr<BrianQCBase> brianqc_integrator_;
#endif
#ifdef USING_gauxc
    /// Integrator object for GauXC based integration
    std::shared_ptr<GauXCBase> gauxc_integrator_;
#endif
    // Caches collocation grids
    std::unordered_map<size_t, std::map<std::string, SharedMatrix>> cache_map_;
    int cache_map_deriv_;

    // GRAC data
    bool grac_initialized_;

    // VV10 dispersion, return vv10_nlc energy
    void prepare_vv10_cache(DFTGrid& nlgrid, SharedMatrix D,
                            std::vector<std::map<std::string, SharedVector>>& vv10_cache,
                            std::vector<std::shared_ptr<PointFunctions>>& nl_point_workers, int ansatz = 1);
    double vv10_nlc(SharedMatrix D, SharedMatrix ret);
    SharedMatrix vv10_nlc_gradient(SharedMatrix D);

    /// Set things up
    void common_init();

   public:
    VBase(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    virtual ~VBase();

    static std::shared_ptr<VBase> build_V(std::shared_ptr<BasisSet> primary,
                                          std::shared_ptr<SuperFunctional> functional, Options& options,
                                          const std::string& type = "RV");

    std::shared_ptr<SuperFunctional> functional() const { return functional_; }
    std::vector<std::shared_ptr<PointFunctions>> properties() const { return point_workers_; }
    std::shared_ptr<DFTGrid> grid() const { return grid_; }
    std::shared_ptr<BlockOPoints> get_block(int block);
    size_t nblocks();

    // Creates a collocation cache map based on stride
    size_t collocation_size() const override;
    void build_collocation_cache(size_t memory);
    void clear_collocation_cache() { cache_map_.clear(); }

    // Set the site of the grac shift
    void set_grac_shift(double value);

    virtual SharedMatrix compute_gradient();
    virtual SharedMatrix compute_hessian();

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

    virtual void initialize();
    virtual void finalize();

    void print_header() const override;
};

// => Derived Classes <= //
class SAP : public VBase {
   protected:
   public:
    SAP(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    ~SAP() override;

    void initialize() override;
    void finalize() override;

    std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) override;
    void print_header() const override;
};

class RV : public VBase {
   protected:
   public:
    RV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    ~RV() override;

    void initialize() override;
    void finalize() override;

    // compute_V assuming same orbitals for different spin. Computes V_alpha, not spin-summed V.
    std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) override;
    bool can_compute_Vx() override { return true; };
    /// Compute the orbital derivative of the KS potential, contract against Dx, and
    /// putting the result in ret. ret[i] is Vx where x = Dx[i]. The "true" vector has
    /// 2^-0.5 Dx[i] for each input spin case and returns **half** the α component of the output.
    /// The singlet flag controls whether to assume singlet spin-integration (β components
    /// are the α components) or triplet (β components are -α components)
    void compute_Vx_full(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret, bool singlet);
    /// A convenience function to call compute_Vx_full for singlets.
    /// And no, we can't just make singlet a default argument. Then compute_Vx has different signatures for
    /// different VBase subclasses, so we can't call compute_Vx from VBase, which breaks the hessian code.
    void compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) override { compute_Vx_full(Dx, ret, true); };
    bool can_compute_gradient() override { return true; };
    SharedMatrix compute_gradient() override;
    bool can_compute_hessian() override { return true; };
    std::vector<SharedMatrix> compute_fock_derivatives() override;
    SharedMatrix compute_hessian() override;
};

class UV : public VBase {
   protected:
   public:
    UV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    ~UV() override;

    void initialize() override;
    void finalize() override;

    std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) override;
    bool can_compute_Vx() override { return true; };
    /// Compute the orbital derivative of the KS potential, contract against Dx, and
    /// putting the result in ret. ret[i] is Vx where x = Dx[i].
    /// ret[2n], ret[2n+1] are alpha and beta Vx where x concatenates Dx[2n] (α) and Dx[2n+1] (β).
    void compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) override;
    bool can_compute_gradient() override { return true; };
    SharedMatrix compute_gradient() override;
    bool can_compute_hessian() override { return true; };
    std::vector<SharedMatrix> compute_fock_derivatives() override;
    SharedMatrix compute_hessian() override;
};
}
#endif
