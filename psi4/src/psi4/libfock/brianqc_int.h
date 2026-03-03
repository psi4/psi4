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

#ifndef LIBFOCK_BRIANQC_INT_H
#define LIBFOCK_BRIANQC_INT_H
#ifdef USING_BrianQC
#include "integrator_manager.h"

namespace psi {

// auxiliary structures passing the grid to BrianQC
struct BrianRadialPoint {
    double r;
    double w;
};

struct BrianAngularPoint {
    double x, y, z;
    double w;
};

struct BrianBlock {
    std::vector<BrianRadialPoint> radialPoints;
    std::vector<BrianAngularPoint> angularPoints;
};

class BrianQCBase : public IntegratorManager {
   using IntegratorManager::IntegratorManager;
    
   void initialize() override;
   void build_grid();

   protected:
    /// Integrator object for GauXC based integration
    std::shared_ptr<GauXC::XCIntegrator<Eigen::MatrixXd>> integrator_;
    /// SG0 and SG1... It's an array of 4 for historical reasons.
    /// Array specifies which named grid. Outer vector specifies atom.
    static std::array<std::vector<std::vector<BrianBlock>>, 4> brian_standard_grids;

};

class BrianRV : public BrianQCBase {
    using BrianQCBase::BrianQCBase;

    std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) override;
};

class BrianUV : public BrianQCBase {
    using BrianQCBase::BrianQCBase;

    std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) override;
};
}

#endif
#endif

