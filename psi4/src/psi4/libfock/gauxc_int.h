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

#ifndef LIBFOCK_GAUXC_INT_H
#define LIBFOCK_GAUXC_INT_H
#ifdef USING_gauxc
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#include <eigen3/Eigen/Core>

#include "integrator_manager.h"

namespace psi {

class GauXCBase : public IntegratorManager {
   using IntegratorManager::IntegratorManager;

   public:
    void initialize() override;
    virtual ExchCXX::Spin spin() const = 0;

   protected:
    /// Integrator object for GauXC based integration
    std::shared_ptr<GauXC::XCIntegrator<Eigen::MatrixXd>> integrator_;
    

};

class GauRV : public GauXCBase {
    using GauXCBase::GauXCBase;

    std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) override;
    ExchCXX::Spin spin() const override { return ExchCXX::Spin::Unpolarized; };
};

class GauUV : public GauXCBase {
    using GauXCBase::GauXCBase;

    std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) override;
    ExchCXX::Spin spin() const override { return ExchCXX::Spin::Polarized; };
};
}

#endif
#endif

