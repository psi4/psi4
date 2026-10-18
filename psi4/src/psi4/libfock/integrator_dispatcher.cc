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

#include "integrator_manager.h"
#include "integrator_dispatcher.h"
#include "v.h"

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

void IntegratorDispatcher::initialize() {
  for (const auto& manager: managers_) manager->initialize();
}

void IntegratorDispatcher::print_header() const {
  outfile->Printf("  ==> Numerical Integrators <==\n\n");
  for (const auto& manager: managers_) manager->print_header();
}

void IntegratorDispatcher::set_D(std::vector<SharedMatrix> Dvec) {
  for (const auto& manager: managers_) manager->set_D(Dvec);
}

void IntegratorDispatcher::set_grac_shift(double shift) {
  // This function isn't final, because GRAC support is currently only in Psi's integrator.
  // Don't make any assumptions about a general IntegratorManager interface
  // for GRAC until after that changes. It's too unstable until then.
  for (const auto& manager: managers_) {
    if (auto grac_manager = dynamic_pointer_cast<VBase>(manager)) {
      grac_manager->set_grac_shift(shift);
    }
  }
}

std::map<std::string, double> IntegratorDispatcher::compute_V(std::vector<SharedMatrix> ret) {
  std::map<std::string, double> rmap;
  for (const auto& manager: managers_) {
    try {
      rmap = manager->compute_V(ret);
      if (!rmap.empty()) {
        quad_values_ = std::move(rmap);
        return quad_values_;
      }
    }
    catch (const std::exception& e) {
      outfile->Printf("Error in manager " + manager->get_name() + ":\n");
      outfile->Printf(e.what());
      timer_off(manager->get_timer_name() + ": Form V"); 
    }
    for (auto& vec: ret) { vec->zero(); }
  }
  throw PSIEXCEPTION("No IntegratorManager managed to compute_V.");
}

SharedMatrix IntegratorDispatcher::compute_gradient() {
  bool computable = false;
  SharedMatrix grad;
  for (const auto& manager: managers_) {
    if (manager->can_compute_gradient()) computable = true;
    else continue;
    try {
      grad = manager->compute_gradient();
    }
    catch (const std::exception& e) {
      outfile->Printf("Error in manager " + manager->get_name() + ":\n");
      outfile->Printf(e.what());
      timer_off(manager->get_timer_name() + ": Form gradient"); 
      continue;
    }
    if (grad != nullptr) return grad;
  }
  if (computable)
    throw PSIEXCEPTION("No IntegratorManager managed to compute_gradient.");
  else throw PSIEXCEPTION("No IntegratorManager has compute_gradient defined.");
}

void IntegratorDispatcher::compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
  bool computable = false;
  for (const auto& manager: managers_) {
    if (manager->can_compute_Vx()) computable = true;
    else continue;
    try {
      manager->compute_Vx(Dx, ret);
    }
    catch (const std::exception& e) {
      outfile->Printf("Error in manager " + manager->get_name() + ":\n");
      outfile->Printf(e.what());
      timer_off(manager->get_timer_name() + ": Form Vx"); 
    }
    return;
  }
  for (auto& vec: ret) { vec->zero(); }
  if (computable)
    throw PSIEXCEPTION("No IntegratorManager managed to compute_Vx.");
  else throw PSIEXCEPTION("No IntegratorManager has compute_Vx defined.");
}

void IntegratorDispatcher::compute_Vx_triplet(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
  bool computable = false;
  for (const auto& manager: managers_) {
    if (manager->can_compute_Vx_triplet()) computable = true;
    else continue;
    try {
      manager->compute_Vx_triplet(Dx, ret);
    }
    catch (const std::exception& e) {
      outfile->Printf("Error in manager " + manager->get_name() + ":\n");
      outfile->Printf(e.what());
      timer_off(manager->get_timer_name() + ": Form Vx"); 
    }
    return;
  }
  for (auto& vec: ret) { vec->zero(); }
  if (computable)
    throw PSIEXCEPTION("No IntegratorManager managed to compute_Vx_triplet.");
  else throw PSIEXCEPTION("No IntegratorManager has compute_Vx_triplet defined.");
}

SharedMatrix IntegratorDispatcher::compute_hessian() {
  bool computable = false;
  SharedMatrix hess;
  for (const auto& manager: managers_) {
    if (manager->can_compute_hessian()) computable = true;
    else continue;
    try {
      hess = manager->compute_hessian();
    }
    catch (const std::exception& e) {
      outfile->Printf("Error in manager " + manager->get_name() + ":\n");
      outfile->Printf(e.what());
      timer_off(manager->get_timer_name() + ": Form hessian"); 
      continue;
    }
    if (hess != nullptr) return hess;
  }
  if (computable)
    throw PSIEXCEPTION("No IntegratorManager managed to compute_hessian.");
  else throw PSIEXCEPTION("No IntegratorManager has compute_hessian.");
}

std::vector<SharedMatrix> IntegratorDispatcher::compute_fock_derivatives() {
  bool computable = false;
  std::vector<SharedMatrix> derivs;
  for (const auto& manager: managers_) {
    if (manager->can_compute_fock_derivatives()) computable = true;
    else continue;
    try {
      derivs = manager->compute_fock_derivatives();
    }
    catch (const std::exception& e) {
      outfile->Printf("Error in manager " + manager->get_name() + ":\n");
      outfile->Printf(e.what());
      timer_off(manager->get_timer_name() + ": Form Fx"); 
      continue;
    }
    if (!derivs.empty()) return derivs;
  }
  if (computable)
    throw PSIEXCEPTION("No IntegratorManager managed to compute_fock_derivatives.");
  else throw PSIEXCEPTION("No IntegratorManager has compute_fock_derivatives.");
}

size_t IntegratorDispatcher::collocation_size() const {
  size_t size = 0;
  // n.b.: this assumes that only one manager collocates. If all collocate, sum.
  for (const auto& manager: managers_) {
    size = std::max(size, manager->collocation_size());
  }
  return size;
}

void IntegratorDispatcher::build_collocation_cache(size_t doubles_per_integrator) {
  for (const auto& manager: managers_) {
    if (auto collocating_integrator = dynamic_pointer_cast<VBase>(manager)) {
      collocating_integrator->build_collocation_cache(doubles_per_integrator);
    } 
  }
}

void IntegratorDispatcher::clear_collocation_cache() {
  for (const auto& manager: managers_) {
    if (auto collocating_integrator = dynamic_pointer_cast<VBase>(manager)) {
      collocating_integrator->clear_collocation_cache();
    } 
  }
}

std::shared_ptr<SuperFunctional> IntegratorDispatcher::functional() {
  std::shared_ptr<SuperFunctional> first_functional;
  for (const auto& manager: managers_) {
    auto current_functional = manager->functional();
    if (!first_functional) first_functional = current_functional;
    else if (first_functional != current_functional) throw PSIEXCEPTION("IntegratorManagers have differnet functionals.");
  }
  return first_functional;
}

std::shared_ptr<VBase> IntegratorDispatcher::psi_manager() const {
  for (const auto& manager: managers_) {
    if (auto vbase = dynamic_pointer_cast<VBase>(manager)) {
      return vbase;
    }
  }
  return nullptr;
}

}

