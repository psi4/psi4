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
#include "psi4/libqt/qt.h"

namespace psi {

void DensityFunctionalIntegratorDispatcher::initialize() {
  std::vector<std::pair<std::string, std::string>> msgs;

  std::erase_if(managers_, [&](const auto& mgr) {
    try {
      mgr->initialize();
      return false;
    } catch (const std::exception& e) {
      const auto& name = mgr->get_name();
      outfile->Printf("Error in manager %s:\n%s", name.c_str(), e.what());
      msgs.emplace_back(name, e.what());
      return true;
    }
  });

  if (managers_.empty())
      throw CompositePsiException("No DensityFunctionalIntegratorManager managed to initialize.", msgs, __FILE__, __LINE__);
}

void DensityFunctionalIntegratorDispatcher::print_header() const {
  outfile->Printf("  ==> Numerical Integrators <==\n\n");
  for (const auto& manager: managers_) manager->print_header();
}

void DensityFunctionalIntegratorDispatcher::set_D(std::vector<SharedMatrix> Dvec) {
  for (const auto& manager: managers_) manager->set_D(Dvec);
}

void DensityFunctionalIntegratorDispatcher::set_grac_shift(double shift) {
  // If the user requests GRAC, that's a change to the low-level numerical integration
  // that most codes don't support. Drop those that can't support it.
  std::vector<std::pair<std::string, std::string>> msgs;
  std::erase_if(managers_, [&](const auto& mgr) {
    if (auto grac_manager = dynamic_pointer_cast<VBase>(mgr)) {
      grac_manager->set_grac_shift(shift);
      return false;
    } else {
      return true;
    }
  });

  if (managers_.empty())
      throw CompositePsiException("No DensityFunctionalIntegratorManager managed to set_grac_shift.", msgs, __FILE__, __LINE__);
}

std::map<std::string, double> DensityFunctionalIntegratorDispatcher::compute_V(std::vector<SharedMatrix> ret) {
  std::map<std::string, double> rmap;
  std::vector<std::pair<std::string, std::string>> msgs;
  for (const auto& manager: managers_) {
    try {
      rmap = manager->compute_V(ret);
      if (!rmap.empty()) {
        quad_values_ = std::move(rmap);
        return quad_values_;
      }
    }
    catch (const std::exception& e) {
      const auto& name = manager->get_name();
      outfile->Printf("Error in manager %s:\n%s", name.c_str(), e.what());
      msgs.emplace_back(name, e.what());
      timer_off(manager->get_timer_name() + ": Form V"); 
    }
    for (auto& vec: ret) { vec->zero(); }
  }
  throw CompositePsiException("No DensityFunctionalIntegratorManager managed to compute_V.", msgs, __FILE__, __LINE__);
}

SharedMatrix DensityFunctionalIntegratorDispatcher::compute_gradient() {
  bool computable = false;
  SharedMatrix grad;
  std::vector<std::pair<std::string, std::string>> msgs;
  for (const auto& manager: managers_) {
    if (manager->can_compute_gradient()) computable = true;
    else continue;
    try {
      grad = manager->compute_gradient();
    }
    catch (const std::exception& e) {
      const auto& name = manager->get_name();
      outfile->Printf("Error in manager %s:\n%s", name.c_str(), e.what());
      msgs.emplace_back(name, e.what());
      timer_off(manager->get_timer_name() + ": Form gradient"); 
      continue;
    }
    if (grad != nullptr) return grad;
  }
  if (computable)
    throw CompositePsiException("No DensityFunctionalIntegratorManager managed to compute_gradient.", msgs, __FILE__, __LINE__);
  else throw PSIEXCEPTION("No DensityFunctionalIntegratorManager has compute_gradient defined.");
}

void DensityFunctionalIntegratorDispatcher::compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
  bool computable = false;
  std::vector<std::pair<std::string, std::string>> msgs;
  for (const auto& manager: managers_) {
    if (manager->can_compute_Vx()) computable = true;
    else continue;
    try {
      manager->compute_Vx(Dx, ret);
      return;
    }
    catch (const std::exception& e) {
      const auto& name = manager->get_name();
      outfile->Printf("Error in manager %s:\n%s", name.c_str(), e.what());
      msgs.emplace_back(name, e.what());
      timer_off(manager->get_timer_name() + ": Form Vx");
      for (auto& vec: ret) { vec->zero(); }
    }
  }
  if (computable)
    throw CompositePsiException("No DensityFunctionalIntegratorManager managed to compute_Vx.", msgs, __FILE__, __LINE__);
  else throw PSIEXCEPTION("No DensityFunctionalIntegratorManager has compute_Vx defined.");
}

void DensityFunctionalIntegratorDispatcher::compute_Vx_triplet(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
  bool computable = false;
  std::vector<std::pair<std::string, std::string>> msgs;
  for (const auto& manager: managers_) {
    if (manager->can_compute_Vx_triplet()) computable = true;
    else continue;
    try {
      manager->compute_Vx_triplet(Dx, ret);
      return;
    }
    catch (const std::exception& e) {
      const auto& name = manager->get_name();
      outfile->Printf("Error in manager %s:\n%s", name.c_str(), e.what());
      msgs.emplace_back(name, e.what());
      timer_off(manager->get_timer_name() + ": Form Vx");
      for (auto& vec: ret) { vec->zero(); }
    }
  }
  if (computable)
    throw CompositePsiException("No DensityFunctionalIntegratorManager managed to compute_Vx_triplet.", msgs, __FILE__, __LINE__);
  else throw PSIEXCEPTION("No DensityFunctionalIntegratorManager has compute_Vx_triplet defined.");
}

SharedMatrix DensityFunctionalIntegratorDispatcher::compute_hessian() {
  bool computable = false;
  SharedMatrix hess;
  std::vector<std::pair<std::string, std::string>> msgs;
  for (const auto& manager: managers_) {
    if (manager->can_compute_hessian()) computable = true;
    else continue;
    try {
      hess = manager->compute_hessian();
    }
    catch (const std::exception& e) {
      const auto& name = manager->get_name();
      outfile->Printf("Error in manager %s:\n%s", name.c_str(), e.what());
      msgs.emplace_back(name, e.what());
      timer_off(manager->get_timer_name() + ": Form hessian"); 
      continue;
    }
    if (hess != nullptr) return hess;
  }
  if (computable)
    throw CompositePsiException("No DensityFunctionalIntegratorManager managed to compute_hessian.", msgs, __FILE__, __LINE__);
  else throw PSIEXCEPTION("No DensityFunctionalIntegratorManager has compute_hessian.");
}

std::vector<SharedMatrix> DensityFunctionalIntegratorDispatcher::compute_fock_derivatives() {
  bool computable = false;
  std::vector<SharedMatrix> derivs;
  std::vector<std::pair<std::string, std::string>> msgs;
  for (const auto& manager: managers_) {
    if (manager->can_compute_fock_derivatives()) computable = true;
    else continue;
    try {
      derivs = manager->compute_fock_derivatives();
    }
    catch (const std::exception& e) {
      const auto& name = manager->get_name();
      outfile->Printf("Error in manager %s:\n%s", name.c_str(), e.what());
      msgs.emplace_back(name, e.what());
      timer_off(manager->get_timer_name() + ": Form Fx"); 
      continue;
    }
    if (!derivs.empty()) return derivs;
  }
  if (computable)
    throw CompositePsiException("No DensityFunctionalIntegratorManager managed to compute_fock_derivatives.", msgs, __FILE__, __LINE__);
  else throw PSIEXCEPTION("No DensityFunctionalIntegratorManager has compute_fock_derivatives.");
}

size_t DensityFunctionalIntegratorDispatcher::collocation_size() const {
  size_t size = 0;
  // n.b.: this assumes that only one manager collocates. If all collocate, sum.
  for (const auto& manager: managers_) {
    size = std::max(size, manager->collocation_size());
  }
  return size;
}

void DensityFunctionalIntegratorDispatcher::build_collocation_cache(size_t doubles_per_integrator) {
  for (const auto& manager: managers_) {
    if (auto collocating_integrator = dynamic_pointer_cast<VBase>(manager)) {
      collocating_integrator->build_collocation_cache(doubles_per_integrator);
    } 
  }
}

void DensityFunctionalIntegratorDispatcher::clear_collocation_cache() {
  for (const auto& manager: managers_) {
    if (auto collocating_integrator = dynamic_pointer_cast<VBase>(manager)) {
      collocating_integrator->clear_collocation_cache();
    } 
  }
}

std::shared_ptr<SuperFunctional> DensityFunctionalIntegratorDispatcher::functional() {
  std::shared_ptr<SuperFunctional> first_functional;
  for (const auto& manager: managers_) {
    auto current_functional = manager->functional();
    if (!first_functional) first_functional = current_functional;
    else if (first_functional != current_functional) throw PSIEXCEPTION("DensityFunctionalIntegratorManagers have different functionals.");
  }
  return first_functional;
}

std::shared_ptr<VBase> DensityFunctionalIntegratorDispatcher::psi_manager() const {
  for (const auto& manager: managers_) {
    if (auto vbase = dynamic_pointer_cast<VBase>(manager)) {
      return vbase;
    }
  }
  return nullptr;
}

}

