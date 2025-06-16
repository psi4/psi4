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

#ifndef _psi_src_bin_psimrcc_updater_h_
#define _psi_src_bin_psimrcc_updater_h_

#include "psi4/liboptions/liboptions.h"

/**
 *  @file updater.h
 *  @ingroup (PSIMRCC)
 */

namespace psi {
namespace psimrcc {

class Hamiltonian;

/**
 *  @class Updater
 *  @brief Containts the procedure for updating the amplitudes
 */
class Updater {
   public:
    Updater(std::shared_ptr<PSIMRCCWfn> wfn, Options &options);
    virtual ~Updater();
    virtual void update(int cycle, Hamiltonian *heff) = 0;
    void zero_internal_amps();
    void zero_t1_internal_amps();
    void zero_internal_delta_amps();

   protected:
    Options &options_;
    std::shared_ptr<PSIMRCCWfn> wfn_;
};

class MkUpdater : public Updater {
   public:
    MkUpdater(std::shared_ptr<PSIMRCCWfn> wfn, Options &options);
    ~MkUpdater() override;
    void update(int cycle, Hamiltonian *heff) override;
};

class BWUpdater : public Updater {
   public:
    BWUpdater(std::shared_ptr<PSIMRCCWfn> wfn, Options &options);
    ~BWUpdater() override;
    void update(int cycle, Hamiltonian *heff) override;
};

}  // namespace psimrcc
}  // namespace psi

#endif  // _psi_src_bin_psimrcc_updater_h_
