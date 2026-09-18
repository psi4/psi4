/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef _psi_src_lib_libpsi4util_memory_ledger_h_
#define _psi_src_lib_libpsi4util_memory_ledger_h_

#include <atomic>
#include <cstddef>
#include <cstdint>

namespace psi {

/*! \brief Process-wide tally of the budget-sized buffers that are resident right now.
 *
 * scf_initialize() divides the SCF memory budget between the JK object and the DFT
 * collocation cache by reading the global memory setting.  On its own that setting
 * describes an empty process, so a driver that starts an SCF while earlier
 * wavefunctions are still alive hands the new SCF the whole budget a second time:
 * SAPT(DFT) holds the dimer and both monomers, and the GRAC procedure holds the
 * neutral wavefunction while the cation SCF runs.  Objects that park a budget-sized
 * buffer therefore keep a MemoryClaim, so the driver can subtract what is already
 * spoken for and divide up only what is genuinely left.
 *
 * A claim allocates nothing and polices nothing -- it only reports.  The unit is
 * doubles, matching JK::set_memory() and VBase::build_collocation_cache().
 */
class MemoryClaim {
   public:
    MemoryClaim() = default;
    ~MemoryClaim() { set(0); }

    /// A claim stands for one particular buffer, so it is no more copyable than the buffer is.
    MemoryClaim(const MemoryClaim&) = delete;
    MemoryClaim& operator=(const MemoryClaim&) = delete;

    /// Report the doubles this holder has resident *now*; pass 0 once the buffer is gone.
    void set(size_t doubles) {
        total().fetch_add(static_cast<std::int64_t>(doubles) - static_cast<std::int64_t>(held_),
                          std::memory_order_relaxed);
        held_ = doubles;
    }

    /// Doubles this claim last reported.
    size_t held() const { return held_; }

    /// Doubles reported by every live claim in the process.
    static size_t committed() {
        const std::int64_t tally = total().load(std::memory_order_relaxed);
        return tally > 0 ? static_cast<size_t>(tally) : 0;
    }

   private:
    size_t held_ = 0;

    static std::atomic<std::int64_t>& total() {
        static std::atomic<std::int64_t> tally{0};
        return tally;
    }
};

}  // namespace psi

#endif  // _psi_src_lib_libpsi4util_memory_ledger_h_
