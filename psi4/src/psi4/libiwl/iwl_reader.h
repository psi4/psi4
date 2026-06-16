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

#ifndef _psi_src_lib_libiwl_iwl_reader_h_
#define _psi_src_lib_libiwl_iwl_reader_h_

#include <memory>
#include "iwl.hpp"  // brings in PSI_API (via psio.hpp) and class IWL

namespace psi {

class PSIO;

/*!
** IWLReader -- sequential, range-for-friendly reader over an IWL integral file.
**
** Replaces the hand-rolled fetch/lastbuf/idx/inbuf/labels/values loop that is
** copy-pasted across ~15 call sites. Example:
**
**     IWLReader eri(psio(), PSIF_SO_TEI);
**     for (const auto& I : eri) {
**         // I.p, I.q, I.r, I.s, I.value
**     }
**
** Notes:
**  - The labels are returned exactly as stored on disk. In particular some
**    writers (psimrcc) encode a sentinel in the sign of p; callers that relied
**    on `std::abs((int)label[0])` must apply std::abs() to `I.p` themselves.
**  - Construction opens the file (PSIO_OPEN_OLD) and primes the first bucket,
**    matching the historical iwl_buf_init(..., readflag=1) behavior.
**  - On destruction the file is kept by default (set_keep(false) to erase),
**    matching the dominant iwl_buf_close(&Buf, 1) usage.
*/
class PSI_API IWLReader {
   public:
    struct Entry {
        int p, q, r, s;
        double value;
    };

    IWLReader(std::shared_ptr<PSIO> psio, int unit);

    void set_keep(bool keep) { buf_.set_keep_flag(keep); }

    //! Fetch the next integral. Returns false once the file is exhausted.
    bool next(Entry& e);

    //! Minimal input iterator so IWLReader works in a range-for.
    class iterator {
       public:
        iterator() = default;  // end sentinel
        explicit iterator(IWLReader* reader) : reader_(reader) { advance(); }

        const Entry& operator*() const { return cur_; }
        const Entry* operator->() const { return &cur_; }
        iterator& operator++() {
            advance();
            return *this;
        }
        bool operator==(const iterator& o) const { return done_ == o.done_; }
        bool operator!=(const iterator& o) const { return done_ != o.done_; }

       private:
        void advance() { done_ = !(reader_ && reader_->next(cur_)); }

        IWLReader* reader_ = nullptr;
        Entry cur_{};
        bool done_ = true;
    };

    iterator begin() { return iterator(this); }
    iterator end() { return iterator(); }

   private:
    IWL buf_;
    int cur_ = 0;  // our cursor within the current bucket
};

}  // namespace psi

#endif
