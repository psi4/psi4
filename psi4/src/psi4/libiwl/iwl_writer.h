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

#ifndef _psi_src_lib_libiwl_iwl_writer_h_
#define _psi_src_lib_libiwl_iwl_writer_h_

#include <memory>
#include "iwl.hpp"  // brings in PSI_API (via psio.hpp) and class IWL

namespace psi {

class PSIO;

/*!
** IWLWriter -- RAII writer over an IWL integral file.
**
** Replaces the iwl_buf_init / write_value... / iwl_buf_flush / iwl_buf_close
** boilerplate scattered across the integral-sort and density-dump code. The
** buffer is flushed (with the last-buffer marker) and closed automatically on
** destruction. Example:
**
**     {
**         IWLWriter out(psio(), PSIF_MO_TPDM, 1.0e-16);
**         for (...) out.write(p, q, r, s, value);
**     }   // flush(lastbuf) + close happen here
**
** Notes:
**  - write() applies the file's cutoff and the 16-bit Label overflow check
**    (it delegates to the shared write_value path), so behavior matches the
**    legacy calls exactly.
**  - `dirac` requests the Mulliken->Dirac index swap (q<->r) some callers used
**    via iwl_buf_wrt_val(..., dirac=1).
**  - The file is kept by default; pass keep=false (or call set_keep(false)) to
**    erase it on close.
*/
class PSI_API IWLWriter {
   public:
    IWLWriter(std::shared_ptr<PSIO> psio, int unit, double cutoff, bool keep = true)
        : buf_(psio.get(), unit, cutoff, /*oldfile*/ 0, /*readflag*/ 0) {
        buf_.set_keep_flag(keep);
    }

    ~IWLWriter() { close(); }

    // Flush the partial bucket (with the last-buffer marker) and close the
    // underlying file now. Idempotent; the destructor calls it if you don't.
    // Use this when the same psio unit must be reopened before the writer
    // would otherwise go out of scope.
    void close() {
        if (!closed_) {
            buf_.flush(/*lastbuf*/ 1);
            buf_.close();
            closed_ = true;
        }
    }

    void set_keep(bool keep) { buf_.set_keep_flag(keep); }

    void write(int p, int q, int r, int s, double value, bool dirac = false) {
        buf_.write_value(p, q, r, s, value, /*printflag*/ 0, std::string("outfile"), dirac ? 1 : 0);
    }

    //! Direct access to the underlying buffer for the rare caller that needs it.
    IWL& buffer() { return buf_; }

   private:
    IWL buf_;
    bool closed_ = false;
};

}  // namespace psi

#endif
