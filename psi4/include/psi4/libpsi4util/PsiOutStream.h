/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef _psi_src_lib_libpsi4util_psioutstream_h_
#define _psi_src_lib_libpsi4util_psioutstream_h_

#include "psi4/pragma.h"
#include <vector>
#include <string>
#include <iostream>
#include <ostream>

namespace psi {

class PSI_API PsiOutStream {
   private:
    std::ostream* stream_;
    bool is_cout_;
    std::vector<char> buffer_;

   public:
    PsiOutStream(std::string fname = "", std::ios_base::openmode mode = std::ostream::trunc);
    ~PsiOutStream();

    void Printf(const char* fmt, ...);
    void Printf(std::string fp);
    void MakeBanner(std::string header);
    void Flush();

    std::ostream* stream() { return stream_; }

    // Incase we want to overload << again
    // template <class T>
    // PsiOutStream& operator<<(T&& x) {
    //     *stream_ << std::forward<T>(x);
    //     return (*this);
    // }

    // PsiOutStream& operator<<(std::ostream& (*oper)(std::ostream&)) {
    //     *stream_ << oper;
    //     return (*this);
    // }
};

}  // namespace psi
#endif
