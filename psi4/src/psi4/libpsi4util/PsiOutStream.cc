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

#include "PsiOutStream.h"

#include "psi4/libpsi4util/exception.h"

#include <cstdio>
#include <fstream>
#include <cstdarg>

namespace psi {

PsiOutStream::PsiOutStream(std::string fname, std::ios_base::openmode mode) {
    if (fname == "") {
        stream_ = &std::cout;
        is_cout_ = true;
    } else {
        std::ofstream* tmpf = new std::ofstream(fname, mode);
        if (!tmpf->is_open()) {
            std::ostringstream oss;
            oss << "PsiOutStream: Failed to open file " << fname << ".";
            throw PSIEXCEPTION(oss.str());
        }

        stream_ = tmpf;
        is_cout_ = false;
    }

    // This is sort of big, but vsnprintf does not appear to work correctly on all OS's.
    // Im looking at YOU CentOS
    buffer_.resize(512000);
}

PsiOutStream::~PsiOutStream() {
    if (!is_cout_) {
        delete stream_;
    }
}

void PsiOutStream::Printf(const char* format, ...) {
    // We can check if the buffer is large enough
    va_list args;
    va_start(args, format);
    int left = vsnprintf(buffer_.data(), buffer_.size(), format, args);

    if (left < 0) {
        // Encoding error?!?
        va_end(args);
        throw PSIEXCEPTION("PsiOutStream: vsnprintf encoding error!");
    } else if (left >= buffer_.size()) {
        // Buffer was too small! Try again
        std::vector<char> tmp_buffer(left + 1);
        left = vsnprintf(tmp_buffer.data(), left + 1, format, args);
        if (left < 0) {
            va_end(args);
            throw PSIEXCEPTION("PsiOutStream: vsnprintf encoding error!");
        }
    }
    // Everything is cool

    va_end(args);
    (*stream_) << buffer_.data() << std::flush;
}
void PsiOutStream::Printf(std::string fp) { (*stream_) << fp << std::flush; }

void PsiOutStream::Flush() { stream_->flush(); }
}  // namespace psi
