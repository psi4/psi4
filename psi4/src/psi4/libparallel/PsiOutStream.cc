/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include <stdarg.h>
#include <stdio.h>
namespace psi{
void Destructor(std::ostream* Stream_){}

void PsiOutStream::Buffer2Stream() {
    if (this->ImSpecial()) {
        (*Stream_) << Buffer_.str();
        this->Flush();
    }
    this->EmptyBuffer();
}

void PsiOutStream::MakeBanner(const std::string& message, const char delimiter, const int width) {
    std::string symbols(width, delimiter);
    (*this) << symbols << std::endl;
    int size = message.size();
    int nspaces = 2;  // Number of spaces between message, on each side
    if (size < width - 2 * (nspaces + 1)) {
        // Divide message in half, extra char goes left
        int lsize = (size - (size % 2)) / 2 + (size % 2);
        int rsize = size - lsize;
        // Number of times to print character
        int lchars = width / 2 - nspaces - lsize;
        int rchars = width / 2 - nspaces - rsize;
        std::string lsym(lchars, delimiter);
        std::string rsym(rchars, delimiter);
        std::string spaces(nspaces, ' ');
        (*this) << lsym << spaces << message << spaces << rsym << std::endl;
    }
    (*this) << symbols << std::endl;
}

PsiOutStream::PsiOutStream(SharedOutStream Stream) {
    if (this->ImSpecial()) {
        Stream_ = (Stream ? Stream : SharedOutStream(&std::cout, Destructor));
    }
}

void PsiOutStream::Printf(const char* format, ...) {
    // We don't know how long the fully expanded string is so lets guess our average print is about
    // a line
    const int guess_limit = 512;
    char* buffer = new char[guess_limit];
    va_list args;
    va_start(args, format);
    int left = vsnprintf(buffer, guess_limit, format, args);

    if (left < 0) {
        // Encoding error?!?
        throw PSIEXCEPTION("PsiOutStream: vsnprintf encoding error!");
    } else if (left >= guess_limit) {
        // Buffer was too small! Try again
        delete[] buffer;
        buffer = new char[left + 1];
        left = vsnprintf(buffer, guess_limit, format, args);
        if (left < 0) {
            throw PSIEXCEPTION("PsiOutStream: vsnprintf encoding error!");
        }
    }
    // Everything is cool

    va_end(args);
    Write2Buffer(buffer);
    delete[] buffer;
}

std::ostream& PsiOutStream::Write2Buffer(StreamManips fp) {
    Buffer_ << fp;
    this->Buffer2Stream();
    return Buffer_;
}
}
