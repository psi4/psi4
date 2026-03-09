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

#include "psi4/libpsi4util/exception.h"

#ifndef _MSC_VER
#include <execinfo.h>
#include <cxxabi.h>
#endif

#include <array>
#include <cstdlib>

namespace psi {

PsiException::PsiException(std::string msg, const char *_file, int _line) noexcept : runtime_error(msg) {
    file_ = _file;
    line_ = _line;

    std::stringstream message;
    message << std::endl << "Fatal Error: " << msg << std::endl;
    message << "Error occurred in file: " << file_ << " on line: " << line_ << std::endl;

// Disable stack trace printing on Windows
#ifndef _MSC_VER

    constexpr size_t stacksize = 5;
    std::array<void *, stacksize> Stack;
    char **strings;
    int size = backtrace(Stack.data(), stacksize);
    int status = -1;

    message << "The most recent " << (size < stacksize ? size : stacksize) << " function calls were:" << std::endl << std::endl;
    strings = backtrace_symbols(Stack.data(), size);

    for (int i = 0; i < size; i++) {
        // This part from https://panthema.net/2008/0901-stacktrace-demangled/
        char *begin_name = nullptr, *begin_offset = nullptr, *end_offset = nullptr;
        for (char *p = strings[i]; *p; ++p) {
            if (*p == '(')
                begin_name = p;
            else if (*p == '+')
                begin_offset = p;
            else if (*p == ')' && begin_offset) {
                end_offset = p;
                break;
            }
        }
        if (begin_name && begin_offset && end_offset && begin_name < begin_offset) {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';
            char *demangledname = abi::__cxa_demangle(begin_name, 0, 0, &status);
            if (status == 0) message << demangledname << std::endl;
            ::free(demangledname);
        }
    }
#endif

    msg_ = message.str();
}

PsiException::PsiException(const PsiException &copy) noexcept
    : runtime_error(copy.msg_), msg_(copy.msg_), file_(strdup(copy.file_)), line_(copy.line_) {}

void PsiException::rewrite_msg(std::string msg) noexcept { msg_ = msg; }

const char *PsiException::what() const noexcept {
    // stringstream sstr;
    // sstr << msg_ << "\n";
    // sstr << location();
    // return sstr.str().c_str();
    return msg_.c_str();
}

const char *PsiException::file() const noexcept { return file_; }

int PsiException::line() const noexcept { return line_; }

std::string PsiException::location() const noexcept {
    std::stringstream sstr;
    sstr << "file: " << file_ << "\n";
    sstr << "line: " << line_;
    return sstr.str();
}

PsiException::~PsiException() noexcept {}

SanityCheckError::SanityCheckError(std::string message, const char *_file, int _line) noexcept
    : PsiException(message, _file, _line) {
    std::stringstream sstr;
    sstr << "sanity check failed! " << message;
    rewrite_msg(sstr.str());
}

SanityCheckError::~SanityCheckError() noexcept {}

SystemError::SystemError(int eno, const char *_file, int _line) noexcept : PsiException("", _file, _line) {
    std::stringstream sstr;
    sstr << "SystemError:  " << strerror(eno);
    rewrite_msg(sstr.str());
}

SystemError::~SystemError() noexcept {}

InputException::InputException(std::string msg, std::string param_name, int value, const char *_file, int _line) noexcept
    : PsiException(msg, _file, _line) {
    write_input_msg<int>(msg, param_name, value);
}

InputException::InputException(std::string msg, std::string param_name, std::string value, const char *_file,
                               int _line) noexcept
    : PsiException(msg, _file, _line) {
    write_input_msg<std::string>(msg, param_name, value);
}

InputException::InputException(std::string msg, std::string param_name, double value, const char *_file,
                               int _line) noexcept
    : PsiException(msg, _file, _line) {
    write_input_msg<double>(msg, param_name, value);
}

InputException::InputException(std::string msg, std::string param_name, const char *_file, int _line) noexcept
    : PsiException(msg, _file, _line) {
    write_input_msg<std::string>(msg, param_name, "in input");
}

template <class T>
void InputException::write_input_msg(std::string msg, std::string param_name, T value) noexcept {
    std::stringstream sstr;
    sstr << msg << "\n";
    sstr << "value " << value << " is incorrect"
         << "\n";
    sstr << "please change " << param_name << " in input";
    rewrite_msg(sstr.str());
}

template <class T>
StepSizeError<T>::StepSizeError(std::string value_name, T max, T actual, const char *_file, int _line) noexcept
    : LimitExceeded<T>(value_name + " step size", max, actual, _file, _line) {}

template <class T>
StepSizeError<T>::~StepSizeError() noexcept {}

template <class T>
ResourceAllocationError<T>::ResourceAllocationError(std::string resource_name, T max, T actual, const char *_file,
                                                    int _line) throw()
    : LimitExceeded<size_t>(resource_name, max, actual, _file, _line) {}

template <class T>
ResourceAllocationError<T>::~ResourceAllocationError() throw() {}

FeatureNotImplemented::FeatureNotImplemented(std::string module_name, std::string feature_name, const char *_file,
                                             int _line) throw()
    : PsiException("psi exception", _file, _line) {
    std::stringstream sstr;
    sstr << feature_name << " not implemented in " << module_name;
    rewrite_msg(sstr.str());
}

FeatureNotImplemented::~FeatureNotImplemented() throw() {}
}
