/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#ifndef _MSC_VER
#include <execinfo.h>
#include <cxxabi.h>
#endif

#include <vector>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include "psi4/libpsi4util/exception.h"

namespace psi {

PsiException::PsiException(std::string msg, const char *_file, int _line) noexcept : runtime_error(msg) {
    file_ = _file;
    line_ = _line;

    std::stringstream message;
    message << std::endl << "Fatal Error: " << msg << std::endl;
    message << "Error occurred in file: " << file_ << " on line: " << line_ << std::endl;

// Disable stack trace printing on Windows
#ifndef _MSC_VER

    std::vector<void *> Stack(5);
    char **strings;
    int size = backtrace(&Stack[0], 5);
    int status = -1;

    message << "The most recent " << (size < 5 ? size : 5) << " function calls were:" << std::endl << std::endl;
    strings = backtrace_symbols(&Stack[0], size);

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

const char *PsiException::location() const noexcept {
    std::stringstream sstr;
    sstr << "file: " << file_ << "\n";
    sstr << "line: " << line_;
    return sstr.str().c_str();
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
MaxIterationsExceeded<T>::MaxIterationsExceeded(std::string routine_name, T max, const char *_file, int _line) noexcept
    : LimitExceeded<T>(routine_name + " iterations", max, max, _file, _line) {}

template <class T>
MaxIterationsExceeded<T>::~MaxIterationsExceeded() noexcept {}

template <class T>
ConvergenceError<T>::ConvergenceError(std::string routine_name, T max, double _desired_accuracy,
                                      double _actual_accuracy, const char *_file, int _line) noexcept
    : MaxIterationsExceeded<T>(routine_name + " iterations", max, _file, _line),
      desired_acc_(_desired_accuracy),
      actual_acc_(_actual_accuracy) {
    std::stringstream sstr;
    sstr << "could not converge " << routine_name << ".  desired " << _desired_accuracy << " but got "
         << _actual_accuracy << "\n";
    sstr << LimitExceeded<T>::description();
    PsiException::rewrite_msg(sstr.str());
}

template <class T>
ConvergenceError<T>::~ConvergenceError() noexcept {}

template <class T>
double ConvergenceError<T>::desired_accuracy() const noexcept {
    return desired_acc_;
}

template <class T>
double ConvergenceError<T>::actual_accuracy() const noexcept {
    return actual_acc_;
}

template <>
ConvergenceError<int>::ConvergenceError(std::string routine_name, int max, double _desired_accuracy,
                                        double _actual_accuracy, const char *_file, int _line) throw()
    : MaxIterationsExceeded<int>(routine_name + " iterations", max, _file, _line),
      desired_acc_(_desired_accuracy),
      actual_acc_(_actual_accuracy) {
    std::stringstream sstr;
    sstr << "could not converge " << routine_name << ".  desired " << _desired_accuracy << " but got "
         << _actual_accuracy << "\n";
    sstr << LimitExceeded<int>::description();
    PsiException::rewrite_msg(sstr.str());
}

template <>
ConvergenceError<int>::~ConvergenceError() throw() {}

template <>
double ConvergenceError<int>::desired_accuracy() const throw() {
    return desired_acc_;
}

template <>
double ConvergenceError<int>::actual_accuracy() const throw() {
    return actual_acc_;
}

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
