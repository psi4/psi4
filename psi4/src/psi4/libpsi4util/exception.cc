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
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */


#include <execinfo.h>
#include <cxxabi.h>
#include <vector>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include "psi4/libpsi4util/exception.h"
using namespace std;

namespace psi {

PsiException::PsiException(string msg,
                           const char *_file,
                           int _line) throw()
        : runtime_error(msg)
{
    file_ = _file;
    line_ = _line;

    std::stringstream message;
    message << std::endl << "Fatal Error: " << msg << std::endl;
    message << "Error occurred in file: " << file_ << " on line: " << line_ << std::endl;

    std::vector<void *> Stack(5);
    char **strings;
    int size = backtrace(&Stack[0], 5);
    int status = -1;

    message << "The most recent " << (size < 5 ? size : 5) << " function calls were:" << std::endl << std::endl;
    strings = backtrace_symbols(&Stack[0], size);

    for (int i = 0; i < size; i++) {
        //This part from https://panthema.net/2008/0901-stacktrace-demangled/
        char *begin_name = NULL, *begin_offset = NULL, *end_offset = NULL;
        for (char *p = strings[i]; *p; ++p) {
            if (*p == '(') begin_name = p;
            else if (*p == '+')begin_offset = p;
            else if (*p == ')' && begin_offset) {
                end_offset = p;
                break;
            }
        }
        if (begin_name && begin_offset && end_offset
            && begin_name < begin_offset) {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';
            char *demangledname =
                    abi::__cxa_demangle(begin_name, 0, 0, &status);
            if (status == 0)
                message << demangledname << std::endl;
            ::free(demangledname);
        }
    }

    msg_ = message.str();
}

PsiException::PsiException(const PsiException& copy) throw()
        : runtime_error(copy.msg_), msg_(copy.msg_), file_(strdup(copy.file_)),
          line_(copy.line_)
{
}

void
PsiException::rewrite_msg(string msg) throw()
{
    msg_ = msg;
}

const char *
PsiException::what() const throw()
{
    //stringstream sstr;
    //sstr << msg_ << "\n";
    //sstr << location();
    //return sstr.str().c_str();
    return msg_.c_str();
}

const char *
PsiException::file() const throw()
{
    return file_;
}

int
PsiException::line() const throw()
{
    return line_;
}

const char *
PsiException::location() const throw()
{
    stringstream sstr;
    sstr << "file: " << file_ << "\n";
    sstr << "line: " << line_;
    return sstr.str().c_str();
}

PsiException::~PsiException() throw()
{
}

SanityCheckError::SanityCheckError(
        string message,
        const char *_file,
        int _line
) throw()
        : PsiException(message, _file, _line)
{
    stringstream sstr;
    sstr << "sanity check failed! " << message;
    rewrite_msg(sstr.str());
}

SanityCheckError::~SanityCheckError() throw()
{ }

SystemError::SystemError(
        int eno,
        const char *_file,
        int _line
) throw()
        : PsiException("", _file, _line)
{
    stringstream sstr;
    sstr << "SystemError:  " << strerror(eno);
    rewrite_msg(sstr.str());
}

SystemError::~SystemError() throw()
{ }

InputException::InputException(
        string msg,
        string param_name,
        int value,
        const char *_file,
        int _line
) throw() : PsiException(msg, _file, _line)
{
    write_input_msg<int>(msg, param_name, value);
}

InputException::InputException(
        string msg,
        string param_name,
        string value,
        const char *_file,
        int _line
) throw() : PsiException(msg, _file, _line)
{
    write_input_msg<string>(msg, param_name, value);
}

InputException::InputException(
        string msg,
        string param_name,
        double value,
        const char *_file,
        int _line
) throw() : PsiException(msg, _file, _line)
{
    write_input_msg<double>(msg, param_name, value);
}

InputException::InputException(
        string msg,
        string param_name,
        const char *_file,
        int _line
) throw() : PsiException(msg, _file, _line)
{
    write_input_msg<string>(msg, param_name, "in input");
}

template<class T>
void InputException::write_input_msg(
        string msg,
        string param_name,
        T value
) throw()
{
    stringstream sstr;
    sstr << msg << "\n";
    sstr << "value " << value << " is incorrect" << "\n";
    sstr << "please change " << param_name << " in input";
    rewrite_msg(sstr.str());
}

template<class T>
StepSizeError<T>::StepSizeError(
        string value_name,
        T max,
        T actual,
        const char *_file,
        int _line) throw()
        : LimitExceeded<T>(value_name + " step size", max, actual, _file, _line)
{
}

template<class T>
StepSizeError<T>::~StepSizeError() throw()
{ }

template<class T>
MaxIterationsExceeded<T>::MaxIterationsExceeded(
        string routine_name,
        T max,
        const char *_file,
        int _line)  throw()
        : LimitExceeded<T>(routine_name + " iterations", max, max, _file, _line)
{
}

template<class T>
MaxIterationsExceeded<T>::~MaxIterationsExceeded() throw()
{ }

template<class T>
ConvergenceError<T>::ConvergenceError(
        string routine_name,
        T max,
        double _desired_accuracy,
        double _actual_accuracy,
        const char *_file,
        int _line) throw()
        : MaxIterationsExceeded<T>(routine_name + " iterations", max, _file, _line), desired_acc_(_desired_accuracy),
          actual_acc_(_actual_accuracy)
{
    stringstream sstr;
    sstr << "could not converge " << routine_name << ".  desired " << _desired_accuracy << " but got " <<
    _actual_accuracy << "\n";
    sstr << LimitExceeded<T>::description();
    PsiException::rewrite_msg(sstr.str());
}

template<class T>
ConvergenceError<T>::~ConvergenceError() throw()
{ }

template<class T>
double
ConvergenceError<T>::desired_accuracy() const throw()
{ return desired_acc_; }

template<class T>
double
ConvergenceError<T>::actual_accuracy() const throw()
{ return actual_acc_; }

template<>
ConvergenceError<int>::ConvergenceError(
        string routine_name,
        int max,
        double _desired_accuracy,
        double _actual_accuracy,
        const char *_file,
        int _line) throw()
        : MaxIterationsExceeded<int>(routine_name + " iterations", max, _file, _line), desired_acc_(_desired_accuracy),
          actual_acc_(_actual_accuracy)
{
    stringstream sstr;
    sstr << "could not converge " << routine_name << ".  desired " << _desired_accuracy << " but got " <<
    _actual_accuracy << "\n";
    sstr << LimitExceeded<int>::description();
    PsiException::rewrite_msg(sstr.str());
}

template<>
ConvergenceError<int>::~ConvergenceError() throw()
{ }

template<>
double
ConvergenceError<int>::desired_accuracy() const throw()
{ return desired_acc_; }

template<>
double
ConvergenceError<int>::actual_accuracy() const throw()
{ return actual_acc_; }

template<class T>
ResourceAllocationError<T>::ResourceAllocationError(
        string resource_name,
        T max,
        T actual,
        const char *_file,
        int _line)  throw()
        : LimitExceeded<size_t>(resource_name, max, actual, _file, _line)
{
}

template<class T>
ResourceAllocationError<T>::~ResourceAllocationError() throw()
{ }

FeatureNotImplemented::FeatureNotImplemented(
        string module_name,
        string feature_name,
        const char *_file,
        int _line
) throw()
        : PsiException("psi exception", _file, _line)
{
    stringstream sstr;
    sstr << feature_name << " not implemented in " << module_name;
    rewrite_msg(sstr.str());
}

FeatureNotImplemented::~FeatureNotImplemented() throw()
{
}

}
