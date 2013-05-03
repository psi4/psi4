/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <libyeti/exception.h>
#include <sstream>
#include <iostream>
#include <cstdlib>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

void
yeti::exception_crash(
    const std::string &excname,
    const std::string &msg,
    const std::string &file,
    int line
)
{
    cerr << excname << endl;
    cerr << msg << endl;
    cerr << file << endl;
    cerr << line << endl;
    abort();
}

YetiException::YetiException(
    const string& msg,
    const string& file,
    int line
) throw() : runtime_error(msg), 
    msg_(msg), 
    file_(file), 
    line_(line)
{
}

void
YetiException::rewrite_msg(const string& msg) throw()
{
    msg_.assign(msg);
}
    
const char* 
YetiException::what() const throw()
{
    stringstream sstr;
    sstr << msg_ << "\n";
    //sstr << location();
    sstr << "file: " << file_ << "\n";
    sstr << "line: " << line_;
    sstr << "      " << endl;

    return sstr.str().c_str();
}

std::string
YetiException::file() const throw()
{
    return file_;
}

int 
YetiException::line() const throw()
{
    return line_;
}

std::string
YetiException::location() const throw()
{
    stringstream sstr;
    sstr << "file: " << file_ << "\n";
    sstr << "line: " << line_;
    return sstr.str();
}

YetiException::~YetiException() throw()
{
}

SanityCheckError::SanityCheckError(
    const string& message,
    const string& file,
    int line
    ) throw() 
  : YetiException(message, file, line)
{
    stringstream sstr;
    sstr << "sanity check failed! " << message;
    rewrite_msg(sstr.str());
}

SanityCheckError::~SanityCheckError() throw() {}

InputException::InputException(
    const string& msg,
    const string& param_name,
    int value,
    const string& file,
    int line
) throw() : YetiException(msg, file, line)
{
    write_input_msg<int>(msg, param_name, value);
}

InputException::InputException(
    const string& msg,
    const string& param_name,
    const string& value,
    const string& file,
    int line
) throw() : YetiException(msg, file, line)
{
    write_input_msg<string>(msg, param_name, value);
}

InputException::InputException(
    const string& msg,
    const string& param_name,
    double value,
    const string& file,
    int line
) throw() : YetiException(msg, file, line)
{
    write_input_msg<double>(msg, param_name, value);
}

InputException::InputException(
    const string& msg,
    const string& param_name,
    const string& file,
    int line
) throw() : YetiException(msg, file, line)
{
    write_input_msg<string>(msg, param_name, "in input");
}

template <
    class T
> void
InputException::write_input_msg(
    const string& msg,
    const string& param_name,
    T value
) throw()
{
    stringstream sstr;
    sstr << msg << "\n";
    sstr << "value " << value << " is incorrect" << "\n";
    sstr << "please change " << param_name << " in input";
    rewrite_msg(sstr.str());
}

ResourceAllocationError::ResourceAllocationError(
    const string& resource_name,
    size_t max,
    size_t actual,
    const string& file,
    int line)  throw()
    : LimitExceeded<size_t>(resource_name, max, actual, file, line)
{
}

ResourceAllocationError::~ResourceAllocationError() throw() {}

FeatureNotImplemented::FeatureNotImplemented(
    const string& module_name,
    const string& feature_name,
    const string& file,
    int line
) throw() 
 : YetiException("exception", file, line)
{
    stringstream sstr;
    sstr << feature_name << " not implemented in " << module_name;
    rewrite_msg(sstr.str());
}

FeatureNotImplemented::~FeatureNotImplemented() throw()
{
}

