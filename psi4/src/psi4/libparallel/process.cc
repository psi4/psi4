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

#include "psi4/psi4-dec.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/extern.h"

#include <regex>

//MKL Header
#ifdef USING_LAPACK_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

// Apple doesn't provide the environ global variable
#if defined(__APPLE__) && !defined(environ)
#   include <crt_externs.h>
#   define environ (*_NSGetEnviron())
#endif
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace psi {

Process::Environment Process::environment;
const std::string empty_;

// Need to split each entry by the first '=', left side is key, right the value
void Process::Environment::initialize()
{
    // Setup defaults
    environment_["PSI_SCRATCH"] = "/tmp/";
    environment_["PSI_DATADIR"] = "";


    // Go through user provided environment overwriting defaults if necessary
    int i = 0;
    if (environ) {
        while (environ[i] != NULL) {
            std::vector <std::string> strs = split(environ[i], "=");
            if (strs.size() > 1) {
                environment_[strs[0]] = strs[1];
            }
            ++i;
        }
    }

    nthread_ = 1;

#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif
}

void Process::Environment::set_n_threads(int nthread)
{
    nthread_ = nthread;
#ifdef _OPENMP
    omp_set_num_threads(nthread_);
#endif
#ifdef USING_LAPACK_MKL
    mkl_set_num_threads(nthread_);
#endif

    // HACK: TODO: CC-pthread codes should ask us how many threads
    // No, this didn't work, back this out for now (and we won't need
    // this once the final solution is in anyway)  --CDS
    // Process::environment.options.set_global_int("NUM_THREADS",nthread);
}

const std::string &Process::Environment::operator()(const std::string &key) const
{
    // Search for the key:
    std::map<std::string, std::string>::const_iterator it = environment_.find(key);

    if (it == environment_.end())
        return empty_;      // Not found return empty std::string.
    else
        return it->second;  // Found, return the value
}

std::string Process::Environment::operator()(const std::string &key)
{
    // Search for the key:
    std::map<std::string, std::string>::const_iterator it = environment_.find(key);

    if (it == environment_.end())
        return std::string(); // Not found return empty std::string.
    else
        return it->second; // Found, return the value
}

const std::string Process::Environment::set(const std::string &key, const std::string &value)
{
    const std::string &old = operator()(key);
    environment_[key] = value;

    // Attempt to set the variable in the system.
#ifdef HAVE_SETENV
    setenv(key.c_str(), value.c_str(), 1);
#elif HAVE_PUTENV
    size_t len = key.length() + value.length() + 2; // 2 = 1 (equals sign) + 1 (null char)
    char *str = new char[len];
    sprintf(str, "%s=%s", key.c_str(), value.c_str());
    // we give up ownership of the memory allocation
    putenv(str);
#else
#pragma error setenv and putenv not available.
#endif

    return std::string();
}

void Process::Environment::set_molecule(const std::shared_ptr <Molecule> &molecule)
{
    molecule_ = molecule;
}

std::shared_ptr <Molecule> Process::Environment::molecule() const
{
    return molecule_;
}

void Process::Environment::set_legacy_molecule(const std::shared_ptr <Molecule> &legacy_molecule)
{
    legacy_molecule_ = legacy_molecule;
}

std::shared_ptr <Molecule> Process::Environment::legacy_molecule() const
{
    return legacy_molecule_;
}

void Process::Environment::set_legacy_wavefunction(const std::shared_ptr <Wavefunction> &legacy_wavefunction)
{
    legacy_wavefunction_ = legacy_wavefunction;
}

std::shared_ptr <Wavefunction> Process::Environment::legacy_wavefunction() const
{
    return legacy_wavefunction_;
}

Process::Environment Process::get_environment()
{
    return environment;
}

unsigned long int Process::Environment::get_memory() const
{ return memory_; }

void Process::Environment::set_memory(unsigned long int m)
{ memory_ = m; }

int Process::Environment::get_n_threads() const
{ return nthread_; }

void die_if_not_converged()
{
    outfile->Printf("Iterations did not converge.");

    if (Process::environment.options.get_bool("DIE_IF_NOT_CONVERGED"))
        throw PSIEXCEPTION("Iterations did not converge.");
    else {
        outfile->Printf("Iterations did not converge.");
    }
}

}
