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

#include "atomic.h"
#include <libmints/mints.h>
#include <libfock/cubature.h>
#include <libfock/points.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <physconst.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

AtomicDensity::AtomicDensity() : 
    print_(1),
    debug_(0),
    bench_(0)
{
}
AtomicDensity::~AtomicDensity() 
{
}
boost::shared_ptr<AtomicDensity> AtomicDensity::build(const std::string& type, boost::shared_ptr<BasisSet> basis, boost::shared_ptr<Matrix> D, Options& options)
{
    // TODO
}

StockholderDensity::StockholderDensity() : AtomicDensity() 
{
}
StockholderDensity::~StockholderDensity()
{
}
void StockholderDensity::print_header() const 
{
    fprintf(outfile,"  ==> Stockholder Atomic Densities <==\n\n");
    fflush(outfile);
}
void StockholderDensity::compute()
{
    // TODO
}

}
