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

/*
 *  dfmp2_energy.h
 *
 *
 *
 */

#ifndef DFMP2ENERGY_H
#define DFMP2ENERGY_H

#include <libmints/wavefunction.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Options;
class PSIO;
class Chkpt;

namespace dfmp2 {

class DFMP2Energy : public Wavefunction {
public:
    DFMP2Energy(Options & options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    virtual ~DFMP2Energy() {}

    double compute_E();
};

}}

#endif
