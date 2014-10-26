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

#ifndef CIWAVE_H
#define CIWAVE_H

// Forward declarations
namespace boost {
template<class T> class shared_ptr;
}

namespace psi {
class Wavefunction;
class Options;
}

namespace psi { namespace detci {

class CIWavefunction : public Wavefunction
{
public:
    CIWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    virtual ~CIWavefunction();

    double compute_energy();
    PsiReturnType cas_update();

private:
    void init();
};

}}

#endif // CIWAVE_H

