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

#ifndef FROZENNO_H
#define FROZENNO_H
#include"psi4-dec.h"
#include<libmints/wavefunction.h>

namespace boost {
template<class T> class shared_ptr;
}
namespace psi{
  class Wavefunction;
}

namespace psi{namespace fnocc{

// base class
class FrozenNO : public Wavefunction {
  public:
    FrozenNO(boost::shared_ptr<Wavefunction>wfn,Options&options);
    ~FrozenNO();

    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
    void ComputeNaturalOrbitals();

  protected:

    // mp2 energy in full basis
    double emp2;
    long int nso,nmo,ndocc,nvirt,nfzc,nfzv,ndoccact,nvirt_no;

    void common_init();
};

class DFFrozenNO : public FrozenNO {
  public:
    DFFrozenNO(boost::shared_ptr<Wavefunction>wfn,Options&options);
    ~DFFrozenNO();

    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
    void ComputeNaturalOrbitals();

  protected:

    void ModifyCa(double*Dab);
    void ModifyCa_occ(double*Dij);
    void BuildFock(long int nQ,double*Qso,double*F);
};

}}

#endif
