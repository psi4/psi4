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

#ifndef PCM_H
#define PCM_H

#include <vector>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <boost/shared_ptr.hpp>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class PCM {
  public:
    enum CalcType {Total, NucAndEle, EleOnly};
    PCM() {};
    PCM(Options &options, boost::shared_ptr<PSIO> psio, int nirrep, boost::shared_ptr<BasisSet> basisset);
    ~PCM();
    double compute_E(SharedMatrix &D, CalcType type = NucAndEle);//Total); this should be the default (once an advanced option is available)
    SharedMatrix compute_V();
    SharedMatrix compute_V_electronic(); // This is needed by the CC code (and maybe the LR-SCF code)

  protected:
    // The number of tesserae in PCMSolver.
    int ntess_;
    // The number of irreducible tesserae in PCMSolver.
    int ntessirr_;
    /// A matrix to hold the charges and {x,y,z} coordinates of the tesserae
    SharedMatrix tess_Zxyz_;
    // A scratch array to hold the electronic potential values at the tesserae
    double * tess_pot_e_;
    // A scratch array to hold the nuclear potential values at the tesserae (unchanging)
    double * tess_pot_n_;
    // A scratch array to hold the total potential values at the tesserae
    double * tess_pot_;
    // A scratch array to hold the electronic charges at the tesserae
    double * tess_charges_e_;
    // A scratch array to hold the nuclear charges at the tesserae
    double * tess_charges_n_;
    // A scratch array to hold the charges at the tesserae
    double * tess_charges_;
    // Calculate energy using total charges and potentials
    double compute_E_total(SharedMatrix &D);
    // Calculate energy separating between charges and potentials
    double compute_E_separate(SharedMatrix &D);
    // Calculate electronic polarization energy (U_ee) only (for CC step)
    double compute_E_electronic(SharedMatrix &D);

    // Current basis set (for puream and nao/nso info)
    boost::shared_ptr<BasisSet> basisset_;

    // The AO->SO transformation matrix, which is used for transforming
    // matrices between pure and Cartesian representations.
    SharedMatrix my_aotoso_;

    // Factory for the electrostatic integrals
    PCMPotentialInt* potential_int_;

    // print level
    int pcm_print_;

};

typedef boost::shared_ptr<psi::PCM> SharedPCM;

} // psi
#endif
