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

#include <liboptions/liboptions.h>
#include <libtrans/integraltransform.h>

#ifdef EXTERN
    #undef EXTERN
    #define EXTERN extern
#else
    #define EXTERN
#endif

//#define INDEX(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))
#define ID(x) ints.DPD_ID(x)

using namespace boost;

namespace psi{
class Wavefunction;

namespace mollerplesset2{
    // Nasty, nasty global variables.
    EXTERN int nirreps, nmo;
    EXTERN shared_ptr<PSIO> psio;
    EXTERN int *aOccOrbsPI, *bOccOrbsPI, *aVirOrbsPI, *bVirOrbsPI;
    EXTERN int *mopi, *clsdpi, *openpi, *frzcpi, *frzvpi;
    EXTERN int numAOcc, numBOcc, numAVir, numBVir;
    EXTERN double eSCF;
    EXTERN double *aOccEvals, *bOccEvals, *aVirEvals, *bVirEvals;
    
    EXTERN double plugin_mp2_unrestricted(boost::shared_ptr<Wavefunction> wfn, Options &options);
    EXTERN double plugin_mp2_restricted(boost::shared_ptr<Wavefunction> wfn, Options &options);
}} // Namespaces
