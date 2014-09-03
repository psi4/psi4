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

#ifndef LIBFRAG_H_
#define LIBFRAG_H_

#include <string>
#include <boost/python.hpp>
//#include "LibFragTypes.h"
#include "MBEFragSet.h"
#include "libpsio/MOFile.h"

namespace psi{
namespace LibFrag {
class GMBE;
class LibFragOptions;

typedef boost::python::str PyStr;
typedef boost::python::list PyList;


/** \brief This class is a mess from a c++ point of view; however my python
 *   is too terrible to do this right (mirror the classes comprising this
 *   class in python).
 *
 *   Basically any function that we may need in python from libfrag has
 *   an extension here.
 *
 */
class LibFragHelper {
   private:
      /** \brief These are the actual fragments, dimers, etc.
       *
       *   For the MBE this looks like:
       *   Systems[0]={Monomers}
       *   Systems[1]={Dimers}
       *   ...
       *   Systems[n-1]={n-mers}
       *
       *   For the GMBE this looks like:
       *   Systems[0]={Monomers and monomer intersections}
       *   Systems[1]={n-mers, n-mer intersections}
       *
       *   Python doesn't need to know what it is running it just
       *   needs to know what to loop over...
       *
       *
       *
       */
      std::vector<MBEFragSet> Systems_;

      ///The actual energy expansion
      boost::shared_ptr<GMBE> Expansion_;

      ///The list of options
      boost::shared_ptr<LibFragOptions> DaOptions_;

      ///Things we need to do MPI on
      ///@{
      ///A vector of the MO files
      std::vector<MOFile> MOFiles_;

      ///A vector of the charge sets we read in
      std::vector<boost::shared_ptr<double[]> > ChargeSets_;
      ///@}

      ///A wrapper for the broadcast/receive operations of Synch
      void BroadCastWrapper(const int i, const int j,std::string& comm,
            std::vector<psi::MOFile>& tempfiles,
            std::vector<boost::shared_ptr<double[]> >& tempChargeSets,
            bool bcast);

   public:

      ///Starts a timer so that we know how long it takes us to run
      LibFragHelper();

      ///Stops the timer
      ~LibFragHelper();

      ///Sets up the class and makes the fragments
      void Fragment_Helper(PyStr& FragMethod, const int N,
            PyStr& EmbedMethod,PyStr& CapMethod,PyStr& BSSEMethod);

      ///Makes "N"-mers (N>1)
      void NMer_Helper(const int N);

      boost::python::list Embed_Helper(const int N, const int x);

      ///Returns a string: cap_symbol <carts>\n for each cap in "N"-mer "NMer"
      std::string Cap_Helper(const int NMer,const int N);

      ///Prints the energies, in a nice, neat table
      void PrintEnergy(PyList& Energies, const int N,PyStr& Name);

      ///Returns the highest n-body approximate energy available
      double CalcEnergy(PyList& Energies,PyStr& Name);

      int GetNNMers(const int i) {
         if (i<Systems_.size()) return Systems_[i].size();
         else return 0;
      }

      int GetNFrags() {
         return GetNNMers(0);
      }

      boost::python::list GetNMerN(const int NMer, const int N);

      boost::python::list GetGhostNMerN(const int NMer, const int N);

      /** \brief Gathers Relevant data from after each calculation
       *
       *
       *   If we are using the MBE and have not severed any bonds
       *   (laziness on the MBE part, severing bonds makes things complicated)
       *   this will store our MO coefficient files.  Then when we run the
       *   N-mers we call WriteMOs(N,x), where x is the x-th N-Mer, and
       *   we use the direct sum as an initial SCF guess.
       *
       *   This also will read in point charges from wavefunction if we
       *   are using APC embedding
       *
       */
      void GatherData();

      ///Constructs a MOFile that is the direct sum of the N fragments in
      ///the "x"-th "N"-mer
      void WriteMOs(const int N, const int x);

      int IsGMBE();

      ///Are we iterating for whatever embedding reason:
      int Iterate(const int itr);

      ///Do we need to run the fragments
      bool RunFrags();

      /** \brief This function is called after a batch is run. It synchronizes
       *         all the processes
       *
       *   \param[in] Comm The communicator we are synching over
       *   \param[in] N The MBE we just completed (determines what data needs
       *                synched)
       */
      void Synchronize(boost::python::str& Comm,const int N,const int itr);

};
}}//End namespaces

#endif /* LIBFRAG_H_ */
