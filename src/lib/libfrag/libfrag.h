/*
 * libfrag.h
 *
 *  Created on: May 13, 2014
 *      Author: richard
 */

#ifndef LIBFRAG_H_
#define LIBFRAG_H_

#include <string>
#include <boost/python.hpp>
#include "FragOptions.h"
#include "LibFragTypes.h"
#include "libpsio/MOFile.h"

namespace LibFrag {
class GMBE;
class BSSEer;



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
       *   Systems[0]={Monomers}
       *   Systems[1]={n-mers, positive intersections}
       *   Systems[2]={negative intersections}
       *
       *   Python doesn't need to know what it is running...
       *
       *
       */
      std::vector<NMerSet> Systems;

      ///The actual energy expansion
      boost::shared_ptr<GMBE> Expansion;

      ///The way we are correcting for BSSE (if we are)
      boost::shared_ptr<BSSEer> BSSEFactory;

      ///The factory for adding caps
      boost::shared_ptr<Capper> CapFactory;

      ///The list of options
      FragOptions DaOptions;

      ///A vector of the MO files
      std::vector<psi::MOFile> MOFiles;
   public:

      ///Sets up the class and makes the fragments
      void Fragment_Helper(boost::python::str& FragMethod, const int N,
            boost::python::str& EmbedMethod,
            boost::python::str& CapMethod,
            boost::python::str& BSSEMethod);

      void NMer_Helper(const int N);

      void Embed_Helper(boost::python::str& EmbedMethod);

      ///Returns a string: cap_symbol <carts>\n for each cap in "N"-mer "NMer"
      std::string Cap_Helper(const int NMer,const int N);

      ///Returns the highest n-body approximate energy available
      double CalcEnergy(boost::python::list& Energies);

      int GetNNMers(const int i) {
         if (i<Systems.size()) return Systems[i].size();
         else return 0;
      }

      int GetNFrags() {
         return GetNNMers(0);
      }
      boost::python::list GetNMerN(const int NMer, const int N);

      boost::python::list GetGhostNMerN(const int NMer, const int N);

      /** \brief Uses the new libpsio/MOFile class, which is an
       *   object version of the MO coefficient file.
       *
       *   If we are using the MBE and have not severed any bonds
       *   (laziness on the MBE part, severing bonds makes things complicated)
       *   this will store our MO coefficient files.  Then when we run the
       *   N-mers we call WriteMOs(N,x), where x is the x-th N-Mer, and
       *   we use the direct sum as an initial SCF guess.
       */
      void ReadMOs();

      ///Constructs a MOFile that is the direct sum of the N fragments in
      ///the "x"-th "N"-mer
      void WriteMOs(const int N, const int x);
      ///Technically need to check for embedding, but it's not ready yet
      int RunFrags();

      int IsGMBE();

      /** \brief This function is called after a batch is run. It synchronizes
       *         all the processes
       *
       *   \param[in] Comm The communicator we are synching over
       *   \param[in] N The MBE we just completed (determines what data needs
       *                synched)
       */
      void Synchronize(boost::python::str& Comm,const int N);
      ~LibFragHelper() {
      }
};
}

#endif /* LIBFRAG_H_ */
