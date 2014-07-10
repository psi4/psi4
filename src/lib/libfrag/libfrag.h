/*
 * libfrag.h
 *
 *  Created on: May 13, 2014
 *      Author: richard
 */

#ifndef LIBFRAG_H_
#define LIBFRAG_H_

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>
#include "FragOptions.h"
namespace psi {
class Molecule;
}

namespace LibFrag {
class MBEFrag;
class GMBE;
class LibFragHelper {
   private:
      typedef boost::shared_ptr<MBEFrag> SharedFrag;
      typedef std::vector<SharedFrag> NMerSet;
      //Agreeing to store monomers, then dimers, etc.
      std::vector<NMerSet> Systems;
      boost::shared_ptr<GMBE> Expansion;
      FragOptions DaOptions;
   public:
      void Fragment_Helper(boost::python::str& FragMethod);
      void NMer_Helper(const int N);
      void Embed_Helper(boost::python::str& EmbedMethod);
      void Cap_Helper(boost::python::str& CapMethod);
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
      ~LibFragHelper() {
      }
};
}

#endif /* LIBFRAG_H_ */
