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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_FRAGSYSGUTS_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_FRAGSYSGUTS_H_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "MBEProp.h"
#include "../LibFragMolecule.h"
#include "../SuperCell.h"
namespace psi{
namespace LibMolecule{
class Fragment;

class ScaleFactors:public LibFrag::MBEProp<double>{
   public:
      ScaleFactors(const int N):LibFrag::MBEProp<double>(N,"Scale Factors"){}
};

class NMers:public std::vector<std::vector<boost::shared_ptr<Fragment> > >{
   private:
      typedef std::vector<std::vector<boost::shared_ptr<Fragment> > > Base_t;
   public:
      ///This is a handy map for associating a given n-mer w/another
      PsiMap<SerialNumber,SerialNumber> SNLookUp_;
      ScaleFactors ScaleFacts_;
      typedef std::vector<boost::shared_ptr<Fragment> >::iterator NMerItr_t;
      typedef std::vector<boost::shared_ptr<Fragment> >::const_iterator cNMerItr_t;
      ///Makes this a container to hold up to N-mers (N=2 is dimers, etc.)
      NMers(const int N):Base_t(N),ScaleFacts_(N){}
      ///Returns the maximum number of unions we took
      int N()const{return Base_t::size();}
      int NNMers(const int N)const{return Base_t::operator[](N-1).size();}
      ///Returns an iterator to the start of the i-mers (i=2 is dimers, etc.)
      NMerItr_t NMerBegin(const int i){return Base_t::operator[](i-1).begin();}
      ///Returns an iterator to the end of the i-mers (i=2 is dimers, etc.)
      NMerItr_t NMerEnd(const int i){return Base_t::operator[](i-1).end();}
      ///Returns a constant iterator to the start of the i-mers (i=2 is dimers, etc.)
      cNMerItr_t NMerBegin(const int i)const{return Base_t::operator[](i-1).begin();}
      ///Returns a constant iterator to the end of the i-mers (i=2 is dimers, etc.)
      cNMerItr_t NMerEnd(const int i)const{return Base_t::operator[](i-1).end();}
      ///Adds the i-mer given by NMer (i=2 is dimers, etc.)
      void AddNMer(const int i, boost::shared_ptr<Fragment>& NMer){
         Base_t::operator[](i-1).push_back(NMer);}
      std::string PrintOut(const int Value)const;

};

class FragSysGuts{
   protected:
      NMers NMers_;
      ///Ultimately calls the next fxn with Set passed as Set1 and Set2
      void MakeNMers(const std::vector<boost::shared_ptr<Fragment> >& Set);
      /** \brief Fills NMers_ such that one monomer comes from Set1 and
       *         rest from Set2
       *
       *   We assume that Set1 is a subset of Set2.  Moreover we assume that
       *   Set2 is {Set1,Possible Other Frags}.  This is also the call that
       *   performs our distance screening, which is why the prior statement
       *   is important, as otherwise there is no telling which fragments
       *   will get thrown out.
       */
      void MakeNMers(const std::vector<boost::shared_ptr<Fragment> >& Set1,
                     const std::vector<boost::shared_ptr<Fragment> >& Set2);
      FragSysGuts(boost::shared_ptr<Molecule> System2Frag,const int N);
      FragSysGuts(boost::shared_ptr<SuperCell> System2Frag,const int N);
};

}}//End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_IMPLEMENTATIONS_FRAGSYSGUTS_H_ */
