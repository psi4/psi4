/*
 * Fragmenter.h
 *
 *  Created on: May 15, 2014
 *      Author: richard
 */

#ifndef FRAGMENTER_H_
#define FRAGMENTER_H_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "FragOptions.h"
namespace psi{
class Molecule;
}
typedef psi::Molecule Mol;
typedef boost::shared_ptr<Mol> SharedMol;
namespace LibFrag{
class MBEFrag;
typedef boost::shared_ptr<MBEFrag> SharedFrag;
typedef std::vector<SharedFrag> NMerSet;

class Fragmenter{
	public:
		virtual void Fragment(SharedMol& Mol2Frag,NMerSet& Fragments)=0;
		virtual ~Fragmenter(){}
};

class UDFragmenter:public Fragmenter{
	void Fragment(SharedMol& Mol2Frag,NMerSet& Fragments);
};
}


#endif /* FRAGMENTER_H_ */
