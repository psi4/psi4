#ifndef libscf_solver_FFactory_h
#define libscf_solver_FFactory_h
/*
 *  FunctionalFactory.h
 *  Definition of class FunctionalFactory for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/24/10.
 *
 */
#include <libscf_solver/functional.h>
#include <libmints/properties.h>
#include <stdlib.h>
#include <string>
using namespace psi;

namespace psi { namespace scf {
/*! \ingroup SCF */
//! Functional Interface definition 
class FunctionalFactory {
public:
	/** Constructor
	* Does nothing, returns object
	*/
	FunctionalFactory() {/* Does nothing */}
	/** Destructor
	* Does nothing
	*/
	~FunctionalFactory() {/* Does nothing */}
	/** getFunctional 
	* @param id string for given functional or alias
	*  ie: 'B3LYP'
	* @return functional object pointer corresponding to id
	*/
	Functional * getFunctional(string id, int block_size);
};

}}
#endif
