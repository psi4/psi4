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
	/** getExchangeFunctional 
	* @param x_id string for given exchange functional or alias
	* @param c_id string for given correlation functional 
	*  ie: 'B88' or 'B3LYP' would both return the B88 Exchange functional
	* @return functional object pointer corresponding to id
	*/
	Functional * getExchangeFunctional(std::string x_id, std::string c_id, int block_size);
	/** getCorrelationFunctional 
	* @param x_id string for given exchange functional or alias
	* @param c_id string for given correlation functional 
	*  ie: 'LYP' or 'B3LYP' would both return the LYP Correlation functional
	* @return functional object pointer corresponding to id
	*/
	Functional * getCorrelationFunctional(std::string x_id, std::string c_id, int block_size);
};

}}
#endif
