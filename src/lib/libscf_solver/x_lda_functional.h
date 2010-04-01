/*
 *  X_LDA_Functional.h
 *  Definition of class X_LDA_Functional for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/24/10.
 *
 */

#include <libmints/properties.h>
#include <libscf_solver/functional.h>
#include <stdlib.h>
#include <string>
using namespace psi;

namespace psi { namespace scf {
/*! \ingroup SCF */
//! X_LDA_Functional definition 
class X_LDA_Functional: public Functional {
public:
	/** Constructor
	* Does nothing, returns object
	*/
	X_LDA_Functional() {/* Does nothing */}
	/** Constructor
	* Does nothing, returns object
	*/
	X_LDA_Functional(int) {/* Does nothing */}
	/** Destructor
	* Does nothing
	*/
	~X_LDA_Functional() {/* Does nothing */}
	/** 
	* @return false, no spin dependence
	*/
	const bool hasSpinDependence() const { return false; }
	/** Does this functional have exact exchange?
	* @return false, no exact exchange
	*/
	const bool hasExactExchange() const {return false;}
	/** Exact exchange coefficient (0.0 if no exact exchange)
	* @return 0.0
	*/
	const double getExactExchangeCoefficient() const {return 0.0;}
	/** Exact Coulomb coefficient (1.0 if usual)
	* @return 1.0
	*/
	const double getExactCoulombCoefficient() const {return 1.0;}
	/** Does this functional need electron density?
	* @return true, needs density
	*/
	const bool needsDensity() const {return true; }
	/** Does this functional need electron density gradients?
	* @return true if so, false otherwise
	*/
	const bool needsDensityGradient() const {return false; }
	/** Does this functional need electron density jacobians?
	* @return true if so, false otherwise
	*/
	const bool needsDensityJacobian() const {return false; }
	/** Does this functional need electron density laplacians?
	* @return true if so, false otherwise
	*/
	const bool needsDensityLaplacian() const {return false; }
	/** Does this functional need kinetic energy density?
	* @return true if so, false otherwise
	*/
	const bool needsKEDensity() const {return false; }
	/** Functional Value at specified point properties 
	* @return functional value
	*/
	double getValue(shared_ptr<Properties> prop);
	/** Functional Gradient w.r.t. spin up denisty
	*  at specified point properties 
	* @return functional gradient
	*/
	double getGradientA(shared_ptr<Properties> prop);
	/** Functional name
	* @return functional or alias name ie: 'B3LYP'
	*/
	string getName();
	/** Functional description
	* @return description of functional or especially
	* substituent functionals of aliases
	*/
	string getDescription();
	/** Functional citation
	* @return citation of functional
	*/
	string getCitation();	
};

}}
