#ifndef LIBSCF_X_B3_FUNCTIONAL
#define LIBSCF_X_B3_FUNCTIONAL
/*
 *  X_B3_Functional.h
 *  Definition of class X_B3_Functional for use in KS-DFT
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
//! X_B3_Functional definition 
class X_B3_Functional: public Functional {
protected:
        //const static double ao_ = 0.2;
        //const static double ax_ = 0.72;
public:
	/** Constructor
	* Allocates required arrays, returns object
	*/
	X_B3_Functional(int); 
	/** Destructor
	* Frees arrays
	*/
	~X_B3_Functional();
	/** Is this functional a GGA?
	* @return true, it's the first GGA 
	*/
	bool isGGA() {return true;}
	/** Does this functional have exact exchange?
	* @return true, finally hybrid
	*/
	bool hasExactExchange() {return true;}
	/** Exact exchange coefficient (0.0 if no exact exchange)
	* @return 0.20 (80% LDA, 20% exact exchange)
	*/
	double getExactExchangeCoefficient() {return 0.20;}
	/** Exact Coulomb coefficient (1.0 if usual)
	* @return 1.0
	*/
	double getExactCoulombCoefficient() {return 1.0;}
	/** Does this functional need electron density?
	* @return true, needs density
	*/
	bool needsDensity() {return true; }
	/** Does this functional need electron density gradients?
	* @return true if so, false otherwise
	*/
	bool needsDensityGradient() {return true; }
	/** Does this functional need electron density jacobians?
	* @return true if so, false otherwise
	*/
	bool needsDensityHessian() {return false; }
	/** Does this functional need electron density laplacians?
	* @return true if so, false otherwise
	*/
	bool needsDensityLaplacian() {return false; }
	/** Does this functional need kinetic energy density?
	* @return true if so, false otherwise
	*/
	bool needsKEDensity() {return false; }
	/** Functional Value at specified point properties
        * for RKS type SCF 
	* 
	*/
	void computeFunctionalRKS(shared_ptr<Properties> prop);
	/** Functional Value at specified point properties
        * for UKS type SCF 
	* 
	*/
	void computeFunctionalUKS(shared_ptr<Properties> prop);
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
#endif
