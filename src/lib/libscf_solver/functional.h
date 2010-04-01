#ifndef libscf_solver_Functional_H
#define libscf_solver_Functional_H
/*
 *  Functional.h
 *  Definition of class Functional for use in KS-DFT
 *
 *  Created by Robert Parrish on 02/24/10.
 *
 */
#include <libmints/properties.h>
#include <libutil/ref.h>
#include <stdlib.h>
#include <string>
using namespace psi;

namespace psi { namespace scf {
class FunctionalFactory;

/*! \ingroup SCF */
//! Functional Interface definition 
class Functional {
public:
	/** Constructor
	* Does nothing, returns NO object
	*/
	Functional() {/* Does nothing */}
	/** Constructor
	* Does nothing, returns object
	*/
	Functional(int) {/* Does nothing */}
	/** Destructor
	* Does nothing
	*/
	~Functional() {/* Does nothing */}
	/** Does this functional depend on spin?
	* @return true if so, false otherwise
	*/
	virtual bool hasSpinDependence() {return false; }
	/** Does this functional have exact exchange?
	* @return true if so, false otherwise
	*/
	virtual bool hasExactExchange() {return false; }
	/** Exact exchange coefficient (0.0 if no exact exchange)
	* @return coefficient of exact exchange to add to F^{KS}
	*/
	virtual double getExactExchangeCoefficient() {return 0.0;}
	/** Exact Coulomb coefficient (1.0 if most of the time)
	* @return coefficient of Coulomb J matrix to add to F^{KS}
	*/
	virtual double getExactCoulombCoefficient() {return 1.0;}
	/** Does this functional need electron density?
	* @return true if so, false otherwise
	*/
	virtual bool needsDensity() {return true; }
	/** Does this functional need electron density gradients?
	* @return true if so, false otherwise
	*/
	virtual bool needsDensityGradient() {return false; }
	/** Does this functional need electron density jacobians?
	* @return true if so, false otherwise
	*/
	virtual bool needsDensityJacobian() {return false; } 
	/** Does this functional need electron density laplacians?
	* @return true if so, false otherwise
	*/
	virtual bool needsDensityLaplacian() {return false; } 
	/** Does this functional need kinetic energy density?
	* @return true if so, false otherwise
	*/
	virtual bool needsKEDensity()  {return false; } 
	/** Does this functional depend on spin?
	* @return true if so, false otherwise
	*/
	virtual double getValue(shared_ptr<Properties> prop)  {return 0.0;}
	/** Does this functional depend on spin?
	* @return true if so, false otherwise
	*/
	virtual double getGradientA(shared_ptr<Properties> prop)  {return 0.0;}
	/** Functional name
	* @return functional or alias name ie: 'B3LYP'
	*/
	virtual string getName() {return string("None");}
	/** Functional description
	* @return description of functional or especially
	* substituent functionals of aliases
	*/
	virtual string getDescription() {return string("None"); }
	/** Functional citation
	* @return citation of functional
	*/
	virtual string getCitation() {return string("None"); }
};
typedef shared_ptr<Functional> SharedFunctional;
}}
#endif
