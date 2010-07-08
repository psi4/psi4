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
#include <stdlib.h>
#include <string>
using namespace psi;

namespace psi { namespace scf {
class FunctionalFactory;

/*! \ingroup SCF */
//! Functional Interface definition 
class Functional {
protected:
        double *value_; //LSDA
        double *gradA_; //LSDA
        double *gradB_; //LSDA
        double *gradAA_; //GGA
        double *gradAB_; //GGA
        double *gradBB_; //GGA
        int block_size_;
        // Assignment of this variable moved to functionalfactory.cc to be standards compliant
        const static double tol;

public:
	/** Constructor
	* Allocates memory
	*/
	Functional(int block_size) {block_size_ = block_size; }
	/** Destructor
	* Cleans up memory
	*/
	~Functional() { }

	/** Is this function a GGA? 
	* @return true if so, false otherwise
	*/
	virtual bool isGGA() {return false; }
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
	/** Does this functional need electron density hessians?
	* @return true if so, false otherwise
	*/
	virtual bool needsDensityHessian() {return false; } 
	/** Does this functional need electron density laplacians?
	* @return true if so, false otherwise
	*/
	virtual bool needsDensityLaplacian() {return false; } 
	/** Does this functional need kinetic energy density?
	* @return true if so, false otherwise
	*/
	virtual bool needsKEDensity()  {return false; } 
        /** Compute the RKS functional value based on the given properties
        * object over npoints (garanteed to be in the Properties object)
        */	
        virtual void computeFunctionalRKS(shared_ptr<Properties>){}
        /** Compute the UKS functional value based on the given properties
        * object over npoints (garanteed to be in the Properties object)
        */	
        virtual void computeFunctionalUKS(shared_ptr<Properties>){}
        /*
	* Functional Value
	*/
	virtual double *getValue()  {return value_;}
	/** 
	* Functional Gradient w.r.t. alpha density
	*/
	virtual double *getGradientA()  {return gradA_;}
	/** 
	* Functional Gradient w.r.t. beta density
	*/
	virtual double *getGradientB()  {return gradB_;}
	/** 
	* Functional Gradient w.r.t. sigma_{AB}
        * Needed for GGAs only!
	*/
	virtual double *getGradientAB()  {return gradAB_;}
	/** 
	* Functional Gradient w.r.t. sigma_{BB}
        * Needed for GGAs only!
	*/
	virtual double *getGradientBB()  {return gradBB_;}
	/** 
	* Functional Gradient w.r.t. sigma_{AA}
        * Needed for GGAs only!
	*/
	virtual double *getGradientAA()  {return gradAA_;}
	/** Functional name
	* @return functional or alias name ie: 'B3LYP'
	*/
	virtual std::string getName() {return std::string("None");}
	/** Functional description
	* @return description of functional or especially
	* substituent functionals of aliases
	*/
	virtual std::string getDescription() {return std::string("None"); }
	/** Functional citation
	* @return citation of functional
	*/
	virtual std::string getCitation() {return std::string("None"); }
};
typedef shared_ptr<Functional> SharedFunctional;
}}
#endif
