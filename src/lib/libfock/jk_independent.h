/**
 * jk_independent.h
 * 
 * Created by Bill March, 2/25/13
 *
 * This is a version of the JK class that takes two classes as template 
 * arguments and uses them to create the J and K matrices. This allows us
 * to specify two different methods for the two matrices, such as CFMM and
 * LinK.
 *
 * Problem with this: how do I instantiate it? I'll have to add to RHF or 
 * where ever this happens.
 */


#ifndef _jk_independent_h
#define _jk_independent_h

#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libfock/apps.h>
#include <libfock/jk.h>

#include <libmints/typedefs.h>
#include <libmints/sieve.h>

namespace psi {
  
  template <class JDriver, class KDriver>
  class JKIndependent : public JK
  {
    
  protected:
    
    JDriver j_driver_;
    KDriver k_driver_;
    
    // true if we need to call j_driver_ for J computation and k_driver_ for K.
    // if false, then we assume j_driver does everything we need it to
    bool do_separately_;
    
    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const { return allow_desymmetrization_; }
    /// Setup integrals, files, etc
    virtual void preiterations();
    /// Compute J/K for current C/D
    virtual void compute_JK();
    /// Delete integrals, files, etc
    virtual void postiterations();
    
    /// Common initialization
    void common_init();
    
  public:
    // => Constructors < = //
    
    JKIndependent(boost::shared_ptr<BasisSet> primary, bool do_separately);
    /// Destructor
    virtual ~JKIndependent();
    
    // => Accessors <= //
    
    /**
     * Print header information regarding JK
     * type on output file
     */
    virtual void print_header() const;
    
  }; // class
  
} // namespace

// definitions
#include "jk_independent_impl.h"

#endif




