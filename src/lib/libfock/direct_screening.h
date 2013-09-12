//
//  direct_screening.h
//  
//
//  Created by William March on 3/6/13.
//
//

#ifndef __psi_libfock_direct_screening__
#define __psi_libfock_direct_screening__

#include <cfloat>

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
  
  class DirectScreening {
    
  protected:
    
    boost::shared_ptr<BasisSet> basis_;
    
    bool do_J_;
    bool do_K_;
    bool do_wK_;
    
    //std::vector<SharedMatrix >& D_;
    std::vector<SharedMatrix> D_;
    
    // the coulomb and exchange matrices, to be filled in and passed back out
    std::vector<SharedMatrix> J_;
    std::vector<SharedMatrix> K_;
    //std::vector<SharedMatrix> wK_;
    
    // ERIs
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri_;
    /// Integral factory (must be retained for Spherical Transforms)
    boost::shared_ptr<IntegralFactory> factory_;
    /// ERI Sieve
    boost::shared_ptr<ERISieve> sieve_;
    
    ///////////////////// functions ////////////////////////////
  
    
  public:
    
    DirectScreening(boost::shared_ptr<BasisSet> basis_in,
         std::vector<SharedMatrix>& density_in);
    
    ~DirectScreening();
    
    void Compute();
    
    // Assuming this always gets called before Compute()
    void Update(const std::vector<SharedMatrix>& D_new);
    
    std::vector<SharedMatrix>& J();
    
    std::vector<SharedMatrix>& K();
    
    void set_do_J(bool do_it);
    
    void set_do_K(bool do_it);
    
    void print_header() const;
    
    
    
  }; // class

}


#endif /* defined(____direct_screening__) */
