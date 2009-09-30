#ifndef _psi_src_lib_libmints_wavefunction_h
#define _psi_src_lib_libmints_wavefunction_h

/*!
    \file libmints/wavefunction.h
    \ingroup MINTS
*/

#include <libmints/factory.h>
#include <libutil/ref.h>
#include <libmints/molecule.h>
#include <libmints/basisset.h>
#include <libpsio/psio.hpp>

#define MAX_IOFF 30000
extern int ioff[MAX_IOFF];

#define MAX_DF 500
extern double df[MAX_DF];

#define MAX_BC 20
extern double bc[MAX_BC][MAX_BC];

#define MAX_FAC 100
extern double fac[MAX_FAC];

#define INDEX2(i, j) ( i >= j ? ioff[i] + j : ioff[j] + i )
#define INDEX4(i, j, k, l) ( INDEX2( INDEX2(i, j), INDEX2(k, l) ) )

namespace psi {

//! Simple wavefunction base class.
class Wavefunction {
protected:

    shared_ptr<BasisSet> basisset_;
    shared_ptr<Molecule> molecule_;

    // PSI file access variables
    shared_ptr<PSIO> psio_;
    shared_ptr<Chkpt> chkpt_;

    MatrixFactory factory_;
    long int memory_;
    unsigned int debug_;
    double energy_threshold_;
    double density_threshold_;

private:
    // Wavefunction() {}
    void common_init();

public:
    /// Set the PSIO object. Note: Wavefunction assumes ownership of the object. DO NOT DELETE!
    Wavefunction(shared_ptr<PSIO> psio);
    Wavefunction(shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    
    virtual ~Wavefunction();

    /// Compute energy. Subclasses override this function to compute its energy.
    virtual double compute_energy();

    /// Initialize internal variables from checkpoint file.
    void init_with_chkpt();

    shared_ptr<Molecule> molecule() { return molecule_; }

    static void initialize_singletons();
};

}

#endif
