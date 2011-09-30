#ifndef __math_test_uhf_h__
#define __math_test_uhf_h__

#include <libpsio/psio.hpp>
#include "hf.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {
class Matrix;
class Vector;
namespace scf {

class UHF : public HF {
protected:
    boost::shared_ptr<Matrix> Dt_, Dtold_;
    boost::shared_ptr<Matrix> Ga_, Gb_, J_, Ka_, Kb_;

    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    virtual double compute_E();

    virtual void form_G();
    virtual void form_F();

    void save_fock();
    bool diis();

    bool test_convergency();
    void save_information();
    void compute_spin_contamination();

    void common_init();

    void save_density_and_energy();

    // Finalize memory/files
    virtual void finalize();

public:
    UHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    UHF(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~UHF();

    virtual bool restricted() const { return false; }
};

}}

#endif
