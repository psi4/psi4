#ifndef __rohf_psi_h__
#define __rohf_psi_h__

#include <vector>
#include <libpsio/psio.hpp>
#include "hf.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi { namespace scf {

class ROHF : public HF {
protected:
    boost::shared_ptr<Matrix> Feff_;
    boost::shared_ptr<Matrix> soFeff_;
    boost::shared_ptr<Matrix> Dt_old_;
    boost::shared_ptr<Matrix> Dt_;
    boost::shared_ptr<Matrix> Ct_;
    boost::shared_ptr<Matrix> Ga_;
    boost::shared_ptr<Matrix> Gb_;
    boost::shared_ptr<Matrix> Ka_;
    boost::shared_ptr<Matrix> Kb_;
    boost::shared_ptr<Matrix> moFa_;
    boost::shared_ptr<Matrix> moFb_;

    void form_initialF();
    void form_initial_C();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();

    void form_G();
    void form_F();

    void save_fock();
    bool diis();

    bool test_convergency();

    void save_information();
    // Finalize memory/files
    virtual void finalize();

    void save_density_and_energy();

    void common_init();
public:
    ROHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    ROHF(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~ROHF();
    virtual bool restricted() const { return false; }
};

}}

#endif
