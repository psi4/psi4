#ifndef __math_test_uhf_h__
#define __math_test_uhf_h__

#include <libpsio/psio.hpp>
#include "hf.h"

using namespace psi;

class UHF : public HF {
protected:
	RefMatrix Fa_, Fb_;
	RefMatrix Da_, Db_, Dt_;
	RefMatrix Ca_, Cb_;
	RefMatrix Ga_, Gb_;
	int use_out_of_core_;
	
	double *p_j_;
	double *p_k_;
	
	void allocate_PK();
	void form_initialF();
	void form_C();
	void form_D();
	double compute_initial_E();
	double compute_E();
	
	void form_G();
	void form_G_from_PK();
	void form_PK();
	void form_F();
	
	bool test_convergency();
	void save_information();
	
    void common_init();
public:
	UHF(psi::PSIO *psio, psi::Chkpt *chkpt = 0);
    UHF(Ref<psi::PSIO> &psio, Ref<psi::Chkpt> &chkpt);
	virtual ~UHF();
	
	double compute_energy();
};

#endif
