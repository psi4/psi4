#ifndef three_index_H
#define three_index_H

#include <psi4-dec.h>
#include <psiconfig.h>

namespace psi {

class PSIO;
class BasisSet;
class Matrix;

enum ThreeStorage { Disk, SemiDirect, Direct };
enum ThreeAlgorithm { FitThenContract, ContractThenFit };
enum ThreeFitting { Cholesky, QR };

class ThreeIndexTensor {

protected:
    shared_ptr<BasisSet> primary_basis_;
    shared_ptr<BasisSet> zero_;
    shared_ptr<PSIO> psio_;
    double schwarz_cutoff_;
    unsigned long int memory_;
    int print_;
    int nthread_;

    int nbf_;
    unsigned long int nshell_pairs_;
    unsigned long int nfun_pairs_;

    double max_global_val_;
    int* schwarz_shells_;
    int* schwarz_funs_;

    long int* schwarz_shells_reverse_;
    long int* schwarz_funs_reverse_;

    double* schwarz_shell_vals_;
    double* schwarz_fun_vals_;

    void form_schwarz_ints();

public:
    ThreeIndexTensor(shared_ptr<PSIO>, shared_ptr<BasisSet>);
    virtual ~ThreeIndexTensor();
    virtual void finalize();

    // In doubles
    void set_memory(unsigned long int mem) { memory_ = mem; }
    void set_nthread(int nthread) { nthread_ = nthread; }
    void set_print(int print) { print_ = print; }

    virtual void print(FILE* out = outfile) = 0;
    void print_python() { print(); };

    void form_schwarz_sieve(double cutoff);
    unsigned long int get_nshell_pairs() const { return nshell_pairs_; }
    unsigned long int get_nfun_pairs() const { return nfun_pairs_; }
    // I_global = arr[2*I_local], J_global = arr[2*I_local + 1]
    // These are only defined up to nshell_pairs_ and nfun_pairs_, respectively
    int* get_schwarz_shells() const { return schwarz_shells_; }
    int* get_schwarz_funs() const { return schwarz_funs_; }
    // Canonical compound indexingi, -1 if not present
    long int* get_schwarz_shells_reverse() const { return schwarz_shells_reverse_; }
    long int* get_schwarz_funs_reverse() const { return schwarz_funs_reverse_; }

    shared_ptr<BasisSet> get_primary_basis() const { return primary_basis_; }
    // number of finished fitting vectors
    virtual int get_nfit() = 0;
    // number of raw fitting vectors
    virtual int get_nraw() = 0;
    // Convenience routines
    virtual void form_Qmn_disk() = 0;
    virtual void form_Qmi_disk(shared_ptr<Matrix> C_act_virt) = 0;
    virtual void form_Qma_disk(shared_ptr<Matrix> C_act_virt) = 0;
    virtual void form_Qii_disk(shared_ptr<Matrix> C1_act_occ, shared_ptr<Matrix> C2_act_occ) = 0;
    virtual void form_Qia_disk(shared_ptr<Matrix> C_act_occ, shared_ptr<Matrix> C_act_virt) = 0;
    virtual void form_Qaa_disk(shared_ptr<Matrix> C1_act_virt, shared_ptr<Matrix> C2_act_virt) = 0;

};

class DFTensor : public ThreeIndexTensor {

protected:
    shared_ptr<BasisSet> auxiliary_basis_;
    shared_ptr<BasisSet> poisson_basis_;

    shared_ptr<Matrix> fitting_metric_;

    bool poisson_;
    int ngaussian_;
    int npoisson_;
    int naux_;
    int nfin_;


public:
    DFTensor(shared_ptr<PSIO>, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> auxiliary, shared_ptr<BasisSet> poisson);
    DFTensor(shared_ptr<PSIO>, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> auxiliary);
    void common_init();
    virtual ~DFTensor();

    static shared_ptr<DFTensor> bootstrap_DFTensor();
    virtual void print(FILE* out);

    virtual int get_nfit();
    virtual int get_nraw();

    bool is_poisson() const { return poisson_; }
    shared_ptr<BasisSet> get_auxiliary_basis() const { return auxiliary_basis_; }
    shared_ptr<BasisSet> get_poisson_basis() const { return poisson_basis_; }

    shared_ptr<Matrix> form_fitting_metric();
    shared_ptr<Matrix> form_cholesky_metric();
    shared_ptr<Matrix> form_qr_metric(double max_cond = 1.0E-10);

    virtual void form_Qmn_disk();
    virtual void form_Qmi_disk(shared_ptr<Matrix> C_act_virt);
    virtual void form_Qma_disk(shared_ptr<Matrix> C_act_virt);
    virtual void form_Qii_disk(shared_ptr<Matrix> C1_act_occ, shared_ptr<Matrix> C2_act_occ);
    virtual void form_Qia_disk(shared_ptr<Matrix> C_act_occ, shared_ptr<Matrix> C_act_virt);
    virtual void form_Qaa_disk(shared_ptr<Matrix> C1_act_virt, shared_ptr<Matrix> C2_act_virt);
    virtual void form_Amn_disk();
    virtual void form_Ami_disk(shared_ptr<Matrix> C_act_virt);
    virtual void form_Ama_disk(shared_ptr<Matrix> C_act_virt);
    virtual void form_Aii_disk(shared_ptr<Matrix> C1_act_occ, shared_ptr<Matrix> C2_act_occ);
    virtual void form_Aia_disk(shared_ptr<Matrix> C_act_occ, shared_ptr<Matrix> C_act_virt);
    virtual void form_Aaa_disk(shared_ptr<Matrix> C1_act_virt, shared_ptr<Matrix> C2_act_virt);
    void disk_tensor(shared_ptr<Matrix> C1, shared_ptr<Matrix> C2, bool, bool, bool, const std::string &);
};

class DFSCFTensor : public DFTensor {

protected:

public:
    DFSCFTensor(shared_ptr<PSIO>, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> auxiliary, shared_ptr<BasisSet> poisson);
    DFSCFTensor(shared_ptr<PSIO>, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> auxiliary);

    virtual ~DFSCFTensor();



};

class CDTensor : public ThreeIndexTensor {

protected:
    double delta_;

public:
    CDTensor(shared_ptr<PSIO>, shared_ptr<BasisSet>, double delta);
    virtual ~CDTensor();

    virtual void print(FILE* out);

    virtual int get_nfit();
    virtual int get_nraw();
    void form_Qia();

    virtual void form_Qmn_disk(){}
    virtual void form_Qmi_disk(shared_ptr<Matrix> C_act_virt){}
    virtual void form_Qma_disk(shared_ptr<Matrix> C_act_virt){}
    virtual void form_Qii_disk(shared_ptr<Matrix> C1_act_occ, shared_ptr<Matrix> C2_act_occ){}
    virtual void form_Qia_disk(shared_ptr<Matrix> C_act_occ, shared_ptr<Matrix> C_act_virt){}
    virtual void form_Qaa_disk(shared_ptr<Matrix> C1_act_virt, shared_ptr<Matrix> C2_act_virt){}
};

}
#endif
