#ifndef CCSD_H
#define CCSD_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// file identifiers - and those not defined in psifiles.h
#include<psifiles.h>
#define PSIF_DCC_ABCI4 265
#define PSIF_DCC_ABCI  266
#define PSIF_DCC_ABCI2 267

// psi headers
#include"psi4-dec.h"
#include<libiwl/iwl.h>
#include<libpsio/psio.hpp>
#include<libpsio/psio.h>
#include<libmints/wavefunction.h>

long int Position(long int i,long int j);

namespace psi{ namespace qci{

class CoupledCluster: public Wavefunction{
  public:

    CoupledCluster(boost::shared_ptr<psi::Wavefunction> reference_wavefunction,Options &options);
    ~CoupledCluster();

    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }

  protected:

    void common_init();
    virtual void finalize();

    /// is t2 on disk or held in main memory?
    bool t2_on_disk;

    /// which cc method?
    bool mp2_only,mp3_only,mp4_only,isccsd;
    int ccmethod;

    /// flag for low-memory triples algorithm
    bool isLowMemory;

    /// define qci/cc/mp tasks
    void DefineTasks();
    void DefineLinearTasks();
    void DefineQuadraticTasks();
    long int ncctasks,nqtasks,nltasks;

    /// task parameters - not used currently
    struct CCTaskParams{
        int mtile,ntile,ktile;
    };
    CCTaskParams*CCParams,*CCSubParams1,*CCSubParams2,*QParams,*LParams;

    /// cc/qci/mp task
    struct CCTask{
        void(psi::qci::CoupledCluster::*func)(CCTaskParams);
        double flopcount;
        char*name;
    };
    CCTask*CCTasklist,*CCSubTasklist1,*CCSubTasklist2,*LTasklist,*QTasklist;

    /// solve qcisd/ccsd equations
    PsiReturnType CCSDIterations();

    /// SCS-MP2 function and variables
    void SCS_MP2();
    double emp2, emp2_os, emp2_ss, emp2_os_fac, emp2_ss_fac;

    /// SCS-CCSD function and variables
    void SCS_CCSD();
    void Local_SCS_CCSD();
    void Local_SCS_MP2();
    double eccsd, eccsd_os, eccsd_ss, eccsd_os_fac, eccsd_ss_fac;

    /// cc or qci (t)
    PsiReturnType triples();
    PsiReturnType lowmemory_triples();
    PsiReturnType local_triples();
    double et;

    /// mp4 triples
    void mp4_triples();
    double emp4_t;
  
    void WriteBanner();

    /// allocate memory
    void AllocateMemory();

    /// some cc/qci diagrams
    void CPU_t1_vmeai(CCTaskParams params);
    void CPU_t1_vmeni(CCTaskParams params);
    void CPU_t1_vmaef(CCTaskParams params);
    void CPU_I1ab(CCTaskParams params);
    void CPU_I1pij_I1ia_lessmem(CCTaskParams params);
    void CPU_I2p_abci_refactored_term2(CCTaskParams params);

    /// linear diagrams for mp4
    void I2iabj_linear(CCTaskParams params);
    void I2iajb_linear(CCTaskParams params);
    void I2ijkl_linear(CCTaskParams params);
    void I2piajk_linear(CCTaskParams params);
    void CPU_t1_vmeni_linear(CCTaskParams params);
    void CPU_t1_vmaef_linear(CCTaskParams params);
    void CPU_I2p_abci_refactored_term1_linear(CCTaskParams params);
    void CPU_t1_vmeai_linear(CCTaskParams params);
    void Vabcd1_linear(CCTaskParams params);
    void Vabcd2_linear(CCTaskParams params);

    /// linear diagrams for mp4
    void I2iabj_quadratic(CCTaskParams params);
    void I2ijkl_quadratic(CCTaskParams params);
    void I2iajb_quadratic(CCTaskParams params);
    void CPU_I1ab_quadratic(CCTaskParams params);
    void CPU_I1pij_I1ia_lessmem_quadratic(CCTaskParams params);

    /// mp2
    void MP2();

    /// mp4(sdq)
    void MP4_SDQ();

    /// components of mp3 and mp4 energies
    double emp3_os,emp3_ss,emp3,emp4_sd_os,emp4_sd_ss,emp4_sd;
    double emp4_q_os,emp4_q_ss,emp4_q;

    /// Update t1
    void UpdateT1(long int iter);
    void UpdateT1_mp4(long int iter);

    /// Update t2
    void UpdateT2(long int iter);
    void UpdateT2_mp4(long int iter);

    /// evaluate energy
    double CheckEnergy();

    /// the n^6 cc/qci diagrams
    void I2ijkl(CCTaskParams params);
    void I2piajk(CCTaskParams params);
    void Vabcd1(CCTaskParams params);
    void Vabcd2(CCTaskParams params);
    void Vabcd(CCTaskParams params);
    void K(CCTaskParams params);
    void TwoJminusK(CCTaskParams params);

    /// DIIS functions
    void DIIS(double*c,long int nvec,long int n);
    void DIISOldVector(long int iter,int diis_iter,int replace_diis_iter);
    double DIISErrorVector(int diis_iter,int replace_diis_iter,int iter);
    void DIISNewAmplitudes(int diis_iter);
    long int maxdiis;
    double*diisvec;

    /// basic parameters
    long int ndoccact,ndocc,nvirt,nso,nmotemp,nmo,nfzc,nfzv,nvirt_no;

    /// available memory
    long int memory;

    /// maximum number of iterations
    long int maxiter;

    /// energy convergence 
    double e_conv;

    /// amplitude convergence
    double r_conv; 

    /// orbital energies
    double *eps;

    /// reference energy
    double escf;

    /// workspace buffers.
    double*integrals,*tempt,*tempv;

    /// t1 and t2 buffers
    double*tb,*t1;

    /// buffers for singles residual and a couple of tiny intermediates
    double *w1,*I1,*I1p;

    /// define tiling
    void DefineTilingCPU();
    long int ovtilesize,lastovtile,lastov2tile,ov2tilesize;
    long int tilesize,lasttile,maxelem;
    long int ntiles,novtiles,nov2tiles;
};
}};

#endif
