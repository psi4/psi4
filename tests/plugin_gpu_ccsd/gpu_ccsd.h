#ifndef GPU_CCSD_H
#define GPU_CCSD_H

#define PSIF_IJAK  250
#define PSIF_ABCI  251
#define PSIF_ABCI2 252
#define PSIF_ABCI3 253
#define PSIF_ABCI4 254
#define PSIF_ABCI5 255
#define PSIF_ABCD1 256
#define PSIF_ABCD2 257
#define PSIF_AKJC2 259
#define PSIF_KLCD  260
#define PSIF_IJKL  261
#define PSIF_OVEC 262
#define PSIF_EVEC 263
#define PSIF_AKJC 264
#define PSIF_KLCD2 265
#define PSIF_T2 266

#include"psi4-dec.h"
namespace boost {
template<class T> class shared_ptr;
}

long Position(int i,int j);

// gpu ccsd class
namespace psi{

class GPUCoupledCluster{
  public:
    GPUCoupledCluster();
    ~GPUCoupledCluster();
    void Initialize(Options &options);
    void WriteBanner(Options &options);
    void WriteIntegrals(double*tei);
    void ReadIntegrals();
    void AllocateMemory();

    PsiReturnType CCSDIterations(Options&options);

    bool t2_on_disk;

    /**
      * easy diagrams to be evaluated on the cpu:
      */
    void CPU_t1_vmeai();
    void CPU_t1_vmeni();
    void CPU_t1_vmaef();
    void CPU_I1ab();
    void CPU_I2p_abci_refactored();
    void CPU_I1pij_I1ia_lessmem();

    /**
      * update amplitudes
      */
    void UpdateT1(int iter);
    void UpdateT2(int iter);
    void TCEPA();

    double CheckEnergy();

    /**
      * DIIS stuff
      */
    void DIIS(double*c,long int nvec,long int n);
    void DIISOldVector(long int iter,int diis_iter,int replace_diis_iter);
    double DIISErrorVector(int diis_iter,int replace_diis_iter,int iter);
    void DIISNewAmplitudes(int diis_iter);

    /**
      *  SCS-MP2 function and variables
      */
    void SCS_MP2();
    double emp2_os,emp2_ss,emp2_os_fac,emp2_ss_fac;

    /**
      *  SCS-CCSD function and variables
      */
    void SCS_CCSD();
    double eccsd_os,eccsd_ss,eccsd_os_fac,eccsd_ss_fac;

    int ndoccact,ndocc,nvirt,nso,nmotemp,nmo,nirreps,*docc,nfzc,nfzv,*fzc,*fzv,*orbs,*sorbs;
    double *oei,*tei,*Fock,*eps;
    boost::shared_ptr<Vector> eps_test;
    double escf,enuc,efzc,emp2,eccsd,et;
    int maxiter;
    double conv;

    // integrals:
    double*E2ijak,*E2abci,*E2ijkl,*E2akjc_2,*E2klcd_1,*Symabcd1,*Symabcd2;
    double*E2ijak3,*E2ijak2,*E2ijakCopy;

    // other buffers:
    double*tempu,*tb,*wb,*w1,*I1,*I1p,*t1;
    double*tempt,*tempv;

    // diis:
    int maxdiis;
    double*diisvec;

    // memory (in bytes)
    long int memory;

    // GPU specific variables and functions:
    // for tiling:
    double left,wasted;
    int ovtilesize,novtiles,lastovtile,lastov2tile,ov2tilesize,nov2tiles;
    int tilesize,ntiles,lasttile;
    // gpu arrays:
    double*gput2,*gput1,*gpuw,*gpuv,*gpuw1,*gputempw;
    // threads and blocks
    int nthreads,nblocks,num;

    void CudaInit();
    void CudaFinalize();
    void DefineTiling();
    void AllocateGPUMemory();

};
};

#endif
