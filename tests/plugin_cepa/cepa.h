#ifndef CEPA_H
#define CEPA_H

// output files
#define PSIF_IJAK  251
#define PSIF_IJAK2 252
#define PSIF_ABCI3 253
#define PSIF_ABCI5 254
#define PSIF_ABCD1 255
#define PSIF_ABCD2 256
#define PSIF_AKJC2 257
#define PSIF_KLCD  258
#define PSIF_IJKL  259
#define PSIF_OVEC  260
#define PSIF_EVEC  261
#define PSIF_R2    262
#define PSIF_TEMP  263
#define PSIF_T2    264

namespace boost {
template<class T> class shared_ptr;
}

long int Position(long int i,long int j);

// gpu cepa class
namespace psi{

class CoupledPair{
  public:
    CoupledPair();
    ~CoupledPair();

    /**
      * Flag to indicate if t2 is stored in core memory or 
      * needs to be read from disk.  Default false.
      */
    bool t2_on_disk;

    /**
      * Define CEPA Tasks.  most diagrams are designated as
      * independent tasks.  some will be tiled out as 
      * separate tasks either so we can do them with limited
      * memory or for distribution among many processors 
     */
    void DefineTasks();
    long int ncepatasks;
    /**
      * CEPA Tasks parameters.  for terms that are tiled, the
      * let the task know which tile it should be working on.
      * NOTE: the tiling is really meant for the parallel
      * code, not this one.
     */
    struct CepaTaskParams{
        int mtile,ntile,ktile;
    };
    CepaTaskParams*CepaParams,*CCSubParams1,*CCSubParams2;
    /**
      * CEPA Tasks. this struct contains a pointer to
      * the task function, and some other stupid info
      * like how many flops the task will take.
     */
    struct CepaTask{
        void(psi::CoupledPair::*func)(CepaTaskParams);
        double flopcount;
    };
    CepaTask*CepaTasklist,*CCSubTasklist1,*CCSubTasklist2;

    /**
      * this function solves the CEPA equations (requires a minimum of 3o^2v^2 memory)
     */
    PsiReturnType CEPAIterations(Options&options);
  
    void WriteBanner(Options&options);

    /**
      * grab parameters, transform/sort integrals
     */
    void Initialize(Options &options);

    /**
      * allocate memory, define tiling of gigantic diragrams
     */
    void AllocateMemory(Options&options);

    /**
      * some CC diagrams.  these were "easy" ones that i used to
      * do on the cpu while the gpu worked.  maybe we'll do that
      * again eventually.
      */
    void CPU_t1_vmeai(CepaTaskParams params);
    void CPU_t1_vmeni(CepaTaskParams params);
    void CPU_t1_vmaef(CepaTaskParams params);

    /**
      * this diagram required ov^3 storage for an
      * intermediate in the Piecuch CPC.  this formulation
      * only requires o^3v storage...at the expense of
      * 4 extra O(N^5) terms.
      */
    void CPU_I2p_abci_refactored_term1(CepaTaskParams params);

    /**
      * Update t1
      */
    void UpdateT1(long int iter);

    /**
      * Update t2
      */
    void UpdateT2(long int iter);
    /**
      * Get the energy for that iteration. If there is a diis extrapolation,
      * the energy is evaluated after that step.
      */
    double CheckEnergy();

    /**
      * what level of cepa? 0,1,2,3.  default 0
      */
    int cepa_level;
    char*cepa_type;

    /**
      * construct an array of pair energies
      */
    void PairEnergy();
    double*pair_energy;

    /**
      * the N^6 CEPA diagrams.
      */
    void I2ijkl(CepaTaskParams params);
    void I2piajk(CepaTaskParams params);
    void Vabcd1(CepaTaskParams params);
    void Vabcd2(CepaTaskParams params);
    void I2iabj(CepaTaskParams params);
    void I2iajb(CepaTaskParams params);

    /**
      * DIIS stuff
      */
    void DIIS(double*c,long int nvec,long int n);
    void DIISOldVector(long int iter,int diis_iter,int replace_diis_iter);
    double DIISErrorVector(int diis_iter,int replace_diis_iter,int iter);
    void DIISNewAmplitudes(int diis_iter);
    long int maxdiis;
    double*diisvec;

    /**
      * basic parameters
      */
    long int ndoccact,ndocc,nvirt,nso,nmotemp,nmo,nirreps,memory;
    int maxiter,*docc,nfzc,nfzv,*fzc,*fzv,*orbs,*sorbs,nvirt_no;
    double conv,*oei,*tei,*Fock,*eps,scale_t;
    boost::shared_ptr<Vector> eps_test;
    double escf,enuc,efzc,emp2,ecepa,et;

    /**
      * workspace buffers.
      */
    double*integrals,*tempt,*tempv;

    /**
      * t1 and t2
      */
    double*tb,*t1;

    /**
      * the singles residual and a couple of tiny intermediates
      */
    double *w1,*I1,*I1p;

    /**
      * define tiling.
      */
    void DefineTilingCPU();
    long int ovtilesize,lastovtile,lastov2tile,ov2tilesize;
    long int tilesize,lasttile,maxelem;
    long int ntiles,novtiles,nov2tiles;

    /**
      *  SCS-MP2 function and variables
      */
    void SCS_MP2();
    double emp2_os,emp2_ss,emp2_os_fac,emp2_ss_fac;

    /**
      *  SCS-CEPA function and variables
      */
    void SCS_CEPA();
    double ecepa_os,ecepa_ss,ecepa_os_fac,ecepa_ss_fac;

};
};

#endif
