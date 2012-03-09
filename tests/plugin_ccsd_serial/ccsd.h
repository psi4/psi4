#ifndef CCSD_H
#define CCSD_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// output files
#define PSIF_IJAK  250
#define PSIF_IJAK2 251
#define PSIF_ABCI  252
#define PSIF_ABCI2 253
#define PSIF_ABCI3 254
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

// psi headers
#include"psi4-dec.h"
#include <libplugin/plugin.h>
#include<boost/shared_ptr.hpp>
#include<liboptions/liboptions.h>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>
#include<libmints/matrix.h>
#include<libmints/vector.h>
#include<libchkpt/chkpt.h>
#include<libiwl/iwl.h>
#include <libpsio/psio.hpp>

// gpu helper class 
#include"gpuhelper.h"

namespace boost {
template<class T> class shared_ptr;
}

long int Position(long int i,long int j);

// gpu ccsd class
namespace psi{

  class GPUHelper;

  /**
    * get the address of an element somewhere in a file.
   */
  psio_address psio_get_address(psio_address start, long int shift);


class CoupledCluster{
  public:
    CoupledCluster();
    ~CoupledCluster();

    /**
      * Number of threads
      */
    int nthreads;

    /**
      * Flag to indicate if t2 is stored in core memory or 
      * needs to be read from disk.  Default false.
      */
    bool t2_on_disk;

    /**
      * Define CCSD Tasks.  most diagrams are designated as
      * independent tasks.  some will be tiled out as 
      * separate tasks either so we can do them with limited
      * memory or for distribution among many processors 
     */
    void DefineTasks();
    long int ncctasks;
    /**
      * CCSD Tasks parameters.  for terms that are tiled, the
      * let the task know which tile it should be working on.
      * NOTE: the tiling is really meant for the parallel
      * code, not this one.
     */
    struct CCTaskParams{
        int mtile,ntile,ktile;
    };
    CCTaskParams*CCParams,*CCSubParams1,*CCSubParams2;
    /**
      * CCSD Tasks. this struct contains a pointer to
      * the task function, and some other stupid info
      * like how many flops the task will take.
     */
    struct CCTask{
        void(psi::CoupledCluster::*func)(CCTaskParams);
        double flopcount;
    };
    CCTask*CCTasklist,*CCSubTasklist1,*CCSubTasklist2;

    /**
      * this function solves the CCSD equations (requires a minimum of 4o^2v^2 memory)
     */
    PsiReturnType CCSDIterations(Options&options);
  
    /**
      * delete temporary files
     */
    void cleanup();

    void WriteBanner();
    void RandomIntegralFiles();

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
    void CPU_t1_vmeai(CCTaskParams params);
    void CPU_t1_vmeni(CCTaskParams params);
    void CPU_t1_vmaef(CCTaskParams params);
    void CPU_I1ab(CCTaskParams params);
    void CPU_I1pij_I1ia_lessmem(CCTaskParams params);

    /**
      * this diagram required ov^3 storage for an
      * intermediate in the Piecuch CPC.  this formulation
      * only requires o^3v storage...at the expense of
      * 4 extra O(N^5) terms.  
      * 
      * 2/9/12: Now reduced to 3 total terms by dumping
      * one in with I2piajk and one in with I2iabj.  The
      * total cost of the original diagram was reduced 
      * from 3o^2v^3 tp 3o^2v^2.
      */
    void CPU_I2p_abci_refactored_term2(CCTaskParams params);

    /**
      * Update t1
      */
    void UpdateT1(long int iter);

    /**
      * Update t2
      */
    void UpdateT2(long int iter);

    /**
      * get a better guess at t2 than mp2
      */
    void TCEPA();

    /**
      * Get the energy for that iteration. If there is a diis extrapolation,
      * the energy is evaluated after that step.
      */
    double CheckEnergy();

    /**
      * the N^6 CC diagrams.
      */
    void I2ijkl(CCTaskParams params);
    void I2piajk(CCTaskParams params);
    void Vabcd1(CCTaskParams params);
    void Vabcd2(CCTaskParams params);
    void I2iabj(CCTaskParams params);
    void I2iajb(CCTaskParams params);

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
    double escf,enuc,efzc,emp2,eccsd,et;

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
     * helper class definied in gpuhelper.h
     */
    boost::shared_ptr<GPUHelper> helper_;

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

};
};

#endif
