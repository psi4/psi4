#ifndef CCSD_H
#define CCSD_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// output files
#define PSIF_IJAK  250
#define PSIF_IJAK2  266
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
#define PSIF_R2 264
#define PSIF_TEMP 265


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
      * master and slave functions
     */
    void worker();
    void master();

    /**
      * this function solves the CCSD equations (requires a minimum of 4o^2v^2 memory)
     */
    PsiReturnType CCSDIterations();
  
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
      * transform integrals
     */
    void Transformation();

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
      * this diagram was required ov^3 storage for an
      * intermediate in the Piecuch CPC.  this formulation
      * only requires o^3v storage...at the expense of
      * 4 extra O(N^5) terms.
      */
    void CPU_I2p_abci_refactored(CCTaskParams params);
    void CPU_I2p_abci_refactored_term1(CCTaskParams params);
    void CPU_I2p_abci_refactored_term2(CCTaskParams params);
    void CPU_I2p_abci_refactored_term3(CCTaskParams params);

    /**
      * Update t1
      */
    void UpdateT1(long int iter);

    /**
      * Update t2 - returns the energy for that iteration
      */
    double UpdateT2(long int iter);

    /**
      * some N^6 CC diagrams.  most of these can be
      * tiled accross multiple nodes. TODO: I2ijkl needs
      * to be parallelized still.
      */
    void I2ijkl(CCTaskParams params);
    void I2piajk(CCTaskParams params);
    void I2piajk_tiled(CCTaskParams params);
    void Vabcd1(CCTaskParams params);
    void Vabcd1_tiled(CCTaskParams params);
    void Vabcd2(CCTaskParams params);
    void Vabcd2_tiled(CCTaskParams params);

    /**
      * one of the really evil N^6 diagrams. forming the
      * intermediate costs you 2o^3v^3.  using it costs
      * o^3v^3.  can be tiled accross multiple nodes.
      */
    long int niabjtasks,niabjtiles,iabjtilesize,lastiabjtile;
    void Distribute_I2iabj(CCTaskParams params);
    void I2iabj_BuildIntermediate(CCTaskParams params);
    void I2iabj_BuildIntermediate1(CCTaskParams params);
    void I2iabj_BuildIntermediate2(CCTaskParams params);
    void I2iabj_UseIntermediate(CCTaskParams params);

    /**
      * the other really evil N^6 diagram. forming the
      * intermediate costs you o^3v^3.  using it costs
      * 2o^3v^3.  can be tiled accross multiple nodes.
      */
    long int niajbtasks,niajbtiles,iajbtilesize,lastiajbtile;
    void Distribute_I2iajb(CCTaskParams params);
    void I2iajb_BuildIntermediate(CCTaskParams params);
    void I2iajb_UseIntermediate1(CCTaskParams params);
    void I2iajb_UseIntermediate2(CCTaskParams params);

    /**
      * older versions of the above o^3v^3 ones.
      */
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
      * integral arrays.  i don't actually use any of these anymore, but
      * my code will probably scream at me if i remove them
      */
    double*E2ijak,*E2abci,*E2ijkl,*E2akjc_2,*E2klcd_1,*Symabcd1,*Symabcd2;
    double*E2ijak3,*E2ijak2,*E2ijakCopy,*integrals;

    /**
      * extra buffers.  i don't wb or tempu anymore, but
      * my code will probably scream at me if i remove them
      */
    double*tempu,*tb,*wb,*w1,*I1,*I1p,*t1;
    double*tempt,*tempv,energy;

    /**
      * define tiling.  right now, tilesizes give us gemms we
      * can do in memory.  TODO: this will control parallelization
      * as well.
      */
    void DefineTilingCPU();
    long int ovtilesize,lastovtile,lastov2tile,ov2tilesize;
    long int tilesize,lasttile,maxelem;
    long int ntiles,novtiles,nov2tiles;

    /**
     * helper class definied in gpuhelper.h
     */
    boost::shared_ptr<GPUHelper> helper_;

};
};

#endif
