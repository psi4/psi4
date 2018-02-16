/*
 *@BEGIN LICENSE
 *
 * GPU-accelerated density-fitted coupled-cluster, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef PLUGIN_CCSD_H
#define PLUGIN_CCSD_H

#include<psi4/libpsi4util/PsiOutStream.h>

#include<pthread.h>
#include<sys/types.h>
#include"psi4/psi4-dec.h"
#include<psi4/libiwl/iwl.h>
#include<psi4/libpsio/psio.hpp>
#include<psi4/libpsio/psio.h>
#include<psi4/libmints/wavefunction.h>
#include<psi4/psifiles.h>

#include<sys/times.h>
#include"gpuhelper.h"

// cuda libraries
#include<cuda.h>
#include<cublas.h>
#include<cuda_runtime.h>

// header for DFCoupledCluster
#include"dfcc_fnocc.h"

typedef struct {
        int id;
} parm;

namespace psi{ namespace fnocc{

// GPU DFCC class
class GPUDFCoupledCluster : public DFCoupledCluster{

  public:
    GPUDFCoupledCluster(std::shared_ptr<psi::Wavefunction>wfn,Options&options);
    ~GPUDFCoupledCluster();

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
    double compute_energy();
    void common_init();
    void pthreadCCResidual(int id);
  protected:

    /// t1-transformed 3-index integrals
    void T1Integrals();
    /// t1-transformed 3-index integrals (from SCF) for Fock build
    void T1Fock();

    /// cc diagrams:
    virtual void CCResidual();
    //virtual void saveCCResidual();

    /// extra storage for gpu (ac|bd) function: 1/2o(o+1)v(v+1)
    double * tempr, ** tempr2;
    void useVabcd1();
    void workinguseVabcd1();

    /// triples
    PsiReturnType triples();

    /// (ac|bd) diagram - this will be the first target to gpuify
    void brokenVabcd1();
    void tiledVabcd1();
    void slowVabcd1();
    void saveVabcd1();
    void workingVabcd1();
    void cpuVabcd1();
    void nearlythereVabcd1();
    virtual void Vabcd1();
    void FinishVabcd1();

    void UpdateT2();

    /// initialize cuda
    void CudaInit();

    /// initialize cuda
    void CudaFinalize();

    /// define tiling based on available gpu memory
    void DefineTiling();

    /// allocate memory for gpu
    void AllocateGPUMemory();

    /// allocate memory for cpu
    virtual void AllocateMemory();

    pthread_t base_thread;

    /// gpu buffers
    double ** gpubuffer;
    //double * gput2;
    //double * gpuv;
    //double * gpur2;
    //double * gputemp;
    //double * gput1;
    //double * gpur1;


    /// GPU-specific variables ... some are obsolete
    int num_gpus;
    double left, wasted;
    long int ovtilesize, novtiles, lastovtile, lastov2tile, ov2tilesize, nov2tiles;
    long int tilesize, ntiles, lasttile;
    long int ncputhreads,ngputhreads, nblocks, num;
    bool gpudone;
    bool cpudone;
    long int last_a;


    std::shared_ptr<GPUHelper> helper_;
    //__global__ void GPUKernel_Iqdb(int a,int v,int nQ,double * in,double * out);

};

    /// GPU kernels
    //__global__ void GPUKernel_Iqdb(int a,int v,int nQ,double * in,double * out);
    //__global__ void GPUKernel_Vm(int a,int v,double * in,double * out);
    //__device__ int  GPUKernel_Position(int i,int j);

}}


#endif
