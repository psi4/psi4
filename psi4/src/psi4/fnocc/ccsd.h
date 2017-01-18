/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef CCSD_H
#define CCSD_H

#include "psi4/psi4-dec.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/psifiles.h"

long int Position(long int i,long int j);

namespace psi{ namespace fnocc{

class CoupledCluster: public Wavefunction{
  public:

    CoupledCluster(std::shared_ptr<Wavefunction> reference_wavefunction,Options &options);
    ~CoupledCluster();

    double compute_energy();

  protected:

    int iter;
    int brueckner_iter;

    void common_init();
    void finalize();

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
        void(psi::fnocc::CoupledCluster::*func)(CCTaskParams);
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
    double eccsd, eccsd_os, eccsd_ss, eccsd_os_fac, eccsd_ss_fac;

    /// cc or qci (t)
    PsiReturnType triples();
    PsiReturnType lowmemory_triples();
    double et;

    /// mp4 triples
    void mp4_triples();
    double emp4_t;

    void WriteBanner();
    void WriteOptions();

    /// allocate memory
    virtual void AllocateMemory();

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
    void DIIS(double*c,long int nvec,long int n,int replace_diis_iter);
    void DIISOldVector(long int iter,int diis_iter,int replace_diis_iter);
    double DIISErrorVector(int diis_iter,int replace_diis_iter,int iter);
    void DIISNewAmplitudes(int diis_iter,int&replace_diis_iter);
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

// DF CC class
class DFCoupledCluster : public CoupledCluster{

  public:
    DFCoupledCluster(SharedWavefunction ref_wfn, Options&options);
    ~DFCoupledCluster();

    double compute_energy();

  protected:
    void finalize();

    /// CCSD iterations
    PsiReturnType CCSDIterations();

    void WriteBanner();

    /// allocate memory
    virtual void AllocateMemory();

    /// update t1 amplitudes
    void UpdateT1();

    /// update t2 amplitudes
    virtual void UpdateT2();

    /// v^4 CC diagram
    virtual void Vabcd1();

    /// workspace buffers.
    double*Abij,*Sbij;

    /// check energy
    virtual double CheckEnergy();

    ///  3-index integrals for density fitting.
    bool ischolesky_;
    long int nQ;
    long int nQ_scf;
    double*Qov,*Qvv,*Qoo;
    void  ThreeIndexIntegrals();

    /// more 3-index stuff for t1-transformed integrals
    double * Ca_L, * Ca_R, **Ca;
    double *Fij, *Fab, *Fia, *Fai;
    SharedMatrix H;

    /// generate t1-transformed 3-index integrals
    virtual void T1Integrals();
    /// generate t1-transformed Fock matrix
    virtual void T1Fock();

    /// evaluate cc diagrams
    virtual void CCResidual();

    /// SCS-MP2 function and variables
    virtual void SCS_MP2();

    /// SCS-CCSD function and variables
    virtual void SCS_CCSD();
};

// coupled pair class
class CoupledPair : public CoupledCluster{

  public:
    CoupledPair(std::shared_ptr<psi::Wavefunction>wfn,Options&options);
    ~CoupledPair();

    double compute_energy();

  protected:

    /// coupled pair iterations
    PsiReturnType CEPAIterations();

    /// free memory
    void finalize();

    /// pair energies
    void PairEnergy();
    double * pair_energy;

    /// what kind of coupled pair method?
    char * cepa_type;
    int cepa_level;

    /// check energy
    double CheckEnergy();

    /// check energy for coupled pair methods that have an energy functional
    double VariationalEnergy();
    double evar;

    /// update t1 amplitudes
    void UpdateT1();

    /// update t2 amplitudes
    void UpdateT2();

    /// scs functions
    void SCS_CEPA();

    /// compute opdm - only valid for cisd, acpf, aqcc, and cepa(0)
    void OPDM();

    /// banner
    void WriteBanner();
};

}}

#endif
