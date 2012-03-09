#include"psi4-dec.h"
#include<libmints/wavefunction.h>
#include<libmints/vector.h>
#include<libpsio/psio.hpp>
#include<sys/times.h>

#include"blas.h"
#include"gpu_ccsd.h"
#include"gpu_kernels.h"

inline void Check_CUDA_Error(FILE*fp,const char *message){
  cudaError_t error = cudaGetLastError();
  if (error!=cudaSuccess) {
     fprintf(fp,"\n  ERROR: %s: %s\n\n", message, cudaGetErrorString(error) );
     fflush(fp);
     exit(-1);
  }
}

using namespace psi;
using namespace boost;

namespace psi{
  void ReadTEIs(double*tei,Options&options);
}

// position in a symmetric packed matrix
long Position(int i,int j){
  if (i<j){
    return ((j*(j+1))>>1)+i;
  }
  return ((i*(i+1))>>1)+j;
}

namespace psi{
GPUCoupledCluster::GPUCoupledCluster()
{}
GPUCoupledCluster::~GPUCoupledCluster()
{}

void GPUCoupledCluster::WriteBanner(Options &options){
  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                        CCSD                         *\n");
  fprintf(outfile, "        *           Coupled-Cluster Singles Doubles           *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                   Eugene DePrince                   *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);
}

/*================================================================
  
  Initialize:
  set essential variables (ndocc...).  Read 1- and 2-electron
  integrals into core.  Sort and write integrals.
  
================================================================*/
void GPUCoupledCluster::Initialize(Options &options){

  // grab the reference wave function and its parameters
  boost::shared_ptr<Wavefunction> ref = Process::environment.reference_wavefunction();

  if (ref.get() !=NULL){
     escf    = Process::environment.globals["SCF TOTAL ENERGY"];
     nirreps = ref->nirrep();
     sorbs   = ref->nsopi();
     orbs    = ref->nmopi();
     docc    = ref->doccpi();
     fzc     = ref->frzcpi();
     fzv     = ref->frzvpi();
  }
  if (nirreps>1){
     //throw PsiException("plugin_gpu_ccsd requires symmetry c1",__FILE__,__LINE__);
  }
  nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
  int full=0;
  for (int h=0; h<nirreps; h++){
      nfzc   += fzc[h];
      nfzv   += fzv[h];
      nso    += sorbs[h];
      full   += orbs[h];
      nmo    += orbs[h]-fzc[h]-fzv[h];
      ndocc  += docc[h];//-fzc[h];
  }
  ndoccact = ndocc - nfzc;
  nvirt  = nmo - ndoccact;

  if (nvirt<ndoccact){
     throw PsiException("plugin_gpu_ccsd requires more virtual orbitals than active doubly occupied orbitals",__FILE__,__LINE__);
  }

  // get paramters from input 
  conv    = options.get_double("R_CONVERGENCE");
  maxiter = options.get_int("MAXITER");
  maxdiis = options.get_int("DIIS_MAX_VECS");

  // memory is from process::environment, but can override that
  memory = Process::environment.get_memory();
  if (options["MEMORY"].has_changed()){
     memory  = options.get_int("MEMORY");
     memory *= (long int)1024*1024;
  }

  // SCS MP2 and CCSD
  emp2_os_fac = options.get_double("MP2_SCALE_OS");
  emp2_ss_fac = options.get_double("MP2_SCALE_SS");
  eccsd_os_fac = options.get_double("CC_SCALE_OS");
  eccsd_ss_fac = options.get_double("CC_SCALE_SS");

  //boost::shared_ptr<Matrix> Ca = ref->Ca();

  nmotemp = full;//Ca->colspi()[0];

  // orbital energies
  /*eps_test = ref->epsilon_a();
  int i;
  double*tmpeps = eps_test->pointer();
  eps = (double*)malloc(nmo*sizeof(double));
  for (i=0; i<nmo; i++) eps[i] = tmpeps[i+nfzc];
  eps_test.reset();*/

  // orbital energies
  eps = (double*)malloc(nmo*sizeof(double));
  int count=0;
  for (int h=0; h<nirreps; h++){
      eps_test = ref->epsilon_a();
      for (int norb = fzc[h]; norb<docc[h]; norb++){
          eps[count++] = eps_test->get(h,norb);
      }
  }
  for (int h=0; h<nirreps; h++){
      eps_test = ref->epsilon_a();
      for (int norb = docc[h]; norb<orbs[h]-fzv[h]; norb++){
          eps[count++] = eps_test->get(h,norb);
      }
  }
  eps_test.reset();


  // so->mo tei transformation (no, we're just reading from disk)
  struct tms total_tmstime;
  const long clk_tck = sysconf(_SC_CLK_TCK);

  double time_start,user_start,sys_start,time_stop,user_stop,sys_stop;

  // sort integrals and write them to disk
  times(&total_tmstime);
  time_start = time(NULL);
  user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;

  int ntri = nmotemp*(nmotemp+1)/2;
  ntri = ntri*(ntri+1)/2;
  tei = (double*)malloc(sizeof(double)*ntri);
  memset((void*)tei,'\0',ntri*sizeof(double));

  // read integrals 
  ReadTEIs(tei,options);

  // sort integrals and write them to disk
  WriteIntegrals(tei);

  times(&total_tmstime);
  time_stop = time(NULL);
  user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;

  fprintf(outfile,"  Time for integral sort:           %6.2lf s (user)\n",user_stop-user_start);
  fprintf(outfile,"                                    %6.2lf s (system)\n",sys_stop-sys_start);
  fprintf(outfile,"                                    %6d s (total)\n",(int)time_stop-(int)time_start);

  // free teis
  free(tei);

  // t2 is in core
  t2_on_disk = false;
}

/*===================================================================

  solve ccsd equations

===================================================================*/
PsiReturnType GPUCoupledCluster::CCSDIterations(Options&options){

  struct tms total_tmstime;
  const long clk_tck = sysconf(_SC_CLK_TCK);
  time_t iter_start,iter_stop,time_start,time_stop;
  double user_start,user_stop,sys_start,sys_stop;

  int count,j,replace_diis_iter,diis_iter,iter;
  int o = ndoccact;
  int v = nvirt;
  int oo1o2 = o*(o+1)/2;
  int vv1o2 = v*(v+1)/2;

  diis_iter=iter=0;
  replace_diis_iter=1;
  double nrm=1.0;
  double Eold=1.0e9;
  eccsd=0.0;

  // get device parameters, set nblocks,nthreads,dimgrid
  CudaInit();
  // TODO: once nvcc can be used with psi, add dimgrid to GPUCoupledCluster class
  dim3 dimgrid (nblocks,num);

  // define tiling for v^4 and ov^3 diagrams on the gpu
  DefineTiling();

  // allocate memory on the gpu
  AllocateGPUMemory();

  cudaStream_t stream;
  cudaEvent_t estart,estop;
  cudaEventCreate(&estart);
  cudaEventCreate(&estop);

  fprintf(outfile,"\n");
  fprintf(outfile,
    "  Begin singles and doubles coupled cluster iterations\n\n");
  fprintf(outfile,
    "   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");
  fflush(outfile);

  boost::shared_ptr<PSIO> psio(new PSIO());

  // start timing the iterations
  times(&total_tmstime);
  time_start = time(NULL);
  user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;

  //TCEPA();
  while(iter<maxiter && nrm>conv){
      iter_start = time(NULL);

      memset((void*)w1,'\0',o*v*sizeof(double));
      memset((void*)wb,'\0',o*o*v*v*sizeof(double));
      if (iter>0){


         stream = NULL;

         // copy amplitudes to the device
         cudaEventRecord(estart,stream);
             cudaMemcpyAsync(gput1,t1,sizeof(double)*o*v,cudaMemcpyHostToDevice,stream);
             cudaMemcpyAsync(gput2,tb,sizeof(double)*o*o*v*v,cudaMemcpyHostToDevice,stream);
         cudaEventRecord(estop,stream);

         while( cudaEventQuery(estop) == cudaErrorNotReady );

         cudaEventRecord(estart,stream);

            //==========================================================
            //
            // I(ij,kl)
            //
            //==========================================================
            cudaMemcpyAsync(gpuw,E2klcd_1,sizeof(double)*o*o*v*v,cudaMemcpyHostToDevice,stream);
            GPUt2Plust1_and_E2klcd3<<<dimgrid,nthreads>>>(o,v,gput2,gput1,gpuw,gpuv);

            // Build and use I2(ij,kl)
            cudaMemcpyAsync(gpuw,E2ijkl,sizeof(double)*o*o*o*o,cudaMemcpyHostToDevice,stream);
            cublasDgemm('n','n',o*o,o*o,v*v,1.0,gput2,o*o,gpuv,v*v,1.0,gpuw,o*o);

            cudaMemcpyAsync(gpuv,E2ijakCopy,sizeof(double)*o*o*o*v,cudaMemcpyHostToDevice,stream);
            cublasDgemm('n','n',o,o*o*o,v,2.0,gput1,o,gpuv,v,1.0,gpuw,o);

            cublasDgemm('n','n',o*o,v*v,o*o,0.5,gpuw,o*o,gput2,o*o,0.0,gputempw,o*o);

            AddPermutedOnGPU<<<dimgrid,nthreads>>>(o,v,gputempw,gpuw);

            //==========================================================
            //
            // I(ia,jk)
            //
            //==========================================================

            // Build and use I2(ia,jk)
            cudaMemcpyAsync(gputempw,E2ijak2,sizeof(double)*o*o*o*v,cudaMemcpyHostToDevice,stream);

            for (j=0; j<novtiles-1; j++){
                cudaMemcpyAsync(gpuv,E2abci+j*v*v*ovtilesize,sizeof(double)*v*v*ovtilesize,cudaMemcpyHostToDevice,stream);
                cublasDgemm('n','n',o*o,ovtilesize,v*v,1.0,gput2,o*o,gpuv,v*v,1.0,gputempw+j*o*o*ovtilesize,o*o);
            }
            cudaMemcpyAsync(gpuv,E2abci+j*v*v*ovtilesize,sizeof(double)*v*v*lastovtile,cudaMemcpyHostToDevice,stream);
            cublasDgemm('n','n',o*o,lastovtile,v*v,1.0,gput2,o*o,gpuv,v*v,1.0,gputempw+j*o*o*ovtilesize,o*o);

            cublasDgemm('n','n',o*o*v,v,o,-1.0,gputempw,o*o*v,gput1,o,0.0,gpuv,o*o*v);

            GPUFill_I2iajk_and_c2Sym1<<<dimgrid,nthreads>>>(o,v,gpuv,gpuw,gput2,gputempw);

            //==========================================================
            //
            // t(ij,ef) v(ef,ab) - scuseria/jansen version
            //
            //==========================================================
            for (j=0; j<ntiles-1; j++){
                cudaMemcpyAsync(gpuv,Symabcd1+j*tilesize*vv1o2,tilesize*vv1o2*sizeof(double),cudaMemcpyHostToDevice,stream);
                cublasDgemm('n','n',oo1o2,tilesize,vv1o2,0.5,gputempw,oo1o2,gpuv,vv1o2,0.0,gput2+j*tilesize*oo1o2,oo1o2);
            }
            j=ntiles-1;
            cudaMemcpyAsync(gpuv,Symabcd1+j*tilesize*vv1o2,lasttile*vv1o2*sizeof(double),cudaMemcpyHostToDevice,stream);
            cublasDgemm('n','n',oo1o2,lasttile,vv1o2,0.5,gputempw,o*(o+1)/2,gpuv,vv1o2,0.0,gput2+j*tilesize*oo1o2,oo1o2);
            GPUSymmAdd1<<<dimgrid,nthreads>>>(o,v,gput2,gpuw);

            cudaMemcpyAsync(gput2,tb,sizeof(double)*o*o*v*v,cudaMemcpyHostToDevice,stream);
            GPUc2Sym2_onefunction<<<dimgrid,nthreads>>>(o,v,gput2,gput1,gputempw);
            for (j=0; j<ntiles-1; j++){
                cudaMemcpyAsync(gpuv,Symabcd2+j*tilesize*vv1o2,tilesize*vv1o2*sizeof(double),cudaMemcpyHostToDevice,stream);
                cublasDgemm('n','n',oo1o2,tilesize,vv1o2,0.5,gputempw,oo1o2,gpuv,vv1o2,0.0,gput2+j*tilesize*oo1o2,oo1o2);
            }
            j=ntiles-1;
            cudaMemcpyAsync(gpuv,Symabcd2+j*tilesize*vv1o2,lasttile*v*(v+1)/2*sizeof(double),cudaMemcpyHostToDevice,stream);
            cublasDgemm('n','n',oo1o2,lasttile,vv1o2,0.5,gputempw,oo1o2,gpuv,vv1o2,0.0,gput2+j*tilesize*oo1o2,oo1o2);

            GPUSymmAdd2<<<dimgrid,nthreads>>>(o,v,gput2,gpuw);
            cudaMemcpyAsync(tempu,gpuw,sizeof(double)*o*o*v*v,cudaMemcpyDeviceToHost,stream);

         // end of gpu section
         cudaEventRecord(estop,stream);

         // evaluate a few diagrams on the cpu while gpu runs
         CPU_t1_vmeni();
         CPU_t1_vmaef();
         CPU_I2p_abci_refactored();

         count=0;
         while( cudaEventQuery(estop) == cudaErrorNotReady )count++;
         //if (count==0)fprintf(outfile,"          Warning: CPU time exceeds GPU vabcd stream time\n");
            
         F_DAXPY(o*o*v*v,1.0,tempu,1,wb,1);

         // reset cublas stream to null stream
         cublasSetKernelStream(NULL);

         // build I2(ia,bj), which contributes to t2
         cublasSetKernelStream(stream);
         cudaEventRecord(estart,stream);
             cudaMemcpyAsync(gpuw,E2klcd_1,sizeof(double)*o*o*v*v,cudaMemcpyHostToDevice,stream);

             cudaMemcpyAsync(gput2,tb,sizeof(double)*o*o*v*v,cudaMemcpyHostToDevice,stream);
             GPUt2Plus2t1<<<dimgrid,nthreads>>>(o,v,gpuv,gput2,gput1,gputempw);

             cudaMemcpyAsync(gput2,E2klcd_1,sizeof(double)*o*o*v*v,cudaMemcpyHostToDevice,stream);
             cublasDgemm('n','n',o*v,o*v,o*v,-0.5,gpuv,o*v,gpuw,o*v,1.0,gput2,o*v);

             GPUv2MinusHalfv2<<<dimgrid,nthreads>>>(o,v,gpuv,gpuw);
             cublasDgemm('n','n',o*v,o*v,o*v,1.0,gputempw,o*v,gpuv,o*v,1.0,gput2,o*v);

             for (j=0; j<nov2tiles-1; j++){
                 cudaMemcpyAsync(gpuv,E2abci+j*v*ov2tilesize,sizeof(double)*v*ov2tilesize,cudaMemcpyHostToDevice,stream);
                 cublasDgemm('n','n',o,ov2tilesize,v,1.0,gput1,o,gpuv,v,0.0,gpuw+j*o*ov2tilesize,o);
             }
             j=nov2tiles-1;
             cudaMemcpyAsync(gpuv,E2abci+j*v*ov2tilesize,sizeof(double)*v*lastov2tile,cudaMemcpyHostToDevice,stream);
             cublasDgemm('n','n',o,lastov2tile,v,1.0,gput1,o,gpuv,v,0.0,gpuw+j*o*ov2tilesize,o);

             GPUPermute_iabj_to_aijb<<<dimgrid,nthreads>>>(o,v,gpuw,gpuv);

             cudaMemcpyAsync(gpuw,E2ijakCopy,sizeof(double)*o*o*o*v,cudaMemcpyHostToDevice,stream);
             cublasDgemm('n','n',o*o*v,v,o,-1.0,gpuw,o*o*v,gput1,o,1.0,gpuv,o*o*v);

             GPUFill_I2iabj<<<dimgrid,nthreads>>>(o,v,gput2,gpuv,gpuw);

             // use I2(ia,bj)
             GPU2t2Minust2<<<dimgrid,nthreads>>>(o,v,gpuv,gputempw);
             cublasDgemm('n','t',o*v,o*v,o*v,1.0,gpuw,o*v,gpuv,o*v,0.0,gput2,o*v);

             GPUFill_t2_I2iajb<<<dimgrid,nthreads>>>(o,v,gput2,gpuw);
             cudaMemcpyAsync(tempu,gpuw,sizeof(double)*o*o*v*v,cudaMemcpyDeviceToHost,stream);
         cudaEventRecord(estop,stream);

         // read integrals and do one small diagram on the cpu
         CPU_I1ab();

         count=0;
         while( cudaEventQuery(estop) == cudaErrorNotReady )count++;
         //if (count==0) fprintf(outfile,"          Warning: CPU time exceeds GPU I(ia,bj) stream time\n");

         psio->open(PSIF_ABCI4,PSIO_OPEN_OLD);
         psio->read_entry(PSIF_ABCI4,"E2abci4",(char*)&E2abci[0],v*v*v*o*sizeof(double));
         psio->close(PSIF_ABCI4,1);

         cublasSetKernelStream(NULL);
         F_DAXPY(o*o*v*v,1.0,tempu,1,wb,1);

         // build I2(ia,jb), which contributes to t2
         cublasSetKernelStream(stream);
         cudaEventRecord(estart,stream);
             cudaMemcpyAsync(gpuw,E2klcd_1,sizeof(double)*o*o*v*v,cudaMemcpyHostToDevice,stream);
             GPUt2Plus2t1_and_E2klcd2<<<dimgrid,nthreads>>>(o,v,gpuv,gputempw,gput1,gpuw,gput2);

             cudaMemcpyAsync(gpuw,E2akjc_2,sizeof(double)*o*o*v*v,cudaMemcpyHostToDevice,stream);
             cublasDgemm('n','n',o*v,o*v,o*v,-0.5,gpuv,o*v,gput2,o*v,1.0,gpuw,o*v);

             for (j=0; j<nov2tiles-1; j++){
                 cudaMemcpyAsync(gpuv,E2abci+j*v*ov2tilesize,sizeof(double)*v*ov2tilesize,cudaMemcpyHostToDevice,stream);
                 cublasDgemm('n','n',o,ov2tilesize,v,1.0,gput1,o,gpuv,v,0.0,gput2+j*o*ov2tilesize,o);

             }
             j=nov2tiles-1;
             cudaMemcpyAsync(gpuv,E2abci+j*v*ov2tilesize,sizeof(double)*v*lastov2tile,cudaMemcpyHostToDevice,stream);
             cublasDgemm('n','n',o,lastov2tile,v,1.0,gput1,o,gpuv,v,0.0,gput2+j*o*ov2tilesize,o);
             GPUFillI2iajb1<<<dimgrid,nthreads>>>(o,v,gput2,gpuw);

             cudaMemcpyAsync(gpuv,E2ijak3,sizeof(double)*o*o*o*v,cudaMemcpyHostToDevice,stream);
             cublasDgemm('n','n',o*o*v,v,o,-1.0,gpuv,o*o*v,gput1,o,0.0,gput2,o*o*v);
             GPUFillI2iajb2<<<dimgrid,nthreads>>>(o,v,gput2,gpuw);

             // Use I2(ia,jb)
             cublasDgemm('n','t',o*v,o*v,o*v,-1.0,gpuw,o*v,gputempw,o*v,0.0,gput2,o*v);
             GPUFill_t2_I2iabj1<<<dimgrid,nthreads>>>(o,v,gput2,gpuv);

             // Use I2(ia,jb) again
             GPUPermute_tikbc<<<dimgrid,nthreads>>>(gputempw,gput2,o,v);

             cublasDgemm('n','t',o*v,o*v,o*v,-1.0,gpuw,o*v,gput2,o*v,0.0,gputempw,o*v);
             GPUFill_t2_I2iabj2<<<dimgrid,nthreads>>>(o,v,gputempw,gpuv);

             cudaMemcpyAsync(tempu,gpuv,sizeof(double)*o*o*v*v,cudaMemcpyDeviceToHost,stream);
         cudaEventRecord(estop,stream);

         // do some small diagrams on the cpu
         CPU_t1_vmeai();
         CPU_I1pij_I1ia_lessmem();

         count=0;
         while( cudaEventQuery(estop) == cudaErrorNotReady )count++;
         //if (count==0) fprintf(outfile,"          Warning: CPU I1pij+I1ia time exceeds GPU I(ia,jb) stream time\n");

         cublasSetKernelStream(NULL);

         F_DAXPY(o*o*v*v,1.0,tempu,1,wb,1);

         // TODO:  find a better place to do this.  refill E2abci
         psio->open(PSIF_ABCI,PSIO_OPEN_OLD);
         psio->read_entry(PSIF_ABCI,"E2abci",(char*)&E2abci[0],v*v*v*o*sizeof(double));
         psio->close(PSIF_ABCI,1);
      }

      // update the amplitudes
      Eold = eccsd;
      UpdateT1(iter);
      UpdateT2(iter);

      // add vector to list for diis
      DIISOldVector(iter,diis_iter,replace_diis_iter);

      // diis error vector and convergence check
      nrm = DIISErrorVector(diis_iter,replace_diis_iter,iter);

      // diis extrapolation
      if (diis_iter>1){
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,o*o*v*v+o*v);
         else                   DIIS(diisvec,maxdiis,o*o*v*v+o*v);
         DIISNewAmplitudes(diis_iter);
      }
      eccsd = CheckEnergy();

      if (diis_iter<=maxdiis) diis_iter++;
      else if (replace_diis_iter<maxdiis) replace_diis_iter++;
      else replace_diis_iter = 1;

      iter_stop = time(NULL);
      fprintf(outfile,"  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,eccsd,eccsd-Eold,nrm,(int)iter_stop-(int)iter_start);
      fflush(outfile);
      iter++;
      if (iter==1){
         emp2 = eccsd;
         SCS_MP2();
      }
  }
  times(&total_tmstime);
  time_stop = time(NULL);
  user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;
  psio.reset();

  if (iter==maxiter){
     throw PsiException("  CCSD iterations did not converge.",__FILE__,__LINE__);
  }

  SCS_CCSD();

  fprintf(outfile,"\n");
  fprintf(outfile,"  CCSD iterations converged!\n");
  fprintf(outfile,"\n");
  if (options.get_bool("SCS_MP2")){
     fprintf(outfile,"        OS SCS-MP2 correlation energy:  %20.12lf\n",emp2_os*emp2_os_fac);
     fprintf(outfile,"        SS SCS-MP2 correlation energy:  %20.12lf\n",emp2_ss*emp2_ss_fac);
     fprintf(outfile,"        SCS-MP2 correlation energy:     %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac);
     fprintf(outfile,"      * SCS-MP2 total energy:           %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac+escf);
     fprintf(outfile,"\n");
  }
  fprintf(outfile,"        OS MP2 correlation energy:      %20.12lf\n",emp2_os);
  fprintf(outfile,"        SS MP2 correlation energy:      %20.12lf\n",emp2_ss);
  fprintf(outfile,"        MP2 correlation energy:         %20.12lf\n",emp2);
  fprintf(outfile,"      * MP2 total energy:               %20.12lf\n",emp2+escf);
  fprintf(outfile,"\n");
  if (options.get_bool("SCS_CCSD")){
     fprintf(outfile,"        OS SCS-CCSD correlation energy: %20.12lf\n",eccsd_os*eccsd_os_fac);
     fprintf(outfile,"        SS SCS-CCSD correlation energy: %20.12lf\n",eccsd_ss*eccsd_ss_fac);
     fprintf(outfile,"        SCS-CCSD correlation energy:    %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
     fprintf(outfile,"      * SCS-CCSD total energy:          %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
     fprintf(outfile,"\n");
  }
  fprintf(outfile,"        OS CCSD correlation energy:     %20.12lf\n",eccsd_os);
  fprintf(outfile,"        SS CCSD correlation energy:     %20.12lf\n",eccsd_ss);
  fprintf(outfile,"        CCSD correlation energy:        %20.12lf\n",eccsd);
  fprintf(outfile,"      * CCSD total energy:              %20.12lf\n",eccsd+escf);
  fprintf(outfile,"\n");
  fprintf(outfile,"  Total time for CCSD iterations: %10.2lf s (user)\n",user_stop-user_start);
  fprintf(outfile,"                                  %10.2lf s (system)\n",sys_stop-sys_start);
  fprintf(outfile,"                                  %10d s (total)\n",(int)time_stop-(int)time_start);
  fprintf(outfile,"\n");
  fprintf(outfile,"  Time per iteration:             %10.2lf s (user)\n",(user_stop-user_start)/(iter-1));
  fprintf(outfile,"                                  %10.2lf s (system)\n",(sys_stop-sys_start)/(iter-1));
  fprintf(outfile,"                                  %10.2lf s (total)\n",((double)time_stop-(double)time_start)/(iter-1));

  fflush(stdout);
  fflush(outfile);

  // free some of the cpu memory before exiting in case we end up doing triples next
  cudaFreeHost(E2ijkl);
  cudaFreeHost(E2ijak2);
  cudaFreeHost(E2ijak3);
  cudaFreeHost(E2ijakCopy);
  cudaFreeHost(E2akjc_2);
  cudaFreeHost(Symabcd1);
  cudaFreeHost(Symabcd2);
  cudaFreeHost(tempu);
  cudaFreeHost(w1);
  cudaFreeHost(I1);
  cudaFreeHost(I1p);

  return Success;
}
void GPUCoupledCluster::SCS_CCSD(){
  long int v = nvirt;
  long int o = ndoccact;
  long int iajb,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  for (long int a=o; a<nmo; a++){
      for (long int b=o; b<nmo; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  osenergy += E2klcd_1[iajb]*(tb[ijab]+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                  ssenergy += E2klcd_1[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                  ssenergy += E2klcd_1[iajb]*(t1[(a-o)*o+i]*t1[(b-o)*o+j]-t1[(b-o)*o+i]*t1[(a-o)*o+j]);
                  ijab++;
              }
          }
      }
  }
  eccsd_os = osenergy;
  eccsd_ss = ssenergy;
}
void GPUCoupledCluster::SCS_MP2(){
  long int v = nvirt;
  long int o = ndoccact;
  long int iajb,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  for (long int a=o; a<nmo; a++){
      for (long int b=o; b<nmo; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){
                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  osenergy += E2klcd_1[iajb]*tb[ijab];
                  ssenergy += E2klcd_1[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                  ijab++;
              }
          }
      }
  }
  emp2_os = osenergy;
  emp2_ss = ssenergy;
}
double GPUCoupledCluster::CheckEnergy(){
  long int v = nvirt;
  long int o = ndoccact;
  long int iajb,jaib,ijab=0;
  double energy = 0.0;
  for (long int a=o; a<nmo; a++){
      for (long int b=o; b<nmo; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){
                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);
                  energy += (2.*E2klcd_1[iajb]-E2klcd_1[jaib])*(tb[ijab]+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                  ijab++;
              }
          }
      }
  }
  return energy;
}
/*===================================================================

  allocate memory on gpu

===================================================================*/
void GPUCoupledCluster::AllocateGPUMemory(){
  int o2 = ndoccact*ndoccact;
  int ov = ndoccact*nvirt;
  int v2 = nvirt*nvirt;
  cudaMalloc((void**)&gpuv,sizeof(double)*((left-wasted))); 
  Check_CUDA_Error(outfile,"gpuv");
  cudaMalloc((void**)&gput2,sizeof(double)*o2*v2);      
  Check_CUDA_Error(outfile,"gput2");
  cudaMalloc((void**)&gpuw,sizeof(double)*o2*v2);       
  Check_CUDA_Error(outfile,"gpuw");
  cudaMalloc((void**)&gputempw,sizeof(double)*o2*v2);   
  Check_CUDA_Error(outfile,"gputempw");
  cudaMalloc((void**)&gput1,sizeof(double)*ov);         
  Check_CUDA_Error(outfile,"gput1");
  cudaMalloc((void**)&gpuw1,sizeof(double)*nmo*nmo);      
  Check_CUDA_Error(outfile,"gpuw1");
}
/*===================================================================

  free gpu memory and mapped cpu memory

===================================================================*/
void GPUCoupledCluster::CudaFinalize(){
  cudaFree(gpuv);
  cudaFree(gpuw);
  cudaFree(gputempw);
  cudaFree(gput2);
  cudaFree(gput1);
  cudaFreeHost(E2ijak);
  cudaFreeHost(E2klcd_1);
  cudaFreeHost(E2abci);
  cudaFreeHost(tb);
  cudaFreeHost(t1);
  cudaThreadExit();
}
/*===================================================================

  determine tiling for vabcd and vabci diagrams

===================================================================*/
void GPUCoupledCluster::DefineTiling(){
  int i,v = nvirt;
  int o = ndoccact;
  int ov2 = o*v*v;
  int ov = o*v;
  int o2 = o*o;
  int h = v;
  wasted = 350*1024*1024/8.; // leave an extra 200 mb on there.
  ntiles = -999;

  // check whether blocking of the vabcd diagram is necessary
  if (left-wasted>v*(v+1)/2*v*(v+1)/2){
     tilesize = v*(v+1)/2;
     ntiles = 1;
  }
  else{
     for (i=2; i<=v*(v+1)/2; i++){
         if (left-wasted>(double)tilesize*v*(v+1)/2/i+1){
            tilesize = v*(v+1)/2/i;
            if (i*tilesize < v*(v+1)/2) tilesize++;
            ntiles = i;
            break;
         }
     }
     if (ntiles==-999){
        fprintf(outfile,"\n  error. Not enough device memory.\n\n");
        fflush(outfile);
        exit(0);
     }
  }
  lasttile = v*(v+1)/2 - (ntiles-1)*tilesize;
  if (tilesize<o2){
     fprintf(outfile,"\n  error. Not enough device memory. maybe.\n\n");
     fflush(outfile);
  }
  if (ntiles>1){
     fprintf(outfile,
       "  v(ab,cd) diagram will be evaluated in %5i blocks.\n\n",ntiles); 
     fflush(outfile);
  }
  nov2tiles=1;
  ov2tilesize=ov2/1;
  if (nov2tiles*ov2tilesize<ov2) ov2tilesize++;
  while(h*ov2tilesize>left-wasted){
     nov2tiles++;
     ov2tilesize = ov2/nov2tiles;
     if (nov2tiles*ov2tilesize<ov2) ov2tilesize++;
  }
  if (nov2tiles>1){
     fprintf(outfile,
       "  v(ab,ci) terms will be evaluated in   %5i blocks.\n\n",nov2tiles); 
     fflush(outfile);
  }
  lastov2tile = ov2 - (nov2tiles-1)*ov2tilesize;

  novtiles=1;
  ovtilesize=ov/1;
  if (novtiles*ovtilesize<ov) ovtilesize++;
  while(h*h*ovtilesize>left-wasted){
     novtiles++;
     novtiles++;
     ovtilesize = ov/novtiles;
     if (novtiles*ovtilesize<ov) ovtilesize++;
  }
  lastovtile = ov - (novtiles-1)*ovtilesize;
}
/*===================================================================

  initialize cublas and get device properties

===================================================================*/
void GPUCoupledCluster::CudaInit(){
  cublasInit();
  struct cudaDeviceProp cudaProp;
  int gpu_id;
  cudaGetDevice(&gpu_id);
  cudaGetDeviceProperties( &cudaProp,gpu_id );
  fprintf(outfile,
    "\n  _________________________________________________________\n");
  fprintf(outfile,"  CUDA device properties:\n");
  fprintf(outfile,"  name:                 %20s\n",cudaProp.name);
  fprintf(outfile,"  major version:        %20d\n",cudaProp.major);
  fprintf(outfile,"  minor version:        %20d\n",cudaProp.minor);
  fprintf(outfile,"  canMapHostMemory:     %20d\n",cudaProp.canMapHostMemory);
  fprintf(outfile,"  totalGlobalMem:       %20lu mb\n",
    cudaProp.totalGlobalMem/(1024*1024));
  fprintf(outfile,"  sharedMemPerBlock:    %20lu\n",cudaProp.sharedMemPerBlock);
  fprintf(outfile,"  clockRate:            %20.3f ghz\n",
    cudaProp.clockRate/1.0e6);
  fprintf(outfile,"  regsPerBlock:         %20d\n",cudaProp.regsPerBlock);
  fprintf(outfile,"  warpSize:             %20d\n",cudaProp.warpSize);
  fprintf(outfile,"  maxThreadsPerBlock:   %20d\n",cudaProp.maxThreadsPerBlock);
  fprintf(outfile,
    "  _________________________________________________________\n\n");
  fflush(outfile);
  // device memory left after some arrays
  int o = ndoccact;
  int v = nvirt;
  left = cudaProp.totalGlobalMem/8. - 3*o*o*v*v - o*v-nmo*nmo;

  nthreads=NUMTHREADS;
  num=1;
  if ((o*o*v*v)%nthreads==0)
     nblocks = (o*o*v*v)/nthreads;
  else
     nblocks = (o*o*v*v+nthreads-(o*o*v*v)%nthreads)/nthreads;
  if (nblocks>MAXBLOCKS){
     num = nblocks/MAXBLOCKS+1;
     nblocks = nblocks/num + 1;
  }
}
/*===================================================================

  read integrals from disk

===================================================================*/
void GPUCoupledCluster::ReadIntegrals(){
  int v = nvirt;
  int o = ndoccact;
  int k;

  boost::shared_ptr<PSIO> psio(new PSIO());

  psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJAK,"E2ijak",(char*)&E2ijak[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_IJAK,1);

  psio->open(PSIF_ABCI,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_ABCI,"E2abci",(char*)&E2abci[0],o*v*v*v*sizeof(double));
  psio->close(PSIF_ABCI,1);

  psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_ABCD1,"E2abcd1",(char*)&Symabcd1[0],v*(v+1)/2*v*(v+1)/2*sizeof(double));
  psio->close(PSIF_ABCD1,1);

  psio->open(PSIF_ABCD2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_ABCD2,"E2abcd2",(char*)&Symabcd2[0],v*(v+1)/2*v*(v+1)/2*sizeof(double));
  psio->close(PSIF_ABCD2,1);

  psio->open(PSIF_IJKL,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJKL,"E2ijkl",(char*)&E2ijkl[0],o*o*o*o*sizeof(double));
  psio->close(PSIF_IJKL,1);

  psio->open(PSIF_AKJC2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_AKJC2,"E2akjc2",(char*)&E2akjc_2[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_AKJC2,1);

  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&E2klcd_1[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);

  int i,j,a,id=0;
  for (i=0; i<o; i++){
  for (a=0; a<v; a++){
  for (j=0; j<o; j++){
  for (k=0; k<o; k++){
      E2ijakCopy[id] = E2ijak[id];
      E2ijak2[id++] = E2ijak[i*o*o*v+k*o*v+j*v+a];
  }}}}
  id=0;
  for (j=0; j<o; j++){
  for (i=0; i<o; i++){
  for (k=0; k<o; k++){
  for (a=0; a<v; a++){
      E2ijak3[id++] = E2ijak[i*o*o*v+j*o*v+k*v+a];
  }}}}

}
/*===================================================================

  sort and write integrals to disk

===================================================================*/
void GPUCoupledCluster::WriteIntegrals(double*tei){
  int i,j,k,l,a,b,c,d,e,m,f;
  double dum;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  double*tmpei;
  long int dim = nvirt*(nvirt+1)/2;
  dim = dim*dim;
  tmpei = (double*)malloc(dim*sizeof(double));
  int count=0;

  // E<ij|ak>
  psio->open(PSIF_IJAK,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  for (j=nfzc; j<ndocc; j++){
      for (i=nfzc; i<ndocc; i++){
          for (k=nfzc; k<ndocc; k++){
              for (a=ndocc; a<nmotemp; a++){
                  dum = tei[Position(Position(i-nfzc,a-nfzc),Position(j-nfzc,k-nfzc))];
                  tmpei[count++] = dum;
                  if (count==dim){
                     psio->write(PSIF_IJAK,"E2ijak",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
                     count=0;
                  }
              }
          }
      }
  }
  if (count>0)
     psio->write(PSIF_IJAK,"E2ijak",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_IJAK,1);

  // E<ab|ci>
  psio->open(PSIF_ABCI,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (i=nfzc; i<ndocc; i++){
      for (a=ndocc; a<nmotemp; a++){
          for (b=ndocc; b<nmotemp; b++){
              for (c=ndocc; c<nmotemp; c++){
                  dum = tei[Position(Position(a-nfzc,c-nfzc),Position(b-nfzc,i-nfzc))];
                  tmpei[count++] = dum;
                  if (count==dim){
                     psio->write(PSIF_ABCI,"E2abci",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
                     count=0;
                  }
              }
          }
      }
  }
  if (count>0)
     psio->write(PSIF_ABCI,"E2abci",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_ABCI,1);

  // E<ab|cd>
  psio->open(PSIF_ABCD1,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (a=ndocc; a<nmotemp; a++){
  for (b=ndocc; b<=a; b++){
  for (c=ndocc; c<nmotemp; c++){
  for (d=ndocc; d<=c; d++){
      if (c==d) dum = tei[Position(Position(a-nfzc,c-nfzc),Position(b-nfzc,d-nfzc))];
      else      dum = tei[Position(Position(a-nfzc,c-nfzc),Position(b-nfzc,d-nfzc))]
                    + tei[Position(Position(a-nfzc,d-nfzc),Position(b-nfzc,c-nfzc))];
      tmpei[count++] = dum;
      if (count==dim){
         psio->write(PSIF_ABCD1,"E2abcd1",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
         count=0;
      }

  }}}}
  if (count>0)
     psio->write(PSIF_ABCD1,"E2abcd1",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_ABCD1,1);


  psio->open(PSIF_ABCD2,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (a=ndocc; a<nmotemp; a++){
  for (b=ndocc; b<=a; b++){
  for (c=ndocc; c<nmotemp; c++){
  for (d=ndocc; d<=c; d++){
      dum = tei[Position(Position(a-nfzc,c-nfzc),Position(b-nfzc,d-nfzc))]
          - tei[Position(Position(a-nfzc,d-nfzc),Position(b-nfzc,c-nfzc))];
      tmpei[count++] = dum;
      if (count==dim){
         psio->write(PSIF_ABCD2,"E2abcd2",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
         count=0;
      }
  }}}}
  if (count>0)
     psio->write(PSIF_ABCD2,"E2abcd2",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_ABCD2,1);


  // E<ij|kl>
  psio->open(PSIF_IJKL,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (k=nfzc; k<ndocc; k++){
      for (l=nfzc; l<ndocc; l++){
          for (i=nfzc; i<ndocc; i++){
              for (j=nfzc; j<ndocc; j++){
                  dum = tei[Position(Position(i-nfzc,k-nfzc),Position(j-nfzc,l-nfzc))];
                  tmpei[count++] = dum;
                  if (count==dim){
                     psio->write(PSIF_IJKL,"E2ijkl",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
                     count=0;
                  }
              }
          }
      }
  }
  if (count>0)
     psio->write(PSIF_IJKL,"E2ijkl",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_IJKL,1);

  // E<ak|jc>
  psio->open(PSIF_AKJC2,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (k=nfzc; k<ndocc; k++){
      for (c=ndocc; c<nmotemp; c++){
          for (j=nfzc; j<ndocc; j++){
              for (a=ndocc; a<nmotemp; a++){
                  dum = tei[Position(Position(a-nfzc,c-nfzc),Position(k-nfzc,j-nfzc))];
                  tmpei[count++] = dum;
                  if (count==dim){
                     psio->write(PSIF_AKJC2,"E2akjc2",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
                     count=0;
                  }
              }
          }
      }
  }
  if (count>0)
     psio->write(PSIF_AKJC2,"E2akjc2",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_AKJC2,1);

  // E<kl|cd>
  psio->open(PSIF_KLCD,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (k=nfzc; k<ndocc; k++){
      for (c=ndocc; c<nmotemp; c++){
          for (l=nfzc; l<ndocc; l++){
              for (d=ndocc; d<nmotemp; d++){
                  dum = tei[Position(Position(k-nfzc,c-nfzc),Position(l-nfzc,d-nfzc))];
                  tmpei[count++] = dum;
                  if (count==dim){
                     psio->write(PSIF_KLCD,"E2klcd",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
                     count=0;
                  }
              }
          }
      }
  }
  if (count>0)
     psio->write(PSIF_KLCD,"E2klcd",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_KLCD,1);

  // these won't get deleted:
  /*psio->open(PSIF_ABCI5,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (a=ndocc; a<nmotemp; a++){
  for (b=ndocc; b<nmotemp; b++){
  for (i=nfzc; i<ndocc; i++){
  for (c=ndocc; c<nmotemp; c++){
      dum = tei[Position(Position(a-nfzc,c-nfzc),Position(b-nfzc,i-nfzc))];
      tmpei[count++] = dum;
      if (count==ndoccact*ndoccact*nvirt*nvirt){
         psio->write(PSIF_ABCI5,"E2abci5",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
         count=0;
      }
      //psio->write(PSIF_ABCI5,"E2abci5",(char*)&dum,sizeof(double),addr,&addr);
  }}}}
  if (count>0)
     psio->write(PSIF_ABCI5,"E2abci5",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_ABCI5,1);*/

  psio->open(PSIF_ABCI2,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (a=ndocc; a<nmotemp; a++){
  for (b=ndocc; b<nmotemp; b++){
  for (e=ndocc; e<nmotemp; e++){
  for (m=nfzc; m<ndocc; m++){
      dum = 2.*tei[Position(Position(a-nfzc,b-nfzc),Position(e-nfzc,m-nfzc))]
          -    tei[Position(Position(a-nfzc,e-nfzc),Position(b-nfzc,m-nfzc))];
      tmpei[count++] = dum;
      if (count==dim){
         psio->write(PSIF_ABCI2,"E2abci2",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
         count=0;
      }
  }}}}
  if (count>0)
     psio->write(PSIF_ABCI2,"E2abci2",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_ABCI2,1);

  psio->open(PSIF_ABCI3,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (a=ndocc; a<nmotemp; a++){
  for (f=ndocc; f<nmotemp; f++){
  for (m=nfzc; m<ndocc; m++){
  for (e=ndocc; e<nmotemp; e++){
      dum = tei[Position(Position(a-nfzc,f-nfzc),Position(e-nfzc,m-nfzc))];
      tmpei[count++] = dum;
      if (count==dim){
         psio->write(PSIF_ABCI3,"E2abci3",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
         count=0;
      }
  }}}}
  if (count>0)
     psio->write(PSIF_ABCI3,"E2abci3",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_ABCI3,1);

  psio->open(PSIF_ABCI4,PSIO_OPEN_NEW);
  addr = PSIO_ZERO;
  count=0;
  for (i=nfzc; i<ndocc; i++){
  for (a=ndocc; a<nmotemp; a++){
  for (b=ndocc; b<nmotemp; b++){
  for (c=ndocc; c<nmotemp; c++){
      dum = tei[Position(Position(a-nfzc,b-nfzc),Position(c-nfzc,i-nfzc))];
      tmpei[count++] = dum;
      if (count==dim){
         psio->write(PSIF_ABCI4,"E2abci4",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
         count=0;
      }
  }}}}
  if (count>0)
     psio->write(PSIF_ABCI4,"E2abci4",(char*)&tmpei[0],count*sizeof(double),addr,&addr);
  psio->close(PSIF_ABCI4,1);

  psio.reset();
  free(tmpei);
}
/*===================================================================

  allocate cpu memory

===================================================================*/
void GPUCoupledCluster::AllocateMemory(){

  int o=ndoccact;
  int v=nvirt;
  int dim = v*(v+1)/2; 
  dim = dim*dim;

  // integrals:

  // o^4
  cudaMallocHost((void**)&E2ijkl,o*o*o*o*sizeof(double));      Check_CUDA_Error(outfile,"cpu E2ijkl");

  // o^3v
  cudaMallocHost((void**)&E2ijak,o*o*o*v*sizeof(double));      Check_CUDA_Error(outfile,"cpu E2ijak");
  cudaMallocHost((void**)&E2ijak2,o*o*o*v*sizeof(double));     Check_CUDA_Error(outfile,"cpu E2ijak2");
  cudaMallocHost((void**)&E2ijak3,o*o*o*v*sizeof(double));     Check_CUDA_Error(outfile,"cpu E2ijak3");
  cudaMallocHost((void**)&E2ijakCopy,o*o*o*v*sizeof(double));  Check_CUDA_Error(outfile,"cpu E2ijakCopy");

  // o^2v^2
  cudaMallocHost((void**)&E2akjc_2,o*o*v*v*sizeof(double));    Check_CUDA_Error(outfile,"cpu E2akjc_2");
  cudaMallocHost((void**)&E2klcd_1,o*o*v*v*sizeof(double));    Check_CUDA_Error(outfile,"cpu E2klcd_1");

  // ov^3
  cudaMallocHost((void**)&E2abci,o*v*v*v*sizeof(double));      Check_CUDA_Error(outfile,"cpu E2abci");

  // v(v+1)v(v+1)/4
  cudaMallocHost((void**)&Symabcd1,dim*sizeof(double));        Check_CUDA_Error(outfile,"cpu Symabcd1");
  cudaMallocHost((void**)&Symabcd2,dim*sizeof(double));        Check_CUDA_Error(outfile,"cpu Symabcd2");

  // extra buffers:
  cudaMallocHost((void**)&tempu,o*o*v*v*sizeof(double));  Check_CUDA_Error(outfile,"cpu tempu");
  cudaMallocHost((void**)&tb,o*o*v*v*sizeof(double));     Check_CUDA_Error(outfile,"cpu tb");
  cudaMallocHost((void**)&w1,o*v*sizeof(double));         Check_CUDA_Error(outfile,"cpu w1");
  cudaMallocHost((void**)&I1,v*v*sizeof(double));         Check_CUDA_Error(outfile,"cpu I1");
  cudaMallocHost((void**)&I1p,v*v*sizeof(double));        Check_CUDA_Error(outfile,"cpu I1p");
  cudaMallocHost((void**)&t1,o*v*sizeof(double));         Check_CUDA_Error(outfile,"cpu t1");

  // these don't need to be pinned and will never be touched by a gpu function
  wb    = (double*)malloc(sizeof(double)*o*o*v*v);
  tempt = (double*)malloc(sizeof(double)*(o*o*v*v+o*v));
  tempv = (double*)malloc(sizeof(double)*(o*o*v*v+o*v));

  memset((void*)E2ijak,'\0',o*o*o*v*sizeof(double));
  memset((void*)E2ijak2,'\0',o*o*o*v*sizeof(double));
  memset((void*)E2ijak3,'\0',o*o*o*v*sizeof(double));
  memset((void*)E2ijakCopy,'\0',o*o*o*v*sizeof(double));
  memset((void*)E2ijkl,'\0',o*o*o*o*sizeof(double));
  memset((void*)E2abci,'\0',o*v*v*v*sizeof(double));
  memset((void*)tempu,'\0',o*o*v*v*sizeof(double));
  memset((void*)wb,'\0',o*o*v*v*sizeof(double));
  memset((void*)tb,'\0',o*o*v*v*sizeof(double));
  memset((void*)E2klcd_1,'\0',o*o*v*v*sizeof(double));
  memset((void*)E2akjc_2,'\0',o*o*v*v*sizeof(double));
  memset((void*)w1,'\0',o*v*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));
  memset((void*)I1,'\0',v*v*sizeof(double));
  memset((void*)I1p,'\0',v*v*sizeof(double));

  // DIIS:
  diisvec    = (double*)malloc(sizeof(double)*(maxdiis+1));
  memset((void*)diisvec,'\0',(maxdiis+1)*sizeof(double));
}

void GPUCoupledCluster::CPU_t1_vmeai(){
  int o=ndoccact;
  int v=nvirt;
  for (int e=0; e<v; e++){
      for (int m=0; m<o; m++){
          for (int i=0; i<o; i++){
              F_DCOPY(v,E2akjc_2+m*o*v*v+e*o*v+i*v,1,tempv+e*o*o*v+m*o*v+i,o);
              F_DAXPY(v,-2.0,E2klcd_1+m*o*v*v+e*o*v+i*v,1,tempv+e*o*o*v+m*o*v+i,o);
          }
      }
  }
  F_DGEMV('n',o*v,o*v,-1.0,tempv,o*v,t1,1,1.0,w1,1);

}

void GPUCoupledCluster::CPU_t1_vmeni(){
  int o=ndoccact;
  int v=nvirt;
  for (int a=0; a<v; a++){
      for (int m=0; m<o; m++){
          for (int n=0; n<o; n++){
              F_DCOPY(v,tb+a*v*o*o+m*o+n,o*o,tempt+a*o*o*v+m*o*v+n*v,1);
              F_DAXPY(v,-2.0,tb+a*o*o+m*o+n,o*o*v,tempt+a*o*o*v+m*o*v+n*v,1);
          }
      }
  }
  F_DGEMM('t','n',o,v,o*o*v,1.0,E2ijak,o*o*v,tempt,o*o*v,1.0,w1,o);
}

void GPUCoupledCluster::CPU_t1_vmaef(){
  int m,e,i,f;
  int o=ndoccact;
  int v=nvirt;
  for (f=0; f<v; f++){
      for (m=0; m<o; m++){
          for (e=0; e<v; e++){
              F_DCOPY(o,tb+e*v*o*o+f*o*o+m,o,tempt+f*o*o*v+m*o*v+e*o,1);
              F_DAXPY(o,-2.0,tb+e*v*o*o+f*o*o+m*o,1,tempt+f*o*o*v+m*o*v+e*o,1);
          }
      }
  }

  int tilesize,lasttile,ntiles=1;
  int ov2 = o*v*v;
  // tile v in chunks of o
  tilesize=v;
  for (i=1; i<=v; i++){
      if (o>=(double)tilesize/i){
         tilesize = v/i;
         if (i*tilesize < v) tilesize++;
         ntiles = i;
         break;
      }
  }
  lasttile = v - (ntiles-1)*tilesize;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_ABCI3,PSIO_OPEN_OLD);
  psio_address addr;
  addr = PSIO_ZERO;
  for (i=0; i<ntiles-1; i++){
      psio->read(PSIF_ABCI3,"E2abci3",(char*)&tempv[0],tilesize*ov2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,tilesize,ov2,-1.0,tempt,o,tempv,ov2,1.0,w1+i*tilesize*o,o);
  }
  i=ntiles-1;
  psio->read(PSIF_ABCI3,"E2abci3",(char*)&tempv[0],lasttile*ov2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o,lasttile,ov2,-1.0,tempt,o,tempv,ov2,1.0,w1+i*tilesize*o,o);
  psio->close(PSIF_ABCI3,1);
  psio.reset();
}

void GPUCoupledCluster::CPU_I1ab(){
  int b,m,n,e,a,id;
  int o = ndoccact;
  int v = nvirt;
  // build I1(a,b)
  for (m=0; m<o; m++){
      for (e=0; e<v; e++){
          for (n=0; n<o; n++){
              F_DCOPY(v,E2klcd_1+m*o*v*v+n*v+e,o*v,tempv+m*o*v*v+e*o*v+n*v,1);
          }
      }
  }
  F_DAXPY(o*o*v*v,-2.0,E2klcd_1,1,tempv,1);

  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (n=0; n<o; n++){
              F_DCOPY(v,tb+e*v*o*o+m*o+n,o*o,tempt+m*o*v*v+e*o*v+n*v,1);
              for (b=0; b<v; b++){
                  tempt[id++] += t1[e*o+m]*t1[b*o+n];
              }
          }
      }
  }
  F_DGEMM('n','t',v,v,o*o*v,1.0,tempv,v,tempt,v,0.0,I1,v);

  // add the singles parts to I1(a,b). n^4
  // TODO: can tile in larger blocks
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_ABCI2,PSIO_OPEN_OLD);
  psio_address addr = PSIO_ZERO;
  for (a=0; a<v; a++){
      psio->read(PSIF_ABCI2,"E2abci2",(char*)&tempt[0],v*v*o*sizeof(double),addr,&addr);
      F_DGEMV('t',o*v,v,1.0,tempt,o*v,t1,1,1.0,I1+a*v,1);
  }
  psio->close(PSIF_ABCI2,1);
  psio.reset();

  int l,k,c;
  for (l=0; l<o; l++){
      for (c=0; c<v; c++){
          for (k=0; k<o; k++){
              F_DCOPY(v,tb+c*o*o+l*o+k,o*o*v,tempt+l*o*v*v+c*o*v+k*v,1);
          }
      }
  }
  // use I1(a,b) for doubles residual:
  F_DGEMM('t','n',v,o*o*v,v,1.0,I1,v,tempt,v,0.0,tempv,v);
  for (int a=0; a<v; a++){
      for (int b=0; b<v; b++){
          for (int i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+a*v*o+i*v+b,v*v*o,wb+a*o*o*v+b*o*o+i*o,1);
              F_DAXPY(o,1.0,tempv+i*v*v*o+b*v*o+a,v,wb+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }

  // use I1(a,b) for singles residual - 1st contribution to w1. (n^3)
  F_DGEMM('n','n',o,v,v,1.0,t1,o,I1,v,1.0,w1,o);
}


// CPU_I2p_abci required ov^3 storage.  by refactorizing, we reduce storage to o^3v, but increase cost by 2o^2v^3
// who cares, this sits on the cpu anyway
// TODO: move terms around like in plugin_ccsd_serial
void GPUCoupledCluster::CPU_I2p_abci_refactored(){
  int a,b,i;
  int o = ndoccact;
  int v = nvirt;
  int ov2 = o*v*v;
  int o2v = o*o*v;

  // tilesize * v <= o^2v^2
  int tilesize,lasttile,ntiles=1;
  tilesize=ov2;
  for (i=1; i<=ov2; i++){
      if (o*o*v*v>=(double)tilesize*v/i){
         tilesize = ov2/i;
         if (i*tilesize < ov2) tilesize++;
         ntiles = i;
         break;
      }
  }
  lasttile = ov2 - (ntiles-1)*tilesize;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_ABCI,PSIO_OPEN_OLD);
  psio_address addr;
  addr = PSIO_ZERO;
  for (i=0; i<ntiles-1; i++){
      psio->read(PSIF_ABCI,"E2abci",(char*)&tempv[0],v*tilesize*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,tilesize,v,1.0,t1,o,tempv,v,0.0,tempt+i*tilesize*o,o);
  }
  i=ntiles-1;
  psio->read(PSIF_ABCI,"E2abci",(char*)&tempv[0],v*lasttile*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o,lasttile,v,1.0,t1,o,tempv,v,0.0,tempt+i*tilesize*o,o);
  psio->close(PSIF_ABCI,1);
  psio.reset();
  
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempt+a*o*v+b*o+i,v*v*o,wb+a*o*o*v+b*o*o+i*o,1);
              F_DAXPY(o,1.0,tempt+i*v*v*o+b*o*v+a*o,1,wb+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }


  // now build and use 2 new intermediates:
  F_DGEMM('n','n',o,o2v,v,-1.0,t1,o,E2akjc_2,v,0.0,tempt,o);
  F_DGEMM('n','n',o2v,v,o,1.0,tempt,o2v,t1,o,0.0,tempv,o2v);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+b*v*o*o+a*o*o+i,o,wb+a*o*o*v+b*o*o+i*o,1);
              F_DAXPY(o,1.0,tempv+a*v*o*o+b*o*o+i*o,1,wb+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
  F_DGEMM('t','t',o2v,o,v,-1.0,E2klcd_1,v,t1,o,0.0,tempt,o2v);
  F_DGEMM('t','n',v,o2v,o,1.0,t1,o,tempt,o,0.0,tempv,v);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+i*v*v*o+b*v+a,v*v,wb+a*o*o*v+b*o*o+i*o,1);
              F_DAXPY(o,1.0,tempv+i*v*v+a*v+b,o*v*v,wb+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
}

void GPUCoupledCluster::CPU_I1pij_I1ia_lessmem(){
  int m,j,e,i,a,b;
  int o = ndoccact;
  int v = nvirt;
  int ov2 = o*v*v;
  double*tempw1;
  tempw1 = (double*)malloc(o*v*sizeof(double));

  // build I1(i,a). n^4
  for (m=0; m<o; m++){
      for (e=0; e<v; e++){
          for (j=0; j<o; j++){
              F_DCOPY(v,tb+e*v*o*o+j*o+m,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
              F_DAXPY(v,-2.0,tb+e*v*o*o+m*o+j,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
          }
      }
  }
  for (i=0; i<o; i++){
      for (a=0; a<v; a++){
          for (m=0; m<o; m++){
              F_DCOPY(v,E2klcd_1+i*o*v*v+m*v+a,o*v,tempv+i*v*v*o+a*v*o+m,o);
              F_DAXPY(v,-2.0,E2klcd_1+i*o*v*v+a*o*v+m*v,1,tempv+i*v*v*o+a*v*o+m,o);
          }
      }
  }
  F_DGEMM('t','n',o*v,1,o*v,-1.0,tempv,o*v,t1,o*v,0.0,I1,o*v);

  // use I1(i,a) -> w1
  F_DGEMM('n','n',o*v,1,o*v,-1.0,tempt,o*v,I1,o*v,0.0,tempw1,o*v);
  for (i=0; i<o; i++){
      F_DAXPY(v,1.0,tempw1+i*v,1,w1+i,o);
  }

  // build I1'(i,j)
  F_DGEMM('t','n',o,o,ov2,-1.0,tempt,ov2,E2klcd_1,ov2,0.0,I1p,o);
  
  // only n^4
  for (i=0; i<o; i++){
      for (j=0; j<o; j++){
          for (m=0; m<o; m++){
              F_DCOPY(v,E2ijak+m*o*o*v+i*o*v+j*v,1,tempv+i*o*o*v+j*o*v+m,o);
              F_DAXPY(v,-2.0,E2ijak+i*o*o*v+m*o*v+j*v,1,tempv+i*o*o*v+j*o*v+m,o);
          }
      }
  }
  F_DGEMM('t','n',o*o,1,o*v,-1.0,tempv,o*v,t1,o*v,1.0,I1p,o*o);

  // use I1'(i,j) for singles residual. (n^3)
  F_DGEMM('n','n',o,v,o,-1.0,I1p,o,t1,o,1.0,w1,o);

  // build I1(i,j)
  F_DGEMM('n','n',o,o,v,1.0,t1,o,I1,v,1.0,I1p,o);
  for (m=0; m<o; m++){
      for (e=0; e<v; e++){
          for (j=0; j<o; j++){
              F_DCOPY(v,tb+e*o*o*v+m*o+j,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
          }
      }
  }
  F_DGEMM('n','t',o,ov2,o,-1.0,I1p,o,tempt,ov2,0.0,tempv,o);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+a*o*o*v+b*o+i,v*o,wb+a*o*o*v+b*o*o+i*o,1);
              F_DAXPY(o,1.0,tempv+b*o*o*v+i*v*o+a*o,1,wb+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
  free(tempw1);
}

/*================================================================

   update amplitudes

================================================================*/
void GPUCoupledCluster::UpdateT1(int iter){
  int v = nvirt;
  int o = ndoccact;
  int i,a;
  double tnew,dia;
  for (a=o; a<nmo; a++){
      for (i=0; i<o; i++){
          dia = -eps[i]+eps[a];
          tnew = - (w1[(a-o)*o+i])/dia;
          w1[(a-o)*o+i] = tnew;
      }
  }
  // error vector for diis is in tempv:
  F_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  F_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  F_DCOPY(o*v,w1,1,t1,1);
}
void GPUCoupledCluster::UpdateT2(int iter){
  int v = nvirt;
  int o = ndoccact;
  double tnew,dijab,da,dab,dabi;
  int iajb,jaib,ijab=-1;
  for (long int a=o; a<nmo; a++){
      da = eps[a];
      for (long int b=o; b<nmo; b++){
          dab = da + eps[b];
          for (long int i=0; i<o; i++){
              dabi = dab - eps[i];
              for (long int j=0; j<o; j++){
                  ijab++;
                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);

                  dijab = dabi-eps[j];

                  tnew = - (E2klcd_1[iajb] + wb[ijab])/dijab;
                  tempt[ijab] = tnew;
                  tempu[ijab] = (2.*E2klcd_1[iajb]-E2klcd_1[jaib])*(tnew+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
              }
          }
      }
  }
  F_DCOPY(o*o*v*v,tb,1,tempv,1);
  F_DAXPY(o*o*v*v,-1.0,tempt,1,tempv,1);
  F_DCOPY(o*o*v*v,tempt,1,tb,1);
}

/*================================================================

   gpu functions

================================================================*/

__global__ void AddPermutedOnGPU(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  out[id] = in[id] + in[b*h*ns*ns+a*ns*ns+j*ns+i];
}
__global__ void GPUFill_I2iajk_and_c2Sym1(int ns,int h,double*gpuv,double*gpuw,double*gput2,double*gputempw){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  gpuw[id] += gpuv[id] + gpuv[b*ns*ns*h+a*ns*ns+j*ns+i];
  if (i<=j && a<=b)
     gputempw[GPUPosition(a,b)*ns*(ns+1)/2+GPUPosition(i,j)] = gput2[id]+gput2[b*ns*ns*h+a*ns*ns+i*ns+j];
}
__global__ void GPUFill_I2iajk_add(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  out[id] += in[id] + in[b*ns*ns*h+a*ns*ns+j*ns+i];
}
__global__ void GPUc2Sym1(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  if (i>j) return;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  if (a>b) return;
  out[GPUPosition(a,b)*ns*(ns+1)/2+GPUPosition(i,j)] = in[id]+in[b*ns*ns*h+a*ns*ns+i*ns+j];
}
__global__ void GPUc2Sym2(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  if (i>j) return;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  if (a>b) return;
  out[GPUPosition(a,b)*ns*(ns+1)/2+GPUPosition(i,j)] = in[id]-in[b*ns*ns*h+a*ns*ns+i*ns+j];
}
__global__ void GPUSymmAdd2(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  short sg;
  short sg2=1;
  if (a>b) sg2 = -1;
  if (i>j) sg  = -sg2;
  else     sg  =  sg2;
  out[id] += sg*in[GPUPosition(a,b)*ns*(ns+1)/2+GPUPosition(i,j)];
}
__global__ void GPUSymmAdd1(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  out[id] += in[GPUPosition(a,b)*ns*(ns+1)/2+GPUPosition(i,j)];
}
__global__ void gpuCopyArray(double*out,double*in,int dim){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=dim) return;
  out[id] = in[id];
}

__global__ void  GPUc2Sym2_onefunction(int ns,int h,double*gput2,double*gput1,double*gputempw){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  if (i>j) return;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  if (a>b) return;
  gputempw[GPUPosition(a,b)*ns*(ns+1)/2+GPUPosition(i,j)] =
         gput2[id] - gput2[b*ns*ns*h+a*ns*ns+i*ns+j]
       + gput1[a*ns+i]*gput1[b*ns+j] - gput1[b*ns+i]*gput1[a*ns+j];
}

__global__ void GPUt2Plust1_and_E2klcd3(int ns,int h,double*t2,double*t1,double*E2klcd1,double*E2klcd3){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return; 
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  t2[id] += t1[a*ns+i]*t1[b*ns+j];

  E2klcd3[j*ns*h*h+i*h*h+b*h+a] = E2klcd1[j*ns*h*h+b*ns*h+i*h+a];
}
__global__ void GPUt2Plust1(int ns,int h,double*t2,double*t1){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return; 
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  t2[id] += t1[a*ns+i]*t1[b*ns+j];
}
__device__ int GPUPosition(int i,int j){
  if (i<j){
    return j*(j+1)/2+i;
  }
  return i*(i+1)/2+j;
}
__global__ void GPUt2Plus2t1(int ns,int h,double*out,double*t2,double*t1,double*newt2){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short a = id%h;
  unsigned short j = (id-a)%(ns*h)/h;
  unsigned short b = (id-a-j*h)%(ns*h*h)/(ns*h);
  unsigned short i = (id-a-j*h-b*ns*h)/(ns*h*h);
  //out[id] = t2[i*h*h*ns+a*h*ns+j*h+b] + 2.*t1[a*ns+i]*t1[b*ns+j];
  out[id]  = newt2[i*h*h*ns+a*h*ns+j*h+b] = t2[b*h*ns*ns+a*ns*ns+j*ns+i];
  out[id] += 2.*t1[a*ns+i]*t1[b*ns+j];
}
__global__ void GPUv2MinusHalfv2(int ns,int h,double*out,double*v2){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short a = id%h;
  unsigned short j = (id-a)%(ns*h)/h;
  unsigned short b = (id-a-j*h)%(ns*h*h)/(ns*h);
  unsigned short i = (id-a-j*h-b*ns*h)/(ns*h*h);
  out[id] = (v2[id]-.5*v2[i*h*h*ns+a*h*ns+j*h+b]);
}
__global__ void GPUPermute_iabj_to_aijb(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short b = (id-j)%(h*ns)/ns;
  unsigned short a = (id-j-b*ns)%(ns*h*h)/(h*ns);
  unsigned short i = (id-j-b*ns-a*h*ns)/(ns*h*h);
  out[a*ns*ns*h+i*ns*h+j*h+b] = in[id];
}
__global__ void GPUFill_I2iabj(int ns,int h,double*in1,double*in2,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short a = id%h;
  unsigned short j = (id-a)%(h*ns)/h;
  unsigned short b = (id-a-j*h)%(ns*h*h)/(h*ns);
  unsigned short i = (id-a-j*h-b*h*ns)/(ns*h*h);
  out[id] = in1[id] + in2[a*ns*ns*h+i*ns*h+j*h+b];
}
__global__ void GPU2t2Minust2(int ns,int h,double*out,double*t2){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short e = id%h;
  unsigned short m = (id-e)%(ns*h)/h;
  unsigned short a = (id-e-m*h)%(ns*h*h)/(ns*h);
  unsigned short i = (id-e-m*h-a*ns*h)/(ns*h*h);
  out[id] = 2.*t2[id]-t2[i*ns*h*h+e*ns*h+m*h+a];
}
__global__ void GPUFill_t2_I2iajb(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  out[id] = in[j*ns*h*h+b*h*ns+i*h+a] + in[i*ns*h*h+a*h*ns+j*h+b];
}
__global__ void GPUt2Plus2t1_and_E2klcd2(int ns,int h,double*out,double*t2,double*t1,double*E2klcd1,double*E2klcd2){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short a = id%h;
  unsigned short j = (id-a)%(ns*h)/h;
  unsigned short b = (id-a-j*h)%(ns*h*h)/(ns*h);
  unsigned short i = (id-a-j*h-b*ns*h)/(ns*h*h);
  out[id] = t2[j*h*h*ns+b*h*ns+i*h+a] + 2.*t1[a*ns+i]*t1[b*ns+j];
  //out[id] = t2[a*h*ns*ns+b*ns*ns+i*ns+j] + 2.*t1[a*ns+i]*t1[b*ns+j];

  E2klcd2[id] = E2klcd1[i*h*h*ns+a*h*ns+j*h+b]; 
}
__global__ void GPUt2Plus2t1_2(int ns,int h,double*out,double*t2,double*t1){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short a = id%h;
  unsigned short j = (id-a)%(ns*h)/h;
  unsigned short b = (id-a-j*h)%(ns*h*h)/(ns*h);
  unsigned short i = (id-a-j*h-b*ns*h)/(ns*h*h);
  out[id] = t2[j*h*h*ns+b*h*ns+i*h+a] + 2.*t1[a*ns+i]*t1[b*ns+j];
}
__global__ void GPUFillI2iajb2(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short a = id%h;
  unsigned short j = (id-a)%(h*ns)/h;
  unsigned short b = (id-a-j*h)%(ns*h*h)/(h*ns);
  unsigned short i = (id-a-j*h-b*h*ns)/(ns*h*h);
  out[id] += in[a*ns*ns*h+i*ns*h+j*h+b];
}
__global__ void GPUFillI2iajb1(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short a = id%h;
  unsigned short j = (id-a)%(h*ns)/h;
  unsigned short b = (id-a-j*h)%(ns*h*h)/(h*ns);
  unsigned short i = (id-a-j*h-b*h*ns)/(ns*h*h);
  out[id] += in[i*ns*h*h+a*h*ns+b*ns+j];
}
__global__ void GPUFill_t2_I2iabj1(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short j = id%ns;
  unsigned short i = (id-j)%(ns*ns)/ns;
  unsigned short b = (id-j-i*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-j-i*ns-b*ns*ns)/(ns*ns*h);
  out[id] = in[j*ns*h*h+b*h*ns+i*h+a] + in[i*ns*h*h+a*h*ns+j*h+b];
}
__global__ void GPUFill_t2_I2iabj2(int ns,int h,double*in,double*out){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short i = id%ns;
  unsigned short j = (id-i)%(ns*ns)/ns;
  unsigned short b = (id-i-j*ns)%(ns*ns*h)/(ns*ns);
  unsigned short a = (id-i-j*ns-b*ns*ns)/(ns*ns*h);
  out[id] += in[j*ns*h*h+b*h*ns+i*h+a] + in[i*ns*h*h+a*h*ns+j*h+b];
}
__global__ void GPUPermute_tikbc(double*in,double*out,int ns,int h){
  int blockid = blockIdx.x*gridDim.y + blockIdx.y;
  int id = blockid*blockDim.x + threadIdx.x;
  if (id>=ns*ns*h*h) return;
  unsigned short d = id%h;
  unsigned short k = (id-d)%(ns*h)/h;
  unsigned short c = (id-d-k*h)%(ns*h*h)/(ns*h);
  unsigned short l = (id-d-k*h-c*ns*h)/(ns*h*h);
  out[id] = in[l*h*h*ns+d*h*ns+k*h+c];
}

/*
 * truncated cepa
 */
void GPUCoupledCluster::TCEPA(){
  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,jaib,ijab=0;
  double energy = 0.0;
  memset((void*)tb,'\0',o*o*v*v*sizeof(double));
  double*pair_energy = (double*)malloc(o*o*sizeof(double));

  double tconv = 1.0;
  double di,dij,dija,dijab;
  while (tconv>1e-6){
      tconv = 0.0;
      for (i=0; i<o; i++){
          di = -eps[i];
          for (j=0; j<o; j++){
              dij = di - eps[j];
              energy=0.0;
              for (a=o; a<rs; a++){
                  for (b=o; b<rs; b++){
                      ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                      iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                      jaib = j*v*v*o+(a-o)*v*o+i*v+(b-o);
                      energy += E2klcd_1[iajb]*(2.0*tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                      //energy += (2.*E2klcd_1[iajb]-E2klcd_1[jaib])*(tb[ijab]+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                  }
              }
              pair_energy[i*o+j] = energy;
              for (a=o; a<rs; a++){
                  dija = dij + eps[a];
                  for (b=o; b<rs; b++){
                      dijab = dija + eps[b];
                      ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                      iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                      double t2  = -E2klcd_1[iajb]/(dijab - pair_energy[i*o+j]);
                      double dum = (tb[ijab]-t2);
                      tconv += dum*dum;
                      tb[ijab] = t2;
                  }
              }
          }
      }
      for (a=0; a<v; a++){
          t1[a*o+i] = sqrt(fabs(tb[a*o*o*v+a*o*o+i*o+i]));
      }
      energy=0.0;
      for (i=0; i<o*o; i++) energy += pair_energy[i];
      printf("%20.12lf\n",energy);fflush(stdout);
      tconv = sqrt(tconv);
  }
  free(pair_energy);
}


}//end of namespace psi
