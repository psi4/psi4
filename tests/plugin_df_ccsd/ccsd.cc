#include"psi4-dec.h"
#include<libmints/wavefunction.h>
#include<libmints/vector.h>
#include<sys/times.h>
#ifdef _OPENMP
    #include<omp.h>
#endif

#include"gpuhelper.h"
#include"blas.h"
#include"ccsd.h"

using namespace psi;

// position in a symmetric packed matrix
long int Position(long int i,long int j){
  if (i<j){
    return ((j*(j+1))>>1)+i;
  }
  return ((i*(i+1))>>1)+j;
}

namespace psi{


  /*!
   ** PSIO_GET_ADDRESS(): Given a starting page/offset and a shift length
   ** (in bytes), return the page/offset of the next position in the file.
   ** \ingroup PSIO
   */

  psio_address psio_get_address(psio_address start, long int shift) {
    psio_address address;
    long int bytes_left;

    bytes_left = PSIO_PAGELEN - start.offset; /* Bytes remaining on fpage */

    if (shift >= bytes_left) { /* Shift to later page */
      address.page = start.page + (shift - bytes_left)/PSIO_PAGELEN+ 1;
      address.offset = shift - bytes_left -(address.page - start.page- 1)
          *PSIO_PAGELEN;
    } else { /* Block starts on current page */
      address.page = start.page;
      address.offset = start.offset + shift;
    }

    return address;
  }

CoupledCluster::CoupledCluster()
{}
CoupledCluster::~CoupledCluster()
{}

void CoupledCluster::WriteBanner(){
  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                       DF-CCSD                       *\n");
  fprintf(outfile, "        *              Integral-Direct CCSD with              *\n");
  fprintf(outfile, "        *                   Density Fitting                   *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                   Eugene DePrince                   *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);
}

/*================================================================
  
  Initialize:
  set basic parameters (ndocc...). integral transformation.
  integral sort.
  
================================================================*/
void CoupledCluster::Initialize(Options &options){
 
  // grab the reference wave function and its parameters
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.reference_wavefunction();

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
     throw PsiException("plugin_ccsd requires symmetry c1",__FILE__,__LINE__);
  }
  nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
  long int full=0;
  for (long int h=0; h<nirreps; h++){
      nfzc   += fzc[h];
      nfzv   += fzv[h];
      nso    += sorbs[h];
      full   += orbs[h];
      nmo    += orbs[h]-fzc[h]-fzv[h];
      ndocc  += docc[h];//-fzc[h];
  }
  ndoccact = ndocc - nfzc;
  nvirt  = nmo - ndoccact;

  // for triples, we use nvirt_no in case we've truncated the virtual space:
  nvirt_no = nvirt;
  scale_t = 1.0;

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

  // initialize gpu helper class
  helper_ = boost::shared_ptr<GPUHelper>(new GPUHelper);

  // get device parameters, allocate gpu memory and pinned cpu memory
  helper_->ndoccact = ndoccact;
  helper_->nvirt    = nvirt;
  helper_->nmo      = nmo;

  helper_->CudaInit(options);

  // reduce available memory by the amount required by the helper class
  memory -= helper_->max_mapped_memory;

  // quit if max_mapped_memory exceeds available memory
  if ((double)memory<0){
     throw PsiException("max_mapped_memory must be less than available memory",__FILE__,__LINE__);
  }

  // quit if number of virtuals is less than number of doubly occupied
  if (nvirt<ndoccact){
     throw PsiException("ndocc must be larger than nvirt",__FILE__,__LINE__);
  }

  // so->mo tei transformation
  struct tms total_tmstime;
  const long clk_tck = sysconf(_SC_CLK_TCK);

  long int i;
  double time_start,user_start,sys_start,time_stop,user_stop,sys_stop;

  // grab 3-index integrals
  times(&total_tmstime);
  time_start = time(NULL);
  user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;

  //DensityFittedIntegrals();
  ThreeIndexIntegrals();

  times(&total_tmstime);
  time_stop = time(NULL);
  user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;
  fprintf(outfile,"\n");
  //fprintf(outfile,"  Time for density fitting:         %6.2lf s (user)\n",user_stop-user_start);
  fprintf(outfile,"  Time for 3-index transform:       %6.2lf s (user)\n",user_stop-user_start);
  fprintf(outfile,"                                    %6.2lf s (system)\n",sys_stop-sys_start);
  fprintf(outfile,"                                    %6d s (total)\n",(int)time_stop-(int)time_start);
  fprintf(outfile,"\n");

  // orbital energies
  eps_test = ref->epsilon_a();
  double*tempeps = eps_test->pointer();
  eps = (double*)malloc(nmo*sizeof(double));
  F_DCOPY(nmo,tempeps+nfzc,1,eps,1);
  eps_test.reset();

}
/*===================================================================

  solve ccsd equations

===================================================================*/
PsiReturnType CoupledCluster::CCSDIterations(Options&options){

  // timer stuff:
  struct tms total_tmstime;
  const long clk_tck = sysconf(_SC_CLK_TCK);
  time_t iter_start,iter_stop,time_start,time_stop;
  double user_start,user_stop,sys_start,sys_stop;

  int iter,sg,sg2,diis_iter,replace_diis_iter;
  long int o = ndoccact;
  long int v = nvirt;
  long int arraysize = o*o*v*v;
  long int ov2 = o*v*v;
  long int oo1o2 = o*(o+1)/2;
  long int vv1o2 = v*(v+1)/2;
  double nrm,Eold,s1,start,end,siter;

  iter=0;
  diis_iter=0;
  replace_diis_iter=1;
  nrm=1.0;
  Eold=1.0e9;
  eccsd=0.0;

  fprintf(outfile,"\n");
  fprintf(outfile,
    "  Begin singles and doubles coupled cluster iterations ... on GPUs!\n\n");
  fprintf(outfile,
    "   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");
  fflush(outfile);

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  // zero residual
  psio->open(PSIF_R2,PSIO_OPEN_NEW);
  memset((void*)tempt,'\0',o*o*v*v*sizeof(double));
  psio->write_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);

  // cc diagrams split up as tasks
  DefineTasks();

  // start timing the iterations
  times(&total_tmstime);
  time_start = time(NULL);
  user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;
  while(iter<maxiter && nrm>conv){
      iter_start = time(NULL);

      // evaluate cc diagrams
      if (iter>0){
         memset((void*)w1,'\0',o*v*sizeof(double));
         for (int i=0; i<ncctasks; i++) {
             double start = 0.0;
             #ifdef _OPENMP
                 start = omp_get_wtime();
             #endif
             (*this.*CCTasklist[i].func)(CCParams[i]);
             double end = 0.0;
             #ifdef _OPENMP
                 end = omp_get_wtime();
             #endif
             //printf("task %5i took %10.2lf s\n",i,end-start);fflush(stdout);
         }
      }

      // update the amplitudes and check the energy
      Eold = eccsd;
      UpdateT1(iter);
      eccsd = UpdateT2(iter);

      // add vector to list for diis
      DIISOldVector(iter,diis_iter,replace_diis_iter);

      // diis error vector and convergence check
      nrm = DIISErrorVector(diis_iter,replace_diis_iter,iter);

      // diis extrapolation
      if (diis_iter>2){
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,arraysize+o*v);
         else                   DIIS(diisvec,maxdiis,arraysize+o*v);
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
  fflush(outfile);

  return Success;
}
double CoupledCluster::CheckEnergy(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  double ta,tnew,dijab,da,dab,dabi;
  long int iajb,jaib,ijab=0;
  double energy = 0.0;
  // df (ia|bj) formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);
                  energy += (2.*integrals[iajb]-integrals[jaib])*(tb[ijab]+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                  ijab++;
              }
          }
      }
  }
  return energy;
}
void CoupledCluster::SCS_MP2(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,jaib,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  // df (ia|bj) formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);
                  osenergy += integrals[iajb]*tb[ijab];
                  ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                  ijab++;
              }
          }
      }
  }
  emp2_os = osenergy;
  emp2_ss = ssenergy;
}
void CoupledCluster::SCS_CCSD(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,jaib,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  // df (ia|bj) formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);
                  osenergy += integrals[iajb]*(tb[ijab]+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                  ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                  ssenergy += integrals[iajb]*(t1[(a-o)*o+i]*t1[(b-o)*o+j]-t1[(b-o)*o+i]*t1[(a-o)*o+j]);
                  ijab++;
              }
          }
      }
  }
  eccsd_os = osenergy;
  eccsd_ss = ssenergy;
}

/*===================================================================

  determine tiling for vabcd and vabci diagrams for the cpu
  this determines the size of blocks of integrals that 
  can be read into cpu memory.

  this function doesn't really do anything now that we're fitting
  all the integrals.  well, it still checks to see if there is 
  enough memory to do the calculation.

===================================================================*/
void CoupledCluster::DefineTilingCPU(){
  long int i,v = nvirt;
  long int o = ndoccact;
  long int ov2 = o*v*v;
  long int ov = o*v;
  long int o2 = o*o;

  // number of doubles in total memory
  long int ndoubles = memory/8L;
  // minus storage for other necessary buffers 
  ndoubles -= 3L*o*o*v*v+5L*o*v+v*v+(o+v);

  if (ndoubles<0){
     throw PsiException("out of memory.",__FILE__,__LINE__);
  }

  ntiles = -999L;
  tilesize = v*(v+1L)/2L;
  ntiles = 1L;

  // check whether blocking of the vabcd diagram is necessary
  if (ndoubles>v*(v+1L)/2L*v*(v+1L)/2L){
     tilesize = v*(v+1L)/2L;
     ntiles = 1L;
  }
  else{
     for (i=2L; i<=v*(v+1L)/2L; i++){
         if (ndoubles>tilesize*v*(v+1L)/2L/i+1L){
            tilesize = v*(v+1L)/2L/i;
            if (i*tilesize < v*(v+1L)/2L) tilesize++;
            ntiles = i;
            break;
         }
     }
     if (ntiles==-999L){
        throw PsiException("out of memory: (ab,cd)",__FILE__,__LINE__);
     }
  }
  lasttile = v*(v+1L)/2L - (ntiles-1L)*tilesize;

  fprintf(outfile,"        v(ab,cd) diagrams will be evaluated in %3li blocks.\n",ntiles); 
  fflush(outfile);

  // ov^3 type 1:
  if (v>ndoubles){
     throw PsiException("out of memory: (ab,ci)",__FILE__,__LINE__);
  }
  nov2tiles=1L;
  ov2tilesize=ov2/1L;
  if (nov2tiles*ov2tilesize<ov2) ov2tilesize++;
  while(v*ov2tilesize>ndoubles){
     nov2tiles++;
     ov2tilesize = ov2/nov2tiles;
     if (nov2tiles*ov2tilesize<ov2) ov2tilesize++;
  }
  lastov2tile = ov2 - (nov2tiles-1L)*ov2tilesize;

  fprintf(outfile,"        v(ab,ci) diagrams will be evaluated in %3li blocks over ov2.\n",nov2tiles); 
  fflush(outfile);

  // ov^3 type 2:
  if (v*v>ndoubles){
     throw PsiException("out of memory: (ab,ci)",__FILE__,__LINE__);
  }
  novtiles=1L;
  ovtilesize=ov/1L;
  if (novtiles*ovtilesize<ov) ovtilesize++;
  while(v*v*ovtilesize>ndoubles){
     novtiles++;
     ovtilesize = ov/novtiles;
     if (novtiles*ovtilesize<ov) ovtilesize++;
  }
  lastovtile = ov - (novtiles-1L)*ovtilesize;
  fprintf(outfile,"        v(ab,ci) diagrams will be evaluated in %3li blocks over ov.\n",novtiles); 
  fflush(outfile);
}

/*===================================================================

  allocate cpu memory

===================================================================*/
void CoupledCluster::AllocateMemory(Options&options){

  long int i,o=ndoccact;
  long int v=nvirt;
  long int dim;

  // subtract out storage for df integrals from the total available memory
  memory -= 8L*nQ*(o*o+o*v+v*v);

  // and storage for all of the other necessary buffers 
  memory -= 8L*(4L*o*o*v*v+5L*o*v+v*v+(o+v));

  if (memory<0){
     fprintf(outfile,"\n");
     fprintf(outfile,"  error: not enough memory.  increase available memory by %7.2lf mb\n",-(double)memory/1024./1024.);
     fprintf(outfile,"\n");
     fflush(outfile);
     throw PsiException("not enough memory.",__FILE__,__LINE__);
  }

  /*========================================================
     ccsd memory requirements:
    
     tb:     o^2v^2
     tempt:  o^2v^2+ov
     tempv:  o^2v^2+ov
     integrals: max(o^2v^2, 2v^3, 2nQ*o*v) (this is a minimum)
     other stuff: 2ov+2v^2+(o+v)
    
     total: 3o^2v^2 + 2v^3 + 4ov + 2v^2 + (o+v)  or 
            4o^2v^2 + 4ov + 2v^2 + (o+v)         or
            3o^2v^2 + 2ovnQ + 4ov + 2v^2 + (o+v)

     compare to the requirements for the (T) part:

            2o^2v^2 + 3v^3*nthreads + o^3v + ov
    
  ========================================================*/

  // for the df version, the dimension of the large buffer:
  dim = 2*v*v*v;
  if (2*nQ*o*v>dim) dim = 2*nQ*o*v;
  if (o*o*v*v>dim)  dim = o*o*v*v;

  /*========================================================

     tiling of the v^4 diagram:

  ========================================================*/

  long int ndoubles = (memory + 8L*o*o*v*v) / 8L;
  ntiles=1L;
  tilesize=v/1L;
  if (ntiles*tilesize<v) tilesize++;
  while(2*v*v*v*tilesize>ndoubles){
     ntiles++;
     tilesize = v/ntiles;
     if (ntiles*tilesize<v) tilesize++;
  }
  lasttile = v - (ntiles-1L)*tilesize;

  fprintf(outfile,"        v(ab,cd) diagram will be evaluated in %3li blocks.\n",ntiles); 
  fprintf(outfile,"\n");
  fflush(outfile);

  if (dim<2*v*v*v*tilesize) dim = 2*v*v*v*tilesize;
 

  maxelem = dim;

  double total_memory = dim+2.*(o*o*v*v+o*v)+1.*o*o*v*v+2.*o*v+2.*v*v;
  total_memory *= 8./1024./1024.;
  double df_memory = nQ*(o*o+o*v+v*v);
  df_memory *= 8./1024./1024.;

  fprintf(outfile,"  Total memory requirements: %9.2lf mb\n",df_memory+total_memory);
  fprintf(outfile,"  3-index integrals:         %9.2lf mb\n",df_memory);
  fprintf(outfile,"  CCSD intermediates:        %9.2lf mb\n",total_memory);

  integrals = (double*)malloc(dim*sizeof(double));
  tempt     = (double*)malloc((o*o*v*v+o*v)*sizeof(double));
  tempv     = (double*)malloc((o*o*v*v+o*v)*sizeof(double));
  tb        = (double*)malloc(o*o*v*v*sizeof(double));
  w1        = (double*)malloc(o*v*sizeof(double));
  t1        = (double*)malloc(o*v*sizeof(double));
  I1        = (double*)malloc(v*v*sizeof(double));
  I1p       = (double*)malloc(v*v*sizeof(double));

  memset((void*)integrals,'\0',dim*sizeof(double));
  memset((void*)tempv,'\0',(o*o*v*v+o*v)*sizeof(double));
  memset((void*)tempt,'\0',(o*o*v*v+o*v)*sizeof(double));
  memset((void*)tb,'\0',o*o*v*v*sizeof(double));
  memset((void*)w1,'\0',o*v*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));
  memset((void*)I1,'\0',v*v*sizeof(double));
  memset((void*)I1p,'\0',v*v*sizeof(double));

  // DIIS:
  diisvec    = (double*)malloc(sizeof(double)*(maxdiis+1));
  memset((void*)diisvec,'\0',(maxdiis+1)*sizeof(double));

}

void CoupledCluster::CPU_t1_vmeai(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int i,m,q;

  // df 2 (me|ai) t(m,e)
  for (i=0; i<o; i++){
      F_DCOPY(v,t1+i,o,tempt+i*v,1);
  }
  F_DGEMV('t',o*v,nQ,2.0,Qov,o*v,tempt,1,0.0,integrals,1);
  F_DGEMV('n',o*v,nQ,1.0,Qov,o*v,integrals,1,0.0,tempv,1);
  for (i=0; i<o; i++){
      F_DAXPY(v,1.0,tempv+i*v,1,w1+i,o);
  }

  // df - (ae|mi) t(e,m)
  F_DGEMM('t','t',nQ*v,o,v,-1.0,Qvv,v,t1,o,0.0,integrals,nQ*v);
  for (m=0; m<o; m++){
      for (q=0; q<nQ; q++){
          for (i=0; i<o; i++){
              integrals[nQ*v*o + m*nQ*o + q*o + i] = Qoo[q*o*o+m*o+i];
          }
      }
  }
  F_DGEMM('n','t',o,v,o*nQ,1.0,integrals+nQ*v*o,o,integrals,v,1.0,w1,o);
}

void CoupledCluster::CPU_t1_vmeni(CCTaskParams params){
  long int i,m,e,n,a,id;
  long int o=ndoccact;
  long int v=nvirt;
  for (m=0,id=0; m<o; m++){
  for (e=0; e<v; e++){
  for (n=0; n<o; n++){
  for (a=0; a<v; a++){
      tempt[id++] = 2.*tb[e*v*o*o+a*o*o+m*o+n]-tb[a*v*o*o+e*o*o+m*o+n];
  }}}}

  // df (me|ni) ( 2 t(mn,ea) - t(mn,ae) )
  F_DGEMM('n','n',o*v,nQ,o*v,-1.0,tempt,o*v,Qov,o*v,0.0,integrals,o*v);
  F_DGEMM('n','t',o,v,nQ*o,1.0,Qoo,o,integrals,v,1.0,w1,o);

}

void CoupledCluster::CPU_t1_vmaef(CCTaskParams params){
  long int q,m,e,i,f,a,id;
  long int o=ndoccact;
  long int v=nvirt;

  for (m=0,id=0; m<o; m++){
  for (e=0; e<v; e++){
  for (f=0; f<v; f++){
  for (i=0; i<o; i++){
      tempt[id++] = 2.*tb[e*v*o*o+f*o*o+m*o+i]-tb[e*v*o*o+f*o*o+i*o+m];
  }}}}

  // df
  F_DGEMM('n','n',o*v,nQ,o*v,1.0,tempt,o*v,Qov,o*v,0.0,integrals,o*v);
  F_DGEMM('n','t',o,v,nQ*v,1.0,integrals,o,Qvv,v,1.0,w1,o);
}

void CoupledCluster::CPU_I1ab(CCTaskParams params){
  long int i,j,l,k,c,d;
  long int o = ndoccact;
  long int v = nvirt;
  long int q,b,m,n,e,a,id=0;
  // build I1(a,b)
  boost::shared_ptr<PSIO> psio(new PSIO());

  // df -2 (me|nb) ( t(eb,mn) + t(e,m)t(b,n) )
  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (a=0; a<v; a++){
              for (n=0; n<o; n++){
                  tempt[id++] = tb[e*v*o*o+a*o*o+m*o+n] + t1[e*o+m]*t1[a*o+n];
              }
          }
      }
  }
  F_DGEMM('t','t',nQ,o*v,o*v,-2.0,Qov,o*v,tempt,o*v,0.0,integrals,nQ);
  for (a=0; a<v; a++){
      for (q=0; q<nQ; q++){
          F_DCOPY(o,integrals+a*nQ*o+q,nQ,tempv+q*o,1);
      }
      F_DCOPY(o*nQ,tempv,1,integrals+a*nQ*o,1);
  }
  F_DGEMM('n','n',v,v,o*nQ,1.0,Qov,v,integrals,o*nQ,0.0,I1,v);

  // df (mb|ne) ( t(eb,mn) + t(e,m)t(b,n) )
  for (n=0,id=0; n<o; n++){
      for (e=0; e<v; e++){
          for (a=0; a<v; a++){
              for (m=0; m<o; m++){
                  tempt[id++] = tb[e*v*o*o+a*o*o+m*o+n] + t1[e*o+m]*t1[a*o+n];
              }
          }
      }
  }
  F_DGEMM('t','t',nQ,o*v,o*v,1.0,Qov,o*v,tempt,o*v,0.0,integrals,nQ);
  for (a=0; a<v; a++){
      for (q=0; q<nQ; q++){
          F_DCOPY(o,integrals+a*nQ*o+q,nQ,tempv+q*o,1);
      }
      F_DCOPY(o*nQ,tempv,1,integrals+a*nQ*o,1);
  }
  F_DGEMM('n','n',v,v,o*nQ,1.0,Qov,v,integrals,o*nQ,1.0,I1,v);
  
  // add the singles parts to I1(a,b). n^4
  for (i=0; i<o; i++){
      F_DCOPY(v,t1+i,o,tempt+i*v,1);
  }

  // df ( 2 (ab|me) - (ie|mj) ) t(e,m),  ov^3 term
  F_DGEMV('t',o*v,nQ,2.0,Qov,o*v,tempt,1,0.0,integrals,1);
  F_DGEMV('n',v*v,nQ,1.0,Qvv,v*v,integrals,1,1.0,I1,1);

  F_DGEMM('t','n',v*nQ,o,v,-1.0,Qvv,v,tempt,v,0.0,integrals,v*nQ);
  for (long int q=0; q<nQ; q++){
      for (m=0; m<o; m++){
          for (a=0; a<v; a++){
              integrals[nQ*o*v+q*o*v+m*v+a] = integrals[m*nQ*v+q*v+a];
          }
      }
  }
  F_DGEMM('n','t',v,v,o*nQ,1.0,Qov,v,integrals+nQ*o*v,v,1.0,I1,v);

  id=0;
  for (l=0; l<o; l++){
  for (c=0; c<v; c++){
  for (k=0; k<o; k++){
  for (d=0; d<v; d++){
      tempt[id++] = tb[d*v*o*o+c*o*o+l*o+k];
  }}}}
  // use I1(a,b) for doubles residual:
  helper_->GPUTiledDGEMM_NoThread('t','n',v,o*o*v,v,1.0,I1,v,tempt,v,0.0,tempv,v,0);
  //F_DGEMM('t','n',v,o*o*v,v,1.0,I1,v,tempt,v,0.0,tempv,v);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++)
  for (b=0; b<v; b++){
  for (i=0; i<o; i++)
  for (j=0; j<o; j++)
      tempt[id++] += tempv[j*v*v*o+a*v*o+i*v+b]+tempv[i*v*v*o+b*v*o+j*v+a];
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);

  // use I1(a,b) for singles residual - 1st contribution to w1. (n^3)
  //helper_->GPUTiledDGEMM('n','n',o,v,v,1.0,t1,o,I1,v,1.0,w1,o);
  F_DGEMM('n','n',o,v,v,1.0,t1,o,I1,v,1.0,w1,o);

  psio.reset();
}

void CoupledCluster::CPU_I2p_abci_refactored_term1(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int q,e,a,b,c,i,j,id=0;
  long int ov2 = o*v*v;
  long int o2v = o*o*v;

  // df 
  F_DGEMM('n','n',o,v*nQ,v,1.0,t1,o,Qvv,v,0.0,integrals,o);
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,integrals,o*v,Qov,o*v,0.0,tempv,o*v);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempt[b*v*o*o+a*o*o+i*o+j] = tempv[j*o*v*v+b*o*v+a*o+i];
              }
          }
      }
  }

  // contribute to residual
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  addr = PSIO_ZERO;
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
  for (b=0; b<v; b++){
      F_DAXPY(o*o,1.0,tempt+b*v*o*o+a*o*o,1,tempv+a*v*o*o+b*o*o,1);
  }}
  for (a=0; a<v; a++){
  for (b=0; b<v; b++){
  for (i=0; i<o; i++){
      F_DAXPY(o,1.0,tempt+a*v*o*o+b*o*o+i,o,tempv+a*v*o*o+b*o*o+i*o,1);
  }}}
  psio->write_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);
  psio.reset();
}
void CoupledCluster::CPU_I2p_abci_refactored_term2(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int a,b,c,i,j,id=0;

  // df -(ae|mj) t(e,i) t(b,m)
  F_DGEMM('n','n',o,nQ*v,v,-1.0,t1,o,Qvv,v,0.0,integrals,o);
  F_DGEMM('n','t',o*o,o*v,nQ,1.0,Qoo,o*o,integrals,o*v,0.0,tempv,o*o);
  F_DGEMM('t','n',o*o*v,v,o,1.0,tempv,o,t1,o,0.0,tempt,o*o*v);
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempv[id++] = tempt[b*o*o*v+a*o*o+i*o+j];
              }
          }
      }
  }

  // contribute to residual
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  F_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+a*v*o*o+b*o*o+i*o,1,tempt+b*v*o*o+a*o*o+i,o);
          }
      }
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);
  psio.reset();
}
void CoupledCluster::CPU_I2p_abci_refactored_term3(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int a,b,c,i,j,id;

  // df -(me|bj) t(e,i) t(a,m)
  F_DGEMM('n','n',o,o*nQ,v,-1.0,t1,o,Qov,v,0.0,integrals,o);
  F_DGEMM('n','t',o*v,o*o,nQ,1.0,Qov,o*v,integrals,o*o,0.0,tempt,o*v);
  F_DGEMM('t','t',v,o*o*v,o,1.0,t1,o,tempt,o*o*v,0.0,tempv,v);

  // contribute to residual
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempt[id++] += tempv[i*v*v*o+j*v*v+b*v+a] + tempv[j*v*v*o+i*v*v+a*v+b];
              }
          }
      }
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);
  psio.reset();
}

void CoupledCluster::CPU_I1pij_I1ia_lessmem(CCTaskParams params){

  long int o = ndoccact;
  long int v = nvirt;
  long int m,j,e,f,i,a,b;//,one=1;
  long int ov2 = o*v*v;
  long int id=0;

  // df 2 (ia|jb), formerly E2klcd 
  // ... just build these because there is a term below where the scaling
  // won't improve without a gigantic intermediate ( o^3 v nQ )
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);


  // build I1(i,a). n^4
  boost::shared_ptr<PSIO> psio(new PSIO());
  F_DCOPY(o*o*v*v,integrals,1,tempv,1);
  for (i=0; i<o; i++){
      for (a=0; a<v; a++){
          for (m=0; m<o; m++){
              for (e=0; e<v; e++){
                  tempv[i*v*v*o+a*v*o+m*v+e] -= 0.5*integrals[i*o*v*v+e*o*v+m*v+a];
              }
          }
      }
  }
  for (i=0; i<o; i++) F_DCOPY(v,t1+i,o,tempt+i*v,1);
  F_DGEMV('t',o*v,o*v,2.0,tempv,o*v,tempt,1,0.0,I1,1);

  // use I1(i,a) -> w1
  id=0;
  for (m=0; m<o; m++){
      for (e=0; e<v; e++){
          for (j=0; j<o; j++){
              for (f=0; f<v; f++){
                  tempt[id++] = 2*tb[e*o*o*v+f*o*o+m*o+j]-tb[e*v*o*o+f*o*o+j*o+m];
              }
          }
      }
  }
  F_DGEMV('n',o*v,o*v,1.0,tempt,o*v,I1,1,0.0,tempv,1);
  for (i=0; i<o; i++){
      F_DAXPY(v,1.0,tempv+i*v,1,w1+i,o);
  }

  // build I1'(i,j)
  helper_->GPUTiledDGEMM_NoThread('t','n',o,o,ov2,1.0,tempt,ov2,integrals,ov2,0.0,I1p,o,0);
  //F_DGEMM('t','n',o,o,ov2,1.0,tempt,ov2,integrals,ov2,0.0,I1p,o);
  
  // df 2 (ij|me) t(m,e)
  for (i=0; i<o; i++){
      F_DCOPY(v,t1+i,o,tempt+i*v,1);
  }
  F_DGEMV('t',o*v,nQ,2.0,Qov,o*v,tempt,1,0.0,integrals,1);
  F_DGEMV('n',o*o,nQ,1.0,Qoo,o*o,integrals,1,1.0,I1p,1);

  // df - (ie|mj) t(e,m)
  F_DGEMM('t','t',o*nQ,o,v,-1.0,Qov,v,t1,o,0.0,integrals,o*nQ);
  for (long int q=0; q<nQ; q++){
      for (m=0; m<o; m++){
          for (i=0; i<o; i++){
              integrals[nQ*o*o+q*o*o+m*o+i] = integrals[m*nQ*o+q*o+i];
          }
      }
  }
  F_DGEMM('n','t',o,o,o*nQ,1.0,Qoo,o,integrals+nQ*o*o,o,1.0,I1p,o);


  // use I1'(i,j) for singles residual. (n^3)
  F_DGEMM('n','n',o,v,o,-1.0,I1p,o,t1,o,1.0,w1,o);

  // build I1(i,j)
  F_DGEMM('n','n',o,o,v,1.0,t1,o,I1,v,1.0,I1p,o);
  for (m=0,id=0; m<o; m++){
  for (e=0; e<v; e++){
  for (j=0; j<o; j++){
  for (f=0; f<v; f++){
      tempt[id++] = tb[e*o*o*v+f*o*o+m*o+j];
  }}}}
  //helper_->GPUTiledDGEMM_NoThread('n','t',o,ov2,o,-1.0,I1p,o,tempt,ov2,0.0,tempv,o,0);
  helper_->GPUTiledDGEMM('n','t',o,ov2,o,-1.0,I1p,o,tempt,ov2,0.0,tempv,o);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++)
  for (b=0; b<v; b++){
  for (i=0; i<o; i++)
  for (j=0; j<o; j++)
      tempt[id++] += tempv[a*o*o*v+j*v*o+b*o+i] + tempv[b*o*o*v+i*v*o+a*o+j];
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);

  psio.reset();
}

/*================================================================

   update amplitudes

================================================================*/
void CoupledCluster::UpdateT1(long int iter){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int id=0;
  double tnew,dia;
  if (iter<1){
     memset((void*)t1,'\0',o*v*sizeof(double));
     memset((void*)w1,'\0',o*v*sizeof(double));
  }
  else{
     for (a=o; a<rs; a++){
         for (i=0; i<o; i++){
             dia = -eps[i]+eps[a];
             tnew = - (w1[(a-o)*o+i])/dia;
             w1[(a-o)*o+i] = tnew;
         }
     }
  }
  // error vector for diis is in tempv:
  F_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  F_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  F_DCOPY(o*v,w1,1,t1,1);
}
double CoupledCluster::UpdateT2(long int iter){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  double ta,tnew,dijab,da,dab,dabi;
  long int iajb,jaib,ijab=0;
  double energy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());


  // df (ia|bj) formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

  // we still have the residual in memory in tempv
  //psio->open(PSIF_R2,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=o; a<rs; a++){
      da = eps[a];
      for (b=o; b<rs; b++){
          dab = da + eps[b];
          for (i=0; i<o; i++){
              dabi = dab - eps[i];
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);

                  dijab = dabi-eps[j];

                  tnew = - (integrals[iajb] + tempv[ijab])/dijab;
                  tempt[ijab] = tnew;
                  energy += (2.*integrals[iajb]-integrals[jaib])*(tnew+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                  ijab++;
              }
          }
      }
  }

  // error vectors for diis are in tempv:
  F_DCOPY(o*o*v*v,tempt,1,tempv,1);
  F_DAXPY(o*o*v*v,-1.0,tb,1,tempv,1);
  F_DCOPY(o*o*v*v,tempt,1,tb,1);

  psio.reset();

  return energy;
}

/*================================================================

   diis functions

================================================================*/
void CoupledCluster::DIIS(double*c,long int nvec,long int n){
  long int i,j,k;
  doublereal sum,dum;
  integer*ipiv,nvar;
  nvar = nvec+1;
  doublereal*A = (doublereal*)malloc(sizeof(doublereal)*nvar*nvar);
  doublereal*B = (doublereal*)malloc(sizeof(doublereal)*nvar);
  memset((void*)A,'\0',nvar*nvar*sizeof(double));
  memset((void*)B,'\0',nvar*sizeof(double));
  B[nvec] = -1.;
  ipiv = (integer*)malloc(nvar*sizeof(integer));

  char*evector=(char*)malloc(1000*sizeof(char));

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_EVEC,PSIO_OPEN_OLD);

  for (i=0; i<nvec; i++){
      sprintf(evector,"evector%li",i+1);
      psio->read_entry(PSIF_EVEC,evector,(char*)&tempt[0],n*sizeof(double));
      for (j=i+1; j<nvec; j++){
          sprintf(evector,"evector%li",j+1);
          psio->read_entry(PSIF_EVEC,evector,(char*)&tempv[0],n*sizeof(double));
          sum = F_DDOT(n,tempt,1,tempv,1);
          A[j*nvar+i] = sum;
          A[i*nvar+j] = sum;
      }
      A[i*nvar+i] = F_DDOT(n,tempt,1,tempt,1);
  }
  j = nvec;
  for (i=0; i<nvar; i++){
      A[j*nvar+i] = -1.0;
      A[i*nvar+j] = -1.0;
  }
  A[nvar*nvar-1] = 0.;
  psio->close(PSIF_EVEC,1);
  free(evector);

  integer nrhs,lda,ldb,info;
  nrhs = 1;
  lda = ldb = nvar;
  info = 0;
  DGESV(nvar,nrhs,A,lda,ipiv,B,ldb,info);
  F_DCOPY(nvec,B,1,c,1);

  free(A);
  free(B);
  free(ipiv);
  psio.reset();
}
void CoupledCluster::DIISOldVector(long int iter,int diis_iter,int replace_diis_iter){
  long int j,o = ndoccact;
  long int arraysize,v = nvirt;
  arraysize=o*o*v*v;

  char*oldvector=(char*)malloc(1000*sizeof(char));

  if (diis_iter<=maxdiis && iter<=maxdiis){
     sprintf(oldvector,"oldvector%i",diis_iter);
  }
  else{
     sprintf(oldvector,"oldvector%i",replace_diis_iter);
  }

  boost::shared_ptr<PSIO> psio(new PSIO());
  if (diis_iter==0)
     psio->open(PSIF_OVEC,PSIO_OPEN_NEW);
  else
     psio->open(PSIF_OVEC,PSIO_OPEN_OLD);

  psio_address addr;
  addr = PSIO_ZERO;
  psio->write(PSIF_OVEC,oldvector,(char*)&tb[0],arraysize*sizeof(double),addr,&addr);
  psio->write(PSIF_OVEC,oldvector,(char*)&t1[0],o*v*sizeof(double),addr,&addr);
  psio->close(PSIF_OVEC,1);
  psio.reset();

  free(oldvector);
}
double CoupledCluster::DIISErrorVector(int diis_iter,int replace_diis_iter,int iter){
  double nrm;
  long int i,j,o = ndoccact;
  long int arraysize,v = nvirt;
  arraysize=o*o*v*v;

  char*evector   = (char*)malloc(1000*sizeof(char));
  if (diis_iter<=maxdiis && iter<=maxdiis){
     sprintf(evector,"evector%i",diis_iter);
  }
  else{
     sprintf(evector,"evector%i",replace_diis_iter);
  }

  boost::shared_ptr<PSIO> psio(new PSIO());
  if (diis_iter==0)
     psio->open(PSIF_EVEC,PSIO_OPEN_NEW);
  else
     psio->open(PSIF_EVEC,PSIO_OPEN_OLD);

  nrm = F_DNRM2(arraysize+o*v,tempv,1);
  psio->write_entry(PSIF_EVEC,evector,(char*)&tempv[0],(arraysize+o*v)*sizeof(double));

  psio->close(PSIF_EVEC,1);
  psio.reset();

  free(evector);

  // return convergence
  return nrm;
}
void CoupledCluster::DIISNewAmplitudes(int diis_iter){
  long int o = ndoccact;
  long int arraysize,v = nvirt;
  arraysize=o*o*v*v;

  char*oldvector;
  oldvector=(char*)malloc(1000*sizeof(char));

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_OVEC,PSIO_OPEN_OLD);

  psio_address addr;
  memset((void*)tb,'\0',arraysize*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));

  int max = diis_iter;
  if (max > maxdiis) max = maxdiis;

  for (long int j=1; j<=max; j++){
      addr = PSIO_ZERO;
      sprintf(oldvector,"oldvector%li",j);
      psio->read(PSIF_OVEC,oldvector,(char*)&tempt[0],arraysize*sizeof(double),addr,&addr);
      F_DAXPY(arraysize,diisvec[j-1],tempt,1,tb,1);
      psio->read(PSIF_OVEC,oldvector,(char*)&tempt[0],o*v*sizeof(double),addr,&addr);
      F_DAXPY(o*v,diisvec[j-1],tempt,1,t1,1);
  }
  psio->close(PSIF_OVEC,1);
  free(oldvector);
  psio.reset();
}

/**
 *  Build and use I2ijkl
 */
void CoupledCluster::I2ijkl(CCTaskParams params){
  long int id,i,j,a,b,o,v,k,l;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());

  // df (ia|jb) integrals, formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

  for (j=0; j<o; j++){
      for (i=0; i<o; i++){
          for (b=0; b<v; b++){
              F_DCOPY(v,integrals+j*o*v*v+b*o*v+i*v,1,tempv+j*o*v*v+i*v*v+b*v,1);
          }
      }
  }

  // df (ie|jl) t(e,k)
  F_DGEMM('n','n',o,nQ*o,v,2.0,t1,o,Qov,v,0.0,integrals,o);
  F_DGEMM('n','t',o*o,o*o,nQ,1.0,Qoo,o*o,integrals,o*o,0.0,tempt,o*o);
  for (i=0,id=0; i<o; i++){
      for (j=0; j<o; j++){
          for (k=0; k<o; k++){
              for (l=0; l<o; l++){
                  integrals[id++] = tempt[i*o*o*o+k*o*o+j*o+l];
              }
          }
      }
  }
  
  // df (ik|jl) integrals, formerly E2ijkl
  F_DGEMM('n','t',o*o,o*o,nQ,1.0,Qoo,o*o,Qoo,o*o,0.0,tempt,o*o);
  for (i=0; i<o; i++){
      for (j=0; j<o; j++){
          for (k=0; k<o; k++){
              for (l=0; l<o; l++){
                  integrals[i*o*o*o+k*o*o+j*o+l] += tempt[i*o*o*o+j*o*o+k*o+l];
              }
          }
      }
  }

  F_DCOPY(o*o*v*v,tb,1,tempt,1);
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempt[id++] += t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }
  // (ie|jf) c(ef,kl)
  helper_->GPUTiledDGEMM('n','n',o*o,o*o,v*v,1.0,tempt,o*o,tempv,v*v,1.0,integrals,o*o);

  // I2(ij,kl) -> R(ij,ab)
  helper_->GPUTiledDGEMM('n','n',o*o,v*v,o*o,0.5,integrals,o*o,tempt,o*o,0.0,tempv,o*o);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  F_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+b*v*o*o+a*o*o+i,o,tempt+a*v*o*o+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);
  psio.reset();

}
/**
 *  Build and use I2'iajk
 */
void CoupledCluster::I2piajk(CCTaskParams params){
  long int id,k,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;


  // df (ij|ak), formerly E2ijak2
  F_DGEMM('n','t',o*o,o*v,nQ,1.0,Qoo,o*o,Qov,o*v,0.0,tempt,o*o);
  for (i=0; i<o; i++){
      for (j=0; j<o; j++){
          for (k=0; k<o; k++){
              for (a=0; a<v; a++){
                  tempv[j*o*o*v+a*o*o+k*o+i] = tempt[i*o*o*v+a*o*o+j*o+k];
              }
          }
      }
  }

  F_DCOPY(o*o*v*v,tb,1,tempt,1);
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempt[id++] += t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }

  // two options for df (ie|af) c(ef,jk)
  for (i=0; i<o; i++){
      for (long int q=0; q<nQ; q++){
          F_DCOPY(v,Qov+q*o*v+i*v,1,integrals+v*v*v+q*v,1);
      }
      F_DGEMM('n','t',v,v*v,nQ,1.0,integrals+v*v*v,v,Qvv,v*v,0.0,integrals,v);
      for (a=0; a<v; a++){
          for (b=0; b<v; b++){
              F_DCOPY(v,integrals+a*v*v+b,v,integrals+v*v*v+a*v*v+b*v,1);
          }
      }
      helper_->GPUTiledDGEMM('n','n',o*o,v,v*v,1.0,tempt,o*o,integrals+v*v*v,v*v,1.0,tempv+i*o*o*v,o*o);
  }
  // alternate (more expensive) form for df (ie|af) c(ef,jk)
  /*for (i=0; i<o; i++){
      for (long int q=0; q<nQ; q++){
          for (long int e=0; e<v; e++){
              integrals[o*o*v*nQ+q*v+e] = Qov[q*o*v+i*v+e];
          }
      }
      F_DGEMM('n','n',o*o*v,nQ,v,1.0,tempt,o*o*v,integrals+o*o*v*nQ,v,0.0,integrals,o*o*v);
      F_DGEMM('n','t',o*o,v,v*nQ,1.0,integrals,o*o,Qvv,v,1.0,tempv+i*o*o*v,o*o);
  }*/

  /*addr = PSIO_ZERO;
  psio->open(PSIF_ABCI,PSIO_OPEN_OLD);
  for (j=0; j<novtiles-1; j++){
      psio->read(PSIF_ABCI,"E2abci",(char*)&integrals[0],ovtilesize*v*v*sizeof(double),addr,&addr);
      helper_->GPUTiledDGEMM('n','n',o*o,ovtilesize,v*v,1.0,tempt,o*o,integrals,v*v,1.0,tempv+j*o*o*ovtilesize,o*o);
  }
  j=novtiles-1;
  psio->read(PSIF_ABCI,"E2abci",(char*)&integrals[0],lastovtile*v*v*sizeof(double),addr,&addr);
  helper_->GPUTiledDGEMM('n','n',o*o,lastovtile,v*v,1.0,tempt,o*o,integrals,v*v,1.0,tempv+j*o*o*ovtilesize,o*o);
  psio->close(PSIF_ABCI,1);*/

  helper_->GPUTiledDGEMM_NoThread('n','n',o*o*v,v,o,-1.0,tempv,o*o*v,t1,o,0.0,tempt,o*o*v,0);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  F_DAXPY(o*o*v*v,1.0,tempt,1,tempv,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempt+b*v*o*o+a*o*o+i,o,tempv+a*v*o*o+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);
  psio.reset();
}
/**
 *  Use Vabcd1
 */
void CoupledCluster::Vabcd1(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tb[id++] += t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  double time1,t2,start1,start2,end1,end2,t3,start3,end3;
  time1=t2=start1=start2=end1=end2=0.0;

  psio->open(PSIF_QVVT,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_QVVT,"qvv_transpose",(char*)&Qvv[0],nQ*v*v*sizeof(double));
  psio->close(PSIF_QVVT,1);

  for (long int atile=0; atile<ntiles-1; atile++){
      #ifdef _OPENMP
          start1 = omp_get_wtime();
      #endif
      helper_->GPUTiledDGEMM('t','n',v*v,v*tilesize,nQ,1.0,Qvv,nQ,Qvv+atile*tilesize*v*nQ,nQ,0.0,integrals,v*v);
      #ifdef _OPENMP
          end1 = omp_get_wtime();
      #endif
      #ifdef _OPENMP
          start2 = omp_get_wtime();
      #endif
      for (a=0; a<tilesize; a++){
          for (b=0; b<v; b++){
              for (long int c=0; c<v; c++){
                  F_DCOPY(v,integrals+a*v*v*v+c*v*v+b*v,1,integrals+a*v*v*v+b*v*v+c*v+tilesize*v*v*v,1);
                  //for (long int d=0; d<v; d++){
                  //    integrals[a*v*v*v+b*v*v+c*v+d+tilesize*v*v*v] = integrals[a*v*v*v+c*v*v+b*v+d];
                  //}
              }
          }
          //F_DCOPY(v*v*v,integrals+tilesize*v*v*v,1,integrals+a*v*v*v,1);
      }
      #ifdef _OPENMP
          end2 = omp_get_wtime();
      #endif
      #ifdef _OPENMP
          start3 = omp_get_wtime();
      #endif
      helper_->GPUTiledDGEMM('n','n',o*o,v*tilesize,v*v,1.0,tb,o*o,integrals+tilesize*v*v*v,v*v,1.0,tempv+atile*tilesize*o*o*v,o*o);
      #ifdef _OPENMP
          end3 = omp_get_wtime();
      #endif
      time1 += end1-start1;
      t2 += end2-start2;
      t3 += end3-start3;
  }
  long int atile = ntiles-1;
  #ifdef _OPENMP
      start1 = omp_get_wtime();
  #endif
  helper_->GPUTiledDGEMM('t','n',v*v,v*lasttile,nQ,1.0,Qvv,nQ,Qvv+atile*tilesize*v*nQ,nQ,0.0,integrals,v*v);
  #ifdef _OPENMP
      end1 = omp_get_wtime();
  #endif
  #ifdef _OPENMP
      start2 = omp_get_wtime();
  #endif
  for (a=0; a<lasttile; a++){
      for (b=0; b<v; b++){
          for (long int c=0; c<v; c++){
              //for (long int d=0; d<v; d++){
              //    integrals[a*v*v*v+b*v*v+c*v+d+tilesize*v*v*v] = integrals[a*v*v*v+c*v*v+b*v+d];
              //}
              F_DCOPY(v,integrals+a*v*v*v+c*v*v+b*v,1,integrals+a*v*v*v+b*v*v+c*v+tilesize*v*v*v,1);
          }
      }
      //F_DCOPY(v*v*v,integrals+tilesize*v*v*v,1,integrals+a*v*v*v,1);
  }
  #ifdef _OPENMP
      end2 = omp_get_wtime();
  #endif
  #ifdef _OPENMP
      start3 = omp_get_wtime();
  #endif
  helper_->GPUTiledDGEMM('n','n',o*o,v*lasttile,v*v,1.0,tb,o*o,integrals+tilesize*v*v*v,v*v,1.0,tempv+atile*tilesize*o*o*v,o*o);
  #ifdef _OPENMP
      end3 = omp_get_wtime();
  #endif
  time1 += end1-start1;
  t2 += end2-start2;
  t3 += end3-start3;



  psio->open(PSIF_QVV,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_QVV,"qvv",(char*)&Qvv[0],nQ*v*v*sizeof(double));
  psio->close(PSIF_QVV,1);

  /*for (a=0; a<v; a++){
             #ifdef _OPENMP
                 start3 = omp_get_wtime();
             #endif
      F_DCOPY(nQ*v,Qvv+a,v,tempt,1);
             #ifdef _OPENMP
                 end3 = omp_get_wtime();
             #endif
             #ifdef _OPENMP
                 start1 = omp_get_wtime();
             #endif
      helper_->GPUTiledDGEMM('n','t',v*v,v,nQ,1.0,Qvv,v*v,tempt,v,0.0,integrals,v*v);
             #ifdef _OPENMP
                 end1 = omp_get_wtime();
             #endif
             #ifdef _OPENMP
                 start2 = omp_get_wtime();
             #endif
      helper_->GPUTiledDGEMM('n','t',o*o,v,v*v,1.0,tb,o*o,integrals,v,1.0,tempv+a*o*o*v,o*o);
             #ifdef _OPENMP
                 end2 = omp_get_wtime();
             #endif
          time1 += end1-start1;
          t2 += end2-start2;
          t3 += end3-start3;
  }*/
//printf("%7.2lf %7.2lf %7.2lf\n",time1,t2,t3);
  //for (a=0; a<v; a++){
  //    F_DCOPY(nQ*v,Qvv+a,v,tempt,1);
  //    F_DGEMM('n','n',o*o*v,nQ,v,1.0,tb,o*o*v,tempt,v,0.0,integrals,o*o*v);
  //    F_DGEMM('n','t',o*o,v,v*nQ,1.0,integrals,o*o,Qvv,v,0.0,tempv+a*o*o*v,o*o);
  //}
  psio->write_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);
  psio.reset();

  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tb[id++] -= t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }

  // using SJS packing is very very slow.
  /*F_DCOPY(o*o*v*v,tb,1,tempt,1);
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempt[id++] += t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }
  for (i=0; i<o; i++){
      for (j=i; j<o; j++){
          for (a=0; a<v; a++){
              for (b=0; b<v; b++){
                  //tempv[Position(a,b)*o*(o+1)/2+Position(i,j)] =
                  tempv[a*v*o*(o+1)/2+b*o*(o+1)/2+Position(i,j)] =
                     0.5*(tempt[a*o*o*v+b*o*o+i*o+j]+tempt[b*o*o*v+a*o*o+i*o+j]);
              }
              //tempv[Position(a,a)*o*(o+1)/2+Position(i,j)] =
              //   tempt[a*o*o*v+a*o*o+i*o+j];
          }
      }
  }
  double*tempq = (double*)malloc(v*v*nQ*sizeof(double));
  for (a=0; a<v*v; a++){
      for (long int q=0; q<nQ; q++){
          tempq[a*nQ+q] = Qvv[q*v*v+a];
      }
  }
  double t1,t2,start1,start2,end1,end2;
  t1=t2=start1=start2=end1=end2=0.0;
  long int count = 0;
  long int otri = o*(o+1)/2;
  for (b=0; b<v; b++){
      for (a=0; a<=b; a++){
          //for (long int d=0; d<v; d++){
          //    for (long int c=0; c<v; c++){
          //        double acbd = F_DDOT(nQ,tempq+a*v*nQ+c*nQ,1,tempq+b*v*nQ+d*nQ,1);
          //        double adbc = F_DDOT(nQ,tempq+a*v*nQ+d*nQ,1,tempq+b*v*nQ+c*nQ,1);
          //        integrals[count++]   = acbd+adbc;
          //    }
          //}
             #ifdef _OPENMP
                 start1 = omp_get_wtime();
             #endif
          //helper_->GPUTiledDGEMM('t','n',v,v,nQ,2.0,tempq+a*v*nQ,nQ,tempq+b*v*nQ,nQ,0.0,integrals,v);
          F_DGEMM('t','n',v,v,nQ,2.0,tempq+a*v*nQ,nQ,tempq+b*v*nQ,nQ,0.0,integrals,v);
             #ifdef _OPENMP
                 end1 = omp_get_wtime();
             #endif
          //F_DGEMM('t','n',v,v,nQ,1.0,tempq+b*v*nQ,nQ,tempq+a*v*nQ,nQ,1.0,integrals,v);
             #ifdef _OPENMP
                 start2 = omp_get_wtime();
             #endif
          F_DGEMV('n',otri,v*v,1.0,tempv,otri,integrals,1,0.0,tempt+otri*Position(a,b),1);
          //F_DGEMM('n','t',1,o*(o+1)/2,v*v,1.0,integrals,1,tempv,o*(o+1)/2,0.0,tempt+Position(a,b)*o*(o+1)/2,1);
          //helper_->GPUTiledDGEMM('n','t',1,o*(o+1)/2,v*v,1.0,integrals,1,tempv,o*(o+1)/2,0.0,tempt+Position(a,b)*o*(o+1)/2,1);
             #ifdef _OPENMP
                 end2 = omp_get_wtime();
             #endif
          t1 += end1-start1;
          t2 += end2-start2;
      }
  }
printf("%7.2lf %7.2lf\n",t1,t2);
  free(tempq);*/
  /*psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
  addr = PSIO_ZERO;
  for (j=0; j<ntiles-1; j++){
      psio->read(PSIF_ABCD1,"E2abcd1",(char*)&integrals[0],tilesize*v*(v+1)/2*sizeof(double),addr,&addr);
      helper_->GPUTiledDGEMM('n','n',o*(o+1)/2,tilesize,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  }
  j=ntiles-1;
  psio->read(PSIF_ABCD1,"E2abcd1",(char*)&integrals[0],lasttile*v*(v+1)/2*sizeof(double),addr,&addr);
  helper_->GPUTiledDGEMM('n','n',o*(o+1)/2,lasttile,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  psio->close(PSIF_ABCD1,1);*/

  // contribute to residual
  /*psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempv[a*o*o*v+b*o*o+i*o+j] += .5*tempt[Position(a,b)*o*(o+1)/2+Position(i,j)];
              }
          }
      }
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);
  psio.reset();*/

}

/**
 *  Use Vabcd2
 */
void CoupledCluster::Vabcd2(CCTaskParams params){
return;
  long int id,i,j,a,b,o,v;
  int sg,sg2;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  F_DCOPY(o*o*v*v,tb,1,tempt,1);
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempt[id++] += t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }
  for (i=0; i<o; i++){
      for (j=i; j<o; j++){
          for (a=0; a<v; a++){
              for (b=a; b<v; b++){
                  tempv[Position(a,b)*o*(o+1)/2+Position(i,j)] =
                  //tempv[a*o*(o+1)/2*v+b*o*(o+1)/2+Position(i,j)] =
                    tempt[a*o*o*v+b*o*o+i*o+j]-tempt[b*o*o*v+a*o*o+i*o+j];
              }
          }
      }
  }
  // using SJS packing is very very slow.
  double*tempq = (double*)malloc(v*v*nQ*sizeof(double));
  for (a=0; a<v*v; a++){
      for (long int q=0; q<nQ; q++){
          tempq[a*nQ+q] = Qvv[q*v*v+a];
      }
  }
  for (b=0; b<v; b++){
      for (a=0; a<=b; a++){
          long int count=0;
          for (long int d=0; d<v; d++){
              for (long int c=0; c<=d; c++){
                  double acbd = F_DDOT(nQ,tempq+a*v*nQ+c*nQ,1,tempq+b*v*nQ+d*nQ,1);
                  double adbc = F_DDOT(nQ,tempq+a*v*nQ+d*nQ,1,tempq+b*v*nQ+c*nQ,1);
                  integrals[count++]   = acbd-adbc;
              }
          }
          F_DGEMV('n',o*(o+1)/2,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,1,0.0,tempt+Position(a,b)*o*(o+1)/2,1);
      }
  }
  free(tempq);

  /*psio->open(PSIF_ABCD2,PSIO_OPEN_OLD);
  addr = PSIO_ZERO;
  for (j=0; j<ntiles-1; j++){
      psio->read(PSIF_ABCD2,"E2abcd2",(char*)&integrals[0],tilesize*v*(v+1)/2*sizeof(double),addr,&addr);
      helper_->GPUTiledDGEMM('n','n',o*(o+1)/2,tilesize,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  }
  j = ntiles-1;
  psio->read(PSIF_ABCD2,"E2abcd2",(char*)&integrals[0],lasttile*v*(v+1)/2*sizeof(double),addr,&addr);
  helper_->GPUTiledDGEMM('n','n',o*(o+1)/2,lasttile,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  psio->close(PSIF_ABCD2,1);*/

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          if (a>b) sg2 = -1;
          else     sg2 = 1;
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  if (i>j) sg = -1;
                  else     sg = 1;
                  tempv[a*o*o*v+b*o*o+i*o+j] += .5*sg2*sg*tempt[Position(a,b)*o*(o+1)/2+Position(i,j)];
              }
          }
      }
  }
  //psio->write_entry(PSIF_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);
  psio.reset();
}
/**
 *  Build and use I2iabj
 */
void CoupledCluster::I2iabj(CCTaskParams params){
  long int id,i,j,k,a,b,o,v,mytile;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  for (i=0,id=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  tempt[id++] = tb[b*v*o*o+a*o*o+j*o+i]+2.*t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }
  
  // df (ia|jb) integrals, formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,tempv,o*v);

  // df (ib|me) ( t(mj,ae) + 2t(m,a)t(j,e) )
  F_DGEMM('t','t',nQ,o*v,o*v,-0.5,Qov,o*v,tempt,o*v,0.0,integrals,nQ);
  F_DGEMM('n','n',o*v,o*v,nQ,1.0,Qov,o*v,integrals,nQ,0.0,tempt,o*v);
  for (i=0,id=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  tempv[id++] += tempt[j*o*v*v+a*o*v+i*v+b];
              }
          }
      }
  }

  // o^2v^3 contribution to intermediate
  // df (ib|mj) t(a,m), used to be E2ijak
  F_DGEMM('n','t',o*o,o*v,nQ,1.0,Qoo,o*o,Qov,o*v,0.0,tempt,o*o);
  for (i=0; i<o; i++){
      for (j=0; j<o; j++){
          for (k=0; k<o; k++){
              for (a=0; a<v; a++){
                  integrals[j*o*o*v+i*o*v+k*v+a] = tempt[i*o*o*v+a*o*o+j*o+k];
              }
          }
      }
  }

  helper_->GPUTiledDGEMM_NoThread('n','n',o*o*v,v,o,-1.0,integrals,o*o*v,t1,o,0.0,tempt,o*o*v,0);
  for (i=0,id=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  tempv[id++] += tempt[a*o*o*v+i*o*v+j*v+b];
              }
          }
      }
  }
  // contribute to intermediate
  psio->open(PSIF_TEMP,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_TEMP,"temporary",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_TEMP,1);

  // df (ia|jb) integrals, formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,tempt,o*v);

  F_DCOPY(o*o*v*v,tempt,1,tempv,1);
  for (i=0,id=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  tempv[id++] -= 0.5*tempt[i*v*v*o+a*v*o+j*v+b];
              }
          }
      }
  }
  //memset((void*)tempt,'\0',o*o*v*v*sizeof(double));
  for (i=0,id=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  integrals[i*v*v*o+a*v*o+j*v+b] = tb[b*v*o*o+a*o*o+j*o+i];
              }
          }
      }
  }
  helper_->GPUTiledDGEMM('n','n',o*v,o*v,o*v,1.0,integrals,o*v,tempv,o*v,0.0,tempt,o*v);

  // o^2v^3 piece of intermediate
  // df
  F_DGEMM('n','n',o,nQ*v,v,1.0,t1,o,Qvv,v,0.0,integrals,o);
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,integrals,o*v,Qov,o*v,0.0,tempv,o*v);
  for (i=0; i<o; i++){
      for (a=0; a<v; a++){
          for (b=0; b<v; b++){
              for (j=0; j<o; j++){
                  tempt[i*o*v*v+b*o*v+j*v+a] += tempv[i*o*v*v+b*o*v+a*o+j];
              }
          }
      }
  }

  // contribute to intermediate
  psio->open(PSIF_TEMP,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_TEMP,"temporary",(char*)&tempv[0],o*o*v*v*sizeof(double));
  F_DAXPY(o*o*v*v,1.0,tempt,1,tempv,1);
  psio->close(PSIF_TEMP,0);

  // use I2iabj
  for (j=0,id=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (a=0; a<v; a++){
                  tempt[id++] = 2*tb[a*o*o*v+b*o*o+i*o+j]-tb[b*o*o*v+a*o*o+i*o+j];
              }
          }
      }
  }

  //memset((void*)integrals,'\0',o*o*v*v*sizeof(double));
  helper_->GPUTiledDGEMM('n','n',o*v,o*v,o*v,1.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  // if we KNOW this is the first diagram, we don't need to read in the old
  // residual.
  //psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  //tempt[id++] += integrals[j*o*v*v+b*v*o+i*v+a] + integrals[i*o*v*v+a*v*o+j*v+b];
                  tempt[id++] = integrals[j*o*v*v+b*v*o+i*v+a] + integrals[i*o*v*v+a*v*o+j*v+b];
              }
          }
      }
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);

  psio.reset();

}
/**
 *  Build and use I2iajb
 */
void CoupledCluster::I2iajb(CCTaskParams params){
  long int id,i,j,a,b,o,v,k,c;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  // df (ia|jb) integrals, formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,tempt,o*v);

  for (i=0,id=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  integrals[id] = tb[b*o*o*v+a*o*o+j*o+i] + 2.*t1[a*o+i]*t1[b*o+j];
                  tempv[id++] = tempt[i*v*v*o+a*v*o+j*v+b];
              }
          }
      }
  }

  helper_->GPUTiledDGEMM('n','n',o*v,o*v,o*v,-0.5,integrals,o*v,tempv,o*v,0.0,tempt,o*v);

  // df (ij|ab) integrals, formerly E2akjc2
  F_DGEMM('n','t',v*v,o*o,nQ,1.0,Qvv,v*v,Qoo,o*o,0.0,tempv,v*v);
  for (k=0; k<o; k++){
      for (c=0; c<v; c++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  tempt[k*o*v*v+c*o*v+j*v+a] += tempv[k*o*v*v+j*v*v+a*v+c];
              }
          }
      }
  }

  // df (ie|ab) t(e,j)
  F_DGEMM('n','n',o,nQ*o,v,1.0,t1,o,Qov,v,0.0,integrals,o);
  F_DGEMM('n','t',o*o,v*v,nQ,1.0,integrals,o*o,Qvv,v*v,0.0,tempv,o*o);
  for (i=0,id=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  tempt[id++] += tempv[a*o*o*v+b*o*o+i*o+j];
              }
          }
      }
  }

  // df (ij|mb) integrals, formerly E2ijak2
  // for (ij|mb) t(a,m)
  F_DGEMM('n','t',o*o,o*v,nQ,1.0,Qoo,o*o,Qov,o*v,0.0,tempv,o*o);
  for (i=0; i<o; i++){
      for (j=0; j<o; j++){
          for (k=0; k<o; k++){
              for (a=0; a<v; a++){
                  integrals[j*o*o*v+a*o*o+k*o+i] = tempv[i*o*o*v+a*o*o+j*o+k];
              }
          }
      }
  }

  // TODO: this was a problem with cuda 3.2 vs 4.0
  helper_->GPUTiledDGEMM_NoThread('t','n',o*o*v,v,o,-1.0,integrals,o,t1,o,0.0,tempv,o*o*v,0);
  //F_DGEMM('t','n',o*o*v,v,o,-1.0,integrals,o,t1,o,0.0,tempv,o*o*v);
  for (i=0,id=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (a=0; a<v; a++){
                  //tempt[id++] += tempv[a*o*o*v+i*o*v+j*v+b];
                  tempt[id++] += tempv[a*o*o*v+i*o*v+b*o+j];
              }
          }
      }
  }

  // contribute to intermediate
  //psio->open(PSIF_TEMP,PSIO_OPEN_NEW);
  //psio->write_entry(PSIF_TEMP,"temporary",(char*)&tempt[0],o*o*v*v*sizeof(double));
  //psio->close(PSIF_TEMP,1);

  // use I2iajb
  for (j=0,id=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (a=0; a<v; a++){
                  integrals[id++] = tb[b*v*o*o+a*o*o+j*o+i];
              }
          }
      }
  }
  //psio->open(PSIF_TEMP,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_TEMP,"temporary",(char*)&tempt[0],o*o*v*v*sizeof(double));
  //psio->close(PSIF_TEMP,1);

  helper_->GPUTiledDGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,integrals,o*v,0.0,tempv,o*v);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  integrals[id++] += tempv[j*o*v*v+b*v*o+i*v+a] + tempv[i*o*v*v+a*v*o+j*v+b];
              }
          }
      }
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);

  // use I2iajb
  for (j=0,id=0; j<o; j++){
      for (a=0; a<v; a++){
          for (i=0; i<o; i++){
              for (b=0; b<v; b++){
                  tempv[id++] = tb[b*v*o*o+a*o*o+j*o+i];
              }
          }
      }
  }

  helper_->GPUTiledDGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,tempv,o*v,0.0,integrals,o*v);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              for (i=0; i<o; i++){
                  tempt[id++] += integrals[j*o*v*v+b*v*o+i*v+a] + integrals[i*o*v*v+a*v*o+j*v+b];
              }
          }
      }
  }
  psio->write_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_R2,1);

  psio.reset();
}

/**
 *  Tasks:
 */
void CoupledCluster::DefineTasks(){
  CCTasklist = new CCTask[1000];
  CCParams   = new CCTaskParams[1000];
  long int o = ndoccact;
  long int v = nvirt;

  ncctasks=0;

  CCTasklist[ncctasks++].func = &psi::CoupledCluster::I2iabj;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::I2iajb;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::I2ijkl;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::I2piajk;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::CPU_t1_vmeni;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::CPU_t1_vmaef;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::CPU_I2p_abci_refactored_term1;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::CPU_I2p_abci_refactored_term2;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::CPU_I2p_abci_refactored_term3;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::CPU_I1ab;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::CPU_t1_vmeai;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::CPU_I1pij_I1ia_lessmem;
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::Vabcd1;

  // this is the last diagram that contributes to doubles residual,
  // so we can keep it in memory rather than writing and rereading
  CCTasklist[ncctasks++].func = &psi::CoupledCluster::Vabcd2;
}


} // end of namespace psi
