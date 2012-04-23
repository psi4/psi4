#include<libmints/wavefunction.h>
#include<libmints/vector.h>
#include<libmints/matrix.h>
#include<libpsio/psio.hpp>
#include<sys/times.h>
#ifdef _OPENMP
    #include<omp.h>
#endif

#include"blas.h"
#include"cepa.h"

using namespace psi;

namespace psi{
  void OutOfCoreSort(int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,bool islocal);
};

// position in a symmetric packed matrix
long int Position(long int i,long int j){
  if (i<j){
    return ((j*(j+1))>>1)+i;
  }
  return ((i*(i+1))>>1)+j;
}

namespace psi{

CoupledPair::CoupledPair(boost::shared_ptr<psi::Wavefunction> wfn,Options&options)
{
  wfn_ = wfn;
  options_ = options;
}
CoupledPair::~CoupledPair()
{}

void CoupledPair::WriteBanner(Options&options){
  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  if (options.get_str("CEPA_LEVEL")=="CEPA0"){
     fprintf(outfile, "        *                       CEPA(0)                       *\n");
     fprintf(outfile, "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options.get_str("CEPA_LEVEL")=="CEPA1"){
     fprintf(outfile, "        *                       CEPA(1)                       *\n");
     fprintf(outfile, "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options.get_str("CEPA_LEVEL")=="CEPA2"){
     fprintf(outfile, "        *                       CEPA(2)                       *\n");
     fprintf(outfile, "        *        Coupled Electron Pair Approximation          *\n");
  }
  if (options.get_str("CEPA_LEVEL")=="CEPA3"){
     fprintf(outfile, "        *                       CEPA(3)                       *\n");
     fprintf(outfile, "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options.get_str("CEPA_LEVEL")=="ACPF"){
     fprintf(outfile, "        *                        ACPF                         *\n");
     fprintf(outfile, "        *          Averaged Coupled Pair Functional           *\n");
  }
  else if (options.get_str("CEPA_LEVEL")=="AQCC"){
     fprintf(outfile, "        *                        AQCC                         *\n");
     fprintf(outfile, "        *         Averaged Quadratic Coupled Cluster          *\n");
  }
  else if (options.get_str("CEPA_LEVEL")=="CISD"){
     fprintf(outfile, "        *                        CISD                         *\n");
     fprintf(outfile, "        *      Singles Doubles Configuration Interaction      *\n");
  }


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
void CoupledPair::Initialize(Options &options){
 
  if (wfn_.get() !=NULL){
     escf    = Process::environment.globals["SCF TOTAL ENERGY"];
     nirreps = wfn_->nirrep();
     sorbs   = wfn_->nsopi();
     orbs    = wfn_->nmopi();
     docc    = wfn_->doccpi();
     fzc     = wfn_->frzcpi();
     fzv     = wfn_->frzvpi();
  }
  if (nirreps>1){
     //throw PsiException("plugin_cepa requires symmetry c1",__FILE__,__LINE__);
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

  // SCS MP2 and CEPA
  emp2_os_fac = options.get_double("MP2_SCALE_OS");
  emp2_ss_fac = options.get_double("MP2_SCALE_SS");
  ecepa_os_fac = options.get_double("CEPA_SCALE_OS");
  ecepa_ss_fac = options.get_double("CEPA_SCALE_SS");

  // which cepa level? 0,1,2,3
  // also, -1 = cisd
  // also, -2 = acpf
  // also, -3 = aqcc
  std::string cepa = options.get_str("CEPA_LEVEL");
  if (cepa == "CEPA0") cepa_level = 0;
  if (cepa == "CEPA1") cepa_level = 1;
  if (cepa == "CEPA2") cepa_level = 2;
  if (cepa == "CEPA3") cepa_level = 3;
  if (cepa == "CISD") cepa_level = -1;
  if (cepa == "ACPF") cepa_level = -2;
  if (cepa == "AQCC") cepa_level = -3;
  cepa_type = (char*)malloc(100*sizeof(char));
  if (cepa_level == 0)       sprintf(cepa_type,"CEPA(0)");
  else if (cepa_level == 1)  sprintf(cepa_type,"CEPA(1)");
  else if (cepa_level == 2)  sprintf(cepa_type,"CEPA(2)");
  else if (cepa_level == 3)  sprintf(cepa_type,"CEPA(3)");
  else if (cepa_level == -1) sprintf(cepa_type,"CISD");
  else if (cepa_level == -2) sprintf(cepa_type,"ACPF");
  else if (cepa_level == -3) sprintf(cepa_type,"AQCC");

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

  // sort integrals and write them to disk
  times(&total_tmstime);
  time_start = time(NULL);
  user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;

  OutOfCoreSort(nfzc,nfzv,nmotemp,ndoccact,nvirt,wfn_->isCIM());

  times(&total_tmstime);
  time_stop = time(NULL);
  user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;

  fprintf(outfile,"  Time for integral sort:           %6.2lf s (user)\n",user_stop-user_start);
  fprintf(outfile,"                                    %6.2lf s (system)\n",sys_stop-sys_start);
  fprintf(outfile,"                                    %6d s (total)\n",(int)time_stop-(int)time_start);

  // orbital energies
  eps = (double*)malloc(nmo*sizeof(double));
  int count=0;
  for (int h=0; h<nirreps; h++){
      eps_test = wfn_->epsilon_a();
      for (int norb = fzc[h]; norb<docc[h]; norb++){
          eps[count++] = eps_test->get(h,norb);
      }
  }
  for (int h=0; h<nirreps; h++){
      eps_test = wfn_->epsilon_a();
      for (int norb = docc[h]; norb<orbs[h]-fzv[h]; norb++){
          eps[count++] = eps_test->get(h,norb);
      }
  }
  eps_test.reset();

  // by default, t2 will be held in core
  t2_on_disk = false;

}

/*===================================================================

  solve cepa equations

===================================================================*/
PsiReturnType CoupledPair::CEPAIterations(Options&options){

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
  ecepa=0.0;

  fprintf(outfile,"\n");
  fprintf(outfile,
    "  Begin %s iterations\n\n",cepa_type);
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

  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_NEW);
     psio->write_entry(PSIF_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }
  pair_energy = (double*)malloc(o*o*sizeof(double));

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
         for (int i=0; i<ncepatasks; i++) {
             (*this.*CepaTasklist[i].func)(CepaParams[i]);
         }
      }

      // update the amplitudes and check the energy
      Eold = ecepa;
      PairEnergy();
      if (!options.get_bool("CEPA_NO_SINGLES")){
         UpdateT1(iter);
      }
      UpdateT2(iter);

      // add vector to list for diis
      DIISOldVector(iter,diis_iter,replace_diis_iter);

      // diis error vector and convergence check
      nrm = DIISErrorVector(diis_iter,replace_diis_iter,iter);

      // diis extrapolation
      if (diis_iter>1){
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,arraysize+o*v);
         else                   DIIS(diisvec,maxdiis,arraysize+o*v);
         DIISNewAmplitudes(diis_iter);
         // if cepa_no_singles, zero t1 just to be safe after extrapolation
         if (options.get_bool("CEPA_NO_SINGLES")){
            memset((void*)t1,'\0',o*v*sizeof(double));
         }
      }
      ecepa = CheckEnergy();

      if (diis_iter<=maxdiis) diis_iter++;
      else if (replace_diis_iter<maxdiis) replace_diis_iter++;
      else replace_diis_iter = 1;

      iter_stop = time(NULL);
      fprintf(outfile,"  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,ecepa,ecepa-Eold,nrm,(int)iter_stop-(int)iter_start);
      fflush(outfile);
      iter++;
      if (iter==1) emp2 = ecepa;
      if (iter==1) SCS_MP2();
  }
  times(&total_tmstime);
  time_stop = time(NULL);
  user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;
  psio.reset();

  if (iter==maxiter){
     throw PsiException("  CEPA iterations did not converge.",__FILE__,__LINE__);
  }

  // is this a cim-cepa computation?
  if (wfn_->isCIM()){
     Local_SCS_CEPA();
     ecepa = ecepa_os/ecepa_os_fac + ecepa_ss/ecepa_ss_fac;
  }
  else{
     SCS_CEPA();
  }

  fprintf(outfile,"\n");
  fprintf(outfile,"  %s iterations converged!\n",cepa_type);
  fprintf(outfile,"\n");
  fprintf(outfile,"        OS SCS-MP2 correlation energy:     %20.12lf\n",emp2_os);
  fprintf(outfile,"        SS SCS-MP2 correlation energy:     %20.12lf\n",emp2_ss);
  fprintf(outfile,"        SCS-MP2 correlation energy:        %20.12lf\n",emp2_os+emp2_ss);
  fprintf(outfile,"      * SCS-MP2 total energy:              %20.12lf\n",emp2_os+emp2_ss+escf);
  fprintf(outfile,"\n");
  fprintf(outfile,"        OS MP2 correlation energy:         %20.12lf\n",emp2_os/emp2_os_fac);
  fprintf(outfile,"        SS MP2 correlation energy:         %20.12lf\n",emp2_ss/emp2_ss_fac);
  fprintf(outfile,"        MP2 correlation energy:            %20.12lf\n",emp2);
  fprintf(outfile,"      * MP2 total energy:                  %20.12lf\n",emp2+escf);
  fprintf(outfile,"\n");
  if (cepa_level>=0){
     if (options.get_bool("SCS_CEPA")){
        fprintf(outfile,"        OS SCS-%s correlation energy: %20.12lf\n",cepa_type,ecepa_os);
        fprintf(outfile,"        SS SCS-%s correlation energy: %20.12lf\n",cepa_type,ecepa_ss);
        fprintf(outfile,"        SCS-%s correlation energy:    %20.12lf\n",cepa_type,ecepa_os+ecepa_ss);
        fprintf(outfile,"      * SCS-%s total energy:          %20.12lf\n",cepa_type,ecepa_os+ecepa_ss+escf);
        fprintf(outfile,"\n");
     }
     fprintf(outfile,"        OS %s correlation energy:     %20.12lf\n",cepa_type,ecepa_os/ecepa_os_fac);
     fprintf(outfile,"        SS %s correlation energy:     %20.12lf\n",cepa_type,ecepa_ss/ecepa_ss_fac);
     fprintf(outfile,"        %s correlation energy:        %20.12lf\n",cepa_type,ecepa);
     fprintf(outfile,"      * %s total energy:              %20.12lf\n",cepa_type,ecepa+escf);
  }else{
     if (options.get_bool("SCS_CEPA")){
        fprintf(outfile,"        OS SCS-%s correlation energy:    %20.12lf\n",cepa_type,ecepa_os);
        fprintf(outfile,"        SS SCS-%s correlation energy:    %20.12lf\n",cepa_type,ecepa_ss);
        fprintf(outfile,"        SCS-%s correlation energy:       %20.12lf\n",cepa_type,ecepa_os+ecepa_ss);
        fprintf(outfile,"      * SCS-%s total energy:             %20.12lf\n",cepa_type,ecepa_os+ecepa_ss+escf);
        fprintf(outfile,"\n");
     }
     fprintf(outfile,"        OS %s correlation energy:        %20.12lf\n",cepa_type,ecepa_os/ecepa_os_fac);
     fprintf(outfile,"        SS %s correlation energy:        %20.12lf\n",cepa_type,ecepa_ss/ecepa_ss_fac);
     fprintf(outfile,"        %s correlation energy:           %20.12lf\n",cepa_type,ecepa);
     fprintf(outfile,"      * %s total energy:                 %20.12lf\n",cepa_type,ecepa+escf);
  }
  fprintf(outfile,"\n");
  fprintf(outfile,"  Total time for %s iterations: %10.2lf s (user)\n",cepa_type,user_stop-user_start);
  fprintf(outfile,"                                  %10.2lf s (system)\n",sys_stop-sys_start);
  fprintf(outfile,"                                  %10d s (total)\n",(int)time_stop-(int)time_start);
  fprintf(outfile,"\n");
  fprintf(outfile,"  Time per iteration:             %10.2lf s (user)\n",(user_stop-user_start)/(iter-1));
  fprintf(outfile,"                                  %10.2lf s (system)\n",(sys_stop-sys_start)/(iter-1));
  fprintf(outfile,"                                  %10.2lf s (total)\n",((double)time_stop-(double)time_start)/(iter-1));
  fflush(outfile);

  return Success;
}

/*===================================================================

  determine tiling for vabcd and vabci diagrams for the cpu
  this determines the size of blocks of integrals that 
  can be read into cpu memory.

===================================================================*/
void CoupledPair::DefineTilingCPU(){
  long int i,v = nvirt;
  long int o = ndoccact;
  long int ov2 = o*v*v;
  long int ov = o*v;
  long int o2 = o*o;

  // number of doubles in total memory
  long int ndoubles = memory/8L;
  // minus storage for other necessary buffers 
  ndoubles -= o*o*v*v+2L*(o*o*v*v+o*v)+2L*o*v+2*v*v+(o+v);
  if (t2_on_disk){
     ndoubles += o*o*v*v;
  }else{
     fprintf(outfile,"\n");
     fprintf(outfile,"  Define tiling:\n");
     fprintf(outfile,"\n");
  }

  // if not enough space, check to see if keeping t2 on disk will help
  if (ndoubles<o*o*v*v){
     if (t2_on_disk)
        throw PsiException("out of memory: no amount of tiling can fix this!",__FILE__,__LINE__);
     else{
        tilesize = ov2tilesize = ovtilesize = 0;
        return;
     }
  }

  ntiles = -999L;
  tilesize = v*(v+1L)/2L;
  ntiles = 1L;

  // tiling for vabcd diagram
  ntiles=1L;
  tilesize=v*(v+1L)/2L/1L;
  if (ntiles*tilesize<v*(v+1L)/2L) tilesize++;
  while(v*(v+1L)/2L*tilesize>ndoubles){
     ntiles++;
     tilesize = v*(v+1L)/2L/ntiles;
     if (ntiles*tilesize<v*(v+1L)/2L) tilesize++;
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
void CoupledPair::AllocateMemory(Options&options){

  long int i,o=ndoccact;
  long int v=nvirt;
  long int dim;

  fprintf(outfile,"\n");
  fprintf(outfile,"  available memory =                        %9.2lf mb\n",Process::environment.get_memory()/1024./1024.);
  fprintf(outfile,"  minimum memory requirements for %s =    %9.2lf mb\n",cepa_type,
         8./1024./1024.*(o*o*v*v+2.*(o*o*v*v+o*v)+2.*o*v+2.*v*v+o+v));

  // define tiling for v^4 and ov^3 diagrams according to how much memory is available
  DefineTilingCPU();

  dim = 0;
  if (tilesize*v*(v+1)/2 > dim) dim = tilesize*v*(v+1)/2;
  if (ovtilesize*v*v > dim)     dim = ovtilesize*v*v;
  if (ov2tilesize*v > dim)      dim = ov2tilesize*v;

  // if integrals buffer isn't at least o^2v^2, try tiling again assuming t2 is on disk.
  if (dim<o*o*v*v){
     fprintf(outfile,"\n");
     fprintf(outfile,"  Warning: cannot accomodate T2 in core. T2 will be stored on disk.\n");
     fprintf(outfile,"\n");
     fflush(outfile);
     t2_on_disk = true;
     DefineTilingCPU();
     dim = 0;
     if (tilesize*v*(v+1)/2 > dim) dim = tilesize*v*(v+1)/2;
     if (ovtilesize*v*v > dim)     dim = ovtilesize*v*v;
     if (ov2tilesize*v > dim)      dim = ov2tilesize*v;

     if (dim<o*o*v*v){
        throw PsiException("out of memory: general buffer cannot accomodate T2",__FILE__,__LINE__);
     }

     fprintf(outfile,"\n");
     fprintf(outfile,"  Increase memory by %7.2lf mb to hold T2 in core.\n",o*o*v*v*8L/1024./1024.);
     fprintf(outfile,"\n");
  }

  maxelem = dim;

  double total_memory = 1.*dim+2.*(o*o*v*v+o*v)+1.*o*o*v*v+2.*o*v+2.*v*v;
  if (t2_on_disk) total_memory = 1.*dim+2.*(o*o*v*v+o*v)+2.*o*v+2.*v*v;
  total_memory *= 8./1024./1024.;

  fprintf(outfile,"\n");
  fprintf(outfile,"  Allocate cpu memory (%9.2lf mb).....",total_memory);
  integrals = (double*)malloc(dim*sizeof(double));
  tempt     = (double*)malloc((o*o*v*v+o*v)*sizeof(double));
  tempv     = (double*)malloc((o*o*v*v+o*v)*sizeof(double));

  if (!t2_on_disk) tb = (double*)malloc(o*o*v*v*sizeof(double));
  w1        = (double*)malloc(o*v*sizeof(double));
  t1        = (double*)malloc(o*v*sizeof(double));
  I1        = (double*)malloc(v*v*sizeof(double));
  I1p       = (double*)malloc(v*v*sizeof(double));
  fprintf(outfile,"done.\n");

  fprintf(outfile,"  Initialize cpu memory..................");
  memset((void*)integrals,'\0',dim*sizeof(double));
  memset((void*)tempv,'\0',(o*o*v*v+o*v)*sizeof(double));
  memset((void*)tempt,'\0',(o*o*v*v+o*v)*sizeof(double));
  if (!t2_on_disk) memset((void*)tb,'\0',o*o*v*v*sizeof(double));
  memset((void*)w1,'\0',o*v*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));
  memset((void*)I1,'\0',v*v*sizeof(double));
  memset((void*)I1p,'\0',v*v*sizeof(double));
  fprintf(outfile,"done.\n");

  // DIIS:
  diisvec    = (double*)malloc(sizeof(double)*(maxdiis+1));
  memset((void*)diisvec,'\0',(maxdiis+1)*sizeof(double));
}

void CoupledPair::CPU_t1_vmeai(CepaTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int i,a,m,e,id,one=1;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_AKJC2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_AKJC2,"E2akjc2",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_AKJC2,1);

  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  F_DAXPY(o*o*v*v,-2.0,integrals,1,tempv,1);
  
  for (i=0; i<o; i++){
      F_DCOPY(v,t1+i,o,tempt+i*v,1);
  }
  F_DGEMV('n',o*v,o*v,-1.0,tempv,o*v,tempt,1,0.0,integrals,1);
  for (a=0; a<v; a++){
      F_DAXPY(o,1.0,integrals+a,v,w1+a*o,1);
  }

  psio.reset();
}

void CoupledPair::CPU_t1_vmeni(CepaTaskParams params){
  long int m,e,n,a,id;
  long int o=ndoccact;
  long int v=nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = tempv;
  }

  for (a=0,id=0; a<v; a++){
      for (m=0; m<o; m++){
          for (n=0; n<o; n++){
              for (e=0; e<v; e++){
                  tempt[id++] = 2.*tb[e*v*o*o+a*o*o+m*o+n]-tb[a*v*o*o+e*o*o+m*o+n];
              }
          }
      }
  }
  psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJAK,"E2ijak",(char*)&tempv[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_IJAK,1);
  F_DGEMM('t','n',o,v,o*o*v,-1.0,tempv,o*o*v,tempt,o*o*v,1.0,w1,o);
  psio.reset();
}

void CoupledPair::CPU_t1_vmaef(CepaTaskParams params){
  long int m,e,i,f,a,id;
  long int o=ndoccact;
  long int v=nvirt;

  boost::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = tempv;
  }

  for (f=0,id=0; f<v; f++){
      for (m=0; m<o; m++){
          for (e=0; e<v; e++){
              for (i=0; i<o; i++){
                  tempt[id++] = 2.*tb[e*v*o*o+f*o*o+m*o+i]-tb[e*v*o*o+f*o*o+i*o+m];
              }
          }
      }
  }

  long int tilesize,lasttile,ntiles=1;
  long int ov2 = o*v*v;
  // tile v in chunks of o

  ntiles=1L;
  tilesize=v/1L;
  if (ntiles*tilesize<v) tilesize++;
  while(tilesize*ov2>maxelem){
     ntiles++;
     tilesize = v/ntiles;
     if (ntiles*tilesize<ov2) tilesize++;
  }
  lasttile = v - (ntiles-1L)*tilesize;

  psio->open(PSIF_ABCI3,PSIO_OPEN_OLD);
  psio_address addr;
  addr = PSIO_ZERO;

  for (i=0; i<ntiles-1; i++){
      psio->read(PSIF_ABCI3,"E2abci3",(char*)&integrals[0],tilesize*ov2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,tilesize,ov2,1.0,tempt,o,integrals,ov2,1.0,w1+i*tilesize*o,o);
  }
  i=ntiles-1;
  psio->read(PSIF_ABCI3,"E2abci3",(char*)&integrals[0],lasttile*ov2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o,lasttile,ov2,1.0,tempt,o,integrals,ov2,1.0,w1+i*tilesize*o,o);
  psio->close(PSIF_ABCI3,1);
  psio.reset();

}

// a refactored version of I2p(ab,ci) that avoids ov^3 storage
void CoupledPair::CPU_I2p_abci_refactored_term1(CepaTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int a,b,c,i,j,id=0;
  long int ov2 = o*v*v;
  long int o2v = o*o*v;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_ABCI5,PSIO_OPEN_OLD);
  psio_address addr;
  addr = PSIO_ZERO;

  for (i=0; i<nov2tiles-1; i++){
      psio->read(PSIF_ABCI5,"E2abci5",(char*)&integrals[0],v*ov2tilesize*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,ov2tilesize,v,1.0,t1,o,integrals,v,0.0,tempt+i*ov2tilesize*o,o);
  }
  i=nov2tiles-1;
  psio->read(PSIF_ABCI5,"E2abci5",(char*)&integrals[0],v*lastov2tile*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o,lastov2tile,v,1.0,t1,o,integrals,v,0.0,tempt+i*ov2tilesize*o,o);
  psio->close(PSIF_ABCI5,1);

  // contribute to residual
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

/*================================================================

   update amplitudes

================================================================*/
void CoupledPair::UpdateT1(long int iter){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int id=0;
  double tnew,dia,energy;

  double fac = 1.0;
  if (cepa_level == 0)       fac = 0.0;
  else if (cepa_level == -1) fac = 1.0;
  else if (cepa_level == -2) fac = 1.0/o;
  else if (cepa_level == -3) fac = 1.0-(2.0*o-2.0)*(2.0*o-3.0) / (2.0*o*(2.0*o-1.0));
  energy = ecepa * fac;

  if (iter<1){
     memset((void*)t1,'\0',o*v*sizeof(double));
     memset((void*)w1,'\0',o*v*sizeof(double));
  }
  else{
     for (i=0; i<o; i++){

         if (cepa_level == 1){
            energy = 0.0;
            for (long int k=0; k<o; k++){
                energy += (pair_energy[i*o+k]);
            }
         }else if (cepa_level == 2){
            energy = pair_energy[i*o+i];
         }else if (cepa_level == 3){
            energy = -pair_energy[i*o+i];
            for (long int k=0; k<o; k++){
                energy += 2.0*(pair_energy[i*o+k]);
            }
         }

         for (a=o; a<rs; a++){
             dia = -eps[i]+eps[a];
             tnew = - (w1[(a-o)*o+i])/(dia-energy);
             w1[(a-o)*o+i] = tnew;
         }
     }
  }
  // error vector for diis is in tempv:
  F_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  F_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  F_DCOPY(o*v,w1,1,t1,1);
}
void CoupledPair::Local_SCS_CEPA(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,jaib,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);

  SharedMatrix Rii = wfn_->CIMTransformationMatrix();
  double**Rii_pointer = Rii->pointer();
  // transform E2klcd back from quasi-canonical basis
  for (i=0; i<o; i++){
      for (a=0; a<v; a++){
          for (j=0; j<o; j++){
              for (b=0; b<v; b++){
                  double dum = 0.0;
                  for (int ip=0; ip<o; ip++){
                      dum += tempt[ip*o*v*v+a*o*v+j*v+b]*Rii_pointer[ip][i];
                  }
                  integrals[i*o*v*v+a*o*v+j*v+b] = dum;
              }
          }
      }
  }


  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = tempv;
  }

  // transform t2 back from quasi-canonical basis
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tb[a*o*o*v+b*o*o+i*o+j] += t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  double dum = 0.0;
                  for (int ip=0; ip<o; ip++){
                      dum += tb[a*o*o*v+b*o*o+ip*o+j]*Rii_pointer[ip][i];
                  }
                  tempt[a*o*o*v+b*o*o+i*o+j] = dum;
              }
          }
      }
  }
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tb[a*o*o*v+b*o*o+i*o+j] -= t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }


  SharedVector factor = wfn_->CIMOrbitalFactors();
  double*factor_pointer = factor->pointer();
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);
                  osenergy += integrals[iajb]*(tempt[ijab])*factor_pointer[i];
                  ssenergy += integrals[iajb]*(tempt[ijab]-tempt[(b-o)*o*o*v+(a-o)*o*o+i*o+j])*factor_pointer[i];
                  ijab++;
              }
          }
      }
  }
  ecepa_os = ecepa_os_fac * osenergy;
  ecepa_ss = ecepa_ss_fac * ssenergy;

  psio.reset();
}
void CoupledPair::SCS_CEPA(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,jaib,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = tempv;
  }
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);
                  osenergy += integrals[iajb]*(tb[ijab]);
                  ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                  ijab++;
              }
          }
      }
  }
  ecepa_os = ecepa_os_fac * osenergy;
  ecepa_ss = ecepa_ss_fac * ssenergy;

  psio.reset();
}
void CoupledPair::SCS_MP2(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,jaib,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = tempv;
  }
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
  emp2_os = emp2_os_fac * osenergy;
  emp2_ss = emp2_ss_fac * ssenergy;

  psio.reset();
}
void CoupledPair::PairEnergy(){

  if (cepa_level<1) return;

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  double ta,tnew,dijab,da,dab,dabi;
  long int iajb,jaib,ijab=0;
  double energy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = tempt;
  }
  for (i=0; i<o; i++){
      for (j=0; j<o; j++){
          energy=0.0;
          for (a=o; a<rs; a++){
              for (b=o; b<rs; b++){

                  ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  energy += integrals[iajb]*(2.0*tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
              }
          }
          pair_energy[i*o+j] = energy;
      }
  }

  psio.reset();
}
double CoupledPair::CheckEnergy(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  double ta,tnew,dijab,da,dab,dabi;
  long int iajb,jaib,ijab=0;
  double energy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = tempv;
  }
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);
                  energy += (2.*integrals[iajb]-integrals[jaib])*(tb[ijab]);
                  ijab++;
              }
          }
      }
  }

  psio.reset();

  return energy;
}
void CoupledPair::UpdateT2(long int iter){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  double ta,tnew,dijab,di,dij,dija,energy;
  long int iajb,jaib,ijab=0;
  //double energy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  // we still have the residual in memory in tempv
  //psio->open(PSIF_R2,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));

  double fac = 1.0;
  if (cepa_level == 0)       fac = 0.0;
  else if (cepa_level == -1) fac = 1.0;
  else if (cepa_level == -2) fac = 1.0/o;
  else if (cepa_level == -3) fac = 1.0-(2.0*o-2.0)*(2.0*o-3.0) / (2.0*o*(2.0*o-1.0));
  energy = ecepa * fac;

  for (i=0; i<o; i++){
      di = - eps[i];
      for (j=0; j<o; j++){
          dij = di-eps[j];


          if (cepa_level == 1){
             energy = 0.0;
             for (long int k=0; k<o; k++){
                 energy += 0.5*(pair_energy[i*o+k]+pair_energy[j*o+k]);
             }
          }else if (cepa_level == 2){
             energy = pair_energy[i*o+j];
          }else if (cepa_level == 3){
             energy = -pair_energy[i*o+j];
             for (long int k=0; k<o; k++){
                 energy += (pair_energy[i*o+k]+pair_energy[j*o+k]);
             }
          }
          


          for (a=o; a<rs; a++){
              dija = dij + eps[a];
              for (b=o; b<rs; b++){
                  dijab = dija + eps[b];


                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);


                  ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                  tnew = - (integrals[iajb] + tempv[ijab])/(dijab-energy);
                  tempt[ijab] = tnew;
              }
          }
      }
  }

  // error vectors for diis are in tempv:
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempv,1);
  }
  F_DAXPY(o*o*v*v,-1.0,tempt,1,tempv,1);
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->write_entry(PSIF_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tempt,1,tb,1);
  }

  psio.reset();
}

/*================================================================

   diis functions

================================================================*/
void CoupledPair::DIIS(double*c,long int nvec,long int n){
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
void CoupledPair::DIISOldVector(long int iter,int diis_iter,int replace_diis_iter){
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

  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = integrals;
  }

  psio->write(PSIF_OVEC,oldvector,(char*)&tb[0],arraysize*sizeof(double),addr,&addr);
  psio->write(PSIF_OVEC,oldvector,(char*)&t1[0],o*v*sizeof(double),addr,&addr);
  psio->close(PSIF_OVEC,1);
  psio.reset();

  free(oldvector);
}
double CoupledPair::DIISErrorVector(int diis_iter,int replace_diis_iter,int iter){
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
void CoupledPair::DIISNewAmplitudes(int diis_iter){
  long int o = ndoccact;
  long int arraysize,v = nvirt;
  arraysize=o*o*v*v;

  char*oldvector;
  oldvector=(char*)malloc(1000*sizeof(char));

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_OVEC,PSIO_OPEN_OLD);

  psio_address addr;

  if (t2_on_disk){
     tb = integrals;
  }

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

  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_NEW);
     psio->write_entry(PSIF_T2,"t2",(char*)&tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }

  psio.reset();
}

/**
 *  Build and use I2ijkl
 */
void CoupledPair::I2ijkl(CepaTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempt,1);
  }

  psio->open(PSIF_IJKL,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJKL,"E2ijkl",(char*)&integrals[0],o*o*o*o*sizeof(double));
  psio->close(PSIF_IJKL,1);

  F_DGEMM('n','n',o*o,v*v,o*o,0.5,integrals,o*o,tempt,o*o,0.0,tempv,o*o);

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
void CoupledPair::I2piajk(CepaTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  psio->open(PSIF_IJAK2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJAK2,"E2ijak2",(char*)&tempv[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_IJAK2,1);

  F_DGEMM('n','n',o*o*v,v,o,-1.0,tempv,o*o*v,t1,o,0.0,tempt,o*o*v);

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
void CoupledPair::Vabcd1(CepaTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempt,1);
  }
  for (i=0; i<o; i++){
      for (j=i; j<o; j++){
          for (a=0; a<v; a++){
              for (b=a+1; b<v; b++){
                  tempv[Position(a,b)*o*(o+1)/2+Position(i,j)] =
                     tempt[a*o*o*v+b*o*o+i*o+j]+tempt[b*o*o*v+a*o*o+i*o+j];
              }
              tempv[Position(a,a)*o*(o+1)/2+Position(i,j)] =
                 tempt[a*o*o*v+a*o*o+i*o+j];
          }
      }
  }
  psio->open(PSIF_ABCD1,PSIO_OPEN_OLD);
  addr = PSIO_ZERO;
  for (j=0; j<ntiles-1; j++){
      psio->read(PSIF_ABCD1,"E2abcd1",(char*)&integrals[0],tilesize*v*(v+1)/2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o*(o+1)/2,tilesize,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  }
  j=ntiles-1;
  psio->read(PSIF_ABCD1,"E2abcd1",(char*)&integrals[0],lasttile*v*(v+1)/2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o*(o+1)/2,lasttile,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  psio->close(PSIF_ABCD1,1);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
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
  psio.reset();

}

/**
 *  Use Vabcd2
 */
void CoupledPair::Vabcd2(CepaTaskParams params){
  long int id,i,j,a,b,o,v;
  int sg,sg2;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempt,1);
  }
  for (i=0; i<o; i++){
      for (j=i; j<o; j++){
          for (a=0; a<v; a++){
              for (b=a; b<v; b++){
                  tempv[Position(a,b)*o*(o+1)/2+Position(i,j)] =
                    tempt[a*o*o*v+b*o*o+i*o+j]-tempt[b*o*o*v+a*o*o+i*o+j];
              }
          }
      }
  }
  psio->open(PSIF_ABCD2,PSIO_OPEN_OLD);
  addr = PSIO_ZERO;
  for (j=0; j<ntiles-1; j++){
      psio->read(PSIF_ABCD2,"E2abcd2",(char*)&integrals[0],tilesize*v*(v+1)/2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o*(o+1)/2,tilesize,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  }
  j = ntiles-1;
  psio->read(PSIF_ABCD2,"E2abcd2",(char*)&integrals[0],lasttile*v*(v+1)/2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o*(o+1)/2,lasttile,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  psio->close(PSIF_ABCD2,1);

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
void CoupledPair::I2iabj(CepaTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  F_DCOPY(o*o*v*v,integrals,1,tempv,1);

  // use I2iabj
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = integrals;
  }
  for (j=0,id=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (a=0; a<v; a++){
                  tempt[id++] = 2*tb[a*o*o*v+b*o*o+i*o+j]-tb[b*o*o*v+a*o*o+i*o+j];
              }
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,1.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);

  // contribute to residual
  psio->open(PSIF_R2,PSIO_OPEN_OLD);
  // if we KNOW this is the first diagram, we don't need to read in the old
  // residual.
  //psio->read_entry(PSIF_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
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
void CoupledPair::I2iajb(CepaTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  psio->open(PSIF_AKJC2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_AKJC2,"E2akjc2",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_AKJC2,1);

  // use I2iajb
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = tempv;
  }

  for (j=0,id=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (a=0; a<v; a++){
                  integrals[id++] = tb[b*v*o*o+a*o*o+j*o+i];
              }
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,integrals,o*v,0.0,tempv,o*v);

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
  if (t2_on_disk){
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
     tb = integrals;
  }

  for (j=0,id=0; j<o; j++){
      for (a=0; a<v; a++){
          for (i=0; i<o; i++){
              for (b=0; b<v; b++){
                  tempv[id++] = tb[b*v*o*o+a*o*o+j*o+i];
              }
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,tempv,o*v,0.0,integrals,o*v);

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
void CoupledPair::DefineTasks(){
  CepaTasklist = new CepaTask[1000];
  CepaParams   = new CepaTaskParams[1000];
  long int o = ndoccact;
  long int v = nvirt;

  ncepatasks=0;

  CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::I2iabj;
  CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::I2iajb;
  CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::I2ijkl;
  CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::I2piajk;
  // these diagrams are necessary only if we have single excitations
  if (!options_.get_bool("CEPA_NO_SINGLES")){
     CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::CPU_t1_vmeni;
     CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::CPU_t1_vmaef;
     CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::CPU_I2p_abci_refactored_term1;
     CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::CPU_t1_vmeai;
  }
  CepaTasklist[ncepatasks++].func  = &psi::CoupledPair::Vabcd1;

  // this is the last diagram that contributes to doubles residual,
  // so we can keep it in memory rather than writing and rereading
  CepaTasklist[ncepatasks++].func        = &psi::CoupledPair::Vabcd2;
}


} // end of namespace psi
