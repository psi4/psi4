/*
 *@BEGIN LICENSE
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

#define PSIF_CIM 273 // TODO: move to psifiles.h
#include"psi4-dec.h"
#include<libmints/vector.h>
#include<libmints/matrix.h>
#include<libmints/wavefunction.h>
#include<libqt/qt.h>
#include<sys/times.h>
#include<libciomr/libciomr.h>
#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
    #define omp_get_max_threads() 1
#endif

#include"blas.h"
#include"ccsd.h"
#include<libmints/basisset.h>
#include<libmints/basisset_parser.h>
#include<lib3index/3index.h>

using namespace psi;
using namespace boost;

// position in a symmetric packed matrix
long int Position(long int i,long int j){
  if (i<j){
    return ((j*(j+1))>>1)+i;
  }
  return ((i*(i+1))>>1)+j;
}

namespace psi{ namespace fnocc{

// diagrams for mp3 and mp4
void DefineLinearTasks();
void DefineQuadraticTasks();

// sort
void SortIntegrals(int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,Options&options,bool iscim);
void Sort_OV3_LowMemory(long int memory,long int o,long int v,bool islocal);

// coupled cluster constructor
CoupledCluster::CoupledCluster(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options):
        Wavefunction(options, _default_psio_lib_)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

CoupledCluster::~CoupledCluster()
{
}

/**
  * initialize.  set variables and options_.
  */
void CoupledCluster::common_init() {

  mp2_only = options_.get_bool("RUN_MP2");
  mp4_only = options_.get_bool("RUN_MP4");
  mp3_only = options_.get_bool("RUN_MP3");
  isccsd   = options_.get_bool("RUN_CCSD");

  escf    = reference_wavefunction_->reference_energy();
  doccpi_ = reference_wavefunction_->doccpi();
  soccpi_ = reference_wavefunction_->soccpi();
  frzcpi_ = reference_wavefunction_->frzcpi();
  frzvpi_ = reference_wavefunction_->frzvpi();
  nmopi_  = reference_wavefunction_->nmopi();

  Da_ = SharedMatrix(reference_wavefunction_->Da());
  Ca_ = SharedMatrix(reference_wavefunction_->Ca());
  Fa_ = SharedMatrix(reference_wavefunction_->Fa());
  epsilon_a_= boost::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
  epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
  nalpha_ = reference_wavefunction_->nalpha();
  nbeta_  = reference_wavefunction_->nbeta();

  nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
  for (int h=0; h<nirrep_; h++){
      nfzc   += frzcpi_[h];
      nfzv   += frzvpi_[h];
      nso    += nsopi_[h];
      nmo    += nmopi_[h]-frzcpi_[h]-frzvpi_[h];
      ndocc  += doccpi_[h];
  }
  if (reference_wavefunction_->isCIM()){
     ndoccact = reference_wavefunction_->CIMActiveOccupied();
     nvirt    = reference_wavefunction_->CIMActiveVirtual();
     nfzc     = ndocc - ndoccact;
     nmo      = ndoccact + nvirt;
     nfzv     = nmopi_[0] - ndocc - nvirt;
  }else{
     ndoccact = ndocc - nfzc;
     nvirt    = nmo - ndoccact;
  }

  // for triples, we use nvirt_no in case we've truncated the virtual space:
  nvirt_no = nvirt;

  // get paramters from input 
  e_conv   = options_.get_double("E_CONVERGENCE");
  r_conv   = options_.get_double("R_CONVERGENCE");
  maxiter = options_.get_int("MAXITER");
  maxdiis = options_.get_int("DIIS_MAX_VECS");

  // memory is from process::environment
  memory = Process::environment.get_memory();

  // SCS MP2 and CCSD
  emp2_os_fac = options_.get_double("MP2_OS_SCALE");
  emp2_ss_fac = options_.get_double("MP2_SS_SCALE");
  eccsd_os_fac = options_.get_double("CC_OS_SCALE");
  eccsd_ss_fac = options_.get_double("CC_SS_SCALE");

  // quit if number of virtuals is less than number of doubly occupied
  if (nvirt<ndoccact){
     throw PsiException("ndocc must be less than nvirt",__FILE__,__LINE__);
  }

  // by default, use fast (t)
  isLowMemory = false;

  // by default, t2 will be held in core
  t2_on_disk = false;

  // for df-bccd
  brueckner_iter = 0;
}

void CoupledCluster::finalize() {
  if (!mp2_only) {
      for (int i = 0; i < ncctasks; i++){
          free(CCTasklist[i].name);
      }
  }
  // there is something weird with chkpt_ ... reset it
  chkpt_.reset();
}

double CoupledCluster::compute_energy() {
  PsiReturnType status = Success;

  if (options_.get_bool("RUN_MP2")) {
      tstart();
      WriteBanner();
      AllocateMemory();
      MP2();
      Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
      Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;
      Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2;
      Process::environment.globals["MP2 TOTAL ENERGY"] = emp2 + escf;
      Process::environment.globals["CURRENT ENERGY"] = emp2 + escf;
      Process::environment.globals["CURRENT CORRELATION ENERGY"] = emp2;
      tstop();
      return emp2 + escf;
  }

  // integral sort
  tstart();
  SortIntegrals(nfzc,nfzv,nmo+nfzc+nfzv,ndoccact,nvirt,options_,reference_wavefunction_->isCIM());
  tstop();

  // MP4(SDQ)
  tstart();
  WriteBanner();
  AllocateMemory();
  MP4_SDQ();

  // CCSD or QCISD
  if (!options_.get_bool("RUN_MP3") && !options_.get_bool("RUN_MP4")) {
     status = CCSDIterations();

     // ccsd energy
     if (isccsd) {
        Process::environment.globals["CCSD CORRELATION ENERGY"] = eccsd;
        Process::environment.globals["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
        Process::environment.globals["CCSD SAME-SPIN CORRELATION ENERGY"] = eccsd_ss;
        Process::environment.globals["CCSD TOTAL ENERGY"] = eccsd + escf;
     }else{
        Process::environment.globals["QCISD CORRELATION ENERGY"] = eccsd;
        Process::environment.globals["QCISD OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
        Process::environment.globals["QCISD SAME-SPIN CORRELATION ENERGY"] = eccsd_ss;
        Process::environment.globals["QCISD TOTAL ENERGY"] = eccsd + escf;
     }
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = eccsd;
     Process::environment.globals["CURRENT ENERGY"] = eccsd + escf;
  }else if (options_.get_bool("RUN_MP3")){
     Process::environment.globals["CURRENT ENERGY"] = emp2 + emp3 + escf;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = emp2 + emp3;
  }else{
     Process::environment.globals["CURRENT ENERGY"] = emp2 + emp3 + emp4_sd + emp4_q + escf;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = emp2 + emp3 + emp4_sd + emp4_q;
  }
  tstop();

  // mp2 energy
  Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
  Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;
  Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2;
  Process::environment.globals["MP2 TOTAL ENERGY"] = emp2 + escf;

  if ( !options_.get_bool("RUN_MP2") ) {
      // mp3 energy
      Process::environment.globals["MP3 CORRELATION ENERGY"] = emp2 + emp3;
      Process::environment.globals["MP3 TOTAL ENERGY"] = emp2 + emp3 + escf;

      // mp2.5 energy
      Process::environment.globals["MP2.5 CORRELATION ENERGY"] = emp2 + 0.5*emp3 ;
      Process::environment.globals["MP2.5 TOTAL ENERGY"] = emp2 + 0.5*emp3 + escf;

      // mp4 energy
      if ( !options_.get_bool("RUN_MP3") ) {
          Process::environment.globals["MP4(SDQ) TOTAL ENERGY"] = emp2 + emp3 + emp4_sd + emp4_q + escf;
          Process::environment.globals["MP4(SDQ) CORRELATION ENERGY"] = emp2 + emp3 + emp4_sd + emp4_q;
      }

  }

  // free some memory before triples 
  free(integrals);
  free(w1);
  free(I1);
  free(I1p);
  free(diisvec);
  free(tempt);
  free(tempv);

  if (options_.get_bool("COMPUTE_TRIPLES")){

     if (isccsd) ccmethod = 0;
     else        ccmethod = 1;

     long int o = ndoccact;
     long int v = nvirt;

     // if low memory or CIM, need to generate one last set of integrals
     if (reference_wavefunction_->isCIM() || isLowMemory){
        Sort_OV3_LowMemory(memory - 8L*(long int)(!t2_on_disk)*o*o*v*v,o,v,reference_wavefunction_->isCIM());
     }

     // now there should be space for t2
     if (t2_on_disk){
        tb = (double*)malloc(o*o*v*v*sizeof(double));
     }

     bool do_cc = !options_.get_bool("RUN_MP4") && !options_.get_bool("RUN_MP3");
     bool do_mp = options_.get_bool("COMPUTE_MP4_TRIPLES");

     tstart();
     // triples
     if (reference_wavefunction_->isCIM()){
        if (do_cc) {
           status = local_triples();
        }
        if (do_mp) {
           ccmethod = 2;
           status = local_triples();
        }
     }
     else{
        if (isLowMemory){
           if (do_cc) {
              status = lowmemory_triples();
           }
           if (do_mp) {
              ccmethod = 2;
              status = lowmemory_triples();
           }
        }else{
           if (do_cc) {
              status = triples();
           }
           if (do_mp) {
              ccmethod = 2;
              status = triples();
           }
        }
     }
     if (status == Failure){
        throw PsiException(
           "Whoops, the (T) correction died.",__FILE__,__LINE__);
     }
     tstop();

     // ccsd(t) energy
     if (do_cc) {
        Process::environment.globals["(T) CORRECTION ENERGY"] = et;
        Process::environment.globals["CURRENT CORRELATION ENERGY"] = eccsd + et;
        Process::environment.globals["CURRENT ENERGY"] = eccsd + et + escf;
        if (isccsd) {
           Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = eccsd + et;
           Process::environment.globals["CCSD(T) TOTAL ENERGY"] = eccsd + et + escf;
        }else {
           Process::environment.globals["QCISD(T) CORRELATION ENERGY"] = eccsd + et;
           Process::environment.globals["QCISD(T) TOTAL ENERGY"] = eccsd + et + escf;
        }
     }

     if (do_mp) {
        // mp4 triples:
        Process::environment.globals["MP4(T) CORRECTION ENERGY"] = emp4_t;
        Process::environment.globals["MP4(SDTQ) CORRELATION ENERGY"] = emp2+emp3+emp4_sd+emp4_q+emp4_t;
        Process::environment.globals["MP4(SDTQ) TOTAL ENERGY"] = emp2+emp3+emp4_sd+emp4_q+emp4_t+escf;
        Process::environment.globals["MP4 CORRELATION ENERGY"] = emp2+emp3+emp4_sd+emp4_q+emp4_t;
        Process::environment.globals["MP4 TOTAL ENERGY"] = emp2+emp3+emp4_sd+emp4_q+emp4_t+escf;
        if (!do_cc){
           Process::environment.globals["CURRENT CORRELATION ENERGY"] = emp2+emp3+emp4_sd+emp4_q+emp4_t;
           Process::environment.globals["CURRENT ENERGY"] = emp2+emp3+emp4_sd+emp4_q+emp4_t+escf;
        }
     }

     // if we allocated t2 just for triples, free it
     if (t2_on_disk){
        free(tb);
     }
  }

  // free remaining memory
  if (!t2_on_disk){
     free(tb);
  }
  free(t1);

  finalize();

  return Process::environment.globals["CURRENT ENERGY"];
}

void CoupledCluster::WriteBanner(){
  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile,     "        *****************************************************\n");
  fprintf(outfile,     "        *                                                   *\n");
  if (isccsd)
      fprintf(outfile, "        *                       CCSD                        *\n");
  else if (mp2_only)
      fprintf(outfile, "        *                        MP2                        *\n");
  else if (mp4_only)
      fprintf(outfile, "        *                        MP4                        *\n");
  else if (mp3_only)
      fprintf(outfile, "        *                        MP3                        *\n");
  else
      fprintf(outfile, "        *                       QCISD                       *\n");
  fprintf(outfile,     "        *                  Eugene DePrince                  *\n");
  fprintf(outfile,     "        *                                                   *\n");
  fprintf(outfile,     "        *****************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);
}

/*===================================================================

  solve cc/qci equations

===================================================================*/
PsiReturnType CoupledCluster::CCSDIterations() {

  long int o = ndoccact;
  long int v = nvirt;

  int iter              = 0;
  int diis_iter         = 0;
  int replace_diis_iter = 1;
  double nrm            = 1.0;
  double Eold           = 1.0e9;
  eccsd                 = 0.0;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  fprintf(outfile,"\n");
  if (isccsd) {
     fprintf(outfile,
       "  Begin singles and doubles coupled cluster iterations\n\n");
  }else {
     fprintf(outfile,
       "  Begin singles and doubles quadratic ci iterations\n\n");
  }
  fprintf(outfile,
    "   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");
  fflush(outfile);

  // zero residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_NEW);
  memset((void*)tempt,'\0',o*o*v*v*sizeof(double));
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  if (t2_on_disk || options_.get_bool("NAT_ORBS")){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->write_entry(PSIF_DCC_T2,"t1",(char*)&tempt[0],o*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }

  // start timing the iterations
  struct tms total_tmstime;
  const long clk_tck = sysconf(_SC_CLK_TCK);
  times(&total_tmstime);

  time_t time_start = time(NULL);
  double user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  double sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;

  bool timer = options_.get_bool("CC_TIMINGS");

  double s1,e1;
  while(iter < maxiter){
      time_t iter_start = time(NULL);

      // evaluate cc/qci diagrams
      memset((void*)w1,'\0',o*v*sizeof(double));
      if (timer) fprintf(outfile,"\n");
      for (int i = 0; i < ncctasks; i++) {
          if (timer) s1 = omp_get_wtime();
          (*this.*CCTasklist[i].func)(CCParams[i]);
          if (timer) fprintf(outfile,"        %s ... %6.2lf s\n",CCTasklist[i].name,omp_get_wtime()-s1);
      }
      if (timer) fprintf(outfile,"\n");

      // update the amplitudes
      Eold = eccsd;
      UpdateT1(iter);
      UpdateT2(iter);

      // add vector to list for diis
      DIISOldVector(iter,diis_iter,replace_diis_iter);

      // diis error vector and convergence check
      nrm = DIISErrorVector(diis_iter,replace_diis_iter,iter);

      // diis extrapolation
      if (diis_iter > 1){
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,o*o*v*v+o*v);
         else                   DIIS(diisvec,maxdiis,o*o*v*v+o*v);
         DIISNewAmplitudes(diis_iter,replace_diis_iter);
      }
      eccsd = CheckEnergy();

      if (diis_iter < maxdiis ) {
         replace_diis_iter++;
      }else {
          double min = 1.0e9;
          for (int j = 1; j <= (diis_iter < maxdiis ? diis_iter : maxdiis); j++) {
              if ( fabs( diisvec[j-1] ) < min ) {
                  min = fabs( diisvec[j-1] );
                  replace_diis_iter = j;
              }
          }
      }

      if (diis_iter <= maxdiis) diis_iter++;
      //else if (replace_diis_iter < maxdiis) replace_diis_iter++;
      //else    replace_diis_iter = 1;

      time_t iter_stop = time(NULL);
      fprintf(outfile,"  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,eccsd,eccsd-Eold,nrm,(int)iter_stop-(int)iter_start);
      fflush(outfile);
      iter++;

      // energy and amplitude convergence check
      if (fabs(eccsd - Eold) < e_conv && nrm < r_conv) break;
  }

  // stop timing iterations
  times(&total_tmstime);
  time_t time_stop = time(NULL);
  double user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  double sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;

  if ( iter==maxiter ){
     if (isccsd)
        throw PsiException("  CCSD iterations did not converge.",__FILE__,__LINE__);
     else
        throw PsiException("  QCISD iterations did not converge.",__FILE__,__LINE__);
  }

  if (reference_wavefunction_->isCIM()){
     Local_SCS_CCSD();
     Local_SCS_MP2();
  }
  else{
     SCS_CCSD();
  }

  fprintf(outfile,"\n");
  if (isccsd)
     fprintf(outfile,"  CCSD iterations converged!\n");
  else
     fprintf(outfile,"  QCISD iterations converged!\n");
  fprintf(outfile,"\n");

  // T1 and D1 diagnostics:

  double t1diag = F_DNRM2(o*v,t1,1) / sqrt(2.0 * o);
  fprintf(outfile,"        T1 diagnostic:                   %20.12lf\n",t1diag);
  boost::shared_ptr<Matrix>T (new Matrix(o,o));
  boost::shared_ptr<Matrix>eigvec (new Matrix(o,o));
  boost::shared_ptr<Vector>eigval (new Vector(o));
  double ** Tp = T->pointer();
  for (int i = 0; i < o; i++) {
      for (int j = 0; j < o; j++) {
          double dum = 0.0;
          for (int a = 0; a < v; a++) {
              dum += t1[a*o+i] * t1[a*o+j];
          }
          Tp[i][j] = dum;
      }
  }
  T->diagonalize(eigvec,eigval,descending);
  fprintf(outfile,"        D1 diagnostic:                   %20.12lf\n",sqrt(eigval->pointer()[0]));
  fprintf(outfile,"\n");

  // delta mp2 correction for fno computations:
  if (options_.get_bool("NAT_ORBS")){
      double delta_emp2 = Process::environment.globals["MP2 CORRELATION ENERGY"] - emp2;
      double delta_emp2_os = Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] - emp2_os;
      double delta_emp2_ss = Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] - emp2_ss;

      emp2 += delta_emp2;
      emp2_os += delta_emp2_os;
      emp2_ss += delta_emp2_ss;

      eccsd += delta_emp2;
      eccsd_os += delta_emp2_os;
      eccsd_ss += delta_emp2_ss;

      fprintf(outfile,"        OS MP2 FNO correction:           %20.12lf\n",delta_emp2_os);
      fprintf(outfile,"        SS MP2 FNO correction:           %20.12lf\n",delta_emp2_ss);
      fprintf(outfile,"        MP2 FNO correction:              %20.12lf\n",delta_emp2);
      fprintf(outfile,"\n");
  }

  if (options_.get_bool("SCS_MP2")){
      fprintf(outfile,"        OS SCS-MP2 correlation energy:   %20.12lf\n",emp2_os*emp2_os_fac);
      fprintf(outfile,"        SS SCS-MP2 correlation energy:   %20.12lf\n",emp2_ss*emp2_ss_fac);
      fprintf(outfile,"        SCS-MP2 correlation energy:      %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac);
      fprintf(outfile,"      * SCS-MP2 total energy:            %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac+escf);
      fprintf(outfile,"\n");
  }
  fprintf(outfile,"        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
  fprintf(outfile,"        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
  fprintf(outfile,"        MP2 correlation energy:          %20.12lf\n",emp2);
  fprintf(outfile,"      * MP2 total energy:                %20.12lf\n",emp2+escf);
  fprintf(outfile,"\n");
  fprintf(outfile,"        OS MP2.5 correlation energy:     %20.12lf\n",emp2_os+0.5*emp3_os);
  fprintf(outfile,"        SS MP2.5 correlation energy:     %20.12lf\n",emp2_ss+0.5*emp3_ss);
  fprintf(outfile,"        MP2.5 correlation energy:        %20.12lf\n",emp2+0.5*emp3);
  fprintf(outfile,"      * MP2.5 total energy:              %20.12lf\n",emp2+0.5*emp3+escf);
  fprintf(outfile,"\n");
  fprintf(outfile,"        OS MP3 correlation energy:       %20.12lf\n",emp2_os+emp3_os);
  fprintf(outfile,"        SS MP3 correlation energy:       %20.12lf\n",emp2_ss+emp3_ss);
  fprintf(outfile,"        MP3 correlation energy:          %20.12lf\n",emp2+emp3);
  fprintf(outfile,"      * MP3 total energy:                %20.12lf\n",emp2+emp3+escf);
  fprintf(outfile,"\n");
  fprintf(outfile,"        OS MP4(SDQ) correlation energy:  %20.12lf\n",emp2_os+emp3_os+emp4_sd_os+emp4_q_os);
  fprintf(outfile,"        SS MP4(SDQ) correlation energy:  %20.12lf\n",emp2_ss+emp3_ss+emp4_sd_ss+emp4_q_os);
  fprintf(outfile,"        MP4(SDQ) correlation energy:     %20.12lf\n",emp2+emp3+emp4_sd+emp4_q);
  fprintf(outfile,"      * MP4(SDQ) total energy:           %20.12lf\n",emp2+emp3+emp4_sd+emp4_q+escf);
  fprintf(outfile,"\n");

  if (isccsd) {
     if (options_.get_bool("SCS_CCSD")){
        fprintf(outfile,"        OS SCS-CCSD correlation energy:  %20.12lf\n",eccsd_os*eccsd_os_fac);
        fprintf(outfile,"        SS SCS-CCSD correlation energy:  %20.12lf\n",eccsd_ss*eccsd_ss_fac);
        fprintf(outfile,"        SCS-CCSD correlation energy:     %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
        fprintf(outfile,"      * SCS-CCSD total energy:           %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
        fprintf(outfile,"\n");
     }
     fprintf(outfile,"        OS CCSD correlation energy:      %20.12lf\n",eccsd_os);
     fprintf(outfile,"        SS CCSD correlation energy:      %20.12lf\n",eccsd_ss);
     fprintf(outfile,"        CCSD correlation energy:         %20.12lf\n",eccsd);
     fprintf(outfile,"      * CCSD total energy:               %20.12lf\n",eccsd+escf);
     fprintf(outfile,"\n");
     fprintf(outfile,"  Total time for CCSD iterations:  %10.2lf s (user)\n",user_stop-user_start);
  }else{
     if (options_.get_bool("SCS_CCSD")){
        fprintf(outfile,"        OS SCS-QCISD correlation energy: %20.12lf\n",eccsd_os*eccsd_os_fac);
        fprintf(outfile,"        SS SCS-QCISD correlation energy: %20.12lf\n",eccsd_ss*eccsd_ss_fac);
        fprintf(outfile,"        SCS-QCISD correlation energy:    %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
        fprintf(outfile,"      * SCS-QCISD total energy:          %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
        fprintf(outfile,"\n");
     }
     fprintf(outfile,"        OS QCISD correlation energy:     %20.12lf\n",eccsd_os);
     fprintf(outfile,"        SS QCISD correlation energy:     %20.12lf\n",eccsd_ss);
     fprintf(outfile,"        QCISD correlation energy:        %20.12lf\n",eccsd);
     fprintf(outfile,"      * QCISD total energy:              %20.12lf\n",eccsd+escf);
     fprintf(outfile,"\n");
     fprintf(outfile,"  Total time for QCISD iterations: %10.2lf s (user)\n",user_stop-user_start);
  }

  fprintf(outfile,"                                   %10.2lf s (system)\n",sys_stop-sys_start);
  fprintf(outfile,"                                   %10d s (total)\n",(int)time_stop-(int)time_start);
  fprintf(outfile,"\n");
  fprintf(outfile,"  Time per iteration:              %10.2lf s (user)\n",(user_stop-user_start)/iter);
  fprintf(outfile,"                                   %10.2lf s (system)\n",(sys_stop-sys_start)/iter);
  fprintf(outfile,"                                   %10.2lf s (total)\n",((double)time_stop-(double)time_start)/(iter-1));
  fflush(outfile);

  return Success;
}

/*===================================================================

  determine tiling for vabcd and vabci diagrams for the cpu
  this determines the size of blocks of integrals that 
  can be read into cpu memory.

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
  long int oovv = o*o*v*v;
  ndoubles -= o*o*v*v+2L*(oovv+o*v)+2L*o*v+2*v*v+(o+v);
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
  long int fulltile = v*(v+1L)/2L;
  ntiles=1L;
  tilesize=fulltile/1L;
  if ( ntiles*tilesize < fulltile ) tilesize++;
  while( fulltile*tilesize > ndoubles ){
     ntiles++;
     tilesize = fulltile/ntiles;
     if ( ntiles*tilesize < fulltile ) tilesize++;
  }
  lasttile = fulltile - (ntiles-1L)*tilesize;

  fprintf(outfile,"        v(ab,cd) diagrams will be evaluated in %3li blocks.\n",ntiles); 
  fflush(outfile);

  // ov^3 type 1:
  if (v>ndoubles){
     throw PsiException("out of memory: (ab,ci)",__FILE__,__LINE__);
  }
  nov2tiles=1L;
  ov2tilesize=ov2/1L;
  if ( nov2tiles*ov2tilesize < ov2 ) ov2tilesize++;
  while( v*ov2tilesize > ndoubles ){
     nov2tiles++;
     ov2tilesize = ov2/nov2tiles;
     if ( nov2tiles*ov2tilesize < ov2 ) ov2tilesize++;
  }
  lastov2tile = ov2 - (nov2tiles-1L)*ov2tilesize;

  fprintf(outfile,"        v(ab,ci) diagrams will be evaluated in %3li blocks over ov2.\n",nov2tiles); 
  fflush(outfile);

  // ov^3 type 2:
  if (v*v > ndoubles){
     throw PsiException("out of memory: (ab,ci)",__FILE__,__LINE__);
  }
  novtiles=1L;
  ovtilesize=ov/1L;
  if ( novtiles*ovtilesize < ov ) ovtilesize++;
  while( v*v*ovtilesize > ndoubles ){
     novtiles++;
     ovtilesize = ov/novtiles;
     if ( novtiles*ovtilesize < ov ) ovtilesize++;
  }
  lastovtile = ov - (novtiles-1L)*ovtilesize;
  fprintf(outfile,"        v(ab,ci) diagrams will be evaluated in %3li blocks over ov.\n",novtiles); 
  fflush(outfile);
}

/*===================================================================

  - allocate cpu memory
  - define tiling for diagrams that don't fit in core

===================================================================*/
void CoupledCluster::AllocateMemory() {

  long int nthreads = omp_get_max_threads();
  long int o=ndoccact;
  long int v=nvirt;
  if (!options_.get_bool("RUN_MP2")) {
    fprintf(outfile,"\n");
    fprintf(outfile,"  available memory =                         %9.2lf mb\n",memory/1024./1024.);
    if (isccsd){
       fprintf(outfile,"  minimum memory requirements for CCSD =     %9.2lf mb\n",
           8./1024./1024.*(o*o*v*v+2.*(o*o*v*v+o*v)+2.*o*v+2.*v*v+o+v));
    }else{
       fprintf(outfile,"  minimum memory requirements for QCISD =    %9.2lf mb\n",
           8./1024./1024.*(o*o*v*v+2.*(o*o*v*v+o*v)+2.*o*v+2.*v*v+o+v));
    }
    if (options_.get_bool("COMPUTE_TRIPLES") || options_.get_bool("COMPUTE_MP4_TRIPLES")){
       double tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+3L*v*v*v*nthreads);
       if (tempmem > memory) {
          fprintf(outfile,"\n  <<< warning! >>> switched to low-memory (t) algorithm\n\n");
       }
       if (tempmem > memory || options_.get_bool("TRIPLES_LOW_MEMORY")){
          isLowMemory = true;
          tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+5L*o*o*o*nthreads);
       }
       if (isccsd)
          fprintf(outfile,"  memory requirements for CCSD(T) =          %9.2lf mb\n",tempmem/1024./1024.);
       else
          fprintf(outfile,"  memory requirements for QCISD(T) =         %9.2lf mb\n",tempmem/1024./1024.);
    }

  }
  // orbital energies:
  if (!reference_wavefunction_->isCIM()){
      int count=0;
      eps = (double*)malloc((ndoccact+nvirt)*sizeof(double));
      boost::shared_ptr<Vector> eps_test = reference_wavefunction_->epsilon_a();
      for (int h=0; h<nirrep_; h++){
          for (int norb = frzcpi_[h]; norb<doccpi_[h]; norb++){
              eps[count++] = eps_test->get(h,norb);
          }
      }
      for (int h=0; h<nirrep_; h++){
          for (int norb = doccpi_[h]; norb<nmopi_[h]-frzvpi_[h]; norb++){
              eps[count++] = eps_test->get(h,norb);
          }
      }
  }else{
     // orbital energies in qt ordering:
     long int count = 0;
     eps = (double*)malloc((ndoccact+nvirt)*sizeof(double));
     boost::shared_ptr<Vector> eps_test = reference_wavefunction_->CIMOrbitalEnergies();
     for (int i = 0; i < ndoccact + nvirt; i++){
         eps[i] = eps_test->get(0,i+nfzc);
     }
     eps_test.reset();
  }
  if (options_.get_bool("RUN_MP2")) return;

  // define tiling for v^4 and ov^3 diagrams according to how much memory is available
  DefineTilingCPU();

  long int dim = 0;
  int fulltile = v*(v+1)/2;
  if (tilesize*fulltile > dim) dim = tilesize*fulltile;
  if (ovtilesize*v*v > dim)    dim = ovtilesize*v*v;
  if (ov2tilesize*v > dim)     dim = ov2tilesize*v;

  // if integrals buffer isn't at least o^2v^2, try tiling again assuming t2 is on disk.
  if (dim<o*o*v*v){
     fprintf(outfile,"\n");
     fprintf(outfile,"  Warning: cannot accomodate T2 in core. T2 will be stored on disk.\n");
     fprintf(outfile,"\n");
     fflush(outfile);
     t2_on_disk = true;
     DefineTilingCPU();
     dim = 0;
     if (tilesize*fulltile > dim) dim = tilesize*fulltile;
     if (ovtilesize*v*v > dim)    dim = ovtilesize*v*v;
     if (ov2tilesize*v > dim)     dim = ov2tilesize*v;

     if (dim<o*o*v*v){
        throw PsiException("out of memory: general buffer cannot accomodate T2",__FILE__,__LINE__);
     }

     fprintf(outfile,"\n");
     fprintf(outfile,"  Increase memory by %7.2lf mb to hold T2 in core.\n",o*o*v*v*8L/1024./1024.);
     fprintf(outfile,"\n");
  }

  maxelem = dim;

  long int oovv = o*o*v*v;
  double total_memory = 1.*dim+2.*(oovv+o*v)+1.*o*o*v*v+2.*o*v+2.*v*v;
  if (t2_on_disk) total_memory = 1.*dim+2.*(oovv+o*v)+2.*o*v+2.*v*v;
  total_memory *= 8./1024./1024.;

  fprintf(outfile,"\n");
  fprintf(outfile,"  Allocate cpu memory (%9.2lf mb).....",total_memory);

  integrals = (double*)malloc(dim*sizeof(double));
  tempt     = (double*)malloc((oovv+o*v)*sizeof(double));
  tempv     = (double*)malloc((oovv+o*v)*sizeof(double));

  if (!t2_on_disk) tb = (double*)malloc(o*o*v*v*sizeof(double));
  w1        = (double*)malloc(o*v*sizeof(double));
  t1        = (double*)malloc(o*v*sizeof(double));
  I1        = (double*)malloc(v*v*sizeof(double));
  I1p       = (double*)malloc(v*v*sizeof(double));
  fprintf(outfile,"done.\n");

  fprintf(outfile,"  Initialize cpu memory..................");
  memset((void*)integrals,'\0',dim*sizeof(double));
  memset((void*)tempv,'\0',(oovv+o*v)*sizeof(double));
  memset((void*)tempt,'\0',(oovv+o*v)*sizeof(double));
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

/*===================================================================

  CCSD diagrams

===================================================================*/

/**
  * t1 <-- (ma|ei)
  */
void CoupledCluster::CPU_t1_vmeai(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int i,a,m,e,id,one=1;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAB,"E2ijab",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IJAB,1);

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
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

/**
  * t1 <-- (mn|ei)
  */
void CoupledCluster::CPU_t1_vmeni(CCTaskParams params){
  long int m,e,n,a,id;
  long int o=ndoccact;
  long int v=nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (a=0,id=0; a<v; a++){
      for (m=0; m<o; m++){
          for (n=0; n<o; n++){
              F_DCOPY(v,tb+a*v*o*o+m*o+n,o*o,tempt+a*o*o*v+m*o*v+n*v,1);
              F_DAXPY(v,-2.0,tb+a*o*o+m*o+n,o*o*v,tempt+a*o*o*v+m*o*v+n*v,1);
          }
      }
  }
  psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&tempv[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_DCC_IJAK,1);
  F_DGEMM('t','n',o,v,o*o*v,1.0,tempv,o*o*v,tempt,o*o*v,1.0,w1,o);
  psio.reset();
}

/**
  *  t1 <-- (me|af)
  */
void CoupledCluster::CPU_t1_vmaef(CCTaskParams params){
  long int m,e,i,f,a,id;
  long int o=ndoccact;
  long int v=nvirt;

  boost::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (f=0,id=0; f<v; f++){
      for (m=0; m<o; m++){
          for (e=0; e<v; e++){
              F_DCOPY(o,tb+e*v*o*o+f*o*o+m*o,1,tempt+f*o*o*v+m*o*v+e*o,1);
              F_DAXPY(o,-0.5,tb+e*v*o*o+f*o*o+m,o,tempt+f*o*o*v+m*o*v+e*o,1);
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

  psio->open(PSIF_DCC_ABCI3,PSIO_OPEN_OLD);
  psio_address addr;
  addr = PSIO_ZERO;

  for (i=0; i<ntiles-1; i++){
      psio->read(PSIF_DCC_ABCI3,"E2abci3",(char*)&integrals[0],tilesize*ov2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,tilesize,ov2,2.0,tempt,o,integrals,ov2,1.0,w1+i*tilesize*o,o);
  }
  i=ntiles-1;
  psio->read(PSIF_DCC_ABCI3,"E2abci3",(char*)&integrals[0],lasttile*ov2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o,lasttile,ov2,2.0,tempt,o,integrals,ov2,1.0,w1+i*tilesize*o,o);
  psio->close(PSIF_DCC_ABCI3,1);
  psio.reset();
}

/**
  * Build and use I(a,b)
  */
void CoupledCluster::CPU_I1ab(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int b,m,n,e,a,id=0;
  // build I1(a,b)
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }
 
  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (n=0; n<o; n++){
              F_DCOPY(v,tb+e*v*o*o+m*o+n,o*o,tempt+m*o*v*v+e*o*v+n*v,1);
              if (isccsd) {
                 for (b=0; b<v; b++){
                     tempt[id++] += t1[e*o+m]*t1[b*o+n];
                 }
              }
          }
      }
  }
  F_DCOPY(o*o*v*v,integrals,1,tempv,1);
  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (n=0; n<o; n++){
              F_DAXPY(v,-0.5,integrals+m*o*v*v+n*v+e,o*v,tempv+m*o*v*v+e*o*v+n*v,1);
          }
      }
  }
  F_DGEMM('n','t',v,v,o*o*v,-2.0,tempv,v,tempt,v,0.0,I1,v);

  // add the singles parts to I1(a,b). n^4
  long int i,j,l,k,c,d;
  if (isccsd) {
     double sum=0.;
     psio->open(PSIF_DCC_ABCI2,PSIO_OPEN_OLD);
     psio_address addr;
     addr = PSIO_ZERO;

     for (i=0; i<o; i++){
         F_DCOPY(v,t1+i,o,tempt+i*v,1);
     }

     // try tiling dgemv bc ov^3 might be too large.
     long int v2tilesize,nv2tiles,lastv2tile;
     nv2tiles=1L;
     v2tilesize=v*v/1L;
     if (nv2tiles*v2tilesize<v*v) v2tilesize++;
     while(v2tilesize*o*v>maxelem){
        nv2tiles++;
        v2tilesize = v*v/nv2tiles;
        if (nv2tiles*v2tilesize<v*v) v2tilesize++;
     }
     lastv2tile = v*v - (nv2tiles-1L)*v2tilesize;

     for (i=0; i<nv2tiles-1; i++){
         psio->read(PSIF_DCC_ABCI2,"E2abci2",(char*)&integrals[0],v2tilesize*v*o*sizeof(double),addr,&addr);
         F_DGEMV('t',o*v,v2tilesize,-1.0,integrals,o*v,tempt,1,1.0,I1+i*v2tilesize,1);
     }
     i=nv2tiles-1;
     psio->read(PSIF_DCC_ABCI2,"E2abci2",(char*)&integrals[0],lastv2tile*v*o*sizeof(double),addr,&addr);
     F_DGEMV('t',o*v,lastv2tile,-1.0,integrals,o*v,tempt,1,1.0,I1+i*v2tilesize,1);


     psio->close(PSIF_DCC_ABCI2,1);
  }

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (l=0,id=0; l<o; l++){
      for (c=0; c<v; c++){
          for (k=0; k<o; k++){
              F_DCOPY(v,tb+c*o*o+l*o+k,v*o*o,tempt+l*o*v*v+c*o*v+k*v,1);
          }
      }
  }
  // use I1(a,b) for doubles residual:
  F_DGEMM('t','n',v,o*o*v,v,1.0,I1,v,tempt,v,0.0,tempv,v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+a*v*o+i*v+b,v*v*o,tempt+a*o*o*v+b*o*o+i*o,1);
              F_DAXPY(o,1.0,tempv+i*v*v*o+b*v*o+a,v,tempt+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  // use I1(a,b) for singles residual - 1st contribution to w1. (n^3)
  //F_DGEMM('n','n',o,v,v,1.0,t1,o,I1,v,1.0,w1,o);
  F_DGEMM('n','n',o,v,v,1.0,t1,o,I1,v,1.0,w1,o);

  psio.reset();
}

// a refactored version of I2p(ab,ci) that avoids ov^3 storage.
// it turns out that most of the resulting term can be dumped in 
// with other terms.  all we're left with is this 2o^3v^2 term
// and another o^3v^2 term that was added to I2piajk
void CoupledCluster::CPU_I2p_abci_refactored_term2(CCTaskParams params){
  long int o = ndoccact;
  long int v = nvirt;
  long int a,b,c,i,j,id=0;
  long int ov2 = o*v*v;
  long int o2v = o*o*v;

  boost::shared_ptr<PSIO> psio(new PSIO());

  // now build and use intermediate:
  psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAB,"E2ijab",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IJAB,1);
  F_DGEMM('n','n',o,o2v,v,-1.0,t1,o,tempv,v,0.0,tempt,o);
  F_DGEMM('n','n',o2v,v,o,1.0,tempt,o2v,t1,o,0.0,tempv,o2v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  F_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+a*v*o*o+b*o*o+i*o,1,tempt+b*v*o*o+a*o*o+i,o);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  psio.reset();
}

/**
  * Build and use I(i,j), I'(i,j), and I(i,a)
  */
void CoupledCluster::CPU_I1pij_I1ia_lessmem(CCTaskParams params){

  long int o = ndoccact;
  long int v = nvirt;
  long int m,j,e,f,i,a,b;//,one=1;
  long int ov2 = o*v*v;
  long int id=0;

  // build I1(i,a). n^4
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  F_DCOPY(o*o*v*v,integrals,1,tempv,1);
  for (i=0; i<o; i++){
      for (a=0; a<v; a++){
          for (m=0; m<o; m++){
              F_DAXPY(v,-0.5,integrals+i*o*v*v+m*v+a,o*v,tempv+i*v*v*o+a*v*o+m*v,1);
          }
      }
  }
  for (i=0; i<o; i++) F_DCOPY(v,t1+i,o,tempt+i*v,1);
  F_DGEMV('t',o*v,o*v,2.0,tempv,o*v,tempt,1,0.0,I1,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  // use I1(i,a) -> w1
  id=0;
  memset((void*)tempt,'\0',o*o*v*v);
  for (m=0; m<o; m++){
      for (e=0; e<v; e++){
          for (j=0; j<o; j++){
              F_DCOPY(v,tb+e*o*o*v+m*o+j,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
              F_DAXPY(v,-0.5,tb+e*o*o*v+j*o+m,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
          }
      }
  }
  F_DGEMV('n',o*v,o*v,2.0,tempt,o*v,I1,1,0.0,tempv,1);
  for (i=0; i<o; i++){
      F_DAXPY(v,1.0,tempv+i*v,1,w1+i,o);
  }

  // build I1'(i,j)
  F_DGEMM('t','n',o,o,ov2,2.0,tempt,ov2,integrals,ov2,0.0,I1p,o);
  
  // only n^4
  if (isccsd) {
     psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&tempt[0],o*o*o*v*sizeof(double));
     psio->close(PSIF_DCC_IJAK,1);
     id=0;
     for (i=0; i<o; i++){
         for (j=0; j<o; j++){
             for (e=0; e<v; e++){
                 F_DCOPY(o,tempt+i*o*v+j*v+e,o*o*v,tempv+i*o*o*v+j*o*v+e*o,1);
                 F_DAXPY(o,-2.0,tempt+i*o*o*v+j*v+e,o*v,tempv+i*o*o*v+j*o*v+e*o,1);
             }
         }
     }
     F_DGEMV('t',o*v,o*o,-1.0,tempv,o*v,t1,1,1.0,I1p,1);
  }

  // use I1'(i,j) for singles residual. (n^3)
  //F_DGEMM('n','n',o,v,o,-1.0,I1p,o,t1,o,1.0,w1,o);
  F_DGEMM('n','n',o,v,o,-1.0,I1p,o,t1,o,1.0,w1,o);

  // build I1(i,j)
  if (isccsd) {
     //F_DGEMM('n','n',o,o,v,1.0,t1,o,I1,v,1.0,I1p,o);
     F_DGEMM('n','n',o,o,v,1.0,t1,o,I1,v,1.0,I1p,o);
  }

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (j=0; j<o; j++){
              F_DCOPY(v,tb+e*o*o*v+m*o+j,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
          }
      }
  }
  F_DGEMM('n','t',o,ov2,o,-1.0,I1p,o,tempt,ov2,0.0,tempv,o);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+a*o*o*v+b*o+i,v*o,tempt+a*o*o*v+b*o*o+i*o,1);
              F_DAXPY(o,1.0,tempv+b*o*o*v+i*v*o+a*o,1,tempt+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  psio.reset();
}

/**
 *  Build and use I2ijkl
 */
void CoupledCluster::I2ijkl(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempt,1);
  }

  if (isccsd) {
     for (a=0,id=0; a<v; a++){
         for (b=0; b<v; b++){
             for (i=0; i<o; i++){
                 for (j=0; j<o; j++){
                     tempt[id++] += t1[a*o+i]*t1[b*o+j];
                 }
             }
         }
     }
  }
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  for (j=0; j<o; j++){
      for (i=0; i<o; i++){
          for (b=0; b<v; b++){
              F_DCOPY(v,integrals+j*o*v*v+b*o*v+i*v,1,tempv+j*o*v*v+i*v*v+b*v,1);
          }
      }
  }
  psio->open(PSIF_DCC_IJKL,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJKL,"E2ijkl",(char*)&integrals[0],o*o*o*o*sizeof(double));

  psio->close(PSIF_DCC_IJKL,1);
  F_DGEMM('n','n',o*o,o*o,v*v,1.0,tempt,o*o,tempv,v*v,1.0,integrals,o*o);
  if (isccsd) {
     psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&tempv[0],o*o*o*v*sizeof(double));
     psio->close(PSIF_DCC_IJAK,1);
     F_DGEMM('n','n',o,o*o*o,v,2.0,t1,o,tempv,v,1.0,integrals,o);
  }
  F_DGEMM('n','n',o*o,v*v,o*o,0.5,integrals,o*o,tempt,o*o,0.0,tempv,o*o);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  F_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+b*v*o*o+a*o*o+i,o,tempt+a*v*o*o+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();

}
/**
 *  Build and use I2'iajk.
 *  This contains one of the terms that came out of refactorizing I2'abci
 *  (it used to be I2p_abci_refactored_term3)
 */
void CoupledCluster::I2piajk(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  if (isccsd) {
     if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
     }else{
        F_DCOPY(o*o*v*v,tb,1,tempt,1);
     }

     for (a=0,id=0; a<v; a++){
         for (b=0; b<v; b++){
             for (i=0; i<o; i++){
                 for (j=0; j<o; j++){
                     tempt[id++] += t1[a*o+i]*t1[b*o+j];
                 }
             }
         }
     }
  }
  psio->open(PSIF_DCC_IJAK2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAK2,"E2ijak2",(char*)&tempv[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_DCC_IJAK2,1);

  if (isccsd) {
     addr = PSIO_ZERO;
     psio->open(PSIF_DCC_ABCI,PSIO_OPEN_OLD);
     for (j=0; j<novtiles-1; j++){
         psio->read(PSIF_DCC_ABCI,"E2abci",(char*)&integrals[0],ovtilesize*v*v*sizeof(double),addr,&addr);
         F_DGEMM('n','n',o*o,ovtilesize,v*v,1.0,tempt,o*o,integrals,v*v,1.0,tempv+j*o*o*ovtilesize,o*o);
     }
     j=novtiles-1;
     psio->read(PSIF_DCC_ABCI,"E2abci",(char*)&integrals[0],lastovtile*v*v*sizeof(double),addr,&addr);
     F_DGEMM('n','n',o*o,lastovtile,v*v,1.0,tempt,o*o,integrals,v*v,1.0,tempv+j*o*o*ovtilesize,o*o);
     psio->close(PSIF_DCC_ABCI,1);

     // this used to be part of I2p(ab,ci) ... see notes ...
     psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_IAJB,1);
     F_DGEMM('t','t',o*o*v,o,v,1.0,integrals,v,t1,o,0.0,tempt,o*o*v);
     for (j=0; j<o; j++){
         for (a=0; a<v; a++){
             for (i=0; i<o; i++){
                 F_DAXPY(o,1.0,tempt+i*o*o*v+a*o+j,o*v,tempv+j*o*o*v+a*o*o+i*o,1);
             }
         }
     }
  }

  // use intermediate
  F_DGEMM('n','n',o*o*v,v,o,-1.0,tempv,o*o*v,t1,o,0.0,tempt,o*o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  F_DAXPY(o*o*v*v,1.0,tempt,1,tempv,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempt+b*v*o*o+a*o*o+i,o,tempv+a*v*o*o+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
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
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempt,1);
  }
  if (isccsd){
     for (a=0,id=0; a<v; a++){
         for (b=0; b<v; b++){
             for (i=0; i<o; i++){
                 for (j=0; j<o; j++){
                     tempt[id++] += t1[a*o+i]*t1[b*o+j];
                 }
             }
         }
     }
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
  psio->open(PSIF_DCC_ABCD1,PSIO_OPEN_OLD);
  addr = PSIO_ZERO;
  for (j=0; j<ntiles-1; j++){
      psio->read(PSIF_DCC_ABCD1,"E2abcd1",(char*)&integrals[0],tilesize*v*(v+1)/2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o*(o+1)/2,tilesize,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  }
  j=ntiles-1;
  psio->read(PSIF_DCC_ABCD1,"E2abcd1",(char*)&integrals[0],lasttile*v*(v+1)/2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o*(o+1)/2,lasttile,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  psio->close(PSIF_DCC_ABCD1,1);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tempv[a*o*o*v+b*o*o+i*o+j] += .5*tempt[Position(a,b)*o*(o+1)/2+Position(i,j)];
              }
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();

}

/**
 *  Use Vabcd2
 */
void CoupledCluster::Vabcd2(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  int sg,sg2;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempt,1);
  }
  if (isccsd) {
     for (a=0,id=0; a<v; a++){
         for (b=0; b<v; b++){
             for (i=0; i<o; i++){
                 for (j=0; j<o; j++){
                     tempt[id++] += t1[a*o+i]*t1[b*o+j];
                 }
             }
         }
     }
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
  psio->open(PSIF_DCC_ABCD2,PSIO_OPEN_OLD);
  addr = PSIO_ZERO;
  for (j=0; j<ntiles-1; j++){
      psio->read(PSIF_DCC_ABCD2,"E2abcd2",(char*)&integrals[0],tilesize*v*(v+1)/2*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o*(o+1)/2,tilesize,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  }
  j = ntiles-1;
  psio->read(PSIF_DCC_ABCD2,"E2abcd2",(char*)&integrals[0],lasttile*v*(v+1)/2*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o*(o+1)/2,lasttile,v*(v+1)/2,1.0,tempv,o*(o+1)/2,integrals,v*(v+1)/2,0.0,tempt+j*tilesize*o*(o+1)/2,o*(o+1)/2);
  psio->close(PSIF_DCC_ABCD2,1);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          if (a>b) sg2 = -1;
          else     sg2 = 1;
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  if (i>j) sg = -sg2;
                  else     sg = sg2;
                  tempv[a*o*o*v+b*o*o+i*o+j] += .5*sg*tempt[Position(a,b)*o*(o+1)/2+Position(i,j)];
              }
          }
      }
  }
  //psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
  psio.reset();
}

/**
 *  K from the SJS paper
 */
void CoupledCluster::K(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAB,"E2ijab",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IJAB,1);

  // o^2v^3 work
  if (isccsd) {
      addr = PSIO_ZERO;
      psio->open(PSIF_DCC_ABCI3,PSIO_OPEN_OLD);

      for (j=0; j<nov2tiles-1; j++){
          psio->read(PSIF_DCC_ABCI3,"E2abci3",(char*)&integrals[0],ov2tilesize*v*sizeof(double),addr,&addr);
          F_DGEMM('n','n',o,ov2tilesize,v,1.0,t1,o,integrals,v,0.0,tempv+j*o*ov2tilesize,o);
      }
      j=nov2tiles-1;
      psio->read(PSIF_DCC_ABCI3,"E2abci3",(char*)&integrals[0],lastov2tile*v*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,lastov2tile,v,1.0,t1,o,integrals,v,0.0,tempv+j*o*ov2tilesize,o);
      psio->close(PSIF_DCC_ABCI3,1);

      for (i=0; i<o; i++){
          for (b=0; b<v; b++){
              for (j=0; j<o; j++){
                  F_DAXPY(v,1.0,tempv+b*o*o+i*o+j,o*o*v,tempt+i*o*v*v+b*o*v+j*v,1);
              }
          }
      }

      // stick o^3v^2 work on first tile
      psio->open(PSIF_DCC_IJAK2,PSIO_OPEN_OLD);
      psio->read_entry(PSIF_DCC_IJAK2,"E2ijak2",(char*)&integrals[0],o*o*o*v*sizeof(double));
      psio->close(PSIF_DCC_IJAK2,1);
      // TODO: this was a problem with cuda 3.2 vs 4.0
      F_DGEMM('t','n',o*o*v,v,o,-1.0,integrals,o,t1,o,0.0,tempv,o*o*v);

      for (i=0; i<o; i++){
          for (b=0; b<v; b++){
              for (j=0; j<o; j++){
                  F_DAXPY(v,1.0,tempv+i*o*v+b*o+j,o*o*v,tempt+i*o*v*v+b*o*v+j*v,1);
              }
          }
      }
  }

  // before adding o^3v^3 term, write this part of I2(ia,jb) to disk:
  // written as ... ibja (i think..)
  psio->open(PSIF_DCC_TEMP,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_DCC_TEMP,"temporary_K",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_TEMP,1);

  // o^3v^3 part
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }
  for (i=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              F_DCOPY(v,tb+b*o*o*v+j*o+i,o*o,integrals+i*o*v*v+b*o*v+j*v,1);
              if ( isccsd ) {
                  for (a=0; a<v; a++){
                      integrals[i*o*v*v+b*o*v+j*v+a] += 2.0*t1[a*o+i]*t1[b*o+j];
                  }
              }
          }
      }
  }
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  for (i=0; i<o; i++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              F_DCOPY(v,tempt+i*v*v*o+j*v+b,o*v,tempv+i*o*v*v+b*o*v+j*v,1);
          }
      }
  }
  psio->open(PSIF_DCC_TEMP,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_TEMP,"temporary_K",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_TEMP,1);
  F_DGEMM('n','n',o*v,o*v,o*v,-0.5,integrals,o*v,tempv,o*v,1.0,tempt,o*v);

  // use I2iajb

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = integrals;
  }

  for (j=0; j<o; j++){
      for (a=0; a<v; a++){
          for (i=0; i<o; i++){
              F_DCOPY(v,tb+a*o*o+j*o+i,o*o*v,tempv+j*o*v*v+a*o*v+i*v,1);
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,-1.0,tempt,o*v,tempv,o*v,0.0,integrals,o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0,id=0; a<v; a++){
      for (b=0; b<v; b++){
          for (j=0; j<o; j++){
              F_DCOPY(o,    integrals+j*o*v*v+b*v*o+a,v,tempt+a*o*o*v+b*o*o+j*o,1);
              F_DAXPY(o,1.0,integrals+a*v*o+j*v+b,o*v*v,tempt+a*o*o*v+b*o*o+j*o,1);
              F_DAXPY(o,0.5,integrals+j*o*v*v+a*v*o+b,v,tempt+a*o*o*v+b*o*o+j*o,1);
              F_DAXPY(o,0.5,integrals+b*v*o+j*v+a,o*v*v,tempt+a*o*o*v+b*o*o+j*o,1);
          }
      }
  }

  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  psio.reset();
}
/**
 *  2J - K
 */
void CoupledCluster::TwoJminusK(CCTaskParams params){
  long int id,i,j,a,b,o,v;
  o = ndoccact;
  v = nvirt;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  // o^2v^3 contribution to intermediate
  if ( isccsd ) {
      psio->open(PSIF_DCC_IJAK,PSIO_OPEN_OLD);
      psio->read_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&integrals[0],o*o*o*v*sizeof(double));
      psio->close(PSIF_DCC_IJAK,1);
      F_DGEMM('n','n',o*o*v,v,o,-1.0,integrals,o*o*v,t1,o,0.0,tempt,o*o*v);

      psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
      psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempv[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_IAJB,1);
      for (i=0; i<o; i++){
          for (b=0; b<v; b++){
              for (j=0; j<o; j++){
                  F_DAXPY(v,1.0,tempt+i*o*v+j*v+b,o*o*v,tempv+i*o*v*v+b*o*v+j*v,1);
              }
          }
      }
  }else {
      psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
      psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempv[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_IAJB,1);
  }
  // contribute to intermediate
  psio->open(PSIF_DCC_TEMP,PSIO_OPEN_OLD);
  psio->write_entry(PSIF_DCC_TEMP,"temporary_J",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_TEMP,1);

  // the only o^3v^3 part of 2J-K
  // 1/2 ( 2(ib|me) - (ie|mb) ) ( t(ae,jm) - 1/2t(ea,jm) - t(e,j)t(a,m) )

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  F_DCOPY(o*o*v*v,tempt,1,integrals,1);
  for (i = 0,id=0; i < o; i++){
  for (b = 0; b < v; b++){
  for (int m = 0; m < o; m++){
  for (int e = 0; e < v; e++){
      integrals[id++] -= 0.5 * tempt[i*o*v*v+e*o*v+m*v+b];
  }}}}
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempt;
  }
  if (isccsd) {
  for (i = 0,id=0; i < o; i++){
  for (b = 0; b < v; b++){
  for (int m = 0; m < o; m++){
  for (int e = 0; e < v; e++){
      tempv[id++]    = tb[e*o*o*v+b*o*o+m*o+i] - 0.5 * tb[b*o*o*v+e*o*o+m*o+i] - t1[b*o+m]*t1[e*o+i];
  }}}}
  }else{
  for (i = 0,id=0; i < o; i++){
  for (b = 0; b < v; b++){
  for (int m = 0; m < o; m++){
  for (int e = 0; e < v; e++){
      tempv[id++]    = tb[e*o*o*v+b*o*o+m*o+i] - 0.5 * tb[b*o*o*v+e*o*o+m*o+i];
  }}}}
  }

  F_DGEMM('t','n',o*v,o*v,o*v,1.0,tempv,o*v,integrals,o*v,0.0,tempt,o*v);

  // o^2v^3 piece of intermediate ... also, this is identical to the contribution to the
  // residual that used to be found in I2p_abci_refactored_term1
  addr = PSIO_ZERO;
  psio->open(PSIF_DCC_ABCI,PSIO_OPEN_OLD);

  for (j=0; j<nov2tiles-1; j++){
      psio->read(PSIF_DCC_ABCI,"E2abci",(char*)&integrals[0],ov2tilesize*v*sizeof(double),addr,&addr);
      F_DGEMM('n','n',o,ov2tilesize,v,1.0,t1,o,integrals,v,0.0,tempv+j*o*ov2tilesize,o);
  }
  j=nov2tiles-1;
  psio->read(PSIF_DCC_ABCI,"E2abci",(char*)&integrals[0],lastov2tile*v*sizeof(double),addr,&addr);
  F_DGEMM('n','n',o,lastov2tile,v,1.0,t1,o,integrals,v,0.0,tempv+j*o*ov2tilesize,o);
  psio->close(PSIF_DCC_ABCI,1);

  if ( isccsd ) {
      for (i=0; i<o; i++){
          for (a=0; a<v; a++){
              for (b=0; b<v; b++){
                  F_DAXPY(o,1.0,tempv+i*o*v*v+a*o*v+b*o,1,tempt+i*o*v*v+b*o*v+a,v);
              }
          }
      }
  }

  // contribute to residual from I2p_abci_refactored_term1 ... if we know
  // this is the first diagram, we don't need to read in the old residual.
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,tempv+i*v*v*o+b*o*v+a*o,1,integrals+a*v*o*o+b*o*o+i*o,1);
              F_DAXPY(o,1.0,tempv+i+a*o*v+b*o,v*v*o,integrals+a*v*o*o+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  // contribute to intermediate
  psio->open(PSIF_DCC_TEMP,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_TEMP,"temporary_J",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_TEMP,1);
  F_DAXPY(o*o*v*v,1.0,tempt,1,tempv,1);

  // term from K stored as ibja
  psio->open(PSIF_DCC_TEMP,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_TEMP,"temporary_K",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_TEMP,1);
  // contribute K pieces to intermediate
  for (j=0,id=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (a=0; a<v; a++){
                  tempv[id++] -= 0.5 * tempt[j*o*v*v+b*o*v+i*v+a];
  }}}}

  // use I2iabj
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = integrals;
  }
  for (j=0; j<o; j++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DCOPY(v,tb+b*o*o+i*o+j,o*o*v,tempt+j*o*v*v+b*o*v+i*v,1);
              F_DAXPY(v,-0.5,tb+b*o*o*v+i*o+j,o*o,tempt+j*o*v*v+b*o*v+i*v,1);
          }
      }
  }

  F_DGEMM('n','n',o*v,o*v,o*v,2.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              F_DAXPY(o,1.0,integrals+b*v*o+i*v+a,o*v*v,tempt+a*o*o*v+b*o*o+i*o,1);
              F_DAXPY(o,1.0,integrals+i*o*v*v+a*v*o+b,v,tempt+a*o*o*v+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  psio.reset();
}

/*================================================================

   update amplitudes

================================================================*/
void CoupledCluster::UpdateT2(long int iter){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  // we still have the residual in memory in tempv
  //psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  //psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  #pragma omp parallel for schedule (static)
  for (long int a=o; a<rs; a++){
      double da = eps[a];
      for (long int b=o; b<rs; b++){
          double dab = da + eps[b];
          for (long int i=0; i<o; i++){
              double dabi = dab - eps[i];
              for (long int j=0; j<o; j++){

                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;

                  double dijab = dabi-eps[j];
                  double tnew  = - (integrals[iajb] + tempv[ijab])/dijab;
                  tempt[ijab]  = tnew;

              }
          }
      }
  }

  // error vectors for diis are in tempv:
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempv,1);
  }
  F_DAXPY(o*o*v*v,-1.0,tempt,1,tempv,1);
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tempt,1,tb,1);
  }
  psio.reset();
}
void CoupledCluster::UpdateT1(long int iter){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  #pragma omp parallel for schedule (static)
  for (long int a=o; a<rs; a++){
      for (long int i=0; i<o; i++){
          double dia    = -eps[i]+eps[a];
          double tnew   = - (w1[(a-o)*o+i])/dia;
          w1[(a-o)*o+i] = tnew;
      }
  }
  // error vector for diis is in tempv:
  F_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  F_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  F_DCOPY(o*v,w1,1,t1,1);
}


/*================================================================

   SCS functions.  get energy in terms of spin components

================================================================*/
void DFCoupledCluster::Local_SCS_CCSD(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  double ssenergy = 0.0;
  double osenergy = 0.0;

  boost::shared_ptr<PSIO> psio(new PSIO());

  // df (ia|bj) formerly E2klcd
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,tempt,o*v);

  SharedMatrix Rii = reference_wavefunction_->CIMTransformationMatrix();
  // the dimension of the Rii transformation matrix:
  int nocc = Rii->colspi()[0];
  double**Rii_pointer = Rii->pointer();
  // transform E2iajb back from quasi-canonical basis
  F_DGEMM('N','T',v*v*o,o,o,1.0,tempt,v*v*o,&Rii_pointer[0][0],nocc,0.0,integrals,v*v*o);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }
  // transform t2 back from quasi-canonical basis
  for (int a = 0; a < v; a++){
      for (int b = 0; b < v; b++){
          for (int i = 0; i < o; i++){
              for (int j = 0; j < o; j++){
                  tb[a*o*o*v+b*o*o+i*o+j] += t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }
  F_DGEMM('N','N',o,v*v*o,o,1.0,&Rii_pointer[0][0],nocc,tb,o,0.0,tempt,o);
  // resort t(abji) -> t(abij)
  for (int ab = 0; ab < v*v; ab++){
      for (int i = 0; i < o; i++){
          for (int j = 0; j < o; j++){
              I1[i*o+j] = tempt[ab*o*o+j*o+i];
          }
      }
      F_DCOPY(o*o,I1,1,tempt+ab*o*o,1);
  }
  for (int a = 0; a < v; a++){
      for (int b = 0; b < v; b++){
          for (int i = 0; i < o; i++){
              for (int j = 0; j < o; j++){
                  tb[a*o*o*v+b*o*o+i*o+j] -= t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }

  // energy
  SharedVector factor = reference_wavefunction_->CIMOrbitalFactors();
  double*factor_pointer = factor->pointer();
  for (int b = o; b < rs; b++){
      for (int a = o; a < rs; a++){
          for (int i = 0; i < o; i++){
              for (int j = 0; j < o; j++){
                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  long int ijab = (b-o)*o*o*v+(a-o)*o*o+i*o+j;
                  osenergy += integrals[iajb]*(tempt[ijab])*factor_pointer[i];
                  ssenergy += integrals[iajb]*(tempt[ijab]-tempt[(a-o)*o*o*v+(b-o)*o*o+i*o+j])*factor_pointer[i];
              }
          }
      }
  }
  eccsd_os = osenergy;
  eccsd_ss = ssenergy;
  eccsd = eccsd_os + eccsd_ss;

  psio.reset();
}
void CoupledCluster::Local_SCS_CCSD(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  SharedMatrix Rii = reference_wavefunction_->CIMTransformationMatrix();
  // the dimension of the Rii transformation matrix:
  int nocc = Rii->colspi()[0];
  double**Rii_pointer = Rii->pointer();
  // transform E2iajb back from quasi-canonical basis
  F_DGEMM('N','T',v*v*o,o,o,1.0,tempt,v*v*o,&Rii_pointer[0][0],nocc,0.0,integrals,v*v*o);

  double fac = isccsd ? 1.0 : 0.0;
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }
  // transform t2 back from quasi-canonical basis
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tb[a*o*o*v+b*o*o+i*o+j] += fac*t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }
  F_DGEMM('N','N',o,v*v*o,o,1.0,&Rii_pointer[0][0],nocc,tb,o,0.0,tempt,o);
  // resort t(abji) -> t(abij)
  for (int ab = 0; ab < v*v; ab++){
      for (i = 0; i < o; i++){
          for (j = 0; j < o; j++){
              I1[i*o+j] = tempt[ab*o*o+j*o+i];
          }
      }
      F_DCOPY(o*o,I1,1,tempt+ab*o*o,1);
  }
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  tb[a*o*o*v+b*o*o+i*o+j] -= fac*t1[a*o+i]*t1[b*o+j];
              }
          }
      }
  }

  // energy
  SharedVector factor = reference_wavefunction_->CIMOrbitalFactors();
  double*factor_pointer = factor->pointer();
  for (b=o; b<rs; b++){
      for (a=o; a<rs; a++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){
                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  osenergy += integrals[iajb]*(tempt[ijab])*factor_pointer[i];
                  ssenergy += integrals[iajb]*(tempt[ijab]-tempt[(a-o)*o*o*v+(b-o)*o*o+i*o+j])*factor_pointer[i];
                  ijab++;
              }
          }
      }
  }
  eccsd_os = osenergy;
  eccsd_ss = ssenergy;
  eccsd = eccsd_os + eccsd_ss;

  psio.reset();
}
void CoupledCluster::Local_SCS_MP2(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  SharedMatrix Rii = reference_wavefunction_->CIMTransformationMatrix();
  // the dimension of the Rii transformation matrix:
  int nocc = Rii->colspi()[0];
  double**Rii_pointer = Rii->pointer();
  // transform E2iajb back from quasi-canonical basis
  F_DGEMM('N','T',v*v*o,o,o,1.0,tempt,v*v*o,&Rii_pointer[0][0],nocc,0.0,integrals,v*v*o);

  psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_T2,"first",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_T2,1);

  // transform t2 back from quasi-canonical basis
  F_DGEMM('N','N',o,v*v*o,o,1.0,Rii_pointer[0],nocc,tempv,o,0.0,tempt,o);

  // resort t(abji) -> t(abij)
  for (int ab = 0; ab < v*v; ab++){
      for (i = 0; i < o; i++){
          for (j = 0; j < o; j++){
              I1[i*o+j] = tempt[ab*o*o+j*o+i];
          }
      }
      F_DCOPY(o*o,I1,1,tempt+ab*o*o,1);
  }

  SharedVector factor = reference_wavefunction_->CIMOrbitalFactors();
  double*factor_pointer = factor->pointer();
  for (b=o; b<rs; b++){
      for (a=o; a<rs; a++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  osenergy += integrals[iajb]*(tempt[ijab])*factor_pointer[i];
                  ssenergy += integrals[iajb]*(tempt[ijab]-tempt[(a-o)*o*o*v+(b-o)*o*o+i*o+j])*factor_pointer[i];
                  ijab++;
              }
          }
      }
  }
  emp2_os = osenergy;
  emp2_ss = ssenergy;
  emp2 = emp2_ss + emp2_os;

  psio.reset();
}
void CoupledCluster::SCS_CCSD(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }
  double fac = isccsd ? 1.0 : 0.0;
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  osenergy += integrals[iajb]*(tb[ijab]+fac*t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                  ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                  ssenergy += integrals[iajb]*fac*(t1[(a-o)*o+i]*t1[(b-o)*o+j]-t1[(b-o)*o+i]*t1[(a-o)*o+j]);
                  ijab++;
              }
          }
      }
  }
  eccsd_os = osenergy;
  eccsd_ss = ssenergy;

  psio.reset();
}
void CoupledCluster::SCS_MP2(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  osenergy += integrals[iajb]*tb[ijab];
                  ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                  ijab++;
              }
          }
      }
  }
  emp2_os = osenergy;
  emp2_ss = ssenergy;
  emp2    = emp2_os + emp2_ss;

  psio.reset();
}
double CoupledCluster::CheckEnergy(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  double ta,tnew,dijab,da,dab,dabi;
  long int iajb,jaib,ijab=0;
  double energy = 0.0;
  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }
  double fac = isccsd ? 1.0 : 0.0;
  for (a=o; a<rs; a++){
      for (b=o; b<rs; b++){
          for (i=0; i<o; i++){
              for (j=0; j<o; j++){

                  iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  jaib = iajb + (i-j)*v*(1-v*o);
                  energy += (2.*integrals[iajb]-integrals[jaib])*(tb[ijab]+fac*t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                  ijab++;
              }
          }
      }
  }

  psio.reset();

  return energy;
}

/**
 *  function pointers to CC tasks:
 */
void CoupledCluster::DefineTasks(){
  CCTasklist = new CCTask[1000];
  CCParams   = new CCTaskParams[1000];

  ncctasks=0;

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::K;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"K                      ");

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::TwoJminusK;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"2J-K                   ");

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::I2ijkl;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"I(ij,kl)               ");

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::I2piajk;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"I'(ia,jk)              ");

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::CPU_t1_vmeni;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"t1 <-- (mn|ei)         ");

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::CPU_t1_vmaef;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"t1 <-- (me|af)         ");

  if (isccsd) {
     CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::CPU_I2p_abci_refactored_term2;
     CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
     sprintf(CCTasklist[ncctasks++].name,"I'(ab,ci)              ");
  }

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::CPU_I1ab;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"I(a,b)                 ");

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::CPU_t1_vmeai;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"t1 <-- (ma|ei)         ");

  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::CPU_I1pij_I1ia_lessmem;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"I'(i,j), I(i,j), I(i,a)");

  // mo basis, sjs packing
  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::Vabcd1;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"t2 <-- (ac|bd)+        ");
  // this is the last diagram that contributes to doubles residual,
  // so we can keep it in memory rather than writing and rereading
  CCTasklist[ncctasks].func  = &psi::fnocc::CoupledCluster::Vabcd2;
  CCTasklist[ncctasks].name  = (char*)malloc(100*sizeof(char));
  sprintf(CCTasklist[ncctasks++].name,"t2 <-- (ac|bd)-        ");
}
 
void CoupledCluster::MP4_SDQ(){
  boost::shared_ptr<PSIO> psio(new PSIO());
  int o = ndoccact;
  int v = nvirt;

  // cc diagrams split up as tasks 
  // (define here so i don't get free errors when doing only mp4)
  DefineTasks();

  // define tasks for mp3
  DefineLinearTasks();

  // define tasks for mp4
  DefineQuadraticTasks();

  // 1st-order amplitudes
  fprintf(outfile,"\n");
  if (mp3_only) {
      fprintf(outfile,"  ==>   MP3   <==\n");
  }else {
      fprintf(outfile,"  ==> MP4(SDQ) <==\n");
  }
  fprintf(outfile,"\n");
  fprintf(outfile,"        1st-order doubles amplitudes...............done.\n");
  fprintf(outfile,"        MP2 energy.................................");
  UpdateT2_mp4(0);
  fprintf(outfile,"done.\n");

  // <1|V|1> for mp3
  fprintf(outfile,"        MP3 energy.................................");
  memset((void*)w1,'\0',o*v*sizeof(double));
  for (int i=0; i<nltasks; i++) {
      (*this.*LTasklist[i].func)(LParams[i]);
  }
  // mp3 energy and 2nd-order doubles amplitudes
  UpdateT2_mp4(1);
  fprintf(outfile,"done.\n");

  if (!mp3_only) {
      fprintf(outfile,"        2nd-order singles and doubles amplitudes...");

      // 2nd-order singles amplitudes
      UpdateT1_mp4(1);
      fprintf(outfile,"done.\n");

      // V|2> for S and D parts of mp4
      psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
      psio->write_entry(PSIF_DCC_T2,"second",(char*)&tempt[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_T2,1);
      if (t2_on_disk) {
          psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
          psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
          psio->close(PSIF_DCC_T2,1);
      }else {
          F_DCOPY(o*o*v*v,tempt,1,tb,1);
      }

      fprintf(outfile,"        MP4(SD)....................................");
      memset((void*)w1,'\0',o*v*sizeof(double));
      for (int i=0; i<nltasks; i++) {
          (*this.*LTasklist[i].func)(LParams[i]);
      }
      //F_DCOPY(o*o*v*v,tb,1,tempt,1);
      if (t2_on_disk) {
          tb = tempt;
      }
      psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
      psio->read_entry(PSIF_DCC_T2,"first",(char*)&tb[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_T2,1);
      UpdateT2_mp4(2);
      fprintf(outfile,"done.\n");

      // quadruples: evaluate ccd residual using first-order doubles amplitudes
      fprintf(outfile,"        MP4(Q).....................................");
      for (int i=0; i<nqtasks; i++) {
          (*this.*QTasklist[i].func)(QParams[i]);
      }
      UpdateT2_mp4(3);
      fprintf(outfile,"done.\n");
  }

  if (mp4_only) {
      if ( options_.get_bool("NAT_ORBS") ) {
          double delta_emp2 = Process::environment.globals["MP2 CORRELATION ENERGY"] - emp2;
          double delta_emp2_os = Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] - emp2_os;
          double delta_emp2_ss = Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] - emp2_ss;
 
          emp2 += delta_emp2;
          emp2_os += delta_emp2_os;
          emp2_ss += delta_emp2_ss;

          fprintf(outfile,"\n");
          fprintf(outfile,"        OS MP2 FNO correction:          %20.12lf\n",delta_emp2_os);
          fprintf(outfile,"        SS MP2 FNO correction:          %20.12lf\n",delta_emp2_ss);
          fprintf(outfile,"        MP2 FNO correction:             %20.12lf\n",delta_emp2);
      }
 
      fprintf(outfile,"\n");
      fprintf(outfile,"        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
      fprintf(outfile,"        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
      fprintf(outfile,"        MP2 correlation energy:          %20.12lf\n",emp2);
      fprintf(outfile,"      * MP2 total energy:                %20.12lf\n",emp2+escf);
      fprintf(outfile,"\n");
      fprintf(outfile,"        OS MP2.5 correlation energy:     %20.12lf\n",emp2_os/emp2_os_fac+0.5*emp3_os);
      fprintf(outfile,"        SS MP2.5 correlation energy:     %20.12lf\n",emp2_ss/emp2_ss_fac+0.5*emp3_ss);
      fprintf(outfile,"        MP2.5 correlation energy:        %20.12lf\n",emp2+0.5*emp3);
      fprintf(outfile,"      * MP2.5 total energy:              %20.12lf\n",emp2+0.5*emp3+escf);
      fprintf(outfile,"\n");
      fprintf(outfile,"        OS MP3 correlation energy:       %20.12lf\n",emp2_os/emp2_os_fac+emp3_os);
      fprintf(outfile,"        SS MP3 correlation energy:       %20.12lf\n",emp2_ss/emp2_ss_fac+emp3_ss);
      fprintf(outfile,"        MP3 correlation energy:          %20.12lf\n",emp2+emp3);
      fprintf(outfile,"      * MP3 total energy:                %20.12lf\n",emp2+emp3+escf);
      fprintf(outfile,"\n");
      fprintf(outfile,"        OS MP4(SDQ) correlation energy:  %20.12lf\n",emp2_os/emp2_os_fac+emp3_os+emp4_sd_os+emp4_q_os);
      fprintf(outfile,"        SS MP4(SDQ) correlation energy:  %20.12lf\n",emp2_ss/emp2_ss_fac+emp3_ss+emp4_sd_ss+emp4_q_os);
      fprintf(outfile,"        MP4(SDQ) correlation energy:     %20.12lf\n",emp2+emp3+emp4_sd+emp4_q);
      fprintf(outfile,"      * MP4(SDQ) total energy:           %20.12lf\n",emp2+emp3+emp4_sd+emp4_q+escf);
      fprintf(outfile,"\n");
  }else if (mp3_only){
      if ( options_.get_bool("NAT_ORBS") ) {
          double delta_emp2 = Process::environment.globals["MP2 CORRELATION ENERGY"] - emp2;
          double delta_emp2_os = Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] - emp2_os;
          double delta_emp2_ss = Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] - emp2_ss;
 
          emp2 += delta_emp2;
          emp2_os += delta_emp2_os;
          emp2_ss += delta_emp2_ss;

          fprintf(outfile,"\n");
          fprintf(outfile,"        OS MP2 FNO correction:          %20.12lf\n",delta_emp2_os);
          fprintf(outfile,"        SS MP2 FNO correction:          %20.12lf\n",delta_emp2_ss);
          fprintf(outfile,"        MP2 FNO correction:             %20.12lf\n",delta_emp2);
      }
 
      fprintf(outfile,"\n");
      fprintf(outfile,"        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
      fprintf(outfile,"        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
      fprintf(outfile,"        MP2 correlation energy:          %20.12lf\n",emp2);
      fprintf(outfile,"      * MP2 total energy:                %20.12lf\n",emp2+escf);
      fprintf(outfile,"\n");
      fprintf(outfile,"        OS MP2.5 correlation energy:     %20.12lf\n",emp2_os/emp2_os_fac+0.5*emp3_os);
      fprintf(outfile,"        SS MP2.5 correlation energy:     %20.12lf\n",emp2_ss/emp2_ss_fac+0.5*emp3_ss);
      fprintf(outfile,"        MP2.5 correlation energy:        %20.12lf\n",emp2+0.5*emp3);
      fprintf(outfile,"      * MP2.5 total energy:              %20.12lf\n",emp2+0.5*emp3+escf);
      fprintf(outfile,"\n");
      fprintf(outfile,"        OS MP3 correlation energy:       %20.12lf\n",emp2_os/emp2_os_fac+emp3_os);
      fprintf(outfile,"        SS MP3 correlation energy:       %20.12lf\n",emp2_ss/emp2_ss_fac+emp3_ss);
      fprintf(outfile,"        MP3 correlation energy:          %20.12lf\n",emp2+emp3);
      fprintf(outfile,"      * MP3 total energy:                %20.12lf\n",emp2+emp3+escf);
      fprintf(outfile,"\n");
  }else {
     // guess for cc/qci should be |1> + |2>
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"second",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     F_DAXPY(o*o*v*v,1.0,tempt,1,tb,1);
  }
}

// coupled cluster constructor
DFCoupledCluster::DFCoupledCluster(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options):
        CoupledCluster(reference_wavefunction,options)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

DFCoupledCluster::~DFCoupledCluster()
{
}

double DFCoupledCluster::compute_energy() {

  PsiReturnType status = Success;

  //WriteBanner();
  AllocateMemory();
  status = CCSDIterations();

  // free some memory!
  free(Fij);
  free(Fab);
  //free(Fia);
  //free(Fai);
  free(Qmo);
  free(Abij);
  free(Sbij);
  free(integrals);
  free(w1);
  free(I1);
  free(I1p);
  free(diisvec);
  free(tempt);
  free(tempv);

  // tstart in fnocc
  tstop();

  // mp2 energy
  Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2;
  Process::environment.globals["MP2 TOTAL ENERGY"] = emp2 + escf;
  Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
  Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;

  // ccsd energy
  Process::environment.globals["CCSD CORRELATION ENERGY"] = eccsd;
  Process::environment.globals["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
  Process::environment.globals["CCSD SAME-SPIN CORRELATION ENERGY"] = eccsd_ss;
  Process::environment.globals["CCSD TOTAL ENERGY"] = eccsd + escf;
  Process::environment.globals["CURRENT ENERGY"] = eccsd + escf;

  if (options_.get_bool("COMPUTE_TRIPLES")){
      long int o = ndoccact;
      long int v = nvirt;

      if (!isLowMemory) {
          // write (ov|vv) integrals, formerly E2abci, for (t)
          double *tempq = (double*)malloc(v*nQ*sizeof(double));
          // the buffer integrals was at least 2v^3, so these should definitely fit.
          double *Z     = (double*)malloc(v*v*v*sizeof(double));
          double *Z2    = (double*)malloc(v*v*v*sizeof(double));
          boost::shared_ptr<PSIO> psio(new PSIO());
          psio->open(PSIF_DCC_ABCI,PSIO_OPEN_NEW);
          psio_address addr2 = PSIO_ZERO;
          for (long int i=0; i<o; i++){
              #pragma omp parallel for schedule (static)
              for (long int q=0; q<nQ; q++){
                  for (long int b=0; b<v; b++){
                      tempq[q*v+b] = Qov[q*o*v+i*v+b];
                  }
              }      
              F_DGEMM('n','t',v,v*v,nQ,1.0,tempq,v,Qvv,v*v,0.0,&Z[0],v);
              #pragma omp parallel for schedule (static)
              for (long int a=0; a<v; a++){
                  for (long int b=0; b<v; b++){
                      for (long int c=0; c<v; c++){
                          Z2[a*v*v+b*v+c] = Z[a*v*v+c*v+b];
                      }
                  }
              }
              psio->write(PSIF_DCC_ABCI,"E2abci",(char*)&Z2[0],v*v*v*sizeof(double),addr2,&addr2);
          }
          psio->close(PSIF_DCC_ABCI,1);
          free(tempq);
          free(Z);
          free(Z2);
      } else {
          psio_address addr = PSIO_ZERO;
          double * temp1 = (double*)malloc(( nQ*v > o*v*v ? nQ*v : o*v*v)*sizeof(double));
          double * temp2 = (double*)malloc(o*v*v*sizeof(double));
          boost::shared_ptr<PSIO> psio(new PSIO());
          psio->open(PSIF_DCC_ABCI4,PSIO_OPEN_NEW);
          for (long int a = 0; a < v; a++) {
              #pragma omp parallel for schedule (static)
              for (long int q = 0; q < nQ; q++) {
                  for (long int c = 0; c < v; c++) {
                      temp1[q*v+c] = Qvv[q*v*v+a*v+c];
                  }
              }
              F_DGEMM('n','t',o*v,v,nQ,1.0,Qov,o*v,temp1,v,0.0,temp2,o*v);
              #pragma omp parallel for schedule (static)
              for (long int b = 0; b < v; b++) {
                  for (long int i = 0; i < o; i++) {
                      for (long int c = 0; c < v; c++) {
                          temp1[b*o*v+i*v+c] = temp2[c*o*v+i*v+b];
                      }
                  }
              }
              psio->write(PSIF_DCC_ABCI4,"E2abci4",(char*)&temp1[0],o*v*v*sizeof(double),addr,&addr);
          }
          psio->close(PSIF_DCC_ABCI4,1);
          free(temp1);
          free(temp2);
      }
      free(Qvv);

      double * temp1 = (double*)malloc(o*o*v*v*sizeof(double));
      double * temp2 = (double*)malloc(o*o*v*v*sizeof(double));

      // write (oo|ov) integrals, formerly E2ijak, for (t)
      F_DGEMM('n','t',o*o,o*v,nQ,1.0,Qoo,o*o,Qov,o*v,0.0,temp1,o*o);
      for (int i=0; i<o; i++){
          for (int j=0; j<o; j++){
              for (int k=0; k<o; k++){
                  for (int a=0; a<v; a++){
                      temp2[j*o*o*v+i*o*v+k*v+a] = temp1[i*o*o*v+a*o*o+j*o+k];
                  }
              }
          }
      }
      boost::shared_ptr<PSIO> psio(new PSIO());
      psio->open(PSIF_DCC_IJAK,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&temp2[0],o*o*o*v*sizeof(double));
      psio->close(PSIF_DCC_IJAK,1);

      // df (ov|ov) integrals, formerly E2klcd
      F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,temp1,o*v);
      psio->open(PSIF_DCC_IAJB,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&temp1[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_IAJB,1);

      free(Qov);
      free(Qoo);
      free(temp1);
      free(temp2);

      // triples
      tstart();

      ccmethod = 0;
      if (isLowMemory) status = lowmemory_triples();
      else             status = triples();

      if (status == Failure){
         throw PsiException(
            "Whoops, the (T) correction died.",__FILE__,__LINE__);
      }
      tstop();

      // ccsd(t) energy
      Process::environment.globals["(T) CORRECTION ENERGY"] = et;
      Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = eccsd + et;
      Process::environment.globals["CCSD(T) TOTAL ENERGY"] = eccsd + et + escf;
      Process::environment.globals["CURRENT ENERGY"] = eccsd + et + escf;
  }else {
      free(Qoo);
      free(Qov);
      free(Qvv);
  }

  // free remaining memory
  free(Fia);
  free(Fai);
  free(t1);
  free(tb);

  return Process::environment.globals["CURRENT ENERGY"];
}

void DFCoupledCluster::WriteBanner(){
  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                       DF-CCSD                       *\n");
  fprintf(outfile, "        *                 Density-fitted CCSD                 *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                   Eugene DePrince                   *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n\n");
  fflush(outfile);
}

/*===================================================================

  solve ccsd equations

===================================================================*/
PsiReturnType DFCoupledCluster::CCSDIterations() {

  long int o     = ndoccact;
  long int v     = nvirt;

  //int iter              = 0;
  iter                  = 0;
  int diis_iter         = 0;
  int replace_diis_iter = 1;
  double nrm            = 1.0;
  double Eold           = 1.0e9;
  if (brueckner_iter == 0 ) eccsd = 0.0;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  // zero residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_NEW);
  memset((void*)tempt,'\0',o*o*v*v*sizeof(double));
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  // start timing the iterations
  const long clk_tck = sysconf(_SC_CLK_TCK);
  struct tms total_tmstime;
  times(&total_tmstime);

  time_t time_start = time(NULL);
  double user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  double sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;

  bool timer = options_.get_bool("CC_TIMINGS");

  memset((void*)Fij,'\0',o*o*sizeof(double));
  memset((void*)Fia,'\0',o*v*sizeof(double));
  memset((void*)Fai,'\0',o*v*sizeof(double));
  memset((void*)Fab,'\0',v*v*sizeof(double));

  T1Fock();
  T1Integrals();

  fprintf(outfile,"\n");
  fprintf(outfile,"  Begin singles and doubles coupled cluster iterations\n\n");
  fprintf(outfile,"   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");
  fflush(outfile);

  memset((void*)diisvec,'\0',(maxdiis+1)*sizeof(double));
  while(iter < maxiter) {
      time_t iter_start = time(NULL);

      // evaluate cc diagrams
      memset((void*)w1,'\0',o*v*sizeof(double));
      if (iter > 0 || brueckner_iter > 0){
         CCResidual();
      }

      // update the amplitudes
      UpdateT2();
      UpdateT1();

      // add vector to list for diis
      DIISOldVector(iter,diis_iter,replace_diis_iter);

      // diis error vector and convergence check
      nrm = DIISErrorVector(diis_iter,replace_diis_iter,iter);

      // diis extrapolation
      if (diis_iter>2) {
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,o*o*v*v+o*v);
         else                   DIIS(diisvec,maxdiis,o*o*v*v+o*v);
         DIISNewAmplitudes(diis_iter,replace_diis_iter);
      }

      double start;
      if (timer) start = omp_get_wtime();
      T1Fock();
      T1Integrals();
      if (timer) {
          fprintf(outfile,"        T1-transformed integrals                                        %6.2lf\n",omp_get_wtime() - start);
          fprintf(outfile,"\n");
      }

      Eold = eccsd;
      eccsd = CheckEnergy();

      if (diis_iter < maxdiis ) {
         replace_diis_iter++;
      }else {
          double min = 1.0e9;
          for (int j = 1; j <= (diis_iter < maxdiis ? diis_iter : maxdiis); j++) {
              if ( fabs( diisvec[j-1] ) < min ) {
                  min = fabs( diisvec[j-1] );
                  replace_diis_iter = j;
              }
          }
      }

      if (diis_iter<=maxdiis) diis_iter++;
      //if (replace_diis_iter<maxdiis) replace_diis_iter++;
      //else replace_diis_iter = 1;

      time_t iter_stop = time(NULL);
      fprintf(outfile,"  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,eccsd,eccsd-Eold,nrm,(int)iter_stop-(int)iter_start);
      fflush(outfile);
      iter++;
      if (iter==1){
         emp2 = eccsd;
         SCS_MP2();
      }

      // energy and amplitude convergence check
      if (fabs(eccsd - Eold) < e_conv && nrm < r_conv) break;
  }

  times(&total_tmstime);
  time_t time_stop = time(NULL);
  double user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  double sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;

  if (iter==maxiter){
     throw PsiException("  CCSD iterations did not converge.",__FILE__,__LINE__);
  }
  if (reference_wavefunction_->isCIM()){
      Local_SCS_CCSD();
  }else{
      SCS_CCSD();
  }

  fprintf(outfile,"\n");
  fprintf(outfile,"  CCSD iterations converged!\n");
  fprintf(outfile,"\n");

  // T1 and D1 diagnostics:

  double t1diag = F_DNRM2(o*v,t1,1) / sqrt(2.0 * o);
  fprintf(outfile,"        T1 diagnostic:                  %20.12lf\n",t1diag);
  boost::shared_ptr<Matrix>T (new Matrix(o,o));
  boost::shared_ptr<Matrix>eigvec (new Matrix(o,o));
  boost::shared_ptr<Vector>eigval (new Vector(o));
  double ** Tp = T->pointer();
  for (int i = 0; i < o; i++) {
      for (int j = 0; j < o; j++) {
          double dum = 0.0;
          for (int a = 0; a < v; a++) {
              dum += t1[a*o+i] * t1[a*o+j];
          }
          Tp[i][j] = dum;
      }
  }
  T->diagonalize(eigvec,eigval,descending);
  fprintf(outfile,"        D1 diagnostic:                  %20.12lf\n",sqrt(eigval->pointer()[0]));
  fprintf(outfile,"\n");

  // delta mp2 correction for fno computations:
  if (options_.get_bool("NAT_ORBS")){
      double delta_emp2 = Process::environment.globals["MP2 CORRELATION ENERGY"] - emp2;
      double delta_emp2_os = Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] - emp2_os;
      double delta_emp2_ss = Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] - emp2_ss;

      emp2    += delta_emp2;
      emp2_os += delta_emp2_os;
      emp2_ss += delta_emp2_ss;

      eccsd    += delta_emp2;
      eccsd_os += delta_emp2_os;
      eccsd_ss += delta_emp2_ss;

      fprintf(outfile,"        OS MP2 FNO correction:          %20.12lf\n",delta_emp2_os);
      fprintf(outfile,"        SS MP2 FNO correction:          %20.12lf\n",delta_emp2_ss);
      fprintf(outfile,"        MP2 FNO correction:             %20.12lf\n",delta_emp2);
      fprintf(outfile,"\n");
  }

  if (options_.get_bool("SCS_MP2")){
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
  if (options_.get_bool("SCS_CCSD")){
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

  if (options_.get_bool("COMPUTE_TRIPLES")){
      // need to generate non-t1-transformed 3-index integrals
      F_DCOPY(o*v,t1,1,w1,1);
      memset((void*)t1,'\0',o*v*sizeof(double));
      T1Fock();
      T1Integrals();
      F_DCOPY(o*v,w1,1,t1,1);
  }

  return Success;
}

// t1-transformed 3-index fock matrix (using 3-index integrals from SCF)
void DFCoupledCluster::T1Fock(){
    long int o = ndoccact;
    long int v = nvirt;
    long int full = o+v+nfzc+nfzv;

    // Ca_L = C(1-t1^T)
    // Ca_R = C(1+t1)
    double * Catemp = (double*)malloc(nso*full*sizeof(double));
    if ( reference_wavefunction_->isCIM() ) {
        boost::shared_ptr<PSIO> psio (new PSIO());
        psio->open(PSIF_CIM,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_CIM,"C matrix",(char*)&Catemp[0],nso*full*sizeof(double));
        psio->close(PSIF_CIM,1);
        F_DCOPY(nso*full,&Catemp[0],1,Ca_L,1);
        F_DCOPY(nso*full,&Catemp[0],1,Ca_R,1);
    }else {
        F_DCOPY(nso*full,&Ca[0][0],1,Ca_L,1);
        F_DCOPY(nso*full,&Ca[0][0],1,Ca_R,1);
        F_DCOPY(nso*full,&Ca[0][0],1,Catemp,1);
    }

    #pragma omp parallel for schedule (static)
    for (int mu = 0; mu < nso; mu++) {
        for (int a = 0; a < v; a++) {
            double dum = 0.0;
            for (int i = 0; i < o; i++) {
                dum += Catemp[mu*full+i+nfzc] * t1[a*o+i];
            }
            Ca_L[mu*full + a + ndocc] -= dum;
        }
    }
    #pragma omp parallel for schedule (static)
    for (int mu = 0; mu < nso; mu++) {
        for (int i = 0; i < o; i++) {
            double dum = 0.0;
            for (int a = 0; a < v; a++) {
                dum += Catemp[mu*full+a+ndocc] * t1[a*o+i];
            }
            Ca_R[mu*full + i + nfzc] += dum;
        }
    }
    free(Catemp);

    // (Q|rs)
    boost::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO,"Qso SCF",(char*)&integrals[0],nso*nso*nQ_scf*sizeof(double));
    psio->close(PSIF_DCC_QSO,1);
    F_DGEMM('n','n',full,nso*nQ_scf,nso,1.0,Ca_L,full,integrals,nso,0.0,Qmo,full);
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ_scf; q++) {
        for (int mu = 0; mu < nso; mu++) {
            F_DCOPY(full,Qmo+q*nso*full+mu*full,1,integrals+q*nso*full+mu,full);
        }
    }
    F_DGEMM('n','n',full,full*nQ_scf,nso,1.0,Ca_R,full,integrals,nso,0.0,Qmo,full);

    // build Fock matrix

    memset((void*)Fij,'\0',o*o*sizeof(double));
    memset((void*)Fia,'\0',o*v*sizeof(double));
    memset((void*)Fai,'\0',o*v*sizeof(double));
    memset((void*)Fab,'\0',v*v*sizeof(double));

    // transform H
    double ** hp = H->pointer();
    double * h = (double*)malloc(nmo*nmo*sizeof(double));
    for (int mu = 0; mu < nso; mu++) {
        for (int p = 0; p < nmo; p++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso; nu++) {
                dum += Ca_L[nu*full + p + nfzc] * hp[nu][mu];
            }
            integrals[p*nso+mu] = dum;
        }
    }
    for (int p = 0; p < nmo; p++) {
        for (int q = 0; q < nmo; q++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso; nu++) {
                dum += Ca_R[nu*full+q+nfzc] * integrals[p*nso+nu];
            }
            h[p*nmo+q] = dum;
        }
    }

    // build Fock matrix:  sum_k (Q|kk)
    double * temp2 = (double*)malloc(nQ_scf*sizeof(double));
    double * temp3 = (double*)malloc(full*full*sizeof(double));
    for (int q = 0; q < nQ_scf; q++) {
        double dum = 0.0;
        for (int k = 0; k < ndocc; k++) {
            dum += Qmo[q*full*full + k*full + k];
        }
        temp2[q] = 2.0 * dum;
    }

    // use temp and Qvv as storage for Qmo( q, r, k) and Qmo( q, k, s)
    for (int q = 0; q < nQ_scf; q++) {
        for (int r = 0; r < full; r++) {
            for (int k = 0; k < ndocc; k++) {
                integrals[q*full*ndocc+k*full+r] = Qmo[q*full*full+r*full+k];
                Qvv[q*full*ndocc+k*full+r]  = Qmo[q*full*full+k*full+r];
            }
        }
    }
    // I(r,s) = sum q k (q|ks)(q|rk)
    F_DGEMM('n','t',full,full,nQ_scf*ndocc,1.0,Qvv,full,integrals,full,0.0,temp3,full);

    // Fij
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            double dum = h[i*nmo+j] - temp3[(i+nfzc)*full+(j+nfzc)];
            dum += F_DDOT(nQ_scf,temp2,1,Qmo + (i+nfzc)*full + (j+nfzc) , full*full);
            Fij[i*o+j] = dum;
//if (i==j) printf("%20.12lf %20.12lf\n",Fij[i*o+j],eps[i]);
        }
    }
//printf("\n\n");
//exit(0);
    // Fia
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            double dum = h[i*nmo+a+o] - temp3[(i+nfzc)*full+(a+ndocc)];
            dum += F_DDOT(nQ_scf,temp2,1,Qmo + (i+nfzc)*full + (a+ndocc) , full*full);
            Fia[i*v+a] = dum;
        }
    }
    // Fai
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            double dum = h[(a+o)*nmo+i] - temp3[(a+ndocc)*full+(i+nfzc)];
            dum += F_DDOT(nQ_scf,temp2,1,Qmo + (a+ndocc)*full + (i+nfzc) , full*full);
            Fai[a*o+i] = dum;
        }
    }
    // Fab
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            double dum = h[(a+o)*nmo+b+o] - temp3[(a+ndocc)*full+(b+ndocc)];
            dum += F_DDOT(nQ_scf,temp2,1,Qmo + (a+ndocc)*full + (b+ndocc) , full*full);
            Fab[a*v+b] = dum;
        }
    }

    // replace eps
    for (int i = 0; i < o; i++) {
        eps[i] = Fij[i*o+i];
    }
    for (int a = 0; a < v; a++) {
        eps[a+o] = Fab[a*v+a];
    }

    free(h);
    free(temp2);
    free(temp3);
}
void DFCoupledCluster::T1Integrals(){
    long int o = ndoccact;
    long int v = nvirt;
    long int full = o+v+nfzc+nfzv;

    // Ca_L = C(1-t1^T)
    // Ca_R = C(1+t1)
    double * Catemp = (double*)malloc(nso*full*sizeof(double));
    if ( reference_wavefunction_->isCIM() ) {
        boost::shared_ptr<PSIO> psio (new PSIO());
        psio->open(PSIF_CIM,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_CIM,"C matrix",(char*)&Catemp[0],nso*full*sizeof(double));
        psio->close(PSIF_CIM,1);
        F_DCOPY(nso*full,&Catemp[0],1,Ca_L,1);
        F_DCOPY(nso*full,&Catemp[0],1,Ca_R,1);
    }else {
        F_DCOPY(nso*full,&Ca[0][0],1,Ca_L,1);
        F_DCOPY(nso*full,&Ca[0][0],1,Ca_R,1);
        F_DCOPY(nso*full,&Ca[0][0],1,Catemp,1);
    }

    #pragma omp parallel for schedule (static)
    for (int mu = 0; mu < nso; mu++) {
        for (int a = 0; a < v; a++) {
            double dum = 0.0;
            for (int i = 0; i < o; i++) {
                dum += Catemp[mu*full+i+nfzc] * t1[a*o+i];
            }
            Ca_L[mu*full + a + ndocc] -= dum;
        }
    }
    #pragma omp parallel for schedule (static)
    for (int mu = 0; mu < nso; mu++) {
        for (int i = 0; i < o; i++) {
            double dum = 0.0;
            for (int a = 0; a < v; a++) {
                dum += Catemp[mu*full+a+ndocc] * t1[a*o+i];
            }
            Ca_R[mu*full + i + nfzc] += dum;
        }
    }
    free(Catemp);

    // (Q|rs)
    boost::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO,"Qso CC",(char*)&integrals[0],nso*nso*nQ*sizeof(double));
    psio->close(PSIF_DCC_QSO,1);
    F_DGEMM('n','n',full,nso*nQ,nso,1.0,Ca_L,full,integrals,nso,0.0,Qmo,full);
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        for (int mu = 0; mu < nso; mu++) {
            F_DCOPY(full,Qmo+q*nso*full+mu*full,1,integrals+q*nso*full+mu,full);
        }
    }
    F_DGEMM('n','n',full,full*nQ,nso,1.0,Ca_R,full,integrals,nso,0.0,Qmo,full);

    // (Q|vo)
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                Qov[q*o*v+a*o+i] = Qmo[q*full*full+(a+ndocc)*full+(i+nfzc)];
            }
        }
    }
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio->write_entry(PSIF_DCC_QSO,"qvo",(char*)&Qov[0],o*v*nQ*sizeof(double));
    psio->close(PSIF_DCC_QSO,1);
    // (Q|ov)
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        for (int i = 0; i < o; i++) {
            for (int a = 0; a < v; a++) {
                Qov[q*o*v+i*v+a] = Qmo[q*full*full+(i+nfzc)*full+(a+ndocc)];
            }
        }
    }
    // (Q|oo)
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        for (int i = 0; i < o; i++) {
            for (int j = 0; j < o; j++) {
                Qoo[q*o*o+i*o+j] = Qmo[q*full*full+(i+nfzc)*full+(j+nfzc)];
            }
        }
    }
    // (Q|vv)
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < v; a++) {
            for (int b = 0; b < v; b++) {
                Qvv[q*v*v+a*v+b] = Qmo[q*full*full+(a+ndocc)*full+(b+ndocc)];
            }
        }
    }
}

double DFCoupledCluster::CheckEnergy(){

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
void DFCoupledCluster::SCS_MP2(){

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
  emp2    = emp2_os + emp2_ss;
}
void DFCoupledCluster::SCS_CCSD(){

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
  eccsd    = eccsd_os + eccsd_ss;
}

/*===================================================================

  allocate cpu memory

===================================================================*/
void DFCoupledCluster::AllocateMemory() {

  if (nirrep_>1){
     throw PsiException("df_ccsd requires symmetry c1",__FILE__,__LINE__);
  }

  ischolesky_ = ( options_.get_str("DF_BASIS_CC") == "CHOLESKY" );
  nQ          = (int)Process::environment.globals["NAUX (CC)"];
  nQ_scf      = (int)Process::environment.globals["NAUX (SCF)"];

  // orbital energies
  //boost::shared_ptr<Vector> eps_test = reference_wavefunction_->epsilon_a();
  //double * tempeps                   = eps_test->pointer();
  //eps                                = (double*)malloc(nmo*sizeof(double));
  //F_DCOPY(nmo,tempeps+nfzc,1,eps,1);

  if (!reference_wavefunction_->isCIM()){
      int count=0;
      eps = (double*)malloc((ndoccact+nvirt)*sizeof(double));
      boost::shared_ptr<Vector> eps_test = reference_wavefunction_->epsilon_a();
      for (int h=0; h<nirrep_; h++){
          for (int norb = frzcpi_[h]; norb<doccpi_[h]; norb++){
              eps[count++] = eps_test->get(h,norb);
          }
      }
      for (int h=0; h<nirrep_; h++){
          for (int norb = doccpi_[h]; norb<nmopi_[h]-frzvpi_[h]; norb++){
              eps[count++] = eps_test->get(h,norb);
          }
      }
  }else{
     // orbital energies in qt ordering:
     long int count = 0;
     eps = (double*)malloc((ndoccact+nvirt)*sizeof(double));
     boost::shared_ptr<Vector> eps_test = reference_wavefunction_->CIMOrbitalEnergies();
     for (int i = 0; i < ndoccact + nvirt; i++){
         eps[i] = eps_test->get(0,i+nfzc);
     }
     eps_test.reset();
  }

  long int o = ndoccact;
  long int v = nvirt;

  /*========================================================
     ccsd memory requirements:
    
     tb:     o^2v^2
     tempt:  o^2v^2+ov ( actually o(o+1)v(v+1) + ov )
     tempv:  max (o^2v^2+ov , o*v*nQ)
     integrals: max(2v^3,nQ*nso^2, o^2v^2, 2v^3, 2nQ*o*v) (this is a minimum)
     Abij (SJS v^4 result): o(o+1)v/2
     Sbij (SJS v^4 result): o(o+1)v/2
     other stuff: 2ov+2v^2+(o+v)
    
     total: 3o^2v^2 + 2v^3  + o(o+1)v + 4ov  + 2v^2 + (o+v)  or 
            4o^2v^2         + o(o+1)v + 4ov  + 2v^2 + (o+v)  or
            3o^2v^2 + 2ovnQ + o(o+1)v + 4ov  + 2v^2 + (o+v)

     compare to the requirements for the (T) part:

            2o^2v^2 + 3v^3*nthreads + o^3v + ov
    
  ========================================================*/


  // for the df version, the dimension of the large buffer:
  long int nQmax = nQ > nQ_scf ? nQ : nQ_scf;

  long int dim = 2L*v*v*v;
  if (2*nQmax*o*v>dim)   dim = 2*nQmax*o*v;
  if (o*o*v*v>dim)       dim = o*o*v*v;
  if (nQmax*v*v>dim)     dim = nQmax*v*v;
  if (nQmax*nso*nso>dim) dim = nQmax*nso*nso;

  double total_memory = dim+(o*o*v*v+o*v)+(o*(o+1)*v*(v+1)+o*v)+o*o*v*v+2.*o*v+2.*v*v;
  long int max = nvirt*nvirt*nQmax > (nfzv+ndocc+nvirt)*ndocc*nQmax ? nvirt*nvirt*nQmax : (nfzv+ndocc+nvirt)*ndocc*nQmax;
  double df_memory    = nQ*(o*o+o*v)+max;

  total_memory       *= 8./1024./1024.;
  df_memory          *= 8./1024./1024.;

  fprintf(outfile,"  ==> Memory <==\n\n");
  fprintf(outfile,"        Total memory requirements:       %9.2lf mb\n",df_memory+total_memory);
  fprintf(outfile,"        3-index integrals:               %9.2lf mb\n",df_memory);
  fprintf(outfile,"        CCSD intermediates:              %9.2lf mb\n",total_memory);
  fprintf(outfile,"\n");

  if (1.0 * memory / 1024. / 1024. < total_memory + df_memory) {
     fprintf(outfile,"\n");
     fprintf(outfile,"        error: not enough memory for ccsd.  increase available memory by %7.2lf mb\n",total_memory+df_memory-1.0*memory/1024./1024.);
     fprintf(outfile,"\n");
     fflush(outfile);
     throw PsiException("not enough memory (ccsd).",__FILE__,__LINE__);
  }
  if (options_.get_bool("COMPUTE_TRIPLES")) {
      long int nthreads = omp_get_max_threads();
      double tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+3L*v*v*v*nthreads);
      if (tempmem > memory) {
          fprintf(outfile,"\n        <<< warning! >>> switched to low-memory (t) algorithm\n\n");
      }
      if (tempmem > memory || options_.get_bool("TRIPLES_LOW_MEMORY")){
         isLowMemory = true;
         tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+5L*o*o*o*nthreads);
      }
      fprintf(outfile,"        memory requirements for CCSD(T): %9.2lf mb\n\n",tempmem/1024./1024.);
  }

  // allocate some memory for 3-index tensors
  Qoo = (double*)malloc(ndoccact*ndoccact*nQ*sizeof(double));
  Qov = (double*)malloc(ndoccact*nvirt*nQ*sizeof(double));
  // max (v*v*nQ, full*ndocc*nQ)
  Qvv = (double*)malloc(max*sizeof(double));

  integrals = (double*)malloc(dim*sizeof(double));
  tempt     = (double*)malloc((o*(o+1)*v*(v+1)+o*v)*sizeof(double));
  long int tempvdim = o*o*v*v+o*v;
  if ( nQ * o * v > tempvdim) tempvdim = nQ * o * v;
  tempv     = (double*)malloc(tempvdim*sizeof(double));
  Abij      = (double*)malloc(o*(o+1)/2*v*sizeof(double));
  Sbij      = (double*)malloc(o*(o+1)/2*v*sizeof(double));
  tb        = (double*)malloc(o*o*v*v*sizeof(double));
  w1        = (double*)malloc(o*v*sizeof(double));
  t1        = (double*)malloc(o*v*sizeof(double));
  I1        = (double*)malloc(v*v*sizeof(double));
  I1p       = (double*)malloc(v*v*sizeof(double));

  memset((void*)integrals,'\0',dim*sizeof(double));
  memset((void*)tempv,'\0',tempvdim*sizeof(double));
  memset((void*)tempt,'\0',(o*(o+1)*v*(v+1)+o*v)*sizeof(double));
  memset((void*)tb,'\0',o*o*v*v*sizeof(double));
  memset((void*)w1,'\0',o*v*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));
  memset((void*)I1,'\0',v*v*sizeof(double));
  memset((void*)I1p,'\0',v*v*sizeof(double));
  memset((void*)Abij,'\0',o*(o+1)/2*v*sizeof(double));
  memset((void*)Sbij,'\0',o*(o+1)/2*v*sizeof(double));

  // DIIS:
  diisvec    = (double*)malloc(sizeof(double)*(maxdiis+1));
  memset((void*)diisvec,'\0',(maxdiis+1)*sizeof(double));

  // new 3-index stuff for t1-transformed integrals:
  Fij   = (double*)malloc(o*o*sizeof(double));
  Fia   = (double*)malloc(o*v*sizeof(double));
  Fai   = (double*)malloc(o*v*sizeof(double));
  Fab   = (double*)malloc(v*v*sizeof(double));
  Qmo   = (double*)malloc(nso*nso*(nQ > nQ_scf ? nQ : nQ_scf)*sizeof(double));
  Ca_R  = (double*)malloc(nso*(nmo+nfzc+nfzv)*sizeof(double));
  Ca_L  = (double*)malloc(nso*(nmo+nfzc+nfzv)*sizeof(double));
  Ca = reference_wavefunction_->Ca()->pointer();

  // one-electron integrals
  boost::shared_ptr<MintsHelper> mints(new MintsHelper());
  H = mints->so_kinetic();
  H->add(mints->so_potential());
}

/*================================================================

   update amplitudes

================================================================*/
void DFCoupledCluster::UpdateT1(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  #pragma omp parallel for schedule (static)
  for (long int a=o; a<rs; a++){
      for (long int i=0; i<o; i++){
          double dia = -eps[i]+eps[a];
          double tnew =  -w1[(a-o)*o+i]/dia;
          w1[(a-o)*o+i] = tnew + t1[(a-o)*o+i];
      }
  }
  // error vector for diis is in tempv:
  F_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  F_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  F_DCOPY(o*v,w1,1,t1,1);
}
void DFCoupledCluster::UpdateT2(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  boost::shared_ptr<PSIO> psio(new PSIO());

  // df (ai|bj)
  psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_QSO,"qvo",(char*)&tempv[0],nQ*o*v*sizeof(double));
  psio->close(PSIF_DCC_QSO,1);
  F_DGEMM('n','t',o*v,o*v,nQ,1.0,tempv,o*v,tempv,o*v,0.0,integrals,o*v);

  // we still have the residual in memory in tempv
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
 
  #pragma omp parallel for schedule (static)
  for (long int a=o; a<rs; a++){
      double da = eps[a];
      for (long int b=o; b<rs; b++){
          double dab = da + eps[b];
          for (long int i=0; i<o; i++){
              double dabi = dab - eps[i];
              for (long int j=0; j<o; j++){

                  long int iajb = (a-o)*v*o*o+i*v*o+(b-o)*o+j;
                  long int jaib = iajb + (i-j)*v*(1-v*o);
                  long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;

                  double dijab = dabi-eps[j];
                  double tnew  = - (integrals[iajb] + tempv[ijab])/dijab;
                  tempt[ijab]  = tnew + tb[ijab];
              }
          }
      }
  }
  // error vectors for diis are in tempv:
  F_DCOPY(o*o*v*v,tempt,1,tempv,1);
  F_DAXPY(o*o*v*v,-1.0,tb,1,tempv,1);
  F_DCOPY(o*o*v*v,tempt,1,tb,1);
}

/**
 *  Use Vabcd1
 */
void DFCoupledCluster::Vabcd1(){
  long int o = ndoccact;
  long int v = nvirt;
  long int oov = o*o*v;
  long int oo  = o*o;
  long int otri = o*(o+1)/2;
  long int vtri = v*(v+1)/2;

  boost::shared_ptr<PSIO> psio(new PSIO());

  #pragma omp parallel for schedule (static)
  for (long int i=0; i<o; i++){
      for (long int j=i; j<o; j++){
          long int ij = Position(i,j);
          for (long int a=0; a<v; a++){
              for (long int b=a; b<v; b++){
                  tempt[Position(a,b)*otri+ij] =
                     (tb[a*oov+b*oo+i*o+j]+tb[b*oov+a*oo+i*o+j]);
                  tempt[Position(a,b)*otri+ij+vtri*otri] =
                     (tb[a*oov+b*oo+i*o+j]-tb[b*oov+a*oo+i*o+j]);
              }
              tempt[Position(a,a)*otri+ij] = tb[a*oov+a*oo+i*o+j];
          }
      }
  }

  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));

  int nthreads = omp_get_max_threads();

  double * Iqdb = Qmo;
  double * Vp   = Qmo;
  double * Vcdb = integrals;
  double * Vm   = integrals+v*v*v;

  // qvv transpose
  #pragma omp parallel for schedule (static)
  for (int q = 0; q < nQ; q++) {
      F_DCOPY(v*v,Qvv+q*v*v,1,integrals+q,nQ);
  }
  F_DCOPY(nQ*v*v,integrals,1,Qvv,1);

  double time1 = 0.0;
  double time2 = 0.0;
  double time3 = 0.0;
  for (long int a = 0; a < v; a++) {

      int nb = 0;
      // fill Iqdb for b > a
      double start1 = omp_get_wtime();
      //#pragma omp parallel for schedule (static)
      //for (long int b = a; b < v; b++) {
      //    F_DCOPY(nQ*v,Qvv+b*nQ*v,1,Iqdb+(b-a)*nQ*v,1);
      //}
      nb = v-a;

      //F_DGEMM('t','n',v,v*nb,nQ,1.0,Qvv+a*v*nQ,nQ,Iqdb,nQ,0.0,Vcdb,v);
      F_DGEMM('t','n',v,v*nb,nQ,1.0,Qvv+a*v*nQ,nQ,Qvv+a*v*nQ,nQ,0.0,Vcdb,v);

      #pragma omp parallel for schedule (static)
      for (long int b = a; b < v; b++){
          long int cd = 0;
          long int ind1 = (b-a)*vtri;
          long int ind2 = (b-a)*v*v;
          long int v1,v2;
          for (long int c=0; c<v; c++){
              for (long int d=0; d<=c; d++){
                  Vp[ind1+cd] = Vcdb[ind2+d*v+c] + Vcdb[ind2+c*v+d];
                  Vm[ind1+cd] = Vcdb[ind2+d*v+c] - Vcdb[ind2+c*v+d];
                  cd++;
              }
          }
      }
      double end1 = omp_get_wtime();

      double start2 = omp_get_wtime();
      F_DGEMM('n','n',otri,nb,vtri,0.5,tempt,otri,Vp,vtri,0.0,Abij,otri);
      F_DGEMM('n','n',otri,nb,vtri,0.5,tempt+otri*vtri,otri,Vm,vtri,0.0,Sbij,otri);
      double end2 = omp_get_wtime();

      // contribute to residual
      double start3 = omp_get_wtime();
      #pragma omp parallel for schedule (static)
      for (long int b = a; b < v; b++) {
          for (long int i = 0; i < o; i++) {
              for (long int j = 0; j < o; j++) {
                  int sg = ( i > j ) ? 1 : -1;
                  tempv[a*oo*v+b*oo+i*o+j]    +=    Abij[(b-a)*otri+Position(i,j)]
                                              +  sg*Sbij[(b-a)*otri+Position(i,j)];
                  if (a!=b) {
                     tempv[b*oov+a*oo+i*o+j] +=    Abij[(b-a)*otri+Position(i,j)]
                                             -  sg*Sbij[(b-a)*otri+Position(i,j)];
                  }
              }
          }
      }
      double end3 = omp_get_wtime();

      time1 += end1 - start1;
      time2 += end2 - start2;
      time3 += end3 - start3;
  }

  // contribute to residual
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  // qvv un-transpose
  #pragma omp parallel for schedule (static)
  for (int q = 0; q < nQ; q++) {
      F_DCOPY(v*v,Qvv+q,nQ,integrals+q*v*v,1);
  }
  F_DCOPY(nQ*v*v,integrals,1,Qvv,1);
}

void DFCoupledCluster::CCResidual(){
    bool timer = options_.get_bool("CC_TIMINGS");
    long int o = ndoccact;
    long int v = nvirt;

    double start;

    // C2 = -1/2 t(bc,kj) [ (ki|ac) - 1/2 t(ad,li) (kd|lc) ] 
    //      +    t(bc,ki) [ (kj|ac) - 1/2 t(ad,lj) (kd|lc) ] 
    if (timer) start = omp_get_wtime();
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int l = 0; l < o; l++) {
                for (int d = 0; d < v; d++) {
                    tempt[a*o*o*v+i*o*v+l*v+d] = tb[a*o*o*v+d*o*o+l*o+i];
                }
            }
        }
    }
    #pragma omp parallel for schedule (static)
    for (int l = 0; l < o; l++) {
        for (int d = 0; d < v; d++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempv[l*o*v*v+d*o*v+k*v+c] = integrals[k*o*v*v+d*o*v+l*v+c];
                }
            }
        }
    }
    F_DGEMM('n','n',o*v,o*v,o*v,-0.5,tempv,o*v,tempt,o*v,0.0,integrals,o*v);
    F_DGEMM('n','t',v*v,o*o,nQ,1.0,Qvv,v*v,Qoo,o*o,0.0,tempv,v*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    integrals[a*o*o*v+i*o*v+k*v+c] += tempv[k*o*v*v+i*v*v+a*v+c];
                }
            }
        }
    }
    #pragma omp parallel for schedule (static)
    for (int b = 0; b < v; b++) {
        for (int j = 0; j < o; j++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempt[b*o*o*v+j*o*v+k*v+c] = tb[b*o*o*v+c*o*o+k*o+j];
                }
            }
        }
    }
    F_DGEMM('t','n',o*v,o*v,o*v,-1.0,integrals,o*v,tempt,o*v,0.0,tempv,o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[a*o*o*v+b*o*o+i*o+j] = 0.5 * tempv[b*o*o*v+j*o*v+a*o+i] + tempv[b*o*o*v+i*o*v+a*o+j];
                }
            }
        }
    }

    // first contribution to residual 
    boost::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_R2,PSIO_OPEN_NEW);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        fprintf(outfile,"\n");
        fprintf(outfile,"        C2 = -1/2 t(b,c,k,j) [ (ki|ac) - 1/2 t(a,d,l,i) (kd|lc) ]\n");
        fprintf(outfile,"                + t(b,c,k,i) [ (kj|ac) - 1/2 t(a,d,l,j) (kd|lc) ]       %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // D2: 1/2 U(b,c,j,k) [ L(a,i,k,c) + 1/2 U(a,d,i,l) L(l,d,k,c) ] 
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    F_DCOPY(o*o*v*v,integrals,1,tempv,1);
    #pragma omp parallel for schedule (static)
    for (int l = 0; l < o; l++) {
        for (int d = 0; d < v; d++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempv[l*o*v*v+d*o*v+k*v+c] -= 0.5 * integrals[l*o*v*v+c*o*v+k*v+d];
                }
            }
        }
    }
    #pragma omp parallel for schedule (static)
    for (int l = 0; l < o; l++) {
        for (int d = 0; d < v; d++) {
            for (int a = 0; a < v; a++) {
                for (int i = 0; i < o; i++) {
                    tempt[l*o*v*v+d*o*v+a*o+i] = 2.0 * tb[a*o*o*v+d*o*o+i*o+l]-tb[a*o*o*v+d*o*o+l*o+i];
                }
            }
        }
    }
    F_DGEMM('n','t',o*v,o*v,o*v,1.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO,"qvo",(char*)&tempv[0],nQ*o*v*sizeof(double));
    psio->close(PSIF_DCC_QSO,1);
    F_DGEMM('n','t',o*v,o*v,nQ,2.0,Qov,o*v,tempv,o*v,1.0,integrals,o*v);
    F_DGEMM('n','t',o*o,v*v,nQ,-1.0,Qoo,o*o,Qvv,v*v,0.0,tempv,o*o);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    integrals[a*o*o*v+i*o*v+k*v+c] += tempv[a*o*o*v+c*o*o+k*o+i];
                }
            }
        }
    }
    #pragma omp parallel for schedule (static)
    for (int k = 0; k < o; k++) {
        for (int c = 0; c < v; c++) {
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {
                    tempt[k*o*v*v+c*o*v+b*o+j] = 2.0 * tb[b*o*o*v+c*o*o+j*o+k] - tb[b*o*o*v+c*o*o+k*o+j];
                }
            }
        }
    }
    F_DGEMM('n','n',o*v,o*v,o*v,0.5,tempt,o*v,integrals,o*v,0.0,tempv,o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[a*o*o*v+b*o*o+i*o+j] = tempv[a*o*o*v+i*o*v+b*o+j];
                }
            }
        }
    }
    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
    F_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        fprintf(outfile,"        D2 =  1/2 U(b,c,j,k) [ L(a,i,k,c) + 1/2 U(a,d,i,l) L(l,d,k,c) ] %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // E2 a: t(ac,ij) [ F(bc) - U(bd,kl) (ld|kc) ]
    F_DCOPY(o*o*v*v,tb,1,tempt,1);
    #pragma omp parallel for schedule (static)
    for (int b = 0; b < v; b++) {
        for (int d = 0; d < v; d++) {
            for (int k = 0; k < o; k++) {
                F_DAXPY(o,-0.5,tb+b*o*o*v+d*o*o+k,o,tempt+b*o*o*v+d*o*o+k*o,1);
            }
        }
    }
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    #pragma omp parallel for schedule (static)
    for (int c = 0; c < v; c++) {
        for (int d = 0; d < v; d++) {
            for (int k = 0; k < o; k++) {
                for (int l = 0; l < o; l++) {
                    tempv[c*o*o*v+d*o*o+k*o+l] = integrals[l*o*v*v+d*o*v+k*v+c];
                }
            }
        }
    }
    // overwriting Fab here, but it gets rebuilt every iteration anyway.
    F_DGEMM('t','n',v,v,o*o*v,-2.0,tempv,o*o*v,tempt,o*o*v,1.0,Fab,v);
    #pragma omp parallel for schedule (static)
    for (int c = 0; c < v; c++) {
        for (int a = 0; a < v; a++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[c*o*o*v+a*o*o+i*o+j] = tb[a*o*o*v+c*o*o+i*o+j];
                }
            }
        }
    }
    F_DGEMM('n','n',o*o*v,v,v,1.0,tempt,o*o*v,Fab,v,0.0,tempv,o*o*v);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    tempt[a*o*o*v+b*o*o+i*o+j] = tempv[b*o*o*v+a*o*o+i*o+j];
                }
            }
        }
    }
    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
    F_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        fprintf(outfile,"        E2 =      t(a,c,i,j) [ F(b,c) - U(b,d,k,l) (ld|kc) ]            %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // E2 b: -t(a,b,i,k) [ F(kj) - U(c,d,l,j) (kd|lc) ]
    // note that (kd|lc) should still be in integrals buffer
    #pragma omp parallel for schedule (static)
    for (int j = 0; j < o; j++) {
        for (int d = 0; d < v; d++) {
            for (int l = 0; l < o; l++) {
                for (int c = 0; c < v; c++) {
                    tempt[j*o*v*v+d*o*v+l*v+c] = (2.0 * tb[c*o*o*v+d*o*o+l*o+j] - tb[c*o*o*v+d*o*o+j*o+l] );
                }
            }
        }
    }
    // overwriting Fij here, but it gets rebuilt every iteration anyway.
    F_DGEMM('t','n',o,o,o*v*v,1.0,tempt,o*v*v,integrals,o*v*v,1.0,Fij,o);

    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    F_DGEMM('n','n',o,o*v*v,o,-1.0,Fij,o,tb,o,1.0,tempt,o);

    // R2 = R2 + P(ia,jb) R2
    F_DCOPY(o*o*v*v,tempt,1,integrals,1);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    integrals[a*o*o*v+b*o*o+i*o+j] += tempt[b*o*o*v+a*o*o+j*o+i];
                }
            }
        }
    }
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        fprintf(outfile,"                - t(a,b,i,k) [ F(k,j) - U(c,d,l,j) (kd|lc) ]            %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // B2 = t(ab,kl) [ (ki|lj) + t(cd,ij) (kc|ld) ]
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    #pragma omp parallel for schedule (static)
    for (int k = 0; k < o; k++) {
        for (int l = 0; l < o; l++) {
            for (int c = 0; c < v; c++) {
                for (int d = 0; d < v; d++) {
                    tempv[k*o*v*v+l*v*v+c*v+d] = integrals[k*v*v*o+c*o*v+l*v+d];
                }
            }
        }
    }
    F_DGEMM('n','t',o*o,o*o,nQ,1.0,Qoo,o*o,Qoo,o*o,0.0,integrals,o*o);
    #pragma omp parallel for schedule (static)
    for (int k = 0; k < o; k++) {
        for (int i = 0; i < o; i++) {
            for (int l = 0; l < o; l++) {
                for (int j = 0; j < o; j++) {
                    tempt[k*o*o*o+l*o*o+i*o+j] = integrals[k*o*o*o+i*o*o+l*o+j];
                }
            }
        }
    }
    F_DGEMM('n','n',o*o,o*o,v*v,1.0,tb,o*o,tempv,v*v,1.0,tempt,o*o);
    F_DGEMM('n','n',o*o,v*v,o*o,1.0,tempt,o*o,tb,o*o,0.0,integrals,o*o);

    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    F_DAXPY(o*o*v*v,1.0,tempt,1,integrals,1);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);

    if (timer) {
        fprintf(outfile,"        B2 =      t(a,b,k,l) [ (ki|lj) + t(c,d,i,j) (kc|ld) ]           %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // now singles residual:

    // D1: F(ai)
    F_DCOPY(o*v,Fai,1,w1,1);
 
    // A1 (G):  U(c,d,k,l) (ad|kc)
    #pragma omp parallel for schedule (static)
    for (int d = 0; d < v; d++) {
        for (int i = 0; i < o; i++) {
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    tempt[d*o*o*v+i*o*v+k*v+c] = (2.0*tb[c*o*o*v+d*o*o+k*o+i] - tb[c*o*o*v+d*o*o+i*o+k]);
                }
            }
        }
    }
    F_DGEMM('t','n',o*v,nQ,o*v,1.0,tempt,o*v,Qov,o*v,0.0,tempv,o*v);
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < v; a++) {
            F_DCOPY(v,Qvv+q*v*v+a*v,1,integrals+q*v*v+a,v);
        }
    }
    F_DGEMM('n','t',o,v,v*nQ,1.0,tempv,o,integrals,v,1.0,w1,o);

    if (timer) {
        fprintf(outfile,"        A1 =      U(c,d,k,l) (ad|kc)                                    %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // B1 (H): -U(a,c,k,l) (ki|lc)
    F_DGEMM('n','t',o*v,o*o,nQ,1.0,Qov,o*v,Qoo,o*o,0.0,integrals,o*v);
    #pragma omp parallel for schedule (static)
    for (int i = 0; i < o; i++) {
        for (int c = 0; c < v; c++) {
            for (int k = 0; k < o; k++) {
                for (int l = 0; l < o; l++) {
                    tempv[i*o*o*v+c*o*o+k*o+l] = integrals[k*o*o*v+i*o*v+l*v+c];
                }
            }
        }
    }
    F_DCOPY(o*o*v*v,tb,1,tempt,1);
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int c = 0; c < v; c++) {
            for (int k = 0; k < o; k++) {
                F_DAXPY(o,-0.5,tb+a*o*o*v+c*o*o+k,o,tempt+a*o*o*v+c*o*o+k*o,1);
            }
        }
    }
    F_DGEMM('t','n',o,v,o*o*v,-2.0,tempv,o*o*v,tempt,o*o*v,1.0,w1,o);

    if (timer) {
        fprintf(outfile,"        B1 =    - U(a,c,k,l) (ki|lc)                                    %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // C1
    #pragma omp parallel for schedule (static)
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            double dum = 0.0;
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    dum += Fia[k*v+c] * (2.0*tb[a*o*o*v+c*o*o+i*o+k] - tb[a*o*o*v+c*o*o+k*o+i]);
                }
            }
            w1[a*o+i] += dum;
        }
    }
 
    if (timer) {
        fprintf(outfile,"        C1 =      F(k,c) U(a,c,i,k)                                     %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    Vabcd1();
    if (timer) {
        fprintf(outfile,"        A2 =      t(c,d,i,j) (ac|bd)                                    %6.2lf\n",omp_get_wtime()-start);
    }
}

CoupledPair::CoupledPair(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options):
        CoupledCluster(reference_wavefunction,options)
{
    reference_wavefunction_ = reference_wavefunction;
    common_init();

    // which cepa level? 0,1,2,3
    // also, -1 = cisd
    // also, -2 = acpf
    // also, -3 = aqcc
    std::string cepa = options_.get_str("CEPA_LEVEL");

    // set the wavefunction name
    name_ = cepa;

    if (cepa == "CEPA(0)") cepa_level = 0;
    if (cepa == "CEPA(1)") cepa_level = 1;
    if (cepa == "CEPA(2)") cepa_level = 2;
    if (cepa == "CEPA(3)") cepa_level = 3;
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
}

CoupledPair::~CoupledPair()
{
}

void CoupledPair::WriteBanner(){
  fflush(outfile);
  fprintf(outfile,"\n\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  if (options_.get_str("CEPA_LEVEL")=="CEPA(0)"){
     fprintf(outfile, "        *                       CEPA(0)                       *\n");
     fprintf(outfile, "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="CEPA(1)"){
     fprintf(outfile, "        *                       CEPA(1)                       *\n");
     fprintf(outfile, "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="CEPA(2)"){
     fprintf(outfile, "        *                       CEPA(2)                       *\n");
     fprintf(outfile, "        *        Coupled Electron Pair Approximation          *\n");
  }
  if (options_.get_str("CEPA_LEVEL")=="CEPA(3)"){
     fprintf(outfile, "        *                       CEPA(3)                       *\n");
     fprintf(outfile, "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="ACPF"){
     fprintf(outfile, "        *                        ACPF                         *\n");
     fprintf(outfile, "        *          Averaged Coupled Pair Functional           *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="AQCC"){
     fprintf(outfile, "        *                        AQCC                         *\n");
     fprintf(outfile, "        *         Averaged Quadratic Coupled Cluster          *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="CISD"){
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

double CoupledPair::compute_energy() {
  PsiReturnType status = Success;

  // integral sort
  tstart();
  SortIntegrals(nfzc,nfzv,nmo+nfzc+nfzv,ndoccact,nvirt,options_,reference_wavefunction_->isCIM());
  tstop();

  // solve cepa equations
  tstart();
  WriteBanner();
  AllocateMemory();
  status = CEPAIterations();
  tstop();

  // mp2 energy
  Process::environment.globals["MP2 CORRELATION ENERGY"]               = emp2;
  Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
  Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"]     = emp2_ss;
  Process::environment.globals["MP2 TOTAL ENERGY"]                     = emp2 + escf;

  // cepa energy
  char*cepatype = (char*)malloc(100*sizeof(char));
  if (cepa_level == 0){
     Process::environment.globals["CEPA(0) CORRELATION ENERGY"]               = eccsd;
     Process::environment.globals["CEPA(0) OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
     Process::environment.globals["CEPA(0) SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
     Process::environment.globals["CEPA(0) TOTAL ENERGY"]                     = eccsd + escf;
  }
  if (cepa_level == 1){
     Process::environment.globals["CEPA(1) CORRELATION ENERGY"]               = eccsd;
     Process::environment.globals["CEPA(1) OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
     Process::environment.globals["CEPA(1) SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
     Process::environment.globals["CEPA(1) TOTAL ENERGY"]                     = eccsd + escf;
  }
  if (cepa_level == 2){
     Process::environment.globals["CEPA(2) CORRELATION ENERGY"]               = eccsd;
     Process::environment.globals["CEPA(2) OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
     Process::environment.globals["CEPA(2) SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
     Process::environment.globals["CEPA(2) TOTAL ENERGY"]                     = eccsd + escf;
  }
  if (cepa_level == 3){
     Process::environment.globals["CEPA(3) CORRELATION ENERGY"]               = eccsd;
     Process::environment.globals["CEPA(3) OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
     Process::environment.globals["CEPA(3) SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
     Process::environment.globals["CEPA(3) TOTAL ENERGY"]                     = eccsd + escf;
  }
  if (cepa_level == -1){
     Process::environment.globals["CISD CORRELATION ENERGY"]               = eccsd;
     Process::environment.globals["CISD OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
     Process::environment.globals["CISD SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
     Process::environment.globals["CISD TOTAL ENERGY"]                     = eccsd + escf;
  }
  if (cepa_level == -2){
     Process::environment.globals["ACPF CORRELATION ENERGY"]               = eccsd;
     Process::environment.globals["ACPF OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
     Process::environment.globals["ACPF SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
     Process::environment.globals["ACPF TOTAL ENERGY"]                     = eccsd + escf;
  }
  if (cepa_level == -3){
     Process::environment.globals["AQCC CORRELATION ENERGY"]               = eccsd;
     Process::environment.globals["AQCC OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
     Process::environment.globals["AQCC SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
     Process::environment.globals["AQCC TOTAL ENERGY"]                     = eccsd + escf;
  }
  Process::environment.globals["CURRENT ENERGY"] = eccsd + escf;
  Process::environment.globals["CURRENT CORRELATION ENERGY"] = eccsd;

  // build opdm in case we want properties.  don't build if cim.
  if ( cepa_level<=0 && !reference_wavefunction_->isCIM() ) {
      if (options_.get_bool("NAT_ORBS")) {
          //fprintf(outfile,"\n");
          //fprintf(outfile,"\n");
          //fprintf(outfile,"       <<< Warning >>>  %s OPDM will have no correction for FNO truncation.\n",cepa_type);
          //fprintf(outfile,"\n");
      } else {
          OPDM();
      }
  }

  free(cepatype);

  finalize();

  return eccsd+escf;
}

/*===================================================================

  solve cepa equations

===================================================================*/
PsiReturnType CoupledPair::CEPAIterations(){

  long int o = ndoccact;
  long int v = nvirt;

  iter                  = 0;
  int diis_iter         = 0;
  int replace_diis_iter = 1;
  double nrm            = 1.0;
  double Eold           = 1.0e9;
  eccsd                 = 0.0;

  fprintf(outfile,"\n");
  fprintf(outfile,
    "  Begin %s iterations\n\n",cepa_type);
  fprintf(outfile,
    "   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");
  fflush(outfile);

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  // zero residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_NEW);
  memset((void*)tempt,'\0',o*o*v*v*sizeof(double));
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_NEW);
     psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }
  pair_energy = (double*)malloc(o*o*sizeof(double));

  // cepa diagrams split up as tasks
  DefineLinearTasks();

  // start timing the iterations
  struct tms total_tmstime;
  const long clk_tck = sysconf(_SC_CLK_TCK);
  times(&total_tmstime);

  time_t time_start = time(NULL);
  double user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  double sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;
// TODO e_conv

  while(iter < maxiter){
      time_t iter_start = time(NULL);

      // evaluate cepa diagrams
      if (iter>0){
         memset((void*)w1,'\0',o*v*sizeof(double));
         for (int i=0; i<nltasks; i++) {
             (*this.*LTasklist[i].func)(LParams[i]);
         }
      }

      // update the amplitudes and check the energy
      Eold = eccsd;
      PairEnergy();
      if (!options_.get_bool("CEPA_NO_SINGLES")){
          UpdateT1();
      }
      UpdateT2();

      // add vector to list for diis
      DIISOldVector(iter,diis_iter,replace_diis_iter);

      // diis error vector and convergence check
      nrm = DIISErrorVector(diis_iter,replace_diis_iter,iter);

      // diis extrapolation
      if (diis_iter>1){
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,o*o*v*v+o*v);
         else                   DIIS(diisvec,maxdiis,o*o*v*v+o*v);
         DIISNewAmplitudes(diis_iter,replace_diis_iter);
      }
      // if cepa_no_singles, zero t1
      if (options_.get_bool("CEPA_NO_SINGLES")){
         memset((void*)t1,'\0',o*v*sizeof(double));
      }
      eccsd = CheckEnergy();

      if (diis_iter < maxdiis ) {
         replace_diis_iter++;
      }else {
          double min = 1.0e9;
          for (int j = 1; j <= (diis_iter < maxdiis ? diis_iter : maxdiis); j++) {
              if ( fabs( diisvec[j-1] ) < min ) {
                  min = fabs( diisvec[j-1] );
                  replace_diis_iter = j;
              }
          }
      }

      if (diis_iter<=maxdiis) diis_iter++;
      //else if (replace_diis_iter<maxdiis) replace_diis_iter++;
      //else replace_diis_iter = 1;

      time_t iter_stop = time(NULL);
      fprintf(outfile,"  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,eccsd,eccsd-Eold,nrm,(int)iter_stop-(int)iter_start);
      fflush(outfile);
      iter++;
      if (iter==1) emp2 = eccsd;
      if (iter==1) SCS_MP2();

      if (fabs(eccsd - Eold) < e_conv && nrm < r_conv) break;
  }
  times(&total_tmstime);
  time_t time_stop = time(NULL);
  double user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  double sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;
  psio.reset();

  if (iter==maxiter){
     throw PsiException("  CEPA iterations did not converge.",__FILE__,__LINE__);
  }

  // is this a cim-cepa computation?
  if (reference_wavefunction_->isCIM()){
     Local_SCS_CEPA();
     eccsd = eccsd_os + eccsd_ss;
  }
  else{
     SCS_CEPA();
  }

  fprintf(outfile,"\n");
  fprintf(outfile,"  %s iterations converged!\n",cepa_type);
  fprintf(outfile,"\n");

  // delta mp2 correction for fno computations:
  if (options_.get_bool("NAT_ORBS")){
      double delta_emp2 = Process::environment.globals["MP2 CORRELATION ENERGY"] - emp2;
      double delta_emp2_os = Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] - emp2_os;
      double delta_emp2_ss = Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] - emp2_ss;

      emp2 += delta_emp2;
      emp2_os += delta_emp2_os;
      emp2_ss += delta_emp2_ss;

      eccsd += delta_emp2;
      eccsd_os += delta_emp2_os;
      eccsd_ss += delta_emp2_ss;

      fprintf(outfile,"        OS MP2 FNO correction:             %20.12lf\n",delta_emp2_os);
      fprintf(outfile,"        SS MP2 FNO correction:             %20.12lf\n",delta_emp2_ss);
      fprintf(outfile,"        MP2 FNO correction:                %20.12lf\n",delta_emp2);
      fprintf(outfile,"\n");
  }

  fprintf(outfile,"        OS SCS-MP2 correlation energy:     %20.12lf\n",emp2_os*emp2_os_fac);
  fprintf(outfile,"        SS SCS-MP2 correlation energy:     %20.12lf\n",emp2_ss*emp2_ss_fac);
  fprintf(outfile,"        SCS-MP2 correlation energy:        %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac);
  fprintf(outfile,"      * SCS-MP2 total energy:              %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac+escf);
  fprintf(outfile,"\n");
  fprintf(outfile,"        OS MP2 correlation energy:         %20.12lf\n",emp2_os);
  fprintf(outfile,"        SS MP2 correlation energy:         %20.12lf\n",emp2_ss);
  fprintf(outfile,"        MP2 correlation energy:            %20.12lf\n",emp2);
  fprintf(outfile,"      * MP2 total energy:                  %20.12lf\n",emp2+escf);
  fprintf(outfile,"\n");
  if (cepa_level>=0){
     if (options_.get_bool("SCS_CEPA")){
        fprintf(outfile,"        OS SCS-%s correlation energy: %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac);
        fprintf(outfile,"        SS SCS-%s correlation energy: %20.12lf\n",cepa_type,eccsd_ss*eccsd_ss_fac);
        fprintf(outfile,"        SCS-%s correlation energy:    %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
        fprintf(outfile,"      * SCS-%s total energy:          %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
        fprintf(outfile,"\n");
     }
     fprintf(outfile,"        OS %s correlation energy:     %20.12lf\n",cepa_type,eccsd_os);
     fprintf(outfile,"        SS %s correlation energy:     %20.12lf\n",cepa_type,eccsd_ss);
     fprintf(outfile,"        %s correlation energy:        %20.12lf\n",cepa_type,eccsd);
     fprintf(outfile,"      * %s total energy:              %20.12lf\n",cepa_type,eccsd+escf);
  }else{
     if (options_.get_bool("SCS_CEPA")){
        fprintf(outfile,"        OS SCS-%s correlation energy:    %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac);
        fprintf(outfile,"        SS SCS-%s correlation energy:    %20.12lf\n",cepa_type,eccsd_ss*eccsd_ss_fac);
        fprintf(outfile,"        SCS-%s correlation energy:       %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
        fprintf(outfile,"      * SCS-%s total energy:             %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
        fprintf(outfile,"\n");
     }
     fprintf(outfile,"        OS %s correlation energy:        %20.12lf\n",cepa_type,eccsd_os);
     fprintf(outfile,"        SS %s correlation energy:        %20.12lf\n",cepa_type,eccsd_ss);
     fprintf(outfile,"        %s correlation energy:           %20.12lf\n",cepa_type,eccsd);
     fprintf(outfile,"      * %s total energy:                 %20.12lf\n",cepa_type,eccsd+escf);
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

  free(pair_energy);
  return Success;
}

void CoupledPair::PairEnergy(){

  if (cepa_level<1) return;

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempt;
  }

  for (long int i=0; i<o; i++){
      for (long int j=0; j<o; j++){
          double energy=0.0;
          for (long int a=o; a<rs; a++){
              for (long int b=o; b<rs; b++){
                  long int ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  energy += integrals[iajb]*(2.0*tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
              }
          }
          pair_energy[i*o+j] = energy;
      }
  }
}

void CoupledPair::UpdateT2() {

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  double fac = 1.0;
  if (cepa_level == 0)       fac = 0.0;
  else if (cepa_level == -1) fac = 1.0;
  else if (cepa_level == -2) fac = 1.0/o;
  else if (cepa_level == -3) fac = 1.0-(2.0*o-2.0)*(2.0*o-3.0) / (2.0*o*(2.0*o-1.0));
  double energy = eccsd * fac;

  for (long int i=0; i<o; i++){
      double di = - eps[i];
      for (long int j=0; j<o; j++){
          double dij = di-eps[j];

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

          for (long int a=o; a<rs; a++){
              double dija = dij + eps[a];
              for (long int b=o; b<rs; b++){
                  double dijab = dija + eps[b];

                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  long int ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;

                  double tnew = - (integrals[iajb] + tempv[ijab])/(dijab-energy);
                  tempt[ijab] = tnew;
              }
          }
      }
  }

  // error vectors for diis are in tempv:
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tb,1,tempv,1);
  }
  F_DAXPY(o*o*v*v,-1.0,tempt,1,tempv,1);
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     F_DCOPY(o*o*v*v,tempt,1,tb,1);
  }

  psio.reset();
}

void CoupledPair::UpdateT1() {

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  double fac = 1.0;
  if (cepa_level == 0)       fac = 0.0;
  else if (cepa_level == -1) fac = 1.0;
  else if (cepa_level == -2) fac = 1.0/o;
  else if (cepa_level == -3) fac = 1.0-(2.0*o-2.0)*(2.0*o-3.0) / (2.0*o*(2.0*o-1.0));
  double energy = eccsd * fac;

  for (long int i=0; i<o; i++){

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

      for (long int a=o; a<rs; a++){
          double dia = -eps[i]+eps[a];
          double tnew = - (w1[(a-o)*o+i])/(dia-energy);
          w1[(a-o)*o+i] = tnew;
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

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&tempt[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  SharedMatrix Rii = reference_wavefunction_->CIMTransformationMatrix();
  double**Rii_pointer = Rii->pointer();

  // transform E2iajb back from quasi-canonical basis
  for (long int i=0; i<o; i++){
      for (long int a=0; a<v; a++){
          for (long int j=0; j<o; j++){
              for (long int b=0; b<v; b++){
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
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  // transform t2 back from quasi-canonical basis
  for (long int a=0; a<v; a++){
      for (long int b=0; b<v; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){
                  double dum = 0.0;
                  for (int ip=0; ip<o; ip++){
                      dum += tb[a*o*o*v+b*o*o+ip*o+j]*Rii_pointer[ip][i];
                  }
                  tempt[a*o*o*v+b*o*o+i*o+j] = dum;
              }
          }
      }
  }

  SharedVector factor = reference_wavefunction_->CIMOrbitalFactors();
  double*factor_pointer = factor->pointer();

  double ssenergy = 0.0;
  double osenergy = 0.0;
  for (long int a=o; a<rs; a++){
      for (long int b=o; b<rs; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){

                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  long int jaib = iajb + (i-j)*v*(1-v*o);
                  long int ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;

                  osenergy += integrals[iajb]*(tempt[ijab])*factor_pointer[i];
                  ssenergy += integrals[iajb]*(tempt[ijab]-tempt[(b-o)*o*o*v+(a-o)*o*o+i*o+j])*factor_pointer[i];
              }
          }
      }
  }
  eccsd_os = osenergy;
  eccsd_ss = ssenergy;

  psio.reset();
}
void CoupledPair::SCS_CEPA(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  double ssenergy = 0.0;
  double osenergy = 0.0;
  for (long int a=o; a<rs; a++){
      for (long int b=o; b<rs; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){

                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  long int jaib = iajb + (i-j)*v*(1-v*o);
                  long int ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;

                  osenergy += integrals[iajb]*(tb[ijab]);
                  ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
              }
          }
      }
  }
  eccsd_os = osenergy;
  eccsd_ss = ssenergy;

  psio.reset();
}
double CoupledPair::CheckEnergy(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  double energy = 0.0;
  for (long int a=o; a<rs; a++){
      for (long int b=o; b<rs; b++){
          for (long int i=0; i<o; i++){
              for (long int j=0; j<o; j++){

                  long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  long int jaib = iajb + (i-j)*v*(1-v*o);
                  long int ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;

                  energy += (2.*integrals[iajb]-integrals[jaib])*(tb[ijab]);
              }
          }
      }
  }

  psio.reset();

  return energy;
}

void CoupledPair::finalize(){
  free(integrals);
  free(tempt);
  free(tempv);
  if (!t2_on_disk){
     free(tb);
  }
  free(w1);
  free(t1);
  free(I1);
  free(I1p);
  free(diisvec);

  // there is something weird with chkpt_ ... reset it
  chkpt_.reset();
}


}} // end of namespaces
