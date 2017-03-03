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

#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include"psi4/libqt/qt.h"
#include<sys/times.h>
#include "psi4/libciomr/libciomr.h"
#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() 0.0
#endif

#include"blas.h"
#include"ccsd.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/basisset_parser.h"
#include "psi4/lib3index/3index.h"

using namespace psi;


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
void SortIntegrals(int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,Options&options);
void Sort_OV3_LowMemory(long int memory,long int o,long int v);

// coupled cluster constructor
CoupledCluster::CoupledCluster(SharedWavefunction ref_wfn, Options &options):
        Wavefunction(options)
{
    shallow_copy(ref_wfn);
    reference_wavefunction_ = ref_wfn;
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
  epsilon_a_= std::shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
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
  ndoccact = ndocc - nfzc;
  nvirt    = nmo - ndoccact;

  if (ndoccact <= 0) {
      throw PSIEXCEPTION("Number of active orbitals is zero.");
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
  SortIntegrals(nfzc,nfzv,nmo+nfzc+nfzv,ndoccact,nvirt,options_);
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

     // if low memory, need to generate one last set of integrals
     if ( isLowMemory ){
        Sort_OV3_LowMemory(memory - 8L*(long int)(!t2_on_disk)*o*o*v*v,o,v);
     }

     // now there should be space for t2
     if (t2_on_disk){
         tb = (double*)malloc(o*o*v*v*sizeof(double));
         std::shared_ptr<PSIO> psio (new PSIO());
         psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
         psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tb[0],o*o*v*v*sizeof(double));
         psio->close(PSIF_DCC_T2,1);
     }

     bool do_cc = !options_.get_bool("RUN_MP4") && !options_.get_bool("RUN_MP3");
     bool do_mp = options_.get_bool("COMPUTE_MP4_TRIPLES");

     tstart();
     // triples
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

  outfile->Printf("\n\n");
  outfile->Printf(     "        *****************************************************\n");
  outfile->Printf(     "        *                                                   *\n");
  if (isccsd)
      outfile->Printf( "        *                       CCSD                        *\n");
  else if (mp2_only)
      outfile->Printf( "        *                        MP2                        *\n");
  else if (mp4_only)
      outfile->Printf( "        *                        MP4                        *\n");
  else if (mp3_only)
      outfile->Printf( "        *                        MP3                        *\n");
  else
      outfile->Printf( "        *                       QCISD                       *\n");
  outfile->Printf(     "        *                  Eugene DePrince                  *\n");
  outfile->Printf(     "        *                                                   *\n");
  outfile->Printf(     "        *****************************************************\n");
  outfile->Printf("\n\n");

  WriteOptions();
}
void CoupledCluster::WriteOptions(){
    outfile->Printf("\n");
    outfile->Printf("  ==> Input parameters <==\n\n");
    outfile->Printf("        Freeze core orbitals?               %5s\n",nfzc > 0 ? "yes" : "no");
    outfile->Printf("        Use frozen natural orbitals?        %5s\n",options_.get_bool("NAT_ORBS") ? "yes" : "no");
    outfile->Printf("        r_convergence:                  %5.3le\n",r_conv);
    outfile->Printf("        e_convergence:                  %5.3le\n",e_conv);
    outfile->Printf("        Number of DIIS vectors:             %5li\n",maxdiis);
    outfile->Printf("        Number of frozen core orbitals:     %5li\n",nfzc);
    outfile->Printf("        Number of active occupied orbitals: %5li\n",ndoccact);
    outfile->Printf("        Number of active virtual orbitals:  %5li\n",nvirt);
    outfile->Printf("        Number of frozen virtual orbitals:  %5li\n",nfzv);
}

/*===================================================================

  solve cc/qci equations

===================================================================*/
PsiReturnType CoupledCluster::CCSDIterations() {

  long int o = ndoccact;
  long int v = nvirt;

  iter                  = 0;
  int diis_iter         = 0;
  int replace_diis_iter = 1;
  double nrm            = 1.0;
  double Eold           = 1.0e9;
  eccsd                 = 0.0;

  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  outfile->Printf("\n");
  if (isccsd) {
     outfile->Printf(
       "  Begin singles and doubles coupled cluster iterations\n\n");
  }else {
     outfile->Printf(
       "  Begin singles and doubles quadratic ci iterations\n\n");
  }
  outfile->Printf(
    "   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");


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
      if (timer) outfile->Printf("\n");
      for (int i = 0; i < ncctasks; i++) {
          if (timer) s1 = omp_get_wtime();
          (*this.*CCTasklist[i].func)(CCParams[i]);
          if (timer) outfile->Printf("        %s ... %6.2lf s\n",CCTasklist[i].name,omp_get_wtime()-s1);
      }
      if (timer) outfile->Printf("\n");

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
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,o*o*v*v+o*v,replace_diis_iter);
         else                   DIIS(diisvec,maxdiis,o*o*v*v+o*v,replace_diis_iter);
         DIISNewAmplitudes(diis_iter,replace_diis_iter);
      }
      eccsd = CheckEnergy();

      //if (diis_iter < maxdiis ) {
      //   replace_diis_iter++;
      //}else {
      //    double min = 1.0e9;
      //    for (int j = 1; j <= (diis_iter < maxdiis ? diis_iter : maxdiis); j++) {
      //        if ( fabs( diisvec[j-1] ) < min ) {
      //            min = fabs( diisvec[j-1] );
      //            replace_diis_iter = j;
      //        }
      //    }
      //}

      if (diis_iter <= maxdiis) diis_iter++;
      else if (replace_diis_iter < maxdiis) replace_diis_iter++;
      else    replace_diis_iter = 1;

      time_t iter_stop = time(NULL);
      outfile->Printf("  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,eccsd,eccsd-Eold,nrm,(int)iter_stop-(int)iter_start);

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

  SCS_CCSD();

  outfile->Printf("\n");
  if (isccsd)
     outfile->Printf("  CCSD iterations converged!\n");
  else
     outfile->Printf("  QCISD iterations converged!\n");
  outfile->Printf("\n");

  // T1 and D1 diagnostics:

  double t1diag = C_DNRM2(o*v,t1,1) / sqrt(2.0 * o);
  outfile->Printf("        T1 diagnostic:                   %20.12lf\n",t1diag);
  std::shared_ptr<Matrix>T (new Matrix(o,o));
  std::shared_ptr<Matrix>eigvec (new Matrix(o,o));
  std::shared_ptr<Vector>eigval (new Vector(o));
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
  outfile->Printf("        D1 diagnostic:                   %20.12lf\n",sqrt(eigval->pointer()[0]));
  outfile->Printf("\n");

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

      outfile->Printf("        OS MP2 FNO correction:           %20.12lf\n",delta_emp2_os);
      outfile->Printf("        SS MP2 FNO correction:           %20.12lf\n",delta_emp2_ss);
      outfile->Printf("        MP2 FNO correction:              %20.12lf\n",delta_emp2);
      outfile->Printf("\n");
  }

  if (options_.get_bool("SCS_MP2")){
      outfile->Printf("        OS SCS-MP2 correlation energy:   %20.12lf\n",emp2_os*emp2_os_fac);
      outfile->Printf("        SS SCS-MP2 correlation energy:   %20.12lf\n",emp2_ss*emp2_ss_fac);
      outfile->Printf("        SCS-MP2 correlation energy:      %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac);
      outfile->Printf("      * SCS-MP2 total energy:            %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac+escf);
      outfile->Printf("\n");
  }
  outfile->Printf("        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
  outfile->Printf("        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
  outfile->Printf("        MP2 correlation energy:          %20.12lf\n",emp2);
  outfile->Printf("      * MP2 total energy:                %20.12lf\n",emp2+escf);
  outfile->Printf("\n");
  outfile->Printf("        OS MP2.5 correlation energy:     %20.12lf\n",emp2_os+0.5*emp3_os);
  outfile->Printf("        SS MP2.5 correlation energy:     %20.12lf\n",emp2_ss+0.5*emp3_ss);
  outfile->Printf("        MP2.5 correlation energy:        %20.12lf\n",emp2+0.5*emp3);
  outfile->Printf("      * MP2.5 total energy:              %20.12lf\n",emp2+0.5*emp3+escf);
  outfile->Printf("\n");
  outfile->Printf("        OS MP3 correlation energy:       %20.12lf\n",emp2_os+emp3_os);
  outfile->Printf("        SS MP3 correlation energy:       %20.12lf\n",emp2_ss+emp3_ss);
  outfile->Printf("        MP3 correlation energy:          %20.12lf\n",emp2+emp3);
  outfile->Printf("      * MP3 total energy:                %20.12lf\n",emp2+emp3+escf);
  outfile->Printf("\n");
  outfile->Printf("        OS MP4(SDQ) correlation energy:  %20.12lf\n",emp2_os+emp3_os+emp4_sd_os+emp4_q_os);
  outfile->Printf("        SS MP4(SDQ) correlation energy:  %20.12lf\n",emp2_ss+emp3_ss+emp4_sd_ss+emp4_q_os);
  outfile->Printf("        MP4(SDQ) correlation energy:     %20.12lf\n",emp2+emp3+emp4_sd+emp4_q);
  outfile->Printf("      * MP4(SDQ) total energy:           %20.12lf\n",emp2+emp3+emp4_sd+emp4_q+escf);
  outfile->Printf("\n");

  if (isccsd) {
     if (options_.get_bool("SCS_CCSD")){
        outfile->Printf("        OS SCS-CCSD correlation energy:  %20.12lf\n",eccsd_os*eccsd_os_fac);
        outfile->Printf("        SS SCS-CCSD correlation energy:  %20.12lf\n",eccsd_ss*eccsd_ss_fac);
        outfile->Printf("        SCS-CCSD correlation energy:     %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
        outfile->Printf("      * SCS-CCSD total energy:           %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
        outfile->Printf("\n");
     }
     outfile->Printf("        OS CCSD correlation energy:      %20.12lf\n",eccsd_os);
     outfile->Printf("        SS CCSD correlation energy:      %20.12lf\n",eccsd_ss);
     outfile->Printf("        CCSD correlation energy:         %20.12lf\n",eccsd);
     outfile->Printf("      * CCSD total energy:               %20.12lf\n",eccsd+escf);
     outfile->Printf("\n");
     outfile->Printf("  Total time for CCSD iterations:  %10.2lf s (user)\n",user_stop-user_start);
  }else{
     if (options_.get_bool("SCS_CCSD")){
        outfile->Printf("        OS SCS-QCISD correlation energy: %20.12lf\n",eccsd_os*eccsd_os_fac);
        outfile->Printf("        SS SCS-QCISD correlation energy: %20.12lf\n",eccsd_ss*eccsd_ss_fac);
        outfile->Printf("        SCS-QCISD correlation energy:    %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
        outfile->Printf("      * SCS-QCISD total energy:          %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
        outfile->Printf("\n");
     }
     outfile->Printf("        OS QCISD correlation energy:     %20.12lf\n",eccsd_os);
     outfile->Printf("        SS QCISD correlation energy:     %20.12lf\n",eccsd_ss);
     outfile->Printf("        QCISD correlation energy:        %20.12lf\n",eccsd);
     outfile->Printf("      * QCISD total energy:              %20.12lf\n",eccsd+escf);
     outfile->Printf("\n");
     outfile->Printf("  Total time for QCISD iterations: %10.2lf s (user)\n",user_stop-user_start);
  }

  outfile->Printf("                                   %10.2lf s (system)\n",sys_stop-sys_start);
  outfile->Printf("                                   %10d s (total)\n",(int)time_stop-(int)time_start);
  outfile->Printf("\n");
  outfile->Printf("  Time per iteration:              %10.2lf s (user)\n",(user_stop-user_start)/iter);
  outfile->Printf("                                   %10.2lf s (system)\n",(sys_stop-sys_start)/iter);
  outfile->Printf("                                   %10.2lf s (total)\n",((double)time_stop-(double)time_start)/(iter-1));


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
     outfile->Printf("\n");
     outfile->Printf("  ==> Define tiling <==\n");
     outfile->Printf("\n");
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

  outfile->Printf("        v(ab,cd) diagrams will be evaluated in %3li blocks.\n",ntiles);


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

  outfile->Printf("        v(ab,ci) diagrams will be evaluated in %3li blocks over ov2.\n",nov2tiles);


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
  outfile->Printf("        v(ab,ci) diagrams will be evaluated in %3li blocks over ov.\n",novtiles);

}

/*===================================================================

  - allocate cpu memory
  - define tiling for diagrams that don't fit in core

===================================================================*/
void CoupledCluster::AllocateMemory() {

  long int nthreads = Process::environment.get_n_threads();
  long int o=ndoccact;
  long int v=nvirt;
  if (!options_.get_bool("RUN_MP2")) {
    outfile->Printf("\n");
    outfile->Printf("  ==> Memory <==\n\n");
    outfile->Printf("        available memory =                         %9.2lf mb\n",memory/1024./1024.);
    if (isccsd){
       outfile->Printf("        minimum memory requirements for CCSD =     %9.2lf mb\n",
           8./1024./1024.*(o*o*v*v+2.*(o*o*v*v+o*v)+2.*o*v+2.*v*v+o+v));
    }else{
       outfile->Printf("        minimum memory requirements for QCISD =    %9.2lf mb\n",
           8./1024./1024.*(o*o*v*v+2.*(o*o*v*v+o*v)+2.*o*v+2.*v*v+o+v));
    }
    if (options_.get_bool("COMPUTE_TRIPLES") || options_.get_bool("COMPUTE_MP4_TRIPLES")){
       double tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+3L*v*v*v*nthreads);
       if (tempmem > memory) {
          outfile->Printf("\n        <<< warning! >>> switched to low-memory (t) algorithm\n\n");
       }
       if (tempmem > memory || options_.get_bool("TRIPLES_LOW_MEMORY")){
          isLowMemory = true;
          tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+5L*o*o*o*nthreads);
       }
       if (isccsd)
          outfile->Printf("        memory requirements for CCSD(T) =          %9.2lf mb\n",tempmem/1024./1024.);
       else
          outfile->Printf("        memory requirements for QCISD(T) =         %9.2lf mb\n",tempmem/1024./1024.);
    }

  }
  // orbital energies:
  int count=0;
  eps = (double*)malloc((ndoccact+nvirt)*sizeof(double));
  std::shared_ptr<Vector> eps_test = reference_wavefunction_->epsilon_a();
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
     outfile->Printf("\n");
     outfile->Printf("  Warning: cannot accomodate T2 in core. T2 will be stored on disk.\n");
     outfile->Printf("\n");

     t2_on_disk = true;
     DefineTilingCPU();
     dim = 0;
     if (tilesize*fulltile > dim) dim = tilesize*fulltile;
     if (ovtilesize*v*v > dim)    dim = ovtilesize*v*v;
     if (ov2tilesize*v > dim)     dim = ov2tilesize*v;

     if (dim<o*o*v*v){
        throw PsiException("out of memory: general buffer cannot accomodate T2",__FILE__,__LINE__);
     }

     outfile->Printf("\n");
     outfile->Printf("  Increase memory by %7.2lf mb to hold T2 in core.\n",o*o*v*v*8L/1024./1024.);
     outfile->Printf("\n");
  }

  maxelem = dim;

  long int oovv = o*o*v*v;
  double total_memory = 1.*dim+2.*(oovv+o*v)+1.*o*o*v*v+2.*o*v+2.*v*v;
  if (t2_on_disk) total_memory = 1.*dim+2.*(oovv+o*v)+2.*o*v+2.*v*v;
  total_memory *= 8./1024./1024.;

  outfile->Printf("\n");
  outfile->Printf("  Allocate cpu memory (%9.2lf mb).....",total_memory);

  integrals = (double*)malloc(dim*sizeof(double));
  tempt     = (double*)malloc((oovv+o*v)*sizeof(double));
  tempv     = (double*)malloc((oovv+o*v)*sizeof(double));

  if (!t2_on_disk) tb = (double*)malloc(o*o*v*v*sizeof(double));
  w1        = (double*)malloc(o*v*sizeof(double));
  t1        = (double*)malloc(o*v*sizeof(double));
  I1        = (double*)malloc(v*v*sizeof(double));
  I1p       = (double*)malloc(v*v*sizeof(double));
  outfile->Printf("done.\n");

  outfile->Printf("  Initialize cpu memory..................");
  memset((void*)integrals,'\0',dim*sizeof(double));
  memset((void*)tempv,'\0',(oovv+o*v)*sizeof(double));
  memset((void*)tempt,'\0',(oovv+o*v)*sizeof(double));
  if (!t2_on_disk) memset((void*)tb,'\0',o*o*v*v*sizeof(double));
  memset((void*)w1,'\0',o*v*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));
  memset((void*)I1,'\0',v*v*sizeof(double));
  memset((void*)I1p,'\0',v*v*sizeof(double));
  outfile->Printf("done.\n");

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
  std::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAB,"E2ijab",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IJAB,1);

  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  C_DAXPY(o*o*v*v,-2.0,integrals,1,tempv,1);

  for (i=0; i<o; i++){
      C_DCOPY(v,t1+i,o,tempt+i*v,1);
  }
  F_DGEMV('n',o*v,o*v,-1.0,tempv,o*v,tempt,1,0.0,integrals,1);
  for (a=0; a<v; a++){
      C_DAXPY(o,1.0,integrals+a,v,w1+a*o,1);
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
  std::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (a=0,id=0; a<v; a++){
      for (m=0; m<o; m++){
          for (n=0; n<o; n++){
              C_DCOPY(v,tb+a*v*o*o+m*o+n,o*o,tempt+a*o*o*v+m*o*v+n*v,1);
              C_DAXPY(v,-2.0,tb+a*o*o+m*o+n,o*o*v,tempt+a*o*o*v+m*o*v+n*v,1);
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

  std::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     tb = tempv;
  }

  for (f=0,id=0; f<v; f++){
      for (m=0; m<o; m++){
          for (e=0; e<v; e++){
              C_DCOPY(o,tb+e*v*o*o+f*o*o+m*o,1,tempt+f*o*o*v+m*o*v+e*o,1);
              C_DAXPY(o,-0.5,tb+e*v*o*o+f*o*o+m,o,tempt+f*o*o*v+m*o*v+e*o,1);
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
  std::shared_ptr<PSIO> psio(new PSIO());
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
              C_DCOPY(v,tb+e*v*o*o+m*o+n,o*o,tempt+m*o*v*v+e*o*v+n*v,1);
              if (isccsd) {
                 for (b=0; b<v; b++){
                     tempt[id++] += t1[e*o+m]*t1[b*o+n];
                 }
              }
          }
      }
  }
  C_DCOPY(o*o*v*v,integrals,1,tempv,1);
  for (m=0,id=0; m<o; m++){
      for (e=0; e<v; e++){
          for (n=0; n<o; n++){
              C_DAXPY(v,-0.5,integrals+m*o*v*v+n*v+e,o*v,tempv+m*o*v*v+e*o*v+n*v,1);
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
         C_DCOPY(v,t1+i,o,tempt+i*v,1);
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
              C_DCOPY(v,tb+c*o*o+l*o+k,v*o*o,tempt+l*o*v*v+c*o*v+k*v,1);
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
              C_DAXPY(o,1.0,tempv+a*v*o+i*v+b,v*v*o,tempt+a*o*o*v+b*o*o+i*o,1);
              C_DAXPY(o,1.0,tempv+i*v*v*o+b*v*o+a,v,tempt+a*o*o*v+b*o*o+i*o,1);
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

  std::shared_ptr<PSIO> psio(new PSIO());

  // now build and use intermediate:
  psio->open(PSIF_DCC_IJAB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IJAB,"E2ijab",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IJAB,1);
  F_DGEMM('n','n',o,o2v,v,-1.0,t1,o,tempv,v,0.0,tempt,o);
  F_DGEMM('n','n',o2v,v,o,1.0,tempt,o2v,t1,o,0.0,tempv,o2v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
  C_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,tempv+a*v*o*o+b*o*o+i*o,1,tempt+b*v*o*o+a*o*o+i,o);
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
  std::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_IAJB,1);
  C_DCOPY(o*o*v*v,integrals,1,tempv,1);
  for (i=0; i<o; i++){
      for (a=0; a<v; a++){
          for (m=0; m<o; m++){
              C_DAXPY(v,-0.5,integrals+i*o*v*v+m*v+a,o*v,tempv+i*v*v*o+a*v*o+m*v,1);
          }
      }
  }
  for (i=0; i<o; i++) C_DCOPY(v,t1+i,o,tempt+i*v,1);
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
              C_DCOPY(v,tb+e*o*o*v+m*o+j,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
              C_DAXPY(v,-0.5,tb+e*o*o*v+j*o+m,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
          }
      }
  }
  F_DGEMV('n',o*v,o*v,2.0,tempt,o*v,I1,1,0.0,tempv,1);
  for (i=0; i<o; i++){
      C_DAXPY(v,1.0,tempv+i*v,1,w1+i,o);
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
                 C_DCOPY(o,tempt+i*o*v+j*v+e,o*o*v,tempv+i*o*o*v+j*o*v+e*o,1);
                 C_DAXPY(o,-2.0,tempt+i*o*o*v+j*v+e,o*v,tempv+i*o*o*v+j*o*v+e*o,1);
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
              C_DCOPY(v,tb+e*o*o*v+m*o+j,o*o,tempt+m*o*v*v+e*o*v+j*v,1);
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
              C_DAXPY(o,1.0,tempv+a*o*o*v+b*o+i,v*o,tempt+a*o*o*v+b*o*o+i*o,1);
              C_DAXPY(o,1.0,tempv+b*o*o*v+i*v*o+a*o,1,tempt+a*o*o*v+b*o*o+i*o,1);
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
  std::shared_ptr<PSIO> psio(new PSIO());

  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     C_DCOPY(o*o*v*v,tb,1,tempt,1);
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
              C_DCOPY(v,integrals+j*o*v*v+b*o*v+i*v,1,tempv+j*o*v*v+i*v*v+b*v,1);
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
  C_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,tempv+b*v*o*o+a*o*o+i,o,tempt+a*v*o*o+b*o*o+i*o,1);
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
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;

  if (isccsd) {
     if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
     }else{
        C_DCOPY(o*o*v*v,tb,1,tempt,1);
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
                 C_DAXPY(o,1.0,tempt+i*o*o*v+a*o+j,o*v,tempv+j*o*o*v+a*o*o+i*o,1);
             }
         }
     }
  }

  // use intermediate
  F_DGEMM('n','n',o*o*v,v,o,-1.0,tempv,o*o*v,t1,o,0.0,tempt,o*o*v);

  // contribute to residual
  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  C_DAXPY(o*o*v*v,1.0,tempt,1,tempv,1);
  for (a=0; a<v; a++){
      for (b=0; b<v; b++){
          for (i=0; i<o; i++){
              C_DAXPY(o,1.0,tempt+b*v*o*o+a*o*o+i,o,tempv+a*v*o*o+b*o*o+i*o,1);
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
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     C_DCOPY(o*o*v*v,tb,1,tempt,1);
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
  std::shared_ptr<PSIO> psio(new PSIO());
  psio_address addr;
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     C_DCOPY(o*o*v*v,tb,1,tempt,1);
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
  std::shared_ptr<PSIO> psio(new PSIO());
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
                  C_DAXPY(v,1.0,tempv+b*o*o+i*o+j,o*o*v,tempt+i*o*v*v+b*o*v+j*v,1);
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
                  C_DAXPY(v,1.0,tempv+i*o*v+b*o+j,o*o*v,tempt+i*o*v*v+b*o*v+j*v,1);
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
              C_DCOPY(v,tb+b*o*o*v+j*o+i,o*o,integrals+i*o*v*v+b*o*v+j*v,1);
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
              C_DCOPY(v,tempt+i*v*v*o+j*v+b,o*v,tempv+i*o*v*v+b*o*v+j*v,1);
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
              C_DCOPY(v,tb+a*o*o+j*o+i,o*o*v,tempv+j*o*v*v+a*o*v+i*v,1);
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
              C_DCOPY(o,    integrals+j*o*v*v+b*v*o+a,v,tempt+a*o*o*v+b*o*o+j*o,1);
              C_DAXPY(o,1.0,integrals+a*v*o+j*v+b,o*v*v,tempt+a*o*o*v+b*o*o+j*o,1);
              C_DAXPY(o,0.5,integrals+j*o*v*v+a*v*o+b,v,tempt+a*o*o*v+b*o*o+j*o,1);
              C_DAXPY(o,0.5,integrals+b*v*o+j*v+a,o*v*v,tempt+a*o*o*v+b*o*o+j*o,1);
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
  std::shared_ptr<PSIO> psio(new PSIO());
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
                  C_DAXPY(v,1.0,tempt+i*o*v+j*v+b,o*o*v,tempv+i*o*v*v+b*o*v+j*v,1);
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
  C_DCOPY(o*o*v*v,tempt,1,integrals,1);
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
                  C_DAXPY(o,1.0,tempv+i*o*v*v+a*o*v+b*o,1,tempt+i*o*v*v+b*o*v+a,v);
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
              C_DAXPY(o,1.0,tempv+i*v*v*o+b*o*v+a*o,1,integrals+a*v*o*o+b*o*o+i*o,1);
              C_DAXPY(o,1.0,tempv+i+a*o*v+b*o,v*v*o,integrals+a*v*o*o+b*o*o+i*o,1);
          }
      }
  }
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);

  // contribute to intermediate
  psio->open(PSIF_DCC_TEMP,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_TEMP,"temporary_J",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_TEMP,1);
  C_DAXPY(o*o*v*v,1.0,tempt,1,tempv,1);

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
              C_DCOPY(v,tb+b*o*o+i*o+j,o*o*v,tempt+j*o*v*v+b*o*v+i*v,1);
              C_DAXPY(v,-0.5,tb+b*o*o*v+i*o+j,o*o,tempt+j*o*v*v+b*o*v+i*v,1);
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
              C_DAXPY(o,1.0,integrals+b*v*o+i*v+a,o*v*v,tempt+a*o*o*v+b*o*o+i*o,1);
              C_DAXPY(o,1.0,integrals+i*o*v*v+a*v*o+b,v,tempt+a*o*o*v+b*o*o+i*o,1);
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

  std::shared_ptr<PSIO> psio(new PSIO());
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
     C_DCOPY(o*o*v*v,tb,1,tempv,1);
  }
  C_DAXPY(o*o*v*v,-1.0,tempt,1,tempv,1);
  if (t2_on_disk){
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }else{
     C_DCOPY(o*o*v*v,tempt,1,tb,1);
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
  C_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  C_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  C_DCOPY(o*v,w1,1,t1,1);
}


/*================================================================

   SCS functions.  get energy in terms of spin components

================================================================*/
void CoupledCluster::SCS_CCSD(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;
  long int i,j,a,b;
  long int iajb,ijab=0;
  double ssenergy = 0.0;
  double osenergy = 0.0;
  std::shared_ptr<PSIO> psio(new PSIO());
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
  std::shared_ptr<PSIO> psio(new PSIO());
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
  std::shared_ptr<PSIO> psio(new PSIO());
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
  std::shared_ptr<PSIO> psio(new PSIO());
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
  outfile->Printf("\n");
  if (mp3_only) {
      outfile->Printf("  ==>   MP3   <==\n");
  }else {
      outfile->Printf("  ==> MP4(SDQ) <==\n");
  }
  outfile->Printf("\n");
  outfile->Printf("        1st-order doubles amplitudes...............done.\n");
  outfile->Printf("        MP2 energy.................................");
  UpdateT2_mp4(0);
  outfile->Printf("done.\n");

  // <1|V|1> for mp3
  outfile->Printf("        MP3 energy.................................");
  memset((void*)w1,'\0',o*v*sizeof(double));
  for (int i=0; i<nltasks; i++) {
      (*this.*LTasklist[i].func)(LParams[i]);
  }
  // mp3 energy and 2nd-order doubles amplitudes
  UpdateT2_mp4(1);
  outfile->Printf("done.\n");

  if (!mp3_only) {
      outfile->Printf("        2nd-order singles and doubles amplitudes...");

      // 2nd-order singles amplitudes
      UpdateT1_mp4(1);
      outfile->Printf("done.\n");

      // V|2> for S and D parts of mp4
      psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
      psio->write_entry(PSIF_DCC_T2,"second",(char*)&tempt[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_T2,1);
      if (t2_on_disk) {
          psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
          psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
          psio->close(PSIF_DCC_T2,1);
      }else {
          C_DCOPY(o*o*v*v,tempt,1,tb,1);
      }

      outfile->Printf("        MP4(SD)....................................");
      memset((void*)w1,'\0',o*v*sizeof(double));
      for (int i=0; i<nltasks; i++) {
          (*this.*LTasklist[i].func)(LParams[i]);
      }
      //C_DCOPY(o*o*v*v,tb,1,tempt,1);
      if (t2_on_disk) {
          tb = tempt;
      }
      psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
      psio->read_entry(PSIF_DCC_T2,"first",(char*)&tb[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_T2,1);
      UpdateT2_mp4(2);
      outfile->Printf("done.\n");

      // quadruples: evaluate ccd residual using first-order doubles amplitudes
      outfile->Printf("        MP4(Q).....................................");
      for (int i=0; i<nqtasks; i++) {
          (*this.*QTasklist[i].func)(QParams[i]);
      }
      UpdateT2_mp4(3);
      outfile->Printf("done.\n");
  }

  if (mp4_only) {
      if ( options_.get_bool("NAT_ORBS") ) {
          double delta_emp2 = Process::environment.globals["MP2 CORRELATION ENERGY"] - emp2;
          double delta_emp2_os = Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] - emp2_os;
          double delta_emp2_ss = Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] - emp2_ss;

          emp2 += delta_emp2;
          emp2_os += delta_emp2_os;
          emp2_ss += delta_emp2_ss;

          outfile->Printf("\n");
          outfile->Printf("        OS MP2 FNO correction:          %20.12lf\n",delta_emp2_os);
          outfile->Printf("        SS MP2 FNO correction:          %20.12lf\n",delta_emp2_ss);
          outfile->Printf("        MP2 FNO correction:             %20.12lf\n",delta_emp2);
      }

      outfile->Printf("\n");
      outfile->Printf("        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
      outfile->Printf("        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
      outfile->Printf("        MP2 correlation energy:          %20.12lf\n",emp2);
      outfile->Printf("      * MP2 total energy:                %20.12lf\n",emp2+escf);
      outfile->Printf("\n");
      outfile->Printf("        OS MP2.5 correlation energy:     %20.12lf\n",emp2_os/emp2_os_fac+0.5*emp3_os);
      outfile->Printf("        SS MP2.5 correlation energy:     %20.12lf\n",emp2_ss/emp2_ss_fac+0.5*emp3_ss);
      outfile->Printf("        MP2.5 correlation energy:        %20.12lf\n",emp2+0.5*emp3);
      outfile->Printf("      * MP2.5 total energy:              %20.12lf\n",emp2+0.5*emp3+escf);
      outfile->Printf("\n");
      outfile->Printf("        OS MP3 correlation energy:       %20.12lf\n",emp2_os/emp2_os_fac+emp3_os);
      outfile->Printf("        SS MP3 correlation energy:       %20.12lf\n",emp2_ss/emp2_ss_fac+emp3_ss);
      outfile->Printf("        MP3 correlation energy:          %20.12lf\n",emp2+emp3);
      outfile->Printf("      * MP3 total energy:                %20.12lf\n",emp2+emp3+escf);
      outfile->Printf("\n");
      outfile->Printf("        OS MP4(SDQ) correlation energy:  %20.12lf\n",emp2_os/emp2_os_fac+emp3_os+emp4_sd_os+emp4_q_os);
      outfile->Printf("        SS MP4(SDQ) correlation energy:  %20.12lf\n",emp2_ss/emp2_ss_fac+emp3_ss+emp4_sd_ss+emp4_q_os);
      outfile->Printf("        MP4(SDQ) correlation energy:     %20.12lf\n",emp2+emp3+emp4_sd+emp4_q);
      outfile->Printf("      * MP4(SDQ) total energy:           %20.12lf\n",emp2+emp3+emp4_sd+emp4_q+escf);
      outfile->Printf("\n");
  }else if (mp3_only){
      if ( options_.get_bool("NAT_ORBS") ) {
          double delta_emp2 = Process::environment.globals["MP2 CORRELATION ENERGY"] - emp2;
          double delta_emp2_os = Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] - emp2_os;
          double delta_emp2_ss = Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] - emp2_ss;

          emp2 += delta_emp2;
          emp2_os += delta_emp2_os;
          emp2_ss += delta_emp2_ss;

          outfile->Printf("\n");
          outfile->Printf("        OS MP2 FNO correction:          %20.12lf\n",delta_emp2_os);
          outfile->Printf("        SS MP2 FNO correction:          %20.12lf\n",delta_emp2_ss);
          outfile->Printf("        MP2 FNO correction:             %20.12lf\n",delta_emp2);
      }

      outfile->Printf("\n");
      outfile->Printf("        OS MP2 correlation energy:       %20.12lf\n",emp2_os);
      outfile->Printf("        SS MP2 correlation energy:       %20.12lf\n",emp2_ss);
      outfile->Printf("        MP2 correlation energy:          %20.12lf\n",emp2);
      outfile->Printf("      * MP2 total energy:                %20.12lf\n",emp2+escf);
      outfile->Printf("\n");
      outfile->Printf("        OS MP2.5 correlation energy:     %20.12lf\n",emp2_os/emp2_os_fac+0.5*emp3_os);
      outfile->Printf("        SS MP2.5 correlation energy:     %20.12lf\n",emp2_ss/emp2_ss_fac+0.5*emp3_ss);
      outfile->Printf("        MP2.5 correlation energy:        %20.12lf\n",emp2+0.5*emp3);
      outfile->Printf("      * MP2.5 total energy:              %20.12lf\n",emp2+0.5*emp3+escf);
      outfile->Printf("\n");
      outfile->Printf("        OS MP3 correlation energy:       %20.12lf\n",emp2_os/emp2_os_fac+emp3_os);
      outfile->Printf("        SS MP3 correlation energy:       %20.12lf\n",emp2_ss/emp2_ss_fac+emp3_ss);
      outfile->Printf("        MP3 correlation energy:          %20.12lf\n",emp2+emp3);
      outfile->Printf("      * MP3 total energy:                %20.12lf\n",emp2+emp3+escf);
      outfile->Printf("\n");
  }else {
     // guess for cc/qci should be |1> + |2>
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"first",(char*)&tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_DCC_T2,"second",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
     C_DAXPY(o*o*v*v,1.0,tempt,1,tb,1);
  }
}

}} // end of namespaces
