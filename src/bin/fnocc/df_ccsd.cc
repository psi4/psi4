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

namespace psi{ namespace fnocc{

// diagrams for mp3 and mp4
void DefineLinearTasks();
void DefineQuadraticTasks();

// sort
void SortIntegrals(int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,Options&options,bool iscim);
void Sort_OV3_LowMemory(long int memory,long int o,long int v,bool islocal);

void DFCoupledCluster::Local_SCS_MP2(){

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
  F_DGEMM('N','N',o,v*v*o,o,1.0,&Rii_pointer[0][0],nocc,tb,o,0.0,tempt,o);
  // resort t(abji) -> t(abij)
  for (int ab = 0; ab < v*v; ab++){
      for (int i = 0; i < o; i++){
          for (int j = 0; j < o; j++){
              I1[i*o+j] = tempt[ab*o*o+j*o+i];
          }
      }
      C_DCOPY(o*o,I1,1,tempt+ab*o*o,1);
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
  emp2_os = osenergy;
  emp2_ss = ssenergy;
  emp2 = emp2_os + emp2_ss;

  psio.reset();
}
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
      C_DCOPY(o*o,I1,1,tempt+ab*o*o,1);
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

      if (!isLowMemory && !reference_wavefunction_->isCIM() ) {
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

      // now there should be space for t2
      if (t2_on_disk){
          tb = (double*)malloc(o*o*v*v*sizeof(double));
          psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
          psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tb[0],o*o*v*v*sizeof(double));
          psio->close(PSIF_DCC_T2,1);
      }

      tstart();

      ccmethod = 0;
      if (isLowMemory)                           status = lowmemory_triples();
      else if (reference_wavefunction_->isCIM()) status = local_triples();
      else                                       status = triples();

      if (status == Failure){
         throw PsiException(
            "Whoops, the (T) correction died.",__FILE__,__LINE__);
      }
      tstop();

      // if we allocated t2 just for triples, free it
      if (t2_on_disk){
          free(tb);
      }

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
  if (!t2_on_disk) {
      free(tb);
  }

  return Process::environment.globals["CURRENT ENERGY"];
}

void DFCoupledCluster::WriteBanner(){
  fflush(outfile);
  psi::fprintf(outfile,"\n\n");
  psi::fprintf(outfile, "        *******************************************************\n");
  psi::fprintf(outfile, "        *                                                     *\n");
  psi::fprintf(outfile, "        *                       DF-CCSD                       *\n");
  psi::fprintf(outfile, "        *                 Density-fitted CCSD                 *\n");
  psi::fprintf(outfile, "        *                                                     *\n");
  psi::fprintf(outfile, "        *                   Eugene DePrince                   *\n");
  psi::fprintf(outfile, "        *                                                     *\n");
  psi::fprintf(outfile, "        *******************************************************\n");
  psi::fprintf(outfile,"\n\n");
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

  psi::fprintf(outfile,"\n");
  psi::fprintf(outfile,"  Begin singles and doubles coupled cluster iterations\n\n");
  psi::fprintf(outfile,"   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");
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
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,o*o*v*v+o*v,replace_diis_iter);
         else                   DIIS(diisvec,maxdiis,o*o*v*v+o*v,replace_diis_iter);
         DIISNewAmplitudes(diis_iter,replace_diis_iter);
      }

      double start;
      if (timer) start = omp_get_wtime();
      T1Fock();
      T1Integrals();
      if (timer) {
          psi::fprintf(outfile,"        T1-transformed integrals                                        %6.2lf\n",omp_get_wtime() - start);
          psi::fprintf(outfile,"\n");
      }

      Eold = eccsd;
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

      if (diis_iter<=maxdiis) diis_iter++;
      else if (replace_diis_iter<maxdiis) replace_diis_iter++;
      else replace_diis_iter = 1;

      time_t iter_stop = time(NULL);
      psi::fprintf(outfile,"  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,eccsd,eccsd-Eold,nrm,(int)iter_stop-(int)iter_start);
      fflush(outfile);
      iter++;
      if (iter==1){
         emp2 = eccsd;
         if ( reference_wavefunction_->isCIM() ) {
             Local_SCS_MP2();
         }else {
             SCS_MP2();
         }
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

  psi::fprintf(outfile,"\n");
  psi::fprintf(outfile,"  CCSD iterations converged!\n");
  psi::fprintf(outfile,"\n");

  // T1 and D1 diagnostics:

  double t1diag = C_DNRM2(o*v,t1,1) / sqrt(2.0 * o);
  psi::fprintf(outfile,"        T1 diagnostic:                  %20.12lf\n",t1diag);
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
  psi::fprintf(outfile,"        D1 diagnostic:                  %20.12lf\n",sqrt(eigval->pointer()[0]));
  psi::fprintf(outfile,"\n");

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

      psi::fprintf(outfile,"        OS MP2 FNO correction:          %20.12lf\n",delta_emp2_os);
      psi::fprintf(outfile,"        SS MP2 FNO correction:          %20.12lf\n",delta_emp2_ss);
      psi::fprintf(outfile,"        MP2 FNO correction:             %20.12lf\n",delta_emp2);
      psi::fprintf(outfile,"\n");
  }

  if (options_.get_bool("SCS_MP2")){
      psi::fprintf(outfile,"        OS SCS-MP2 correlation energy:  %20.12lf\n",emp2_os*emp2_os_fac);
      psi::fprintf(outfile,"        SS SCS-MP2 correlation energy:  %20.12lf\n",emp2_ss*emp2_ss_fac);
      psi::fprintf(outfile,"        SCS-MP2 correlation energy:     %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac);
      psi::fprintf(outfile,"      * SCS-MP2 total energy:           %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac+escf);
      psi::fprintf(outfile,"\n");
  }
  psi::fprintf(outfile,"        OS MP2 correlation energy:      %20.12lf\n",emp2_os);
  psi::fprintf(outfile,"        SS MP2 correlation energy:      %20.12lf\n",emp2_ss);
  psi::fprintf(outfile,"        MP2 correlation energy:         %20.12lf\n",emp2);
  psi::fprintf(outfile,"      * MP2 total energy:               %20.12lf\n",emp2+escf);
  psi::fprintf(outfile,"\n");
  if (options_.get_bool("SCS_CCSD")){
      psi::fprintf(outfile,"        OS SCS-CCSD correlation energy: %20.12lf\n",eccsd_os*eccsd_os_fac);
      psi::fprintf(outfile,"        SS SCS-CCSD correlation energy: %20.12lf\n",eccsd_ss*eccsd_ss_fac);
      psi::fprintf(outfile,"        SCS-CCSD correlation energy:    %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
      psi::fprintf(outfile,"      * SCS-CCSD total energy:          %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
      psi::fprintf(outfile,"\n");
  }
  psi::fprintf(outfile,"        OS CCSD correlation energy:     %20.12lf\n",eccsd_os);
  psi::fprintf(outfile,"        SS CCSD correlation energy:     %20.12lf\n",eccsd_ss);
  psi::fprintf(outfile,"        CCSD correlation energy:        %20.12lf\n",eccsd);
  psi::fprintf(outfile,"      * CCSD total energy:              %20.12lf\n",eccsd+escf);
  psi::fprintf(outfile,"\n");

  psi::fprintf(outfile,"  Total time for CCSD iterations: %10.2lf s (user)\n",user_stop-user_start);
  psi::fprintf(outfile,"                                  %10.2lf s (system)\n",sys_stop-sys_start);
  psi::fprintf(outfile,"                                  %10d s (total)\n",(int)time_stop-(int)time_start);
  psi::fprintf(outfile,"\n");
  psi::fprintf(outfile,"  Time per iteration:             %10.2lf s (user)\n",(user_stop-user_start)/(iter-1));
  psi::fprintf(outfile,"                                  %10.2lf s (system)\n",(sys_stop-sys_start)/(iter-1));
  psi::fprintf(outfile,"                                  %10.2lf s (total)\n",((double)time_stop-(double)time_start)/(iter-1));
  fflush(outfile);

  if (options_.get_bool("COMPUTE_TRIPLES")){
      // need to generate non-t1-transformed 3-index integrals
      C_DCOPY(o*v,t1,1,w1,1);
      memset((void*)t1,'\0',o*v*sizeof(double));
      T1Fock();
      T1Integrals();
      C_DCOPY(o*v,w1,1,t1,1);
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
        C_DCOPY(nso*full,&Catemp[0],1,Ca_L,1);
        C_DCOPY(nso*full,&Catemp[0],1,Ca_R,1);
    }else {
        C_DCOPY(nso*full,&Ca[0][0],1,Ca_L,1);
        C_DCOPY(nso*full,&Ca[0][0],1,Ca_R,1);
        C_DCOPY(nso*full,&Ca[0][0],1,Catemp,1);
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
    psio_address addr1  = PSIO_ZERO;
    psio_address addr2  = PSIO_ZERO;
    psio_address addroo = PSIO_ZERO;
    psio_address addrov = PSIO_ZERO;
    psio_address addrvo = PSIO_ZERO;
    psio_address addrvv = PSIO_ZERO;

    long int nrows = 1;
    long int rowsize = nQ_scf;
    while ( rowsize*nso*nso > o*o*v*v ) {
        nrows++;
        rowsize = nQ_scf / nrows;
        if (nrows * rowsize < nQ_scf) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ_scf - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso SCF",(char*)&integrals[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,Ca_L,full,integrals,nso,0.0,tempv,full);
        for (int q = 0; q < rowdims[row]; q++) {
            for (int mu = 0; mu < nso; mu++) {
                C_DCOPY(full,tempv+q*nso*full+mu*full,1,integrals+q*nso*full+mu,nso);
            }
        }
        F_DGEMM('n','n',full,full*rowdims[row],nso,1.0,Ca_R,full,integrals,nso,0.0,tempv,full);
        // full Qmo
        psio->write(PSIF_DCC_QSO,"Qmo SCF",(char*)&tempv[0],rowdims[row]*full*full*sizeof(double),addr2,&addr2);
    }
    delete rowdims;

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

    double * temp3 = (double*)malloc(full*full*sizeof(double));

    memset((void*)temp3,'\0',full*full*sizeof(double));
    psio_address addr = PSIO_ZERO;

    nrows = 1;
    rowsize = nQ_scf;
    while ( rowsize*full*full > o*o*v*v ) {
        nrows++;
        rowsize = nQ_scf / nrows;
        if (nrows * rowsize < nQ_scf) rowsize++;
        if (rowsize == 1) break;
    }
    lastrowsize = nQ_scf - (nrows - 1L) * rowsize;
    rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qmo SCF",(char*)&integrals[0],rowdims[row]*full*full*sizeof(double),addr,&addr);
        for (int q = 0; q < rowdims[row]; q++) {
            // sum k (q|rk) (q|ks)
            F_DGEMM('n','n',full,full,ndocc,-1.0,integrals+q*full*full,full,integrals+q*full*full,full,1.0,temp3,full);

            // sum k (q|kk) (q|rs)
            double dum = 0.0;
            for (int k = 0; k < ndocc; k++) {
                dum += integrals[q*full*full+k*full + k];
            }
            F_DAXPY(full*full,2.0 * dum,integrals+q*full*full,1,temp3,1);
        }
    }
    delete rowdims;
    psio->close(PSIF_DCC_QSO,1);

    // Fij
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            Fij[i*o+j] = h[i*nmo+j] + temp3[(i+nfzc)*full+(j+nfzc)];
        }
    }

    // Fia
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            Fia[i*v+a] = h[i*nmo+a+o] + temp3[(i+nfzc)*full+(a+ndocc)];
        }
    }

    // Fai
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            Fai[a*o+i] = h[(a+o)*nmo+i] + temp3[(a+ndocc)*full+(i+nfzc)];
        }
    }

    // Fab
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            Fab[a*v+b] = h[(a+o)*nmo+b+o] + temp3[(a+ndocc)*full+(b+ndocc)];
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
        C_DCOPY(nso*full,&Catemp[0],1,Ca_L,1);
        C_DCOPY(nso*full,&Catemp[0],1,Ca_R,1);
    }else {
        C_DCOPY(nso*full,&Ca[0][0],1,Ca_L,1);
        C_DCOPY(nso*full,&Ca[0][0],1,Ca_R,1);
        C_DCOPY(nso*full,&Ca[0][0],1,Catemp,1);
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
    psio_address addr1  = PSIO_ZERO;
    psio_address addrvo = PSIO_ZERO;
    long int nrows = 1;
    long int rowsize = nQ;
    while ( rowsize*nso*nso > o*o*v*v ) {
        nrows++;
        rowsize = nQ / nrows;
        if (nrows * rowsize < nQ) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso CC",(char*)&integrals[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,Ca_L,full,integrals,nso,0.0,tempv,full);
        for (int q = 0; q < rowdims[row]; q++) {
            for (int mu = 0; mu < nso; mu++) {
                C_DCOPY(full,tempv+q*nso*full+mu*full,1,integrals+q*nso*full+mu,nso);
            }
        }
        F_DGEMM('n','n',full,full*rowdims[row],nso,1.0,Ca_R,full,integrals,nso,0.0,tempv,full);

        // Qoo
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    Qoo[(q+rowdims[0]*row)*o*o+i*o+j] = tempv[q*full*full+(i+nfzc)*full+(j+nfzc)];
                }
            }
        }
        // Qov
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    Qov[(q+rowdims[0]*row)*o*v+i*v+a] = tempv[q*full*full+(i+nfzc)*full+(a+ndocc)];
                }
            }
        }
        // Qvo
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int a = 0; a < v; a++) {
                for (int i = 0; i < o; i++) {
                    integrals[q*o*v+a*o+i] = tempv[q*full*full+(a+ndocc)*full+(i+nfzc)];
                }
            }
        }
        psio->write(PSIF_DCC_QSO,"qvo",(char*)&integrals[0],rowdims[row]*o*v*sizeof(double),addrvo,&addrvo);
        // Qvv
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    Qvv[(q+rowdims[0]*row)*v*v+a*v+b] = tempv[q*full*full+(a+ndocc)*full+(b+ndocc)];
                }
            }
        }
    }
    delete rowdims;
    psio->close(PSIF_DCC_QSO,1);
}

double DFCoupledCluster::CheckEnergy(){

    long int v  = nvirt;
    long int o  = ndoccact;
    long int rs = nmo;
  
    // df (ia|bj) formerly E2klcd
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

    if (t2_on_disk){
        boost::shared_ptr<PSIO> psio (new PSIO());
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
                    long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;
                    long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                    long int jaib = iajb + (i-j)*v*(1-v*o);
                    energy += (2.*integrals[iajb]-integrals[jaib])*(tb[ijab]+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                }
            }
        }
    }
    return energy;
}
void DFCoupledCluster::SCS_MP2(){

    long int v  = nvirt;
    long int o  = ndoccact;
    long int rs = nmo;
  
    double ssenergy = 0.0;
    double osenergy = 0.0;
  
    // df (ia|bj) formerly E2klcd
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

    if (t2_on_disk){
        boost::shared_ptr<PSIO> psio (new PSIO());
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }

    for (long int a = o; a < rs; a++){
        for (long int b = o; b < rs; b++){
            for (long int i = 0; i < o; i++){
                for (long int j = 0; j < o; j++){

                    long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;
                    long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                    long int jaib = iajb + (i-j)*v*(1-v*o);
                    osenergy += integrals[iajb]*tb[ijab];
                    ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                }
            }
        }
    }
    emp2_os = osenergy;
    emp2_ss = ssenergy;
    emp2    = emp2_os + emp2_ss;
}
void DFCoupledCluster::SCS_CCSD(){

    long int v  = nvirt;
    long int o  = ndoccact;
    long int rs = nmo;
  
    double ssenergy = 0.0;
    double osenergy = 0.0;
  
    // df (ia|bj) formerly E2klcd
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

    if (t2_on_disk){
        boost::shared_ptr<PSIO> psio (new PSIO());
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }

    for (long int a = o; a < rs; a++){
        for (long int b = o; b < rs; b++){
            for (long int i = 0; i < o; i++){
                for (long int j = 0; j < o; j++){

                    long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;
                    long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                    long int jaib = iajb + (i-j)*v*(1-v*o);

                    osenergy += integrals[iajb]*(tb[ijab]+t1[(a-o)*o+i]*t1[(b-o)*o+j]);
                    ssenergy += integrals[iajb]*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                    ssenergy += integrals[iajb]*(t1[(a-o)*o+i]*t1[(b-o)*o+j]-t1[(b-o)*o+i]*t1[(a-o)*o+j]);
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

  long int tempvdim = o*o*v*v+o*v;
  if ( nQ * o * v > tempvdim) tempvdim = nQ * o * v;
  if ( nso * nso > tempvdim)  tempvdim = nso * nso;

  double total_memory = dim+tempvdim+(o*(o+1)*v*(v+1)+o*v)+o*o*v*v+2.*o*v+2.*v*v;
  long int max = nvirt*nvirt*nQmax > (nfzv+ndocc+nvirt)*ndocc*nQmax ? nvirt*nvirt*nQmax : (nfzv+ndocc+nvirt)*ndocc*nQmax;
  double df_memory    = nQ*(o*o+o*v)+max;

  total_memory       *= 8./1024./1024.;
  df_memory          *= 8./1024./1024.;


  double available_memory = (double)memory/1024.0/1024.0;
  double size_of_t2       = 8.0*o*o*v*v/1024.0/1024.0;

  if (available_memory < total_memory + df_memory) {

      if ( available_memory > total_memory + df_memory - size_of_t2) {
          psi::fprintf(outfile,"\n");
          psi::fprintf(outfile,"        Warning: cannot accomodate T2 in core. T2 will be stored on disk.\n");
          psi::fprintf(outfile,"\n");
          fflush(outfile);
          t2_on_disk = true;
      } else {
          psi::fprintf(outfile,"\n");
          psi::fprintf(outfile,"        error: not enough memory for ccsd.  increase available memory by %7.2lf mb\n",
                          total_memory + df_memory - size_of_t2 - available_memory);
          psi::fprintf(outfile,"\n");
          fflush(outfile);
          throw PsiException("not enough memory (ccsd).",__FILE__,__LINE__);
      }

  }

  psi::fprintf(outfile,"  ==> Memory <==\n\n");
  psi::fprintf(outfile,"        Total memory requirements:       %9.2lf mb\n",df_memory+total_memory-size_of_t2*t2_on_disk);
  psi::fprintf(outfile,"        3-index integrals:               %9.2lf mb\n",df_memory);
  psi::fprintf(outfile,"        CCSD intermediates:              %9.2lf mb\n",total_memory-size_of_t2*t2_on_disk);
  psi::fprintf(outfile,"\n");

  if (options_.get_bool("COMPUTE_TRIPLES")) {
      long int nthreads = omp_get_max_threads();
      double tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+3L*v*v*v*nthreads);
      if (tempmem > memory) {
          psi::fprintf(outfile,"\n        <<< warning! >>> switched to low-memory (t) algorithm\n\n");
      }
      if (tempmem > memory || options_.get_bool("TRIPLES_LOW_MEMORY")){
         isLowMemory = true;
         tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+5L*o*o*o*nthreads);
      }
      psi::fprintf(outfile,"        memory requirements for CCSD(T): %9.2lf mb\n\n",tempmem/1024./1024.);
  }
  psi::fprintf(outfile,"  ==> Input parameters <==\n\n");
  psi::fprintf(outfile,"        Freeze core orbitals?               %5s\n",nfzc > 0 ? "yes" : "no");
  psi::fprintf(outfile,"        Use frozen natural orbitals?        %5s\n",options_.get_bool("NAT_ORBS") ? "yes" : "no");
  psi::fprintf(outfile,"        r_convergence:                  %5.3le\n",r_conv);
  psi::fprintf(outfile,"        e_convergence:                  %5.3le\n",e_conv);
  psi::fprintf(outfile,"        Number of DIIS vectors:             %5li\n",maxdiis);
  psi::fprintf(outfile,"        Number of frozen core orbitals:     %5li\n",nfzc);
  psi::fprintf(outfile,"        Number of active occupied orbitals: %5li\n",ndoccact);
  psi::fprintf(outfile,"        Number of active virtual orbitals:  %5li\n",nvirt);
  psi::fprintf(outfile,"        Number of frozen virtual orbitals:  %5li\n",nfzv);
  psi::fprintf(outfile,"\n");


  // allocate some memory for 3-index tensors
  Qoo = (double*)malloc(ndoccact*ndoccact*nQ*sizeof(double));
  Qov = (double*)malloc(ndoccact*nvirt*nQ*sizeof(double));
  // max (v*v*nQ, full*ndocc*nQ)
  Qvv = (double*)malloc(max*sizeof(double));

  integrals = (double*)malloc(dim*sizeof(double));
  tempt     = (double*)malloc((o*(o+1)*v*(v+1)+o*v)*sizeof(double));
  tempv     = (double*)malloc(tempvdim*sizeof(double));
  Abij      = (double*)malloc(o*(o+1)/2*v*sizeof(double));
  Sbij      = (double*)malloc(o*(o+1)/2*v*sizeof(double));
  if (!t2_on_disk) {
      tb        = (double*)malloc(o*o*v*v*sizeof(double));
  }
  w1        = (double*)malloc(o*v*sizeof(double));
  t1        = (double*)malloc(o*v*sizeof(double));
  I1        = (double*)malloc(v*v*sizeof(double));
  I1p       = (double*)malloc(v*v*sizeof(double));

  memset((void*)integrals,'\0',dim*sizeof(double));
  memset((void*)tempv,'\0',tempvdim*sizeof(double));
  memset((void*)tempt,'\0',(o*(o+1)*v*(v+1)+o*v)*sizeof(double));
  if (!t2_on_disk) {
      memset((void*)tb,'\0',o*o*v*v*sizeof(double));
  }
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
  Ca_R  = (double*)malloc(nso*(nmo+nfzc+nfzv)*sizeof(double));
  Ca_L  = (double*)malloc(nso*(nmo+nfzc+nfzv)*sizeof(double));
  Ca = reference_wavefunction_->Ca()->pointer();

  // one-electron integrals
  boost::shared_ptr<MintsHelper> mints(new MintsHelper());
  H = mints->so_kinetic();
  H->add(mints->so_potential());

  if (t2_on_disk) {
     boost::shared_ptr<PSIO> psio (new PSIO());
     psio->open(PSIF_DCC_T2,PSIO_OPEN_NEW);
     psio->write_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_DCC_T2,1);
  }
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
  C_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  F_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  C_DCOPY(o*v,w1,1,t1,1);
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
                  tempt[ijab]  = tnew;
              }
          }
      }
  }
    // error vector is just dt
    C_DCOPY(o*o*v*v,tempt,1,tempv,1);

    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        F_DAXPY(o*o*v*v,1.0,tempt,1,integrals,1);
        psio->write_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
    }else {
        F_DAXPY(o*o*v*v,1.0,tempt,1,tb,1);
    }
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

    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }

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
  
    double * Vcdb = integrals;
    double * Vm   = integrals+v*v*v;
    double * Vp   = Vm;
  
    // qvv transpose
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        C_DCOPY(v*v,Qvv+q*v*v,1,integrals+q,nQ);
    }
    C_DCOPY(nQ*v*v,integrals,1,Qvv,1);
  
    double time1 = 0.0;
    double time2 = 0.0;
    double time3 = 0.0;
    for (long int a = 0; a < v; a++) {
  
        double start1 = omp_get_wtime();
        int nb = v-a;
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
                    cd++;
                }
            }
        }
        double end1 = omp_get_wtime();
  
        double start2 = omp_get_wtime();
        F_DGEMM('n','n',otri,nb,vtri,0.5,tempt,otri,Vp,vtri,0.0,Abij,otri);
        #pragma omp parallel for schedule (static)
        for (long int b = a; b < v; b++){
            long int cd = 0;
            long int ind1 = (b-a)*vtri;
            long int ind2 = (b-a)*v*v;
            long int v1,v2;
            for (long int c=0; c<v; c++){
                for (long int d=0; d<=c; d++){
                    Vm[ind1+cd] = Vcdb[ind2+d*v+c] - Vcdb[ind2+c*v+d];
                    cd++;
                }
            }
        }
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
        C_DCOPY(v*v,Qvv+q,nQ,integrals+q*v*v,1);
    }
    C_DCOPY(nQ*v*v,integrals,1,Qvv,1);
}

void DFCoupledCluster::CCResidual(){
    bool timer = options_.get_bool("CC_TIMINGS");
    long int o = ndoccact;
    long int v = nvirt;

    double start;

    boost::shared_ptr<PSIO> psio (new PSIO());

    // C2 = -1/2 t(bc,kj) [ (ki|ac) - 1/2 t(ad,li) (kd|lc) ] 
    //      +    t(bc,ki) [ (kj|ac) - 1/2 t(ad,lj) (kd|lc) ] 
    if (timer) start = omp_get_wtime();
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
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
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
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
    psio->open(PSIF_DCC_R2,PSIO_OPEN_NEW);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);
    if (timer) {
        psi::fprintf(outfile,"\n");
        psi::fprintf(outfile,"        C2 = -1/2 t(b,c,k,j) [ (ki|ac) - 1/2 t(a,d,l,i) (kd|lc) ]\n");
        psi::fprintf(outfile,"                + t(b,c,k,i) [ (kj|ac) - 1/2 t(a,d,l,j) (kd|lc) ]       %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // D2: 1/2 U(b,c,j,k) [ L(a,i,k,c) + 1/2 U(a,d,i,l) L(l,d,k,c) ] 
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
    C_DCOPY(o*o*v*v,integrals,1,tempv,1);
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
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = integrals;
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
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
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
        psi::fprintf(outfile,"        D2 =  1/2 U(b,c,j,k) [ L(a,i,k,c) + 1/2 U(a,d,i,l) L(l,d,k,c) ] %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // E2 a: t(ac,ij) [ F(bc) - U(bd,kl) (ld|kc) ]
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    C_DCOPY(o*o*v*v,tb,1,tempt,1);
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
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
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
        psi::fprintf(outfile,"        E2 =      t(a,c,i,j) [ F(b,c) - U(b,d,k,l) (ld|kc) ]            %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // E2 b: -t(a,b,i,k) [ F(kj) - U(c,d,l,j) (kd|lc) ]
    // note that (kd|lc) should still be in integrals buffer
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
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
    C_DCOPY(o*o*v*v,tempt,1,integrals,1);
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
        psi::fprintf(outfile,"                - t(a,b,i,k) [ F(k,j) - U(c,d,l,j) (kd|lc) ]            %6.2lf\n",omp_get_wtime()-start);
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
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = integrals;
    }
    F_DGEMM('n','n',o*o,o*o,v*v,1.0,tb,o*o,tempv,v*v,1.0,tempt,o*o);
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    F_DGEMM('n','n',o*o,v*v,o*o,1.0,tempt,o*o,tb,o*o,0.0,integrals,o*o);

    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    F_DAXPY(o*o*v*v,1.0,tempt,1,integrals,1);
    psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);

    if (timer) {
        psi::fprintf(outfile,"        B2 =      t(a,b,k,l) [ (ki|lj) + t(c,d,i,j) (kc|ld) ]           %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // now singles residual:

    // D1: F(ai)
    C_DCOPY(o*v,Fai,1,w1,1);
 
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
            C_DCOPY(v,Qvv+q*v*v+a*v,1,integrals+q*v*v+a,v);
        }
    }
    F_DGEMM('n','t',o,v,v*nQ,1.0,tempv,o,integrals,v,1.0,w1,o);

    if (timer) {
        psi::fprintf(outfile,"        A1 =      U(c,d,k,l) (ad|kc)                                    %6.2lf\n",omp_get_wtime()-start);
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
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = integrals;
    }
    C_DCOPY(o*o*v*v,tb,1,tempt,1);
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
        psi::fprintf(outfile,"        B1 =    - U(a,c,k,l) (ki|lc)                                    %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    // C1
    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
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
        psi::fprintf(outfile,"        C1 =      F(k,c) U(a,c,i,k)                                     %6.2lf\n",omp_get_wtime()-start);
        start = omp_get_wtime();
    }

    Vabcd1();
    if (timer) {
        psi::fprintf(outfile,"        A2 =      t(c,d,i,j) (ac|bd)                                    %6.2lf\n",omp_get_wtime()-start);
    }
}


}} // end of namespaces
