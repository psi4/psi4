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
#include "psi4/libmints/mintshelper.h"
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


namespace psi{ namespace fnocc{

// diagrams for mp3 and mp4
void DefineLinearTasks();
void DefineQuadraticTasks();

// coupled cluster constructor
DFCoupledCluster::DFCoupledCluster(SharedWavefunction ref_wfn, Options &options):
        CoupledCluster(ref_wfn, options)
{
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

      if ( !isLowMemory ) {
          // write (ov|vv) integrals, formerly E2abci, for (t)
          double *tempq = (double*)malloc(v*nQ*sizeof(double));
          // the buffer integrals was at least 2v^3, so these should definitely fit.
          double *Z     = (double*)malloc(v*v*v*sizeof(double));
          double *Z2    = (double*)malloc(v*v*v*sizeof(double));
          std::shared_ptr<PSIO> psio(new PSIO());
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
          std::shared_ptr<PSIO> psio(new PSIO());
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
      for (long int i=0; i<o; i++){
          for (long int j=0; j<o; j++){
              for (long int k=0; k<o; k++){
                  for (long int a=0; a<v; a++){
                      temp2[j*o*o*v+i*o*v+k*v+a] = temp1[i*o*o*v+a*o*o+j*o+k];
                  }
              }
          }
      }
      std::shared_ptr<PSIO> psio(new PSIO());
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

  outfile->Printf("\n\n");
  outfile->Printf( "        *******************************************************\n");
  outfile->Printf( "        *                                                     *\n");
  outfile->Printf( "        *                       DF-CCSD                       *\n");
  outfile->Printf( "        *                 Density-fitted CCSD                 *\n");
  outfile->Printf( "        *                                                     *\n");
  outfile->Printf( "        *                   Eugene DePrince                   *\n");
  outfile->Printf( "        *                                                     *\n");
  outfile->Printf( "        *******************************************************\n");
  outfile->Printf("\n\n");

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

  std::shared_ptr<PSIO> psio(new PSIO());
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

  outfile->Printf("\n");
  outfile->Printf("  Begin singles and doubles coupled cluster iterations\n\n");
  outfile->Printf("   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");


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
          outfile->Printf("        T1-transformed integrals                                        %6.2lf\n",omp_get_wtime() - start);
          outfile->Printf("\n");
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
      outfile->Printf("  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,eccsd,eccsd-Eold,nrm,(int)iter_stop-(int)iter_start);

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
  SCS_CCSD();

  outfile->Printf("\n");
  outfile->Printf("  CCSD iterations converged!\n");
  outfile->Printf("\n");

  // T1 and D1 diagnostics:

  double t1diag = C_DNRM2(o*v,t1,1) / sqrt(2.0 * o);
  outfile->Printf("        T1 diagnostic:                  %20.12lf\n",t1diag);
  std::shared_ptr<Matrix>T (new Matrix(o,o));
  std::shared_ptr<Matrix>eigvec (new Matrix(o,o));
  std::shared_ptr<Vector>eigval (new Vector(o));
  double ** Tp = T->pointer();
  for (long int i = 0; i < o; i++) {
      for (long int j = 0; j < o; j++) {
          double dum = 0.0;
          for (long int a = 0; a < v; a++) {
              dum += t1[a*o+i] * t1[a*o+j];
          }
          Tp[i][j] = dum;
      }
  }
  T->diagonalize(eigvec,eigval,descending);
  outfile->Printf("        D1 diagnostic:                  %20.12lf\n",sqrt(eigval->pointer()[0]));
  outfile->Printf("\n");

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

      outfile->Printf("        OS MP2 FNO correction:          %20.12lf\n",delta_emp2_os);
      outfile->Printf("        SS MP2 FNO correction:          %20.12lf\n",delta_emp2_ss);
      outfile->Printf("        MP2 FNO correction:             %20.12lf\n",delta_emp2);
      outfile->Printf("\n");
  }

  if (options_.get_bool("SCS_MP2")){
      outfile->Printf("        OS SCS-MP2 correlation energy:  %20.12lf\n",emp2_os*emp2_os_fac);
      outfile->Printf("        SS SCS-MP2 correlation energy:  %20.12lf\n",emp2_ss*emp2_ss_fac);
      outfile->Printf("        SCS-MP2 correlation energy:     %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac);
      outfile->Printf("      * SCS-MP2 total energy:           %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac+escf);
      outfile->Printf("\n");
  }
  outfile->Printf("        OS MP2 correlation energy:      %20.12lf\n",emp2_os);
  outfile->Printf("        SS MP2 correlation energy:      %20.12lf\n",emp2_ss);
  outfile->Printf("        MP2 correlation energy:         %20.12lf\n",emp2);
  outfile->Printf("      * MP2 total energy:               %20.12lf\n",emp2+escf);
  outfile->Printf("\n");
  if (options_.get_bool("SCS_CCSD")){
      outfile->Printf("        OS SCS-CCSD correlation energy: %20.12lf\n",eccsd_os*eccsd_os_fac);
      outfile->Printf("        SS SCS-CCSD correlation energy: %20.12lf\n",eccsd_ss*eccsd_ss_fac);
      outfile->Printf("        SCS-CCSD correlation energy:    %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
      outfile->Printf("      * SCS-CCSD total energy:          %20.12lf\n",eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
      outfile->Printf("\n");
  }
  outfile->Printf("        OS CCSD correlation energy:     %20.12lf\n",eccsd_os);
  outfile->Printf("        SS CCSD correlation energy:     %20.12lf\n",eccsd_ss);
  outfile->Printf("        CCSD correlation energy:        %20.12lf\n",eccsd);
  outfile->Printf("      * CCSD total energy:              %20.12lf\n",eccsd+escf);
  outfile->Printf("\n");

  outfile->Printf("  Total time for CCSD iterations: %10.2lf s (user)\n",user_stop-user_start);
  outfile->Printf("                                  %10.2lf s (system)\n",sys_stop-sys_start);
  outfile->Printf("                                  %10d s (total)\n",(int)time_stop-(int)time_start);
  outfile->Printf("\n");
  outfile->Printf("  Time per iteration:             %10.2lf s (user)\n",(user_stop-user_start)/(iter-1));
  outfile->Printf("                                  %10.2lf s (system)\n",(sys_stop-sys_start)/(iter-1));
  outfile->Printf("                                  %10.2lf s (total)\n",((double)time_stop-(double)time_start)/(iter-1));


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

double DFCoupledCluster::CheckEnergy(){

    long int v  = nvirt;
    long int o  = ndoccact;

    // df (ia|bj) formerly E2klcd
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

    if (t2_on_disk){
        std::shared_ptr<PSIO> psio (new PSIO());
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    double energy = 0.0;
    for (long int a = 0; a < v; a++){
        for (long int b = 0; b < v; b++){
            for (long int i = 0; i < o; i++){
                for (long int j = 0; j < o; j++){
                    long int ijab = a*v*o*o + b*o*o + i*o + j;
                    long int iajb = i*v*v*o + a*v*o + j*v + b;
                    long int jaib = j*v*v*o + a*v*o + i*v + b;
                    energy += (2.0*integrals[iajb]-integrals[jaib])*(tb[ijab]+t1[a*o+i]*t1[b*o+j]);
                }
            }
        }
    }
    return energy;
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
          outfile->Printf("\n");
          outfile->Printf("        Warning: cannot accomodate T2 in core. T2 will be stored on disk.\n");
          outfile->Printf("\n");

          t2_on_disk = true;
      } else {
          outfile->Printf("\n");
          outfile->Printf("        error: not enough memory for ccsd.  increase available memory by %7.2lf mb\n",
                          total_memory + df_memory - size_of_t2 - available_memory);
          outfile->Printf("\n");

          throw PsiException("not enough memory (ccsd).",__FILE__,__LINE__);
      }

  }

  outfile->Printf("  ==> Memory <==\n\n");
  outfile->Printf("        Total memory available:          %9.2lf mb\n",available_memory);
  outfile->Printf("\n");
  outfile->Printf("        CCSD memory requirements:        %9.2lf mb\n",df_memory+total_memory-size_of_t2*t2_on_disk);
  outfile->Printf("            3-index integrals:           %9.2lf mb\n",df_memory);
  outfile->Printf("            CCSD intermediates:          %9.2lf mb\n",total_memory-size_of_t2*t2_on_disk);

  if (options_.get_bool("COMPUTE_TRIPLES")) {
      int nthreads = Process::environment.get_n_threads();
      double mem_t = 8.*(2L*o*o*v*v+1L*o*o*o*v+o*v+3L*v*v*v*nthreads);
      outfile->Printf("\n");
      outfile->Printf("        (T) part (regular algorithm):    %9.2lf mb\n",
          mem_t/1024./1024.);
      if (mem_t > memory) {
          outfile->Printf("        <<< warning! >>> switched to low-memory (t) algorithm\n\n");
      }
      if (mem_t > memory || options_.get_bool("TRIPLES_LOW_MEMORY")){
         isLowMemory = true;
         mem_t = 8.*(2L*o*o*v*v+o*o*o*v+o*v+5L*o*o*o*nthreads);
         outfile->Printf("        (T) part (low-memory alg.):      %9.2lf mb\n\n",mem_t/1024./1024.);
      }
  }
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
  outfile->Printf("\n");


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
  std::shared_ptr<MintsHelper> mints(new MintsHelper(basisset_, options_, 0));
  H = mints->so_kinetic();
  H->add(mints->so_potential());

  if (t2_on_disk) {
     std::shared_ptr<PSIO> psio (new PSIO());
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
  C_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  C_DCOPY(o*v,w1,1,t1,1);
}
void DFCoupledCluster::UpdateT2(){

  long int v = nvirt;
  long int o = ndoccact;

  std::shared_ptr<PSIO> psio(new PSIO());

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
  for (long int a = 0; a < v; a++){
      double da = eps[a+o];
      for (long int b = 0; b < v; b++){
          double dab = da + eps[b+o];
          for (long int i = 0; i < o; i++){
              double dabi = dab - eps[i];
              for (long int j = 0; j < o; j++){

                  long int iajb = a*v*o*o + i*v*o + b*o + j;
                  long int ijab = a*v*o*o + b*o*o + i*o + j;

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
      C_DAXPY(o*o*v*v,1.0,tempt,1,integrals,1);
      psio->write_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_T2,1);
  }else {
      C_DAXPY(o*o*v*v,1.0,tempt,1,tb,1);
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

    std::shared_ptr<PSIO> psio(new PSIO());

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

    int nthreads = Process::environment.get_n_threads();

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


}} // end of namespaces
