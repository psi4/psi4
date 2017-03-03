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


namespace psi{ namespace fnocc{

// diagrams for mp3 and mp4
void DefineLinearTasks();
void DefineQuadraticTasks();

// sort
void SortIntegrals(int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,Options&options);

CoupledPair::CoupledPair(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options):
        CoupledCluster(reference_wavefunction, options)
{
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

  outfile->Printf("\n\n");
  outfile->Printf( "        *******************************************************\n");
  outfile->Printf( "        *                                                     *\n");
  if (options_.get_str("CEPA_LEVEL")=="CEPA(0)"){
     outfile->Printf( "        *                       CEPA(0)                       *\n");
     outfile->Printf( "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="CEPA(1)"){
     outfile->Printf( "        *                       CEPA(1)                       *\n");
     outfile->Printf( "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="CEPA(2)"){
     outfile->Printf( "        *                       CEPA(2)                       *\n");
     outfile->Printf( "        *        Coupled Electron Pair Approximation          *\n");
  }
  if (options_.get_str("CEPA_LEVEL")=="CEPA(3)"){
     outfile->Printf( "        *                       CEPA(3)                       *\n");
     outfile->Printf( "        *        Coupled Electron Pair Approximation          *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="ACPF"){
     outfile->Printf( "        *                        ACPF                         *\n");
     outfile->Printf( "        *          Averaged Coupled Pair Functional           *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="AQCC"){
     outfile->Printf( "        *                        AQCC                         *\n");
     outfile->Printf( "        *         Averaged Quadratic Coupled Cluster          *\n");
  }
  else if (options_.get_str("CEPA_LEVEL")=="CISD"){
     outfile->Printf( "        *                        CISD                         *\n");
     outfile->Printf( "        *      Singles Doubles Configuration Interaction      *\n");
  }


  outfile->Printf( "        *                                                     *\n");
  outfile->Printf( "        *                   Eugene DePrince                   *\n");
  outfile->Printf( "        *                                                     *\n");
  outfile->Printf( "        *******************************************************\n");
  outfile->Printf("\n\n");

  WriteOptions();
}

double CoupledPair::compute_energy() {
  PsiReturnType status = Success;

  // integral sort
  tstart();
  SortIntegrals(nfzc,nfzv,nmo+nfzc+nfzv,ndoccact,nvirt,options_);
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
     if (options_.get_bool("CEPA_NO_SINGLES")) {
        Process::environment.globals["LCCD CORRELATION ENERGY"]               = eccsd;
        Process::environment.globals["LCCD OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
        Process::environment.globals["LCCD SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
        Process::environment.globals["LCCD TOTAL ENERGY"]                     = eccsd + escf;
     } else {
        Process::environment.globals["LCCSD CORRELATION ENERGY"]               = eccsd;
        Process::environment.globals["LCCSD OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
        Process::environment.globals["LCCSD SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
        Process::environment.globals["LCCSD TOTAL ENERGY"]                     = eccsd + escf;
        Process::environment.globals["CEPA(0) CORRELATION ENERGY"]               = eccsd;
        Process::environment.globals["CEPA(0) OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
        Process::environment.globals["CEPA(0) SAME-SPIN CORRELATION ENERGY"]     = eccsd_ss;
        Process::environment.globals["CEPA(0) TOTAL ENERGY"]                     = eccsd + escf;
     }
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

  // build opdm in case we want properties.
  if ( cepa_level<=0 ) {
      if (options_.get_bool("NAT_ORBS")) {
          //outfile->Printf("\n");
          //outfile->Printf("\n");
          //outfile->Printf("       <<< Warning >>>  %s OPDM will have no correction for FNO truncation.\n",cepa_type);
          //outfile->Printf("\n");
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

  outfile->Printf("\n");
  outfile->Printf(
    "  Begin %s iterations\n\n",cepa_type);
  outfile->Printf(
    "   Iter  DIIS          Energy       d(Energy)          |d(T)|     time\n");


  std::shared_ptr<PSIO> psio(new PSIO());
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
      Eold = ( cepa_level < 1 ) ? evar : eccsd;
      if ( iter==1 ) Eold = eccsd; // use mp2 energy as Eold on 2nd iteration
      if (cepa_level < 1) {
          evar = VariationalEnergy();
          eccsd = evar; // for updates to t2 and t1
      }

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
         if (diis_iter<maxdiis) DIIS(diisvec,diis_iter,o*o*v*v+o*v,replace_diis_iter);
         else                   DIIS(diisvec,maxdiis,o*o*v*v+o*v,replace_diis_iter);
         DIISNewAmplitudes(diis_iter,replace_diis_iter);
      }
      // if cepa_no_singles, zero t1
      if (options_.get_bool("CEPA_NO_SINGLES")){
         memset((void*)t1,'\0',o*v*sizeof(double));
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

      if (diis_iter<=maxdiis) diis_iter++;
      else if (replace_diis_iter<maxdiis) replace_diis_iter++;
      else replace_diis_iter = 1;

      time_t iter_stop = time(NULL);
      double dume = ( cepa_level < 1 ) ? evar : eccsd ;
      if ( iter==0 ) dume = eccsd; // use mp2 energy on first iteration
      outfile->Printf("  %5i   %i %i %15.10f %15.10f %15.10f %8d\n",
            iter,diis_iter-1,replace_diis_iter,dume,dume-Eold,nrm,(int)iter_stop-(int)iter_start);

      iter++;
      if (iter==1) emp2 = eccsd;
      if (iter==1) SCS_MP2();

      if (fabs(dume - Eold) < e_conv && nrm < r_conv) break;
  }
  times(&total_tmstime);
  time_t time_stop = time(NULL);
  double user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  double sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;
  psio.reset();

  if (iter==maxiter){
     throw PsiException("  CEPA iterations did not converge.",__FILE__,__LINE__);
  }

  SCS_CEPA();

  outfile->Printf("\n");
  outfile->Printf("  %s iterations converged!\n",cepa_type);
  outfile->Printf("\n");

  if (cepa_level < 1) {
      outfile->Printf("  %s variational energy: %20.12lf\n",cepa_type,evar);
      outfile->Printf("  %s transition energy:  %20.12lf\n",cepa_type,eccsd);
      outfile->Printf("\n");
      eccsd = evar;
  }

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

      outfile->Printf("        OS MP2 FNO correction:             %20.12lf\n",delta_emp2_os);
      outfile->Printf("        SS MP2 FNO correction:             %20.12lf\n",delta_emp2_ss);
      outfile->Printf("        MP2 FNO correction:                %20.12lf\n",delta_emp2);
      outfile->Printf("\n");
  }

  outfile->Printf("        OS SCS-MP2 correlation energy:     %20.12lf\n",emp2_os*emp2_os_fac);
  outfile->Printf("        SS SCS-MP2 correlation energy:     %20.12lf\n",emp2_ss*emp2_ss_fac);
  outfile->Printf("        SCS-MP2 correlation energy:        %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac);
  outfile->Printf("      * SCS-MP2 total energy:              %20.12lf\n",emp2_os*emp2_os_fac+emp2_ss*emp2_ss_fac+escf);
  outfile->Printf("\n");
  outfile->Printf("        OS MP2 correlation energy:         %20.12lf\n",emp2_os);
  outfile->Printf("        SS MP2 correlation energy:         %20.12lf\n",emp2_ss);
  outfile->Printf("        MP2 correlation energy:            %20.12lf\n",emp2);
  outfile->Printf("      * MP2 total energy:                  %20.12lf\n",emp2+escf);
  outfile->Printf("\n");
  if (cepa_level>=0){
     if (options_.get_bool("SCS_CEPA")){
        outfile->Printf("        OS SCS-%s correlation energy: %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac);
        outfile->Printf("        SS SCS-%s correlation energy: %20.12lf\n",cepa_type,eccsd_ss*eccsd_ss_fac);
        outfile->Printf("        SCS-%s correlation energy:    %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
        outfile->Printf("      * SCS-%s total energy:          %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
        outfile->Printf("\n");
     }
     outfile->Printf("        OS %s correlation energy:     %20.12lf\n",cepa_type,eccsd_os);
     outfile->Printf("        SS %s correlation energy:     %20.12lf\n",cepa_type,eccsd_ss);
     outfile->Printf("        %s correlation energy:        %20.12lf\n",cepa_type,eccsd);
     outfile->Printf("      * %s total energy:              %20.12lf\n",cepa_type,eccsd+escf);
  }else{
     if (options_.get_bool("SCS_CEPA")){
        outfile->Printf("        OS SCS-%s correlation energy:    %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac);
        outfile->Printf("        SS SCS-%s correlation energy:    %20.12lf\n",cepa_type,eccsd_ss*eccsd_ss_fac);
        outfile->Printf("        SCS-%s correlation energy:       %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac);
        outfile->Printf("      * SCS-%s total energy:             %20.12lf\n",cepa_type,eccsd_os*eccsd_os_fac+eccsd_ss*eccsd_ss_fac+escf);
        outfile->Printf("\n");
     }
     outfile->Printf("        OS %s correlation energy:        %20.12lf\n",cepa_type,eccsd_os);
     outfile->Printf("        SS %s correlation energy:        %20.12lf\n",cepa_type,eccsd_ss);
     outfile->Printf("        %s correlation energy:           %20.12lf\n",cepa_type,eccsd);
     outfile->Printf("      * %s total energy:                 %20.12lf\n",cepa_type,eccsd+escf);
  }
  outfile->Printf("\n");
  outfile->Printf("  Total time for %s iterations: %10.2lf s (user)\n",cepa_type,user_stop-user_start);
  outfile->Printf("                                  %10.2lf s (system)\n",sys_stop-sys_start);
  outfile->Printf("                                  %10d s (total)\n",(int)time_stop-(int)time_start);
  outfile->Printf("\n");
  outfile->Printf("  Time per iteration:             %10.2lf s (user)\n",(user_stop-user_start)/(iter-1));
  outfile->Printf("                                  %10.2lf s (system)\n",(sys_stop-sys_start)/(iter-1));
  outfile->Printf("                                  %10.2lf s (total)\n",((double)time_stop-(double)time_start)/(iter-1));


  free(pair_energy);
  return Success;
}

void CoupledPair::PairEnergy(){

  if (cepa_level<1) return;

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

  std::shared_ptr<PSIO> psio(new PSIO());
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

  std::shared_ptr<PSIO> psio(new PSIO());
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
  C_DCOPY(o*v,w1,1,tempv+o*o*v*v,1);
  C_DAXPY(o*v,-1.0,t1,1,tempv+o*o*v*v,1);
  C_DCOPY(o*v,w1,1,t1,1);
}

void CoupledPair::SCS_CEPA(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

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
double CoupledPair::VariationalEnergy(){

    long int v = nvirt;
    long int o = ndoccact;
    long int rs = nmo;

    // (ai|bj)
    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_IAJB,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&integrals[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_IAJB,1);

    if (t2_on_disk){
       psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
       psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempt[0],o*o*v*v*sizeof(double));
       psio->close(PSIF_DCC_T2,1);
       tb = tempt;
    }

    // read residual (no, it is in tempv still)
    //psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    //psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
    //psio->close(PSIF_DCC_R2,1);

    double fac = 1.0;
    if (cepa_level == 0)       fac = 0.0;
    else if (cepa_level == -1) fac = 1.0;
    else if (cepa_level == -2) fac = 1.0/o;
    else if (cepa_level == -3) fac = 1.0-(2.0*o-2.0)*(2.0*o-3.0) / (2.0*o*(2.0*o-1.0));

    double evar = 0.0;
    double nrm = 1.0;
    for (int i=0; i<o; i++){
        for (int j=0; j<o; j++){
            for (int a=o; a<rs; a++){
                for (int b=o; b<rs; b++){
                    long int ijab = (a-o)*o*o*v+(b-o)*o*o+i*o+j;
                    long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                    double dum = 0.5*(tb[ijab]-tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                    nrm += (2.0*dum*dum + tb[ijab]*tb[ijab]) * fac;
                }
            }
        }
    }
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            nrm += 2.0 * t1[a*o+i]*t1[a*o+i] * fac;
        }
    }
    double ed = 0.0;
    double e0 = 0.0;
    for (int i = 0; i < o; i++) {
        double di = - eps[i];
        for (int j = 0; j < o; j++) {
            double dij = di - eps[j];
            for (int a = o; a < rs; a++) {
                double dija = dij + eps[a];
                for (int b = o; b < rs; b++) {
                    double dijab = dija + eps[b];
                    long int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                    long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;
                    ed += (tempv[ijab] + dijab * tb[ijab]) * (2.0 * tb[ijab] - tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                    e0 += integrals[iajb] * (2.0 * tb[ijab] - tb[(b-o)*o*o*v+(a-o)*o*o+i*o+j]);
                }
            }
        }
    }
    evar = ed + 2.0 * e0;
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            double dia = -eps[i]+eps[a+o];
            evar += 2.0 * (dia*t1[a*o+i] + w1[a*o+i])*t1[a*o+i];
        }
    }
    return evar/nrm;
}
double CoupledPair::CheckEnergy(){

  long int v = nvirt;
  long int o = ndoccact;
  long int rs = nmo;

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

}


}} // end of namespaces
