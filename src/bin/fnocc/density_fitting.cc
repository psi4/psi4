#include<lib3index/dftensor.h>
#include<lib3index/cholesky.h>
#include<sys/times.h>
#include<psifiles.h>
#include<../bin/fnocc/ccsd.h>

using namespace psi;

namespace psi{ namespace fnocc{

void DFCoupledCluster::ThreeIndexIntegrals(){

  struct tms total_tmstime;
  const long clk_tck = sysconf(_SC_CLK_TCK);

  long int i;
  double time_start,user_start,sys_start,time_stop,user_stop,sys_stop;

  times(&total_tmstime);
  time_start = time(NULL);
  user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_start  = ((double) total_tmstime.tms_stime)/clk_tck;

  fprintf(outfile,"\n");
  fprintf(outfile,"  Generate 3-index integrals:\n");
  fprintf(outfile,"\n");

  if (!ischolesky_) {
      // 3-index integrals from auxiliary basis
      int nocc = reference_wavefunction_->doccpi()[0];
      int nvir = reference_wavefunction_->nmopi()[0]-reference_wavefunction_->doccpi()[0];
      int aocc = nocc-reference_wavefunction_->frzcpi()[0];
      int avir = nvir-reference_wavefunction_->frzvpi()[0];

      boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
      boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, reference_wavefunction_->molecule(), "DF_BASIS_CC");
      boost::shared_ptr<DFTensor> DF (new DFTensor(reference_wavefunction_->basisset(),auxiliary,reference_wavefunction_->Ca(),nocc,nvir,aocc,avir,options_));
      SharedMatrix tmpQso = DF->Qso();
      nQ = auxiliary->nbf();
      double ** Qp = tmpQso->pointer();

      boost::shared_ptr<PSIO> psio(new PSIO());
      psio->open(PSIF_DCC_QSO,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_QSO,"qso",(char*)&Qp[0][0],nQ*nso*nso*sizeof(double));
      psio->close(PSIF_DCC_QSO,1);
  }else {
      // Cholesky 3-index integrals
      boost::shared_ptr<BasisSet> primary = reference_wavefunction_->basisset();
      boost::shared_ptr<IntegralFactory> integral (new IntegralFactory(primary,primary,primary,primary));
      double tol = options_.get_double("CHOLESKY_TOLERANCE");
      boost::shared_ptr<CholeskyERI> Ch (new CholeskyERI(boost::shared_ptr<TwoBodyAOInt>(integral->eri()),0.0,tol,memory));
      Ch->choleskify();
      int nch  = Ch->Q();
      boost::shared_ptr<Matrix> L = Ch->L();
      double ** Lp = L->pointer();
      nQ = nch;

      boost::shared_ptr<PSIO> psio(new PSIO());
      psio->open(PSIF_DCC_QSO,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_QSO,"qso",(char*)&Lp[0][0],nQ*nso*nso*sizeof(double));
      psio->close(PSIF_DCC_QSO,1);
  }

  times(&total_tmstime);
  time_stop = time(NULL);
  user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_stop  = ((double) total_tmstime.tms_stime)/clk_tck;

  fprintf(outfile,"\n");
  fprintf(outfile,"  Time to generate 3-index integrals: %6.2lf s (user)\n",user_stop-user_start);
  fprintf(outfile,"                                      %6.2lf s (system)\n",sys_stop-sys_start);
  fprintf(outfile,"                                      %6d s (total)\n",(int)time_stop-(int)time_start);
  fprintf(outfile,"\n");


}

}} // end of namespaces
