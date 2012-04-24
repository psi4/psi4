#include<lib3index/dftensor.h>
#include"ccsd.h"
#include"blas.h"

using namespace psi;

namespace psi{

void CoupledCluster::ThreeIndexIntegrals(){

  fprintf(outfile,"\n");
  fprintf(outfile,"  Transform 3-index integrals:\n");
  fprintf(outfile,"\n");

  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.reference_wavefunction();
  int nocc = ref->doccpi()[0];
  int nvir = ref->nmopi()[0]-ref->doccpi()[0];
  int aocc = nocc-ref->frzcpi()[0];
  int avir = nvir-ref->frzvpi()[0];

  boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
  boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, ref->molecule(), "DF_BASIS_CC");

  boost::shared_ptr<DFTensor> DF (new DFTensor(ref->basisset(),auxiliary,ref->Ca(),nocc,nvir,aocc,avir,Process::environment.options));
  SharedMatrix tmpQoo = DF->Qoo();
  SharedMatrix tmpQov = DF->Qov();
  SharedMatrix tmpQvv = DF->Qvv();

  double**Qoo_pointer = tmpQoo->pointer();
  double**Qov_pointer = tmpQov->pointer();
  double**Qvv_pointer = tmpQvv->pointer();

  nQ = auxiliary->nbf();
  int o = aocc;
  int v = avir;

  // copy 3-index tensors into ccsd
  Qoo = (double*)malloc(o*o*nQ*sizeof(double));
  Qov = (double*)malloc(o*v*nQ*sizeof(double));
  Qvv = (double*)malloc(v*v*nQ*sizeof(double));

  F_DCOPY(nQ*o*o,&Qoo_pointer[0][0],1,Qoo,1);
  F_DCOPY(nQ*o*v,&Qov_pointer[0][0],1,Qov,1);

  for (long int q=0; q<nQ; q++){
      F_DCOPY(v*v,&Qvv_pointer[q][0],1,Qvv+q,nQ);
  }

  boost::shared_ptr<PSIO> psio(new PSIO());
  psio->open(PSIF_QVVT,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_QVVT,"qvv_transpose",(char*)&Qvv[0],nQ*v*v*sizeof(double));
  psio->close(PSIF_QVVT,1);

  F_DCOPY(nQ*v*v,&Qvv_pointer[0][0],1,Qvv,1);
  psio->open(PSIF_QVV,PSIO_OPEN_NEW);
  psio->write_entry(PSIF_QVV,"qvv",(char*)&Qvv[0],nQ*v*v*sizeof(double));
  psio->close(PSIF_QVV,1);
}

} // end of namespace psi
