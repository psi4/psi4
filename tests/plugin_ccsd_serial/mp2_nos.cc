#include"psi4-dec.h"
#include"ccsd.h"
#include <libplugin/plugin.h>
#include<boost/shared_ptr.hpp>
#include<liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include"blas.h"
#include"sort.h"

using namespace psi;
using namespace boost;

namespace psi{
PsiReturnType MP2NaturalOrbitals(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options);
void Transformation(boost::shared_ptr<psi::CoupledCluster>ccsd,double*Ca_motono);

PsiReturnType MP2NaturalOrbitals(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options){
  fprintf(outfile,"        *******************************************************\n");
  fprintf(outfile,"        *                                                     *\n");
  fprintf(outfile,"        *          MP2 Natural Orbitals for CCSD(T)           *\n");
  fprintf(outfile,"        *                                                     *\n");
  fprintf(outfile,"        *******************************************************\n");
  fprintf(outfile,"\n");
  fflush(outfile);

  int o = ccsd->ndoccact;
  int v = ccsd->nvirt;
  double *F  = ccsd->eps;

  boost::shared_ptr<PSIO> psio(new PSIO());

  // point to some memory we've already allocated
  double*amps1 = ccsd->tempt;
  double*amps2 = ccsd->tempv;

  // build mp2 amplitudes for mp2 density
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&amps2[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);
  int ijab = 0;
  for (int a=o; a<o+v; a++){
      double da = F[a];
      for (int b=o; b<o+v; b++){
          double dab = da + F[b];
          for (int i=0; i<o; i++){
              double dabi = dab - F[i];
              for (int j=0; j<o; j++){
                  int iajb = i*v*v*o+(a-o)*v*o+j*v+(b-o);
                  double dijab = dabi-F[j];
                  amps1[ijab++] = - amps2[iajb]/dijab;
              }
          }
      }
  }
  ijab = 0;
  for (int a=o; a<o+v; a++){
      for (int b=o; b<o+v; b++){
          for (int i=0; i<o; i++){
              for (int j=0; j<o; j++){
                  int ijba = (b-o)*o*o*v+(a-o)*o*o+i*o+j;
                  amps2[ijab] = 2.0*amps1[ijab++] - amps1[ijba];
              }
          }
      }
  }

  // build ab block of the density:
  double*Dab = (double*)malloc(v*v*sizeof(double));
  F_DGEMM('t','n',v,v,v*o*o,2.0,amps1,v*o*o,amps2,v*o*o,0.0,Dab,v);

  // diagonalize Dab
  double*eigvalDab=(double*)malloc(v*sizeof(double));
  Diagonalize(v,Dab,eigvalDab);
  // reorder transformation matrix:
  double*temp    = (double*)malloc(v*v*sizeof(double));
  for (int i=0; i<v; i++){
      F_DCOPY(v,Dab+(v-1-i)*v,1,temp+i*v,1);
  }

  // establish cutoff for frozen virtuals
  double cutoff = options.get_double("VIRTUAL_CUTOFF");
  int nvirt_no = 0;
  for (int i=0; i<v; i++) if (eigvalDab[i]>cutoff) nvirt_no++;
  ccsd->nvirt_no = nvirt_no;

  fprintf(outfile,"        Cutoff for significant NO occupancy: %5.3le\n",cutoff);
  fprintf(outfile,"\n");
  fprintf(outfile,"        Number of virtual orbitals in original space:  %5i\n",v);
  fprintf(outfile,"        Number of virtual orbitals in truncated space: %5i\n",nvirt_no);
  fprintf(outfile,"\n");

  // transform Fock matrix to MP2 NO basis
  double*newFock = (double*)malloc(v*v*sizeof(double));
  memset((void*)newFock,'\0',v*v*sizeof(double));
  F_DCOPY(v,ccsd->eps+o,1,newFock,v+1);
  F_DGEMM('n','n',v,nvirt_no,v,1.0,newFock,v,temp,v,0.0,Dab,v);
  F_DGEMM('t','n',nvirt_no,nvirt_no,v,1.0,temp,v,Dab,v,0.0,newFock,nvirt_no);

  // diagonalize new Fock matrix for semi-canonical orbitals
  double*neweps = (double*)malloc(nvirt_no*sizeof(double));
  Diagonalize(nvirt_no,newFock,neweps);
  // construct full mo -> no transformation matrix
  F_DGEMM('n','n',v,nvirt_no,nvirt_no,1.0,temp,v,newFock,nvirt_no,0.0,Dab,v);

  // transform (oo|ov), (ov|vv), and (ov|ov) integrals with libtrans
  Transformation(ccsd,Dab);

  // new orbital energies in truncated space
  F_DCOPY(nvirt_no,neweps,1,ccsd->eps+o,1);

  // transform the t1 amplitudes
  F_DGEMM('n','n',o,nvirt_no,v,1.0,ccsd->t1,o,Dab,v,0.0,amps1,o);
  F_DCOPY(o*v,amps1,1,ccsd->t1,1);

  // transform t2(a_mo,b_mo,i,j) -> t2(a_no,b_mo,i,j)
  F_DGEMM('n','n',o*o*v,nvirt_no,v,1.0,ccsd->tb,o*o*v,Dab,v,0.0,amps1,o*o*v);
  // sort t2(a_no,b_mo,i,j) -> t2(b_mo,a_no,j,i)
  // if we go ahead and swap i & j, we won't need to sort again after this
  for (int a=0; a<nvirt_no; a++){
      for (int b=0; b<v; b++){
          for (int i=0; i<o; i++){
              for (int j=0; j<o; j++){
                  amps2[b*nvirt_no*o*o+a*o*o+j*o+i] = amps1[a*v*o*o+b*o*o+i*o+j];
              }
          }
      }
  }
  // transform t2(b_mo,a_no,j,i) -> t2(b_no,a_mo,j,i)
  F_DGEMM('n','n',o*o*nvirt_no,nvirt_no,v,1.0,amps2,o*o*nvirt_no,Dab,v,0.0,ccsd->tb,o*o*nvirt_no);

  // sort integrals required for triples
  OutOfCoreSortTriples(ccsd->nfzc,ccsd->nfzv,ccsd->nmotemp,ccsd->ndoccact,ccsd->nvirt_no);

  // free memory
  free(neweps);
  free(newFock);
  free(temp);
  free(eigvalDab);
  free(Dab);

  return Success;
}
void Transformation(boost::shared_ptr<psi::CoupledCluster>ccsd,double*Ca_motono){
  shared_ptr<PSIO> psio(_default_psio_lib_);
  std::vector<shared_ptr<MOSpace> > spaces;
  shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  spaces.push_back(MOSpace::all);

  boost::shared_ptr<Wavefunction> ref =
       Process::environment.reference_wavefunction();

  boost::shared_ptr<Matrix> Cc = ref->Ca_subset("SO","FROZEN_OCC");
  boost::shared_ptr<Matrix> Ci = ref->Ca_subset("SO","ACTIVE_OCC");
  boost::shared_ptr<Matrix> Ca = ref->Ca_subset("SO","ACTIVE_VIR");
  boost::shared_ptr<Matrix> Cv = ref->Ca_subset("SO","FROZEN_VIR");

  // create shared matrix to hold so->no transformation vectors
  int*rowspi = new int[1];
  int*colspi = new int[1];
  int v = ccsd->nvirt;
  rowspi[0] = ccsd->ndocc+ccsd->nvirt;
  colspi[0] = ccsd->nvirt_no;
  SharedMatrix Ca_sotono = SharedMatrix (new Matrix("SO to NO transformation",1,rowspi,colspi));
  delete[] rowspi;
  delete[] colspi;

  // fill shared matrix with mo->no transformation vectors
  double**Ca_sotono_pointer = Ca_sotono->pointer(0);
  double**Ca_pointer        = Ca->pointer(0);

  // construct so->no transformation matrix
  F_DGEMM('t','n',Ca_sotono->colspi()[0],Ca_sotono->rowspi()[0],v,1.0,Ca_motono,v,&Ca_pointer[0][0],v,0.0,&Ca_sotono_pointer[0][0],Ca_sotono->colspi()[0]);

  // so->no integral transformation
  IntegralTransform ints(Cc,Ci,Ca_sotono,Cv, spaces, IntegralTransform::Restricted,
           IntegralTransform::IWLOnly, IntegralTransform::QTOrder, IntegralTransform::OccAndVir, false);
  dpd_set_default(ints.get_dpd_id());
  int print = 4;
  ints.set_print(print);
  ints.set_keep_dpd_so_ints(1);
  ints.set_keep_iwl_so_ints(1);
  ints.initialize();
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
}

} // end of namespace
