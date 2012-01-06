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

  long int o = ccsd->ndoccact;
  long int v = ccsd->nvirt;
  double *F  = ccsd->eps;

  boost::shared_ptr<PSIO> psio(new PSIO());

  // allocate memory for a couple of buffers
  long int memory = Process::environment.get_memory();
  //memory *= 1024L*1024L;
  memory -= 8L*(o*o*v*v+o*v);
  // how many tiles for the ov^2 transformation?
  long int ntiles,tilesize,lasttile;
  ntiles = 1L;
  tilesize = o/ntiles;
  if (ntiles*tilesize<o) tilesize++;
  while(2*tilesize*v*v*v*8L>memory){
     ntiles++;
     tilesize = o/ntiles;
     if (ntiles*tilesize<o) tilesize++;
  }
  lasttile = o - (ntiles-1L)*tilesize;
  long int*tilesizes=(long int*)malloc(ntiles*sizeof(long int));
  for (int i=0; i<ntiles-1; i++) tilesizes[i] = tilesize;
  tilesizes[ntiles-1] = lasttile;
  
  long int dim = o*o*v*v;
  if (tilesize*v*v*v>o*o*v*v) dim = tilesize*v*v*v;
  double*amps1 = (double*)malloc(dim*sizeof(double));
  double*amps2 = (double*)malloc(dim*sizeof(double));

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
  fprintf(outfile,"        Number of virtual orbitals in original space:  %5li\n",v);
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

  // transform (oo|ov) integrals
  psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJAK,"E2ijak",(char*)&amps1[0],o*o*o*v*sizeof(double));
  F_DGEMM('t','n',nvirt_no,o*o*o,v,1.0,Dab,v,amps1,v,0.0,amps2,nvirt_no);
  psio->write_entry(PSIF_IJAK,"E2ijak",(char*)&amps2[0],o*o*o*nvirt_no*sizeof(double));
  psio->close(PSIF_IJAK,1);

  // transform (ov|ov) integrals: (ia|jb)
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&amps2[0],o*o*v*v*sizeof(double));
  // transform b of (ia|jb)
  F_DGEMM('t','n',nvirt_no,o*o*v,v,1.0,Dab,v,amps2,v,0.0,amps1,nvirt_no);
  // sort to (iajb) -> (ijba) so we can transform a
  for (int i=0; i<o; i++){
      for (int a=0; a<v; a++){
          for (int j=0; j<o; j++){
              for (int b=0; b<nvirt_no; b++){
                  amps2[i*o*v*nvirt_no+j*v*nvirt_no+b*v+a] = 
                  amps1[i*o*v*nvirt_no+a*o*nvirt_no+j*nvirt_no+b];
              }
          }
      }
  }
  // transform b
  F_DGEMM('t','n',nvirt_no,o*o*nvirt_no,v,1.0,Dab,v,amps2,v,0.0,amps1,nvirt_no);
  // resort back to (iajb)
  for (int i=0; i<o; i++){
      for (int a=0; a<nvirt_no; a++){
          for (int j=0; j<o; j++){
              for (int b=0; b<nvirt_no; b++){
                  amps2[i*o*nvirt_no*nvirt_no+a*o*nvirt_no+j*nvirt_no+b] = 
                  amps1[i*o*nvirt_no*nvirt_no+j*nvirt_no*nvirt_no+b*nvirt_no+a];
              }
          }
      }
  }
  psio->write_entry(PSIF_KLCD,"E2klcd",(char*)&amps2[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);


  // transform (ov|vv) integrals:
  psio_address addrread,addrwrite;
  psio->open(PSIF_ABCI,PSIO_OPEN_OLD);
  addrread  = PSIO_ZERO;
  addrwrite = PSIO_ZERO;

  for (int n=0; n<ntiles; n++){
      psio->read(PSIF_ABCI,"E2abci",(char*)&amps2[0],tilesizes[n]*v*v*v*sizeof(double),addrread,&addrread);
      // transform iabc -> iabC
      F_DGEMM('t','n',nvirt_no,tilesizes[n]*v*v,v,1.0,Dab,v,amps2,v,0.0,amps1,nvirt_no);
      // sort iabC -> aiCb
      for (int i=0; i<tilesizes[n]; i++){
          for (int a=0; a<v; a++){
              for (int b=0; b<v; b++){
                  for (int c=0; c<nvirt_no; c++){
                      amps2[a*tilesizes[n]*nvirt_no*v+i*v*nvirt_no+c*v+b] = 
                      amps1[i*v*v*nvirt_no+a*v*nvirt_no+b*nvirt_no+c];
                  }
              }
          }
      }
      // transform aiCb -> aiCB
      F_DGEMM('t','n',nvirt_no,tilesizes[n]*v*nvirt_no,v,1.0,Dab,v,amps2,v,0.0,amps1,nvirt_no);
      // transform aiCB -> AiCB
      F_DGEMM('n','n',tilesizes[n]*nvirt_no*nvirt_no,nvirt_no,v,1.0,amps1,tilesizes[n]*nvirt_no*nvirt_no,Dab,v,0.0,amps2,tilesizes[n]*nvirt_no*nvirt_no);
      // sort AiCB -> iABC
      for (int i=0; i<tilesizes[n]; i++){
          for (int a=0; a<nvirt_no; a++){
              for (int b=0; b<nvirt_no; b++){
                  for (int c=0; c<nvirt_no; c++){
                      amps1[i*nvirt_no*nvirt_no*nvirt_no+a*nvirt_no*nvirt_no+b*nvirt_no+c] =
                      amps2[a*tilesizes[n]*nvirt_no*nvirt_no+i*nvirt_no*nvirt_no+c*nvirt_no+b];
                  }
              }
          }
      }
      psio->write(PSIF_ABCI,"E2abci",(char*)&amps1[0],tilesizes[n]*nvirt_no*nvirt_no*nvirt_no*sizeof(double),addrwrite,&addrwrite);
  }
  psio->close(PSIF_ABCI,1);

  // transform (oo|ov), (ov|vv), and (ov|ov) integrals with libtrans
  //Transformation(ccsd,Dab);
  // sort integrals required for triples
  //OutOfCoreSortTriples(ccsd->nfzc,ccsd->nfzv,ccsd->nmotemp,ccsd->ndoccact,ccsd->nvirt_no);

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

  // compute scaling factor for (T) from MP2 energy in MO and NO spaces
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd",(char*)&amps2[0],o*o*nvirt_no*nvirt_no*sizeof(double));
  psio->close(PSIF_KLCD,1);
  double emp2=0.0;
  for (int a=o; a<o+nvirt_no; a++){
      double da = ccsd->eps[a];
      for (int b=o; b<o+nvirt_no; b++){
          double dab = da + ccsd->eps[b];
          for (int i=0; i<o; i++){
              double dabi = dab - ccsd->eps[i];
              for (int j=0; j<o; j++){
                  int iajb = i*nvirt_no*nvirt_no*o+(a-o)*nvirt_no*o+j*nvirt_no+(b-o);
                  int jaib = j*nvirt_no*nvirt_no*o+(a-o)*nvirt_no*o+i*nvirt_no+(b-o);
                  double dijab = dabi-ccsd->eps[j];
                  double amp = - amps2[iajb]/dijab;
                  emp2 += amp * (2.0*amps2[iajb]-amps2[jaib]);
              }
          }
      }
  }
  ccsd->scale_t = ccsd->emp2/emp2;
  fprintf(outfile,"        MP2 correlation energy in original space:   %20.12lf\n",ccsd->emp2);
  fprintf(outfile,"        MP2 correlation energy in truncated space:  %20.12lf\n",emp2);
  fprintf(outfile,"\n");
  double eccsd=0.0;
  ijab = 0;
  for (int a=o; a<o+nvirt_no; a++){
      for (int b=o; b<o+nvirt_no; b++){
          for (int i=0; i<o; i++){
              for (int j=0; j<o; j++){
                  int iajb = i*nvirt_no*nvirt_no*o+(a-o)*nvirt_no*o+j*nvirt_no+(b-o);
                  int jaib = j*nvirt_no*nvirt_no*o+(a-o)*nvirt_no*o+i*nvirt_no+(b-o);
                  eccsd += (ccsd->tb[ijab++] + ccsd->t1[(a-o)*o+i]*ccsd->t1[(b-o)*o+j]) 
                        * (2.0*amps2[iajb]-amps2[jaib]);
              }
          }
      }
  }
  ccsd->scale_t = ccsd->eccsd/eccsd;
  fprintf(outfile,"        CCSD correlation energy in original space:  %20.12lf\n",ccsd->eccsd);
  fprintf(outfile,"        CCSD correlation energy in truncated space: %20.12lf\n",eccsd);
  fprintf(outfile,"\n");

  // free memory
  free(tilesizes);
  free(amps1);
  free(amps2);
  free(neweps);
  free(newFock);
  free(temp);
  free(eigvalDab);
  free(Dab);

  return Success;
}

void Transformation(boost::shared_ptr<psi::CoupledCluster>ccsd,double*Ca_motono){
  /*shared_ptr<PSIO> psio(_default_psio_lib_);
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
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);*/
}

} // end of namespace
