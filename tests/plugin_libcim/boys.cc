#include"psi4-dec.h"
#include<psifiles.h>
#include"blas.h"
#include"boys.h"
#include<libmints/mints.h>
#include<libmints/mintshelper.h>
#include<libmints/wavefunction.h>
#include<libmints/matrix.h>
#define au_debye 0.393430307

using namespace psi;
using namespace boost;

namespace psi{

/*
 * Boys localization
 * 
 * This code generates local molecular orbitals that satisfy the Boys 
 * criterion for localization: maximize ( sum <i|r|i>^2 ).
 */
Boys::Boys(Options&options){
  fprintf(outfile,"  ==> Boys Localization <==\n");
  fprintf(outfile,"\n");
  options_ = options;
  Localize();
}
Boys::~Boys()
{}

void Boys::Localize(){
  /*
   * grab reference and parameters
   */
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.reference_wavefunction();
  if (ref.get() == NULL){
     throw PsiException("no wavefunction?",__FILE__,__LINE__);
  }
  nirreps = ref->nirrep();
  if (nirreps>1){
     throw PsiException("boys localization requires symmetry c1",__FILE__,__LINE__);
  }
  int*sorbs = ref->nsopi();
  int*orbs  = ref->nmopi();
  int*docc  = ref->doccpi();
  int*fzc   = ref->frzcpi();
  int*fzv   = ref->frzvpi();
  nso = nmo = ndocc = nvirt = nfzc = nfzv = 0;
  for (long int h=0; h<nirreps; h++){
      nfzc   += fzc[h];
      nfzv   += fzv[h];
      nso    += sorbs[h];
      nmo    += orbs[h];
      ndocc  += docc[h];
  }
  ndoccact = ndocc - nfzc;
  nvirt  = nmo - ndocc;

  /*
   * lmo Fock matrix
   */
  Fock = SharedMatrix (new Matrix("LMO Fock Matix",nmo,nmo));
  double**Fock_pointer = Fock->pointer();

  /*
   * mo -> lmo transformation matrix
   */
  Clmo = SharedMatrix (new Matrix("Localized Molecular Orbitals",nso,nmo));
  Clmo_pointer = Clmo->pointer();

  /*
   * grab dipole integrals
   */
  boost::shared_ptr<MintsHelper> mints(new MintsHelper());
  std::vector< boost::shared_ptr<Matrix> > dipole = mints->so_dipole();
  mu = (double***)malloc(3*sizeof(double**));
  mu[0] = dipole[0]->pointer();
  mu[1] = dipole[1]->pointer();
  mu[2] = dipole[2]->pointer();

  /*
   * transform dipole integrals to mo basis
   */
  boost::shared_ptr<Matrix> Ca = ref->Ca();
  double**Ca_pointer = Ca->pointer();
  SoToMo(Ca->rowspi()[0],Ca->colspi()[0],mu[0],Ca_pointer,Clmo_pointer);
  SoToMo(Ca->rowspi()[0],Ca->colspi()[0],mu[1],Ca_pointer,Clmo_pointer);
  SoToMo(Ca->rowspi()[0],Ca->colspi()[0],mu[2],Ca_pointer,Clmo_pointer);

  memset((void*)&Clmo_pointer[0][0],'\0',nmo*nmo*sizeof(double));
  for (int i=0; i<nmo; i++) Clmo_pointer[i][i] = 1.0;

  /*
   * array to shuffle orbitals between iterations
   * always seed srand the same way so we can reproduce results
   */
  reorder = (int*)malloc(nmo*sizeof(int));
  srand(0);

  double fac = 1.0/au_debye;
  fac *= fac;
  double angle,sum = LocalizationSum();
  fprintf(outfile,"\n");
  fprintf(outfile,"  Initial localization sum: %20.12lf\n",sum*fac);
  fprintf(outfile,"\n");

  /*
   * begin localization procedure
   */
  fprintf(outfile,"  ==> Boys Iterations <==\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"  Iter    Localization sum    Convergence\n");
  double convergence = options_.get_double("BOYS_CONVERGENCE");
  int maxiter        = options_.get_int("BOYS_MAXITER");
  int iter = 1;
  conv=1.0;
  while(conv>convergence && iter<maxiter){
      F_DCOPY(nmo*nmo,&Clmo_pointer[0][0],1,&Fock_pointer[0][0],1);
      ShuffleOrbitals();
      for (int i=nfzc; i<ndocc; i++){
          int ii = reorder[i-nfzc];
          for (int j=i+1; j<ndocc; j++){
              int jj = reorder[j-nfzc];
              GetAngle(ii,jj,angle);
              Rotate(ii,jj,angle);
          }
      }
      /*
       * check convergence
       */
      F_DAXPY(nmo*nmo,-1.0,&Clmo_pointer[0][0],1,&Fock_pointer[0][0],1);
      conv = F_DNRM2(nmo*nmo,&Fock_pointer[0][0],1);
      sum  = LocalizationSum();
      fprintf(outfile,"  %4i        %12.4lf   %12.7lf\n",iter++,sum*fac,conv);
      fflush(outfile);
  }
  fprintf(outfile,"\n");
  fprintf(outfile,"  Final localization sum:   %20.12lf\n",sum*fac);
  fprintf(outfile,"\n");
  fflush(outfile);

  /*
   * build Fock matrix in LMO basis
   */
  memset((void*)&Fock_pointer[0][0],'\0',nmo*nmo*sizeof(double));
  boost::shared_ptr<Vector> eps = ref->epsilon_a();
  F_DCOPY(nmo,eps->pointer(),1,&Fock_pointer[0][0],nmo+1);
  F_DGEMM('t','n',nmo,nmo,nmo,1.0,&Clmo_pointer[0][0],nmo,&Fock_pointer[0][0],nmo,0.0,&mu[0][0][0],nmo);
  F_DGEMM('t','t',nmo,nmo,nmo,1.0,&Clmo_pointer[0][0],nmo,&mu[0][0][0],nmo,0.0,&Fock_pointer[0][0],nmo);
  
  /*
   * build so->lmo transformation matrix.
   */
  F_DGEMM('t','n',nmo,nso,nmo,1.0,&Clmo_pointer[0][0],nmo,&Ca_pointer[0][0],nmo,0.0,&mu[0][0][0],nmo);
  F_DCOPY(nso*nmo,&mu[0][0][0],1,&Clmo_pointer[0][0],1);

  //free(mu);
  free(reorder);
}

/*
 * 2x2 rotation.  update dipole integrals and accumulate mo->lmo transformation.
 */
void Boys::Rotate(int i,int j,double angle){
  double c = cos(angle);
  double s = sin(angle);
  // update dipole integrals
  for (int xyz=0; xyz<3; xyz++){
      for (int k=0; k<nmo; k++){
          double muIk =   c * mu[xyz][i][k] + s * mu[xyz][j][k];
          double muJk = - s * mu[xyz][i][k] + c * mu[xyz][j][k];
          mu[xyz][i][k] = muIk;
          mu[xyz][j][k] = muJk;
      }
      for (int k=0; k<nmo; k++){
          double mukI =   c * mu[xyz][k][i] + s * mu[xyz][k][j];
          double mukJ = - s * mu[xyz][k][i] + c * mu[xyz][k][j];
          mu[xyz][k][i] = mukI;
          mu[xyz][k][j] = mukJ;
      }
  }
  // update mo -> lmo transformation matrix
  for (int k=0; k<nmo; k++){
      double CIk =   c * Clmo_pointer[i][k] + s * Clmo_pointer[j][k];
      double CJk = - s * Clmo_pointer[i][k] + c * Clmo_pointer[j][k];
      Clmo_pointer[i][k] = CIk;
      Clmo_pointer[j][k] = CJk;
  }
}

/*
 * define 2x2 rotation to maximize sum of squares of orbital centroids
 */
void Boys::GetAngle(int i,int j,double&angle){
  double a12   = 0.0;
  double b12   = 0.0;
  for (int xyz=0; xyz<3; xyz++){
      double diff = mu[xyz][i][i]-mu[xyz][j][j];
      double muij = mu[xyz][i][j];
      a12 += muij * muij - 0.25 * diff * diff;
      b12 += muij * diff;
  }
  if (fabs(a12)<1e-10 && fabs(b12)<1e-10){
     angle = 0.0;
     return;
  }
  double val1 = -a12/sqrt(a12*a12+b12*b12);
  //double val2 = b12/sqrt(a12*a12+b12*b12);

  // alpha should be <= pi/2
  angle = 0.25 * acos(val1);
  if (b12<0.0) angle *= -1;
  //angle = 0.25*atan2(b12,-a12);
}
double Boys::LocalizationSum(){
  double sum = 0.0;
  for (int xyz=0; xyz<3; xyz++){
      sum += F_DDOT(ndoccact,&mu[xyz][nfzc][nfzc],nmo+1,&mu[xyz][nfzc][nfzc],nmo+1);
  }
  return sum;
}

/*
 * transform a matrix from the so to mo basis
 */
void Boys::SoToMo(int nso,int nmo,double**mat,double**trans,double**temp){
  F_DGEMM('n','n',nmo,nso,nso,1.0,&trans[0][0],nso,&mat[0][0],nso,0.0,&temp[0][0],nmo);
  F_DGEMM('n','t',nmo,nmo,nso,1.0,&temp[0][0],nmo,&trans[0][0],nmo,0.0,&mat[0][0],nmo);
}

/*
 * shuffle orbitals to keep from getting stuck
 */
void Boys::ShuffleOrbitals(){
  int dim = ndoccact;
  for (int i=0; i<dim; i++) reorder[i] = i+nfzc;
  for (int i=0; i<dim; i++){
      int ind1 = (int)(1.*dim*rand()/RAND_MAX);
      int ind2 = (int)(1.*dim*rand()/RAND_MAX);
      if (ind1==dim) ind1--;
      if (ind2==dim) ind2--;
      while(ind1==ind2){
        ind2 = (int)(1.*dim*rand()/RAND_MAX);
        if (ind2==dim) ind2--;
      }
      int temp = reorder[ind1];
      reorder[ind1] = reorder[ind2];
      reorder[ind2] = temp;
  }
}

}
