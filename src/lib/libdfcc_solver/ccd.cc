#include "ccd.h"

#include <time.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libpsio/aiohandler.h>
#include <libmints/mints.h>

using namespace boost;
using namespace psi;

namespace psi { namespace dfcc {

CCD::CCD(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  psio_->open(DFCC_INT_FILE,PSIO_OPEN_NEW);
  df_integrals();
  mo_integrals();
}

CCD::~CCD()
{
  psio_->close(DFCC_INT_FILE,1);
}

double CCD::compute_energy()
{
  time_t start = time(NULL);
  time_t stop;

  if (options_.get_bool("DIIS"))
    diis_ = shared_ptr<DFCCDIIS>(new DFCCDIIS(DFCC_DIIS_FILE,naocc_*naocc_*
      navir_*navir_,options_.get_int("MAX_DIIS_VECS"),psio_));

  tIAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);
  t2IAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);
  vIAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);
  xIAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);

  double emp2;

  term_1();
  apply_denom();
  emp2 = energy();

  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",emp2);
  fprintf(outfile,"  Total DF-MP2 Energy         %18.10lf\n\n",Eref_+emp2);
  fprintf(outfile,"  Iter       Energy (H)          dE (H)             RMS (H)     Time (s)\n");
  fflush(outfile);

  int iter = 1;
  int done = 0;
  double e_new;
  double e_old = emp2;
  double rms;

  do {

    term_1(); // t_ab^ij <- (ia|jb)
    term_2(); // t_ab^ij <- t_ac^ik [2(jb|kc)-(jk|bc)]
    term_3(); // t_ab^ij <- 2 t_ac^ik [2(kc|ld)-(kd|lc)] t_bd^jl
    term_4(); // t_ab^ij <- - t_ac^ij [2(kc|ld)-(kd|lc)] t_bd^kl
    term_5(); // t_ab^ij <- - t_ab^ik [2(kc|ld)-(kd|lc)] t_cd^jl

    iajb_ibja(tIAJB_);

    term_6(); // t_ab^ij <- - t_ac^ik [2(kc|ld)-(kd|lc)] t_bd^lj
    term_7(); // t_ab^ij <- - t_ac^ki (jb|kc)
              // t_ab^ij <- t_ac^ki (kc|ld) t_bd^lj

    iajb_ibja(t2IAJB_);

    term_8(); // t_ab^ij <- t_cb^il (kc|ld) t_da^jk
    term_9(); // t_ab^ij <- - t_bc^ki (jk|ac)

    iajb_ibja(tIAJB_);
    iajb_ibja(t2IAJB_);

    iajb_ijab(tIAJB_);
    iajb_ijab(t2IAJB_);

    term_10(); // t_ab^ij <- t_ab^kl (ik|jl)
    term_11(); // t_ab^ij <- t_cd^ij (ac|bd)
    term_12(); // t_ab^ij <- t_cd^ij (kc|ld) t_ab^kl

    ijab_iajb(tIAJB_);
    ijab_iajb(t2IAJB_);

    symmetrize(); // t_ab^ij <- 0.5*( t_ab^ij + t_ba^ji )
    apply_denom();
    e_new = energy();
    rms = store_error_vecs();

    stop = time(NULL);

    fprintf(outfile,"  %4d %16.8lf %17.9lf %17.9lf %12ld",iter,e_new,
      e_old-e_new,rms,stop-start);
    fflush(outfile);

    if (options_.get_int("MIN_DIIS_VECS") <= iter &&
        options_.get_bool("DIIS")) {
      diis_->get_new_vector(tIAJB_,xIAJB_);
      fprintf(outfile,"  DIIS\n");
    }
    else {
      fprintf(outfile,"\n");
    }
    fflush(outfile);

    iter++;

    if (iter > options_.get_int("MAXITER")) done = 1;
    if (fabs(e_old-e_new) < pow(10.0,-(double) options_.get_int("E_CONVERGE")))
      done = 1;
    if (rms < pow(10.0,-(double) options_.get_int("T_CONVERGE"))) done = 1;

    e_old = e_new;
  }
  while(!done);

  fprintf(outfile,"\n");
  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",e_old);
  fprintf(outfile,"  Total DF-CCD Energy         %18.10lf\n\n",Eref_+e_old);

  free(tIAJB_);
  free(t2IAJB_);
  free(vIAJB_);
  free(xIAJB_);

  return(0.0);
}

void CCD::print_header()
{
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                       DF-CCD                         *\n");
    fprintf(outfile, "\t*               Coupled-Cluster Doubles                *\n");
    fprintf(outfile, "\t*                with all sorts of shit                *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\n");
    CC::print_header();

}

double CCD::energy()
{
  double energy;

  psio_->read_entry(DFCC_INT_FILE,"G OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  energy = C_DDOT(naocc_*navir_*naocc_*(long int) navir_,t2IAJB_,1,vIAJB_,1);

  C_DCOPY(naocc_*navir_*naocc_*(long int) navir_,tIAJB_,1,xIAJB_,1);
  C_DCOPY(naocc_*navir_*naocc_*(long int) navir_,t2IAJB_,1,tIAJB_,1);

  return(energy);
}

double CCD::store_error_vecs()
{
  double rms;

  C_DAXPY(naocc_*navir_*naocc_*(long int) navir_,-1.0,t2IAJB_,1,xIAJB_,1);

  if (options_.get_bool("DIIS"))
    diis_->store_vectors(t2IAJB_,xIAJB_);

  rms = C_DDOT(naocc_*navir_*naocc_*(long int) navir_,xIAJB_,1,xIAJB_,1);
  rms /= naocc_*navir_*naocc_*(double) navir_;

  return(sqrt(rms));
}

void CCD::apply_denom()
{
  for (int i=0,ia=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++,ia++) {
    for (int j=0,jb=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++,jb++) {
      long int iajb = ia*naocc_*(long int) navir_ + (long int) jb;
      t2IAJB_[iajb] /= evals_aoccp_[i] + evals_aoccp_[j] - evals_avirp_[a]
        - evals_avirp_[b];
    }}
  }}
}

void CCD::symmetrize()
{
  for (int ia=0; ia<naocc_*navir_; ia++) {
    for (int jb=0; jb<ia; jb++) {
      long int iajb = ia*naocc_*(long int) navir_ + (long int) jb;
      long int jbia = jb*naocc_*(long int) navir_ + (long int) ia;
      double tval = t2IAJB_[iajb] + t2IAJB_[jbia];
      t2IAJB_[iajb] = 0.5*tval;
      t2IAJB_[jbia] = 0.5*tval;
  }}
}

void CCD::term_1()
{
  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(t2IAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));
}

void CCD::term_2()
{
  psio_->read_entry(DFCC_INT_FILE,"G OVVO Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_3()
{
  psio_->read_entry(DFCC_INT_FILE,"G OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,0.0,xIAJB_,naocc_*navir_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,2.0,
    xIAJB_,naocc_*navir_,tIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_4()
{
  double **xAB = block_matrix(navir_,navir_);

  C_DGEMM('T','N',navir_,navir_,naocc_*navir_*naocc_,1.0,vIAJB_,navir_,
    tIAJB_,navir_,0.0,xAB[0],navir_);

  C_DGEMM('N','N',naocc_*navir_*naocc_,navir_,navir_,-2.0,tIAJB_,navir_,
    xAB[0],navir_,1.0,t2IAJB_,navir_);

  free_block(xAB);
}

void CCD::term_5()
{
  double **xIJ = block_matrix(naocc_,naocc_);

  C_DGEMM('N','T',naocc_,naocc_,navir_*naocc_*navir_,1.0,vIAJB_,
    navir_*naocc_*navir_,tIAJB_,navir_*naocc_*navir_,0.0,xIJ[0],naocc_);

  C_DGEMM('T','N',naocc_,navir_*naocc_*navir_,naocc_,-2.0,xIJ[0],naocc_,
    tIAJB_,navir_*naocc_*navir_,1.0,t2IAJB_,navir_*naocc_*navir_);

  free_block(xIJ);
}

void CCD::term_6()
{
  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,-2.0,
    xIAJB_,naocc_*navir_,tIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}
/* Non-DF Factorization ... Remember to dick with the read in term_8
void CCD::term_7()
{
  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,0.0,xIAJB_,naocc_*navir_);

  C_DAXPY(naocc_*navir_*naocc_*(long int) navir_,-2.0,xIAJB_,1,t2IAJB_,1);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    xIAJB_,naocc_*navir_,tIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}
*/

void CCD::term_7()
{
  double *B_p_OV = init_array(naocc_*navir_*ndf_);
  double *T_p_OV = init_array(naocc_*navir_*ndf_);

  psio_->read_entry(DFCC_INT_FILE,"OV DF Integrals",(char *)
      &(B_p_OV[0]),naocc_*navir_*ndf_*sizeof(double));

  C_DGEMM('N','N',naocc_*navir_,ndf_,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,B_p_OV,ndf_,0.0,T_p_OV,ndf_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,-2.0,T_p_OV,ndf_,
    B_p_OV,ndf_,1.0,t2IAJB_,naocc_*navir_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,T_p_OV,ndf_,
    T_p_OV,ndf_,1.0,t2IAJB_,naocc_*navir_);

  free(B_p_OV);
  free(T_p_OV);
}

void CCD::term_8()
{
  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  iajb_ibja(vIAJB_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,0.0,xIAJB_,naocc_*navir_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,1.0,
    xIAJB_,naocc_*navir_,tIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_9()
{
  psio_->read_entry(DFCC_INT_FILE,"OOVV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,-2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void CCD::term_10()
{
  double *vIJKL = init_array(naocc_*naocc_*naocc_*(long int) naocc_);

  psio_->read_entry(DFCC_INT_FILE,"OOOO Integrals",(char *)
    &(vIJKL[0]),naocc_*naocc_*naocc_*naocc_*sizeof(double));

  C_DGEMM('N','N',naocc_*naocc_,navir_*navir_,naocc_*naocc_,1.0,vIJKL,
    naocc_*naocc_,tIAJB_,navir_*navir_,1.0,t2IAJB_,navir_*navir_);

  free(vIJKL);
}

/* For testing only
void CCD::term_11()
{
  double *vABCD = init_array(navir_*navir_*navir_*(long int) navir_);

  psio_->read_entry(DFCC_INT_FILE,"VVVV Integrals",(char *)
    &(vABCD[0]),navir_*navir_*navir_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*naocc_,navir_*navir_,navir_*navir_,1.0,tIAJB_,
    navir_*navir_,vABCD,navir_*navir_,1.0,t2IAJB_,navir_*navir_);

  free(vABCD);
}
*/
void CCD::term_11()
{
  int occtri = naocc_*(naocc_+1)/2;
  int virtri = navir_*(navir_+1)/2;
  int svirtri = navir_*(navir_-1)/2;

  double **tpIJAB = block_matrix(occtri,virtri);
  double **tmIJAB = block_matrix(occtri,svirtri);

  for(int i=0; i<naocc_; i++) {
  for(int j=0; j<=i; j++) {
    for(int a=0; a<navir_; a++) {
    for(int b=0; b<=a; b++) {
      int ij = INDEX(i,j);
      int ab = INDEX(a-1,b);
      int cd = INDEX(a,b);
      int ijab = i*naocc_*navir_*navir_ + j*navir_*navir_ + a*navir_ + b;
      int jiab = j*naocc_*navir_*navir_ + i*navir_*navir_ + a*navir_ + b;

      if (a != b) {
        tpIJAB[ij][cd] = 0.5*tIAJB_[ijab];
        tpIJAB[ij][cd] += 0.5*tIAJB_[jiab];
        tmIJAB[ij][ab] = 0.5*tIAJB_[ijab];
        tmIJAB[ij][ab] -= 0.5*tIAJB_[jiab];
      }
      else {
        tpIJAB[ij][cd] = 0.25*tIAJB_[ijab];
        tpIJAB[ij][cd] += 0.25*tIAJB_[jiab];
      }
  }}}}

  int blocksize;
  int loopsize;

  if (navir_ % 2 == 0) {
    blocksize = navir_+1;
    loopsize = virtri/blocksize;
  }
  else {
    blocksize = navir_;
    loopsize = virtri/blocksize;
  }

  double **sIJAB = block_matrix(occtri,virtri);
  double **vVVVVp[2];
  vVVVVp[0] = block_matrix(blocksize,virtri);
  vVVVVp[1] = block_matrix(blocksize,virtri);

  psio_address next_VVVVp = PSIO_ZERO;

  shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

  psio_->read(DFCC_INT_FILE,"VVVV+ Integrals",(char *) &(vVVVVp[0][0][0]),
        blocksize*virtri*(ULI) sizeof(double),next_VVVVp,&next_VVVVp);

  for(int a_read=0; a_read<loopsize; a_read++) {
    if (a_read < loopsize-1)
      aio->read(DFCC_INT_FILE,"VVVV+ Integrals",(char *)
        &(vVVVVp[(a_read+1)%2][0][0]),blocksize*virtri*sizeof(double),
        next_VVVVp,&next_VVVVp);

    C_DGEMM('N','T',occtri,blocksize,virtri,1.0,tpIJAB[0],virtri,
      vVVVVp[a_read%2][0],virtri,1.0,&(sIJAB[0][a_read*blocksize]),virtri);

    if (a_read < loopsize-1)
      aio->synchronize();
  }

  free_block(vVVVVp[0]);
  free_block(vVVVVp[1]);
  free_block(tpIJAB);

  if (navir_ % 2 == 0) {
    blocksize = navir_-1;
    loopsize = svirtri/blocksize;
  }
  else {
    blocksize = navir_;
    loopsize = svirtri/blocksize;
  }

  double **aIJAB = block_matrix(occtri,svirtri);
  double **vVVVVm[2];
  vVVVVm[0] = block_matrix(blocksize,svirtri);
  vVVVVm[1] = block_matrix(blocksize,svirtri);

  psio_address next_VVVVm = PSIO_ZERO;

  psio_->read(DFCC_INT_FILE,"VVVV- Integrals",(char *) &(vVVVVm[0][0][0]),
        blocksize*svirtri*(ULI) sizeof(double),next_VVVVm,&next_VVVVm);

  for(int a_read=0; a_read<loopsize; a_read++) {
    if (a_read < loopsize-1)
      aio->read(DFCC_INT_FILE,"VVVV- Integrals",(char *)
        &(vVVVVm[(a_read+1)%2][0][0]),
        blocksize*svirtri*sizeof(double),next_VVVVm,&next_VVVVm);

    C_DGEMM('N','T',occtri,blocksize,svirtri,1.0,tmIJAB[0],svirtri,
      vVVVVm[a_read%2][0],svirtri,1.0,&(aIJAB[0][a_read*blocksize]),svirtri);

    if (a_read < loopsize-1)
      aio->synchronize();
  }

  free_block(vVVVVm[0]);
  free_block(vVVVVm[1]);
  free_block(tmIJAB);

  for(int i=0,ij=0; i<naocc_; i++) {
  for(int j=0; j<naocc_; j++,ij++) {
    int kl = INDEX(i,j);
    for(int a=0,ab=0; a<navir_; a++) {
    for(int b=0; b<navir_; b++,ab++) {
      int cd = INDEX(a,b);
      int ijab = i*naocc_*navir_*navir_ + j*navir_*navir_ + a*navir_ + b;
      vIAJB_[ijab] = sIJAB[kl][cd];
  }}}}

  for(int i=0; i<naocc_; i++) {
  for(int j=0; j<i; j++) {
    int ij = i*naocc_ + j;
    int ji = j*naocc_ + i;
    int kl = INDEX(i,j);
    for(int a=0; a<navir_; a++) {
    for(int b=0; b<a; b++) {
      int ab = a*navir_ + b;
      int ba = b*navir_ + a;
      int cd = INDEX(a-1,b);
      int ijab = i*naocc_*navir_*navir_ + j*navir_*navir_ + a*navir_ + b;
      int jiba = j*naocc_*navir_*navir_ + i*navir_*navir_ + b*navir_ + a;
      int ijba = i*naocc_*navir_*navir_ + j*navir_*navir_ + b*navir_ + a;
      int jiab = j*naocc_*navir_*navir_ + i*navir_*navir_ + a*navir_ + b;
      vIAJB_[ijab] += aIJAB[kl][cd];
      vIAJB_[jiba] += aIJAB[kl][cd];
      vIAJB_[ijba] -= aIJAB[kl][cd];
      vIAJB_[jiab] -= aIJAB[kl][cd];
  }}}}

  for(int i=0; i<naocc_; i++) {
    int ii = i*naocc_ + i;
    int kk = INDEX(i,i);
    for(int a=0; a<navir_; a++) {
    for(int b=0; b<a; b++) {
      int ab = a*navir_ + b;
      int ba = b*navir_ + a;
      int cd = INDEX(a-1,b);
      int iiab = ii*navir_*navir_ + ab;
      int iiba = ii*navir_*navir_ + ba;
      vIAJB_[iiab] += aIJAB[kk][cd];
      vIAJB_[iiba] -= aIJAB[kk][cd];
  }}}

  free_block(sIJAB);
  free_block(aIJAB);

  C_DAXPY(naocc_*navir_*naocc_*(long int) navir_,1.0,vIAJB_,1,t2IAJB_,1);
}

void CCD::term_12()
{
  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  iajb_ijab(vIAJB_);

  double *vIJKL = init_array(naocc_*naocc_*naocc_*(long int) naocc_);

  C_DGEMM('N','T',naocc_*naocc_,naocc_*naocc_,navir_*navir_,1.0,tIAJB_,
    navir_*navir_,vIAJB_,navir_*navir_,0.0,vIJKL,naocc_*naocc_);

  C_DGEMM('N','N',naocc_*naocc_,navir_*navir_,naocc_*naocc_,1.0,vIJKL,
    naocc_*naocc_,tIAJB_,navir_*navir_,1.0,t2IAJB_,navir_*navir_);

  free(vIJKL);

}

}}
