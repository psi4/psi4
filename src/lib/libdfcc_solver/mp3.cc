#include "mp3.h"
#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libqt/qt.h>
#include <psiconfig.h>
#include <time.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libpsio/aiohandler.h>

using namespace boost;
using namespace psi;

namespace psi { namespace dfcc {

MP3::MP3(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  psio_->open(DFCC_INT_FILE,PSIO_OPEN_NEW);
  df_integrals();
  mo_integrals();
}

MP3::~MP3()
{
  psio_->close(DFCC_INT_FILE,1);
}

double MP3::compute_energy()
{
  tIAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);
  t2IAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);
  vIAJB_ = init_array(naocc_*naocc_*navir_*(ULI) navir_);

  double emp2;

  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(tIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));
  apply_denom(tIAJB_);
  emp2 = energy(tIAJB_);

  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  MP2 Correlation Energy      %18.10lf\n",emp2);
  fprintf(outfile,"  Total DF-MP2 Energy         %18.10lf\n\n",Eref_+emp2);

  term_1(); // t_ab^ij <- t_ac^ik [2(jb|kc)-(jk|bc)]

  iajb_ibja(tIAJB_);

  term_2(); // t_ab^ij <- - t_ac^ki (jb|kc)

  iajb_ibja(t2IAJB_);

  term_3(); // t_ab^ij <- - t_bc^ki (jk|ac)

  iajb_ibja(tIAJB_);
  iajb_ibja(t2IAJB_);

  iajb_ijab(tIAJB_);
  iajb_ijab(t2IAJB_);

  term_4(); // t_ab^ij <- t_ab^kl (ik|jl)
  term_5(); // t_ab^ij <- t_cd^ij (ac|bd)

  ijab_iajb(tIAJB_);
  ijab_iajb(t2IAJB_);

  double emp3;
  symmetrize();
  apply_denom(t2IAJB_);
  emp3 = energy(t2IAJB_);

  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  MP3 Correlation Energy      %18.10lf\n",emp3);
  fprintf(outfile,"  Total DF-MP3 Energy         %18.10lf\n\n",Eref_+emp2+emp3);

  free(tIAJB_);
  free(t2IAJB_);
  free(vIAJB_);

  return(Eref_+emp2+emp3);
}

void MP3::print_header()
{
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                       DF-MP3                         *\n");
    fprintf(outfile, "\t*    Third-Order Moller-Plesset Perturbation Theory    *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\n");
    CC::print_header();

}

double MP3::energy(double *tIAJB)
{
  double energy;

  psio_->read_entry(DFCC_INT_FILE,"G OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  energy = C_DDOT(naocc_*navir_*naocc_*(long int) navir_,tIAJB,1,vIAJB_,1);

  return(energy);
}

void MP3::apply_denom(double *tIAJB)
{
  for (int i=0,ia=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++,ia++) {
    for (int j=0,jb=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++,jb++) {
      long int iajb = ia*naocc_*(long int) navir_ + (long int) jb;
      tIAJB[iajb] /= evals_aoccp_[i] + evals_aoccp_[j] - evals_avirp_[a]
        - evals_avirp_[b];
    }}
  }}
}

void MP3::symmetrize()
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

void MP3::term_1()
{
  psio_->read_entry(DFCC_INT_FILE,"G OVVO Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,0.0,t2IAJB_,naocc_*navir_);
}

void MP3::term_2()
{
  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,-2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void MP3::term_3()
{
  psio_->read_entry(DFCC_INT_FILE,"OOVV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,-2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void MP3::term_4()
{
  double *vIJKL = init_array(naocc_*naocc_*naocc_*(long int) naocc_);

  psio_->read_entry(DFCC_INT_FILE,"OOOO Integrals",(char *)
    &(vIJKL[0]),naocc_*naocc_*naocc_*naocc_*sizeof(double));

  C_DGEMM('N','N',naocc_*naocc_,navir_*navir_,naocc_*naocc_,1.0,vIJKL,
    naocc_*naocc_,tIAJB_,navir_*navir_,1.0,t2IAJB_,navir_*navir_);

  free(vIJKL);
}

void MP3::term_5()
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
#if 0
PSMP3::PSMP3(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  psio_->open(DFCC_INT_FILE,PSIO_OPEN_NEW);
  df_integrals();
  mo_integrals();
  dfints_ = shared_ptr<DFTensor>(new DFTensor(psio_, basisset_, ribasis_,
    true));
  dfints_->form_OV_integrals((ULI)(0.9*(double)doubles_), C_aocc_, C_avir_,
    false, fitting_algorithm_, fitting_condition_, schwarz_cutoff_);
}

PSMP3::~PSMP3()
{
  psio_->close(DFCC_INT_FILE,1);
}

void PSMP3::print_header()
{
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                       PS-MP3                         *\n");
    fprintf(outfile, "\t*    Third-Order Moller-Plesset Perturbation Theory    *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\n");
    CC::print_header();
}

double PSMP3::compute_energy()
{
  mem_ = doubles_ - (long int) 3*naocc_*naocc_*navir_*navir_;
  if (mem_ < 0)
    throw PsiException("Not enough memory", __FILE__,
      __LINE__);

  tIAJB_ = init_array((ULI) naocc_*naocc_*navir_*navir_);
  t2IAJB_ = init_array((ULI) naocc_*naocc_*navir_*navir_);
  vIAJB_ = init_array((ULI) naocc_*naocc_*navir_*navir_);

  double emp2;

  psio_->read_entry(DFCC_INT_FILE,"OVOV Integrals",(char *)
    &(tIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));
  apply_denom(tIAJB_);
  emp2 = energy(tIAJB_);

  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  MP2 Correlation Energy      %18.10lf\n",emp2);
  fprintf(outfile,"  Total DF-MP2 Energy         %18.10lf\n\n",Eref_+emp2);

  term_1(); // t_ab^ij <- 2 t_ac^ik (jb|kc)
  term_2(); // t_ab^ij <- - t_ac^ik (jk|bc)

  iajb_ibja(tIAJB_);

  term_3(); // t_ab^ij <- - t_ac^ki (jb|kc)

  iajb_ibja(t2IAJB_);

  term_4(); // t_ab^ij <- - t_bc^ki (jk|ac)

  iajb_ibja(tIAJB_);
  iajb_ibja(t2IAJB_);

  iajb_ijab(tIAJB_);
  iajb_ijab(t2IAJB_);

  term_5(); // t_ab^ij <- t_ab^kl (ik|jl)
  term_6(); // t_ab^ij <- t_cd^ij (ac|bd)

  ijab_iajb(tIAJB_);
  ijab_iajb(t2IAJB_);

  double emp3;
  symmetrize();
  apply_denom(t2IAJB_);
  emp3 = energy(t2IAJB_);

  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  MP3 Correlation Energy      %18.10lf\n",emp3);
  fprintf(outfile,"  Total DF-MP3 Energy         %18.10lf\n\n",Eref_+emp2+emp3);

  free(tIAJB_);
  free(t2IAJB_);
  free(vIAJB_);

  return(Eref_+emp2+emp3);
}

double PSMP3::energy(double *tIAJB)
{
  double energy;

  psio_->read_entry(DFCC_INT_FILE,"G OVOV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  energy = C_DDOT(naocc_*navir_*naocc_*(long int) navir_,tIAJB,1,vIAJB_,1);

  return(energy);
}

void PSMP3::apply_denom(double *tIAJB)
{
  for (int i=0,ia=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++,ia++) {
    for (int j=0,jb=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++,jb++) {
      long int iajb = ia*naocc_*(long int) navir_ + (long int) jb;
      tIAJB[iajb] /= evals_aoccp_[i] + evals_aoccp_[j] - evals_avirp_[a]
        - evals_avirp_[b];
    }}
  }}
}

void PSMP3::symmetrize()
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

void PSMP3::term_1()
{
  int naux = dfints_->naux();

  if ((long int) 2*naocc_*navir_*naux > (mem_)) {
    throw PsiException("Not enough memory", __FILE__,
      __LINE__);
  }

  shared_ptr<TensorChunk> BpIA;

  BpIA = dfints_->get_ov_iterator(mem_-(long int) naocc_*navir_*naux);

  shared_ptr<Matrix> BpIA_chunk = BpIA->chunk();

  double **B_p_OV = BpIA_chunk->pointer();
  double *T_p_OV = init_array(naocc_*navir_*naux);

  int nblock = BpIA->nblock();

  if (nblock != 1) {
    throw PsiException("Not enough memory", __FILE__,
      __LINE__);
  }

  BpIA->read_block(0);

  C_DGEMM('N','N',naocc_*navir_,naux,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,B_p_OV[0],naux,0.0,T_p_OV,naux);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naux,4.0,T_p_OV,naux,
    B_p_OV[0],naux,0.0,t2IAJB_,naocc_*navir_);

  free(T_p_OV);
}

void PSMP3::term_2()
{
  psio_->read_entry(DFCC_INT_FILE,"OOVV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  iajb_ibja(vIAJB_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,-2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void PSMP3::term_3()
{
  int naux = dfints_->naux();

  if ((long int) 2*naocc_*navir_*naux > (mem_)) {
    throw PsiException("Not enough memory", __FILE__,
      __LINE__);
  }

  shared_ptr<TensorChunk> BpIA;

  BpIA = dfints_->get_ov_iterator(mem_-(long int) naocc_*navir_*naux);

  shared_ptr<Matrix> BpIA_chunk = BpIA->chunk();

  double **B_p_OV = BpIA_chunk->pointer();
  double *S_p_OV = init_array(naocc_*navir_*naux);

  int nblock = BpIA->nblock();

  if (nblock != 1) {
    throw PsiException("Not enough memory", __FILE__,
      __LINE__);
  }

  BpIA->read_block(0);

  C_DGEMM('N','N',naocc_*navir_,naux,naocc_*navir_,1.0,
    tIAJB_,naocc_*navir_,B_p_OV[0],naux,0.0,S_p_OV,naux);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naux,-2.0,S_p_OV,naux,
    B_p_OV[0],naux,1.0,t2IAJB_,naocc_*navir_);

  free(S_p_OV);
}

void PSMP3::term_4()
{
  psio_->read_entry(DFCC_INT_FILE,"OOVV Integrals",(char *)
    &(vIAJB_[0]),naocc_*navir_*naocc_*navir_*sizeof(double));

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,naocc_*navir_,-2.0,
    tIAJB_,naocc_*navir_,vIAJB_,naocc_*navir_,1.0,t2IAJB_,naocc_*navir_);
}

void PSMP3::term_5()
{
  double *vIJKL = init_array(naocc_*naocc_*naocc_*(long int) naocc_);

  psio_->read_entry(DFCC_INT_FILE,"OOOO Integrals",(char *)
    &(vIJKL[0]),naocc_*naocc_*naocc_*naocc_*sizeof(double));

  C_DGEMM('N','N',naocc_*naocc_,navir_*navir_,naocc_*naocc_,1.0,vIJKL,
    naocc_*naocc_,tIAJB_,navir_*navir_,1.0,t2IAJB_,navir_*navir_);

  free(vIJKL);
}

void PSMP3::term_6()
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
#endif
}}
