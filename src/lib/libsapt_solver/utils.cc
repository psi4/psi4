#include "sapt.h"
#include "sapt0.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

void SAPT::zero_disk(int file, const char *array, int rows, int columns)
{
  double *zero = init_array(columns);
  psio_address next_PSIF = PSIO_ZERO;

  for (int i=0; i<rows; i++) {
    psio_->write(file,array,(char *) &(zero[0]),sizeof(double)*columns,
      next_PSIF,&next_PSIF);
  }
  free(zero);
}

void SAPT0::read_all(SAPTDFInts *ints)
{
  int nri = ndf_;
  if (ints->dress_) nri += 3;

  ints->B_p_ = block_matrix(nri,ints->ij_length_);

  psio_->read_entry(ints->filenum_,ints->label_,(char *) 
    &(ints->B_p_[0][0]),sizeof(double)*ndf_*ints->ij_length_);

  C_DCOPY(3*ints->ij_length_,&(ints->B_d_[0][0]),1,
    &(ints->B_p_[ndf_][0]),1);
}

void SAPT0::read_block(Iterator *iter, SAPTDFInts *intA)
{
  bool last_block = false;
  if (iter->curr_block == iter->num_blocks) last_block = true;
  bool dress = false;
  if (intA->dress_) dress = true;
  int block_length = iter->block_size[iter->curr_block-1];
  iter->curr_block++;

//  printf("%d %d %d %d %d %d\n",iter->num_blocks,block_length,
//    iter->curr_block,block_length,last_block,dress);

  iter->curr_size = block_length;
  if (last_block && dress) block_length -= 3;

  psio_->read(intA->filenum_,intA->label_,(char *) &(intA->B_p_[0][0]),
    sizeof(double)*block_length*intA->ij_length_,intA->next_DF_,
    &intA->next_DF_);

  if (dress && last_block) {
    if (intA->dress_) {
      C_DCOPY(3*intA->ij_length_,&(intA->B_d_[0][0]),1,
        &(intA->B_p_[block_length][0]),1);
    }
    else {
      memset(&(intA->B_p_[block_length][0]),'\0',sizeof(double)*3*
        intA->ij_length_);
    }
  }
}

void SAPT0::read_block(Iterator *iter, SAPTDFInts *intA, SAPTDFInts *intB)
{
  bool last_block = false;
  if (iter->curr_block == iter->num_blocks) last_block = true;
  bool dress = false;
  if (intA->dress_ || intB->dress_) dress = true;
  int block_length = iter->block_size[iter->curr_block-1];
  iter->curr_block++;

//  printf("%d %d %d %d %d %d\n",iter->num_blocks,block_length,
//    iter->curr_block,block_length,last_block,dress);

  iter->curr_size = block_length;
  if (last_block && dress) block_length -= 3;

  psio_->read(intA->filenum_,intA->label_,(char *) &(intA->B_p_[0][0]),
    sizeof(double)*block_length*intA->ij_length_,intA->next_DF_,
    &intA->next_DF_);

  psio_->read(intB->filenum_,intB->label_,(char *) &(intB->B_p_[0][0]),
    sizeof(double)*block_length*intB->ij_length_,intB->next_DF_,
    &intB->next_DF_);

  if (dress && last_block) {
    if (intA->dress_) {
      C_DCOPY(3*intA->ij_length_,&(intA->B_d_[0][0]),1,
        &(intA->B_p_[block_length][0]),1);
    }
    else {
      memset(&(intA->B_p_[block_length][0]),'\0',sizeof(double)*3*
        intA->ij_length_);
    }
  
    if (intB->dress_) {
      C_DCOPY(3*intB->ij_length_,&(intB->B_d_[0][0]),1,
        &(intB->B_p_[block_length][0]),1);
    }
    else {
      memset(&(intB->B_p_[block_length][0]),'\0',sizeof(double)*3*
        intB->ij_length_);
    }
  }
}

Iterator SAPT0::get_iterator(long int mem, SAPTDFInts *intA)
{
  long int ij_size = intA->ij_length_;
  int max_length = ndf_;
  if (intA->dress_) max_length += 3;
  if (ij_size > mem)
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  int length = mem/ij_size;
  if (length > max_length) length = max_length;

  return(set_iterator(100,intA));
//return(set_iterator(length,intA));
}

Iterator SAPT0::set_iterator(int length, SAPTDFInts *intA)
{
  if (0 >= length)
    throw PsiException("Not enough memory", __FILE__,__LINE__);

  int max_length = ndf_;
  if (intA->dress_) max_length += 3;

  if (length > max_length)
    length = max_length;

  int num = max_length/length;
  int gimp = max_length%length;

  Iterator iter;
  iter.num_blocks = num;
  if (gimp > 3) iter.num_blocks++;
  iter.curr_block = 1;
  iter.block_size = init_int_array(iter.num_blocks);
  iter.curr_size = 0;

  for (int i=0; i<num; i++) iter.block_size[i] = length;
  if (gimp > 3) {
    iter.block_size[num] = gimp;
  }
  else if (gimp) {
    for (int i=0; i<gimp; i++) iter.block_size[i%num]++;
  }

  int max_block = iter.block_size[0];

  intA->B_p_ = block_matrix(max_block,intA->ij_length_);

  return(iter);
}

Iterator SAPT0::get_iterator(long int mem, SAPTDFInts *intA, SAPTDFInts *intB)
{
  int ij_size = intA->ij_length_ + intB->ij_length_;
  int max_length = ndf_;
  if (intA->dress_ || intB->dress_) max_length += 3;
  if (ij_size > mem)
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  int length = mem/ij_size;
  if (length > max_length) length = max_length;

  return(set_iterator(100,intA,intB));
//return(set_iterator(length,intA,intB));
}

Iterator SAPT0::set_iterator(int length, SAPTDFInts *intA, SAPTDFInts *intB)
{
  if (0 >= length)
    throw PsiException("Not enough memory", __FILE__,__LINE__);

  int max_length = ndf_;
  if (intA->dress_ || intB->dress_) max_length += 3;

  if (length > max_length) 
    length = max_length;
  
  int num = max_length/length;
  int gimp = max_length%length;
    
  Iterator iter;
  iter.num_blocks = num;
  if (gimp > 3) iter.num_blocks++;
  iter.curr_block = 1;
  iter.block_size = init_int_array(iter.num_blocks);
  iter.curr_size = 0;

  for (int i=0; i<num; i++) iter.block_size[i] = length;
  if (gimp > 3) {
    iter.block_size[num] = gimp;
  }
  else if (gimp) {
    for (int i=0; i<gimp; i++) iter.block_size[i%num]++;
  }

  int max_block = iter.block_size[0];

  intA->B_p_ = block_matrix(max_block,intA->ij_length_); 
  intB->B_p_ = block_matrix(max_block,intB->ij_length_); 

  return(iter);
}

SAPTDFInts SAPT0::set_A_AA()
{
  double enuc, NA, NB;

  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);

  SAPTDFInts A_AA;

  A_AA.dress_ = true;
  A_AA.active_ = false;

  A_AA.i_length_ = noccA_;
  A_AA.j_length_ = noccA_;
  A_AA.ij_length_ = noccA_*noccA_;
  A_AA.i_start_ = 0;
  A_AA.j_start_ = 0;

  A_AA.B_d_ = block_matrix(3,noccA_*noccA_);

  A_AA.filenum_ = PSIF_SAPT_AA_DF_INTS;
  A_AA.label_ = "AA RI Integrals";

  for (int a=0; a<noccA_; a++){
    int aa = a*noccA_+a;
    A_AA.B_d_[0][aa] = 1.0;
    A_AA.B_d_[2][aa] = enuc;
    for (int ap=0; ap<noccA_; ap++){
      int aap = a*noccA_+ap;
      A_AA.B_d_[1][aap] = NB*vBAA_[a][ap];
    }
  }

  return(A_AA);
}

SAPTDFInts SAPT0::set_B_BB()
{
  double enuc, NA, NB;

  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);

  SAPTDFInts B_BB;

  B_BB.dress_ = true;
  B_BB.active_ = false;

  B_BB.i_length_ = noccB_;
  B_BB.j_length_ = noccB_;
  B_BB.ij_length_ = noccB_*noccB_; 
  B_BB.i_start_ = 0;
  B_BB.j_start_ = 0;
    
  B_BB.B_d_ = block_matrix(3,noccB_*noccB_);

  B_BB.filenum_ = PSIF_SAPT_BB_DF_INTS;
  B_BB.label_ = "BB RI Integrals";
    
  for (int b=0; b<noccB_; b++){
    int bb = b*noccB_+b;
    B_BB.B_d_[1][bb] = 1.0;
    B_BB.B_d_[2][bb] = enuc;
    for (int bp=0; bp<noccB_; bp++){
      int bbp = b*noccB_+bp;
      B_BB.B_d_[0][bbp] = NA*vABB_[b][bp];
    }
  }

  return(B_BB);
}

SAPTDFInts SAPT0::set_A_AR()
{ 
  double enuc, NA, NB;
  
  NA = 1.0 / NA_;
  NB = 1.0 / NB_; 
  enuc = sqrt(enuc_*NA*NB);
  
  SAPTDFInts A_AR;
  
  A_AR.dress_ = true;
  A_AR.active_ = false;
  
  A_AR.i_length_ = noccA_;
  A_AR.j_length_ = nvirA_;
  A_AR.ij_length_ = noccA_*nvirA_;
  A_AR.i_start_ = 0;
  A_AR.j_start_ = 0;

  A_AR.B_d_ = block_matrix(3,noccA_*nvirA_);
  
  A_AR.filenum_ = PSIF_SAPT_AA_DF_INTS;
  A_AR.label_ = "AR RI Integrals";

  for (int a=0; a<noccA_; a++){
    for (int r=0; r<nvirA_; r++){
      int ar = a*nvirA_+r;
      A_AR.B_d_[1][ar] = NB*vBAA_[a][r+noccA_];
    }
  }

  return(A_AR);
}

SAPTDFInts SAPT0::set_B_BS()
{
  double enuc, NA, NB;
  
  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);
  
  SAPTDFInts B_BS;
  
  B_BS.dress_ = true;
  B_BS.active_ = false;
  
  B_BS.i_length_ = noccB_;
  B_BS.j_length_ = nvirB_;
  B_BS.ij_length_ = noccB_*nvirB_;
  B_BS.i_start_ = 0;
  B_BS.j_start_ = 0;
    
  B_BS.B_d_ = block_matrix(3,noccB_*nvirB_);
  
  B_BS.filenum_ = PSIF_SAPT_BB_DF_INTS;
  B_BS.label_ = "BS RI Integrals";
    
  for (int b=0; b<noccB_; b++){
    for (int s=0; s<nvirB_; s++){
      int bs = b*nvirB_+s;
      B_BS.B_d_[0][bs] = NA*vABB_[b][s+noccB_];
    }
  }

  return(B_BS);
}

SAPTDFInts SAPT0::set_A_AB()
{
  double enuc, NA, NB;

  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);

  SAPTDFInts A_AB;

  A_AB.dress_ = true;
  A_AB.active_ = false;

  A_AB.i_length_ = noccA_;
  A_AB.j_length_ = noccB_;
  A_AB.ij_length_ = noccA_*noccB_;
  A_AB.i_start_ = 0;
  A_AB.j_start_ = 0;

  A_AB.B_d_ = block_matrix(3,noccA_*noccB_);

  A_AB.filenum_ = PSIF_SAPT_AB_DF_INTS;
  A_AB.label_ = "AB RI Integrals";

  for (int a=0; a<noccA_; a++){
    for (int b=0; b<noccB_; b++){
      int ab = a*noccB_+b;
      A_AB.B_d_[0][ab] = sAB_[a][b];
      A_AB.B_d_[1][ab] = NB*vBAB_[a][b];
      A_AB.B_d_[2][ab] = enuc*sAB_[a][b];
  }}

  return(A_AB);
}

SAPTDFInts SAPT0::set_B_AB()
{
  double enuc, NA, NB;

  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);

  SAPTDFInts B_AB;

  B_AB.dress_ = true;
  B_AB.active_ = false;

  B_AB.i_length_ = noccA_;
  B_AB.j_length_ = noccB_;
  B_AB.ij_length_ = noccA_*noccB_;
  B_AB.i_start_ = 0;
  B_AB.j_start_ = 0;

  B_AB.B_d_ = block_matrix(3,noccA_*noccB_);

  B_AB.filenum_ = PSIF_SAPT_AB_DF_INTS;
  B_AB.label_ = "AB RI Integrals";

  for (int a=0; a<noccA_; a++){
    for (int b=0; b<noccB_; b++){
      int ab = a*noccB_+b;
      B_AB.B_d_[0][ab] = NA*vAAB_[a][b];
      B_AB.B_d_[1][ab] = sAB_[a][b];
      B_AB.B_d_[2][ab] = enuc*sAB_[a][b];
  }}

  return(B_AB);
}

SAPTDFInts SAPT0::set_C_AA()
{
  SAPTDFInts C_AA;

  C_AA.dress_ = false;
  C_AA.active_ = false;

  C_AA.i_length_ = noccA_;
  C_AA.j_length_ = noccA_;
  C_AA.ij_length_ = noccA_*noccA_;
  C_AA.i_start_ = 0;
  C_AA.j_start_ = 0;

  C_AA.filenum_ = PSIF_SAPT_AA_DF_INTS;
  C_AA.label_ = "AA RI Integrals";

  return(C_AA);
}

SAPTDFInts SAPT0::set_C_AR()
{
  SAPTDFInts C_AR;

  C_AR.dress_ = false;
  C_AR.active_ = false;

  C_AR.i_length_ = noccA_;
  C_AR.j_length_ = nvirA_;
  C_AR.ij_length_ = noccA_*nvirA_;
  C_AR.i_start_ = 0;
  C_AR.j_start_ = 0;

  C_AR.filenum_ = PSIF_SAPT_AA_DF_INTS;
  C_AR.label_ = "AR RI Integrals";

  return(C_AR);
}

SAPTDFInts SAPT0::set_C_RR()
{
  SAPTDFInts C_RR;

  C_RR.dress_ = false;
  C_RR.active_ = false;

  C_RR.i_length_ = nvirA_;
  C_RR.j_length_ = nvirA_;
  C_RR.ij_length_ = nvirA_*nvirA_;
  C_RR.i_start_ = 0;
  C_RR.j_start_ = 0;

  C_RR.filenum_ = PSIF_SAPT_AA_DF_INTS;
  C_RR.label_ = "RR RI Integrals";

  return(C_RR);
}

SAPTDFInts SAPT0::set_C_BB()
{
  SAPTDFInts C_BB;

  C_BB.dress_ = false;
  C_BB.active_ = false;

  C_BB.i_length_ = noccB_;
  C_BB.j_length_ = noccB_;
  C_BB.ij_length_ = noccB_*noccB_;
  C_BB.i_start_ = 0;
  C_BB.j_start_ = 0;

  C_BB.filenum_ = PSIF_SAPT_BB_DF_INTS;
  C_BB.label_ = "BB RI Integrals";

  return(C_BB);
}

SAPTDFInts SAPT0::set_C_BS()
{
  SAPTDFInts C_BS;

  C_BS.dress_ = false;
  C_BS.active_ = false;

  C_BS.i_length_ = noccB_;
  C_BS.j_length_ = nvirB_;
  C_BS.ij_length_ = noccB_*nvirB_;
  C_BS.i_start_ = 0;
  C_BS.j_start_ = 0;

  C_BS.filenum_ = PSIF_SAPT_BB_DF_INTS;
  C_BS.label_ = "BS RI Integrals";

  return(C_BS);
}

SAPTDFInts SAPT0::set_C_SS()
{
  SAPTDFInts C_SS;

  C_SS.dress_ = false;
  C_SS.active_ = false;

  C_SS.i_length_ = nvirB_;
  C_SS.j_length_ = nvirB_;
  C_SS.ij_length_ = nvirB_*nvirB_;
  C_SS.i_start_ = 0;
  C_SS.j_start_ = 0;

  C_SS.filenum_ = PSIF_SAPT_BB_DF_INTS;
  C_SS.label_ = "SS RI Integrals";

  return(C_SS);
}

SAPTDFInts SAPT0::set_A_RB()
{
  double enuc, NA, NB;

  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);

  SAPTDFInts A_RB;

  A_RB.dress_ = true;
  A_RB.active_ = false;

  A_RB.i_length_ = nvirA_;
  A_RB.j_length_ = noccB_;
  A_RB.ij_length_ = nvirA_*noccB_;
  A_RB.i_start_ = 0;
  A_RB.j_start_ = 0;

  A_RB.B_d_ = block_matrix(3,nvirA_*noccB_);

  A_RB.filenum_ = PSIF_SAPT_AB_DF_INTS;
  A_RB.label_ = "RB RI Integrals";

  for (int r=0; r<nvirA_; r++){
    for (int b=0; b<noccB_; b++){
      int rb = r*noccB_+b;
      A_RB.B_d_[0][rb] = sAB_[r+noccA_][b];
      A_RB.B_d_[1][rb] = NB*vBAB_[r+noccA_][b];
      A_RB.B_d_[2][rb] = enuc*sAB_[r+noccA_][b];
  }}

  return(A_RB);
}

SAPTDFInts SAPT0::set_B_RB()
{
  double enuc, NA, NB;

  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);

  SAPTDFInts B_RB;

  B_RB.dress_ = true;
  B_RB.active_ = false;

  B_RB.i_length_ = nvirA_;
  B_RB.j_length_ = noccB_;
  B_RB.ij_length_ = nvirA_*noccB_;
  B_RB.i_start_ = 0;
  B_RB.j_start_ = 0;

  B_RB.B_d_ = block_matrix(3,nvirA_*noccB_);

  B_RB.filenum_ = PSIF_SAPT_AB_DF_INTS;
  B_RB.label_ = "RB RI Integrals";

  for (int r=0; r<nvirA_; r++){
    for (int b=0; b<noccB_; b++){
      int rb = r*noccB_+b;
      B_RB.B_d_[0][rb] = NA*vAAB_[r+noccA_][b];
      B_RB.B_d_[1][rb] = sAB_[r+noccA_][b];
      B_RB.B_d_[2][rb] = enuc*sAB_[r+noccA_][b];
  }}

  return(B_RB);
}

SAPTDFInts SAPT0::set_A_AS()
{
  double enuc, NA, NB;

  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);

  SAPTDFInts A_AS;

  A_AS.dress_ = true;
  A_AS.active_ = false;

  A_AS.i_length_ = noccA_;
  A_AS.j_length_ = nvirB_;
  A_AS.ij_length_ = noccA_*nvirB_;
  A_AS.i_start_ = 0;
  A_AS.j_start_ = 0;

  A_AS.B_d_ = block_matrix(3,noccA_*nvirB_);

  A_AS.filenum_ = PSIF_SAPT_AB_DF_INTS;
  A_AS.label_ = "AS RI Integrals";

  for (int a=0; a<noccA_; a++){
    for (int s=0; s<nvirB_; s++){
      int as = a*nvirB_+s;
      A_AS.B_d_[0][as] = sAB_[a][s+noccB_];
      A_AS.B_d_[1][as] = NB*vBAB_[a][s+noccB_];
      A_AS.B_d_[2][as] = enuc*sAB_[a][s+noccB_];
  }}

  return(A_AS);
}

SAPTDFInts SAPT0::set_B_AS()
{
  double enuc, NA, NB;

  NA = 1.0 / NA_;
  NB = 1.0 / NB_;
  enuc = sqrt(enuc_*NA*NB);

  SAPTDFInts B_AS;

  B_AS.dress_ = true;
  B_AS.active_ = false;

  B_AS.i_length_ = noccA_;
  B_AS.j_length_ = nvirB_;
  B_AS.ij_length_ = noccA_*nvirB_;
  B_AS.i_start_ = 0;
  B_AS.j_start_ = 0;

  B_AS.B_d_ = block_matrix(3,noccA_*nvirB_);

  B_AS.filenum_ = PSIF_SAPT_AB_DF_INTS;
  B_AS.label_ = "AS RI Integrals";

  for (int a=0; a<noccA_; a++){
    for (int s=0; s<nvirB_; s++){
      int as = a*nvirB_+s;
      B_AS.B_d_[0][as] = NA*vAAB_[a][s+noccB_];
      B_AS.B_d_[1][as] = sAB_[a][s+noccB_];
      B_AS.B_d_[2][as] = enuc*sAB_[a][s+noccB_];
  }}

  return(B_AS);
}

}}
