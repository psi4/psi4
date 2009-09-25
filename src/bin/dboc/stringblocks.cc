/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cstring>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include "stringblocks.h"
#include "linalg.h"

using namespace psi::dboc;

StringBlocks::StringBlocks(int nstr, int nstr_per_block) :
  nstr_(nstr), nstr_per_block_(nstr_per_block)
{
  nblocks_ = (nstr_ % nstr_per_block_) ?  nstr_ / nstr_per_block_ + 1 : nstr_ / nstr_per_block_;
  if (nstr_per_block_ > nstr_)
    throw std::runtime_error("StringBlocks::StringBlocks -- nstr less than nstr_per_block");
}

StringBlocks::~StringBlocks()
{
}

int
StringBlocks::nstr_per_block() const
{
  return nstr_per_block_;
}

int
StringBlocks::nblocks() const
{
  return nblocks_;
}

int
StringBlocks::block(int str) const
{
  return str / nstr_per_block_;
}

int
StringBlocks::rel_to_block_begin(int str) const
{
  return str % nstr_per_block_;
}

int
StringBlocks::begin(int block) const
{
  return block*nstr_per_block_;
}

int
StringBlocks::end(int block) const
{
  return block == nblocks_-1 ? nstr_-1 : (block+1)*nstr_per_block_-1;
}

///////////

#define PSIO_INIT if (!psio_state()) { \
    psio_init(); psio_ipv1_config(); \
    need_to_init_psio_ = 1; \
  }

#define PSIO_OPEN(u,n) if (!psio_open_check(u)) { \
    psio_open((u),n); \
    unit_opened_ = 0; \
  }

#define PSIO_CLOSE(u) if (!unit_opened_) \
    psio_close((u),1);

#define PSIO_DONE if (need_to_init_psio_) \
    psio_done();

StringBlockedMatrix::StringBlockedMatrix(const StringBlocks* strblk_bra, const StringBlocks* strblk_ket,
					 const std::string& prefix) :
  prefix_(prefix),
  strblk_bra_(const_cast<StringBlocks*>(strblk_bra)),
  strblk_ket_(const_cast<StringBlocks*>(strblk_ket)),
  buffer_(0),
  current_brablk_(-1), current_ketblk_(-1),
  need_to_init_psio_(0), unit_opened_(1)
{
  buffer_ = create_matrix(strblk_bra->nstr_per_block(), strblk_ket->nstr_per_block());
  blksize_ = (size_t)strblk_bra->nstr_per_block() * strblk_ket->nstr_per_block() * sizeof(FLOAT);

PSIO_INIT
PSIO_OPEN(psio_unit_,PSIO_OPEN_NEW)

}

StringBlockedMatrix::StringBlockedMatrix(const StringBlockedMatrix& A) :
  prefix_(A.prefix_),
  strblk_bra_(A.strblk_bra_),
  strblk_ket_(A.strblk_ket_),
  buffer_(0),
  current_brablk_(-1), current_ketblk_(-1),
  need_to_init_psio_(0), unit_opened_(1)
{
  buffer_ = create_matrix(strblk_bra_->nstr_per_block(), strblk_ket_->nstr_per_block());
  blksize_ = (size_t)strblk_bra_->nstr_per_block() * strblk_ket_->nstr_per_block() * sizeof(FLOAT);

  void* tmp = memcpy(static_cast<void*>(buffer_[0]),static_cast<void*>(A.buffer_[0]),blksize_);

PSIO_INIT
PSIO_OPEN(psio_unit_,PSIO_OPEN_NEW)

}

StringBlockedMatrix::~StringBlockedMatrix()
{
PSIO_CLOSE(psio_unit_)
PSIO_DONE
  delete_matrix(buffer_,strblk_bra_->nstr_per_block(), strblk_ket_->nstr_per_block());
}

FLOAT**
StringBlockedMatrix::buffer()
{
  return buffer_;
}

std::string
StringBlockedMatrix::key(int brablk, int ketblk)
{
#if 0
  std::ostringstream oss;
  oss << prefix_ << "_blk_" << brablk << "_" << ketblk;
  return oss.str();
#endif
  char result[128];
  sprintf(result,"%s_blk_%d_%d",prefix_.c_str(),brablk,ketblk);
  return std::string(result);
}

void
StringBlockedMatrix::write(int brablk, int ketblk)
{

  if (strblk_bra_->nblocks() == 1 && strblk_ket_->nblocks() == 1)
    return;

  psio_address junk;
  const std::string bkey = key(brablk,ketblk);
  const char* bkey_cstr = bkey.c_str();
  int errcod = psio_write(psio_unit_,
			  bkey_cstr,
			  reinterpret_cast<char*>(buffer_[0]),
			  blksize_,
			  PSIO_ZERO,
			  &junk);
}

void
StringBlockedMatrix::read(int brablk, int ketblk)
{

  if (strblk_bra_->nblocks() == 1 && strblk_ket_->nblocks() == 1)
    return;

  if (current_brablk_ != brablk ||
      current_ketblk_ != ketblk) {

#if 0
    std::cout << "StringBlockedMatrix::read -- ("
	      << brablk << "," << ketblk << ")  ("
	      << current_brablk_ << "," << current_ketblk_ << ")" << std::endl;
#endif

    current_brablk_ = brablk;
    current_ketblk_ = ketblk;
    
    psio_address junk;
    const std::string bkey = key(brablk,ketblk);
    const char* bkey_cstr = bkey.c_str();
    int errcod = psio_read(psio_unit_,
			   bkey_cstr,
			   reinterpret_cast<char*>(buffer_[0]),
			   blksize_,
			   PSIO_ZERO,
			   &junk);
  }
}
