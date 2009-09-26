/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#ifndef _psi3_DBOC_stringblocks_h_
#define _psi3_DBOC_stringblocks_h_

#include <string>
#include <psifiles.h>
#include "float.h"

namespace psi { namespace DBOC {

/// Manages logic of arranging strings into blocks of manageble size
class StringBlocks {
 public:
  StringBlocks(int nstr, int nstr_per_block);
  ~StringBlocks();

  /// How many strings per block?
  int nstr_per_block() const;
  /// How many blocks?
  int nblocks() const;
  /// Size of the block
  int size(int block) const;
  /// To which block does this string belong?
  int block(int str) const;
  /// Index within the block
  int rel_to_block_begin(int str) const;
  /// First string in this block
  int begin(int block) const;
  /// Last string in this block
  int end(int block) const;

 private:
  int nstr_;
  int nblocks_;
  int nstr_per_block_;
};

/// Manages matrices in the basis of blocked strings
class StringBlockedMatrix {
 public:
  StringBlockedMatrix(const StringBlocks* strblk_bra, const StringBlocks* strblk_ket, const std::string& prefix);
  /// Makes a copy of A, including deep copy of the buffer
  StringBlockedMatrix(const StringBlockedMatrix& A);
  ~StringBlockedMatrix();

  StringBlocks* strblk_bra() const { return strblk_bra_; }
  StringBlocks* strblk_ket() const { return strblk_ket_; }

  FLOAT** buffer();
  void read(int brablk, int ketblk);
  void write(int brablk, int ketblk);

 private:
  std::string prefix_;
  FLOAT** buffer_;
  StringBlocks* strblk_bra_;
  StringBlocks* strblk_ket_;
  /// block size in bytes
  size_t blksize_;
  int current_brablk_, current_ketblk_;

  // PSIO status
  static const unsigned int psio_unit_ = PSIF_DBOC;
  int need_to_init_psio_;
  int unit_opened_;

  std::string key(int brablk, int ketblk);
};

}} /* namespace psi::DBOC */

#endif
