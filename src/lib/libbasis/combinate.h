/*! \file 
    \ingroup (BASIS)
    \brief Enter brief description of file here 
*/

// I know "combinate" is not a word

#ifndef _psi_src_lib_libbasis_combinate_h_
#define _psi_src_lib_libbasis_combinate_h_

#include <stdexcept>
#include <psitypes.h>

/** Class StatCombData contains static data of combinatorial type -
    factorials, binomial coefficients, etc. */

class StatCombData {

  int imax_;
  PSI_FLOAT* Fact_;        // Factorial
  PSI_FLOAT* Fact2_;       // Double factorial
  PSI_FLOAT** BinomC_;      // Binomial coefficients

  void check_index(int i) const
    {
      if (i < 0 || i > imax_)
	throw std::runtime_error("ERROR: StatCombData::check_index() -- index out of bounds");
    };

  // No default constructor
  StatCombData();
  // No copy constructor
  StatCombData(const StatCombData&);
  // No assignment operator
  StatCombData& operator=(const StatCombData&);

 public:
  StatCombData(int imax);
  ~StatCombData();

  /// Return n!
  PSI_FLOAT fact(int n) const { check_index(n); return Fact_[n]; };
  /// Return n!!
  PSI_FLOAT fact2(int n) const { check_index(n); return Fact2_[n]; };
  /// Return binomial coefficient (n, m), where n>=m
  PSI_FLOAT binomc(int n, int m) const { check_index(n); check_index(m); return BinomC_[n][m]; };
};

#endif
