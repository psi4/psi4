/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#ifndef _psi3_DBOC_dets_h_
#define _psi3_DBOC_dets_h_

namespace psi{ namespace DBOC {

/// Det is a determinant to be sorted according to "indices" Ia and Ib (can be string indices, block indices, etc.)
struct Det {
  int Ia;
  int Ib;
};

typedef std::pair<int,Det> DetI;  // Determinant + index

inline bool detcomp(DetI i, DetI j)
{
  const int Ia = i.second.Ia;
  const int Ib = i.second.Ib;
  const int Ja = j.second.Ia;
  const int Jb = j.second.Ib;
  if (Ia < Ja) return true;
  if (Ia == Ja) return Ib < Jb;
  return false;
}

}} /* namespace psi::DBOC */

#endif
