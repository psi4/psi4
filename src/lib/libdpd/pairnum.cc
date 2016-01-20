#include "dpd.h"

namespace psi {

int DPD::pairnum(string pair)
{
  vector<string> v = dpd_split(pair);

  int left, right;

  if(v.size() == 2) { // "pq"
    for(int i=0; i < moSpaces.size(); i++) {
      if(v[0] == moSpaces[i]) left = i;
      if(v[1] == moSpaces[i]) right = i;
    }
    if(left == right) return left*5; // unrestricted diagonal pair
    else if(left < right) return moSpaces.size()*5 + 2*(left*moSpaces.size()-left*(left+1)/2) + 2*(right-left-1);
    else if(left > right) return moSpaces.size()*5 + 2*(right*moSpaces.size()-right*(right+1)/2) + 2*(left-right-1) + 1;
  }
  else if(v.size() == 4) { // "p>q+" or "p>q-"
    for(int i=0; i < moSpaces.size(); i++) {
      if(v[0] == moSpaces[i]) left = i;
      if(v[2] == moSpaces[i]) right = i;
    }
    if(left != right) { throw; }
    if(v[3] == "+") return left*5 + 1;
    else if(v[3] == "-") return left*5 + 2;
  }
  else if(v.size() == 5) { // "p>=q+" or "p>=q-"
    for(int i=0; i < moSpaces.size(); i++) {
      if(v[0] == moSpaces[i]) left = i;
      if(v[3] == moSpaces[i]) right = i;
    }
    if(left != right) { throw; }
    if(v[4] == "+") return left*5 + 3;
    else if(v[4] == "-") return left*5 + 4;
  }
}

} // namespace psi
