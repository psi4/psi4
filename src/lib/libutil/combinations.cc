#include <cstdlib>
#include <algorithm>

#include "libutil.h"

namespace psi {

/**
 * Generate combinations of 0,1,...,(n-1) taken k at a time
 * @param n
 * @param k
 * @param combinations a vector<vector<int> > that will store all the combinations
 */
void generate_combinations(int n, int k, std::vector<std::vector<int> >& combinations)
{
  if( (n > 0) && (k > 0)){
    std::vector<int> combination;
    bool* a = new bool[n];
    for(int i=0;i<n-k;++i)
      a[i] = false;
    for(int i=n-k;i<n;++i)
      a[i] = true;
    do{
      combination.clear();
      for ( int i = 0 ; i < n ; ++i){
        if(a[i])
          combination.push_back(i);
      }
      combinations.push_back(combination);
    } while (std::next_permutation(a,a+n));
    delete[] a;  
  }
}

}

