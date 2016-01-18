#include <vector>
#include <libmints/mints.h>

using namespace std;

namespace psi { namespace cctransort {

vector<int> pitzer2qt(vector<Dimension> &spaces)
{
  int nirreps = spaces[0].n();

  Dimension total(nirreps);
  for(int h=0; h < nirreps; h++)
    for(int i=0; i < spaces.size(); i++)
      total[h] += spaces[i][h];
  int nmo = total.sum();

  vector<int> order(nmo);
  order.assign(nmo, 0);

  Dimension offset(nirreps);
  offset[0] = 0;
  for(int h=1; h < nirreps; h++)
    offset[h] = offset[h-1] + total[h-1];

  int count = 0;

  for(int j=0; j < spaces.size(); j++)
    for(int h=0; h < nirreps; h++) {
      int this_offset = offset[h];
      for(int k=0; k < j; k++) this_offset += spaces[k][h];
      for(int i=0; i < spaces[j][h]; i++) 
      order[this_offset + i] = count++;
    }

  return order;
}

}} 
