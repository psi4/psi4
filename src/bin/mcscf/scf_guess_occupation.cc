#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <cstdio>

#include <libchkpt/chkpt.hpp>
#include <libmoinfo/libmoinfo.h>

#include "scf.h"

extern FILE* outfile;

using namespace std;

namespace psi{ namespace mcscf{

void SCF::guess_occupation()
{
  if(moinfo_scf->get_guess_occupation()){
    // Assumes the eigenvalues of some Fock operator
    // are in the SBlockVector epsilon
    vector<std::pair<double, int> > evals;

    for(int h = 0; h < nirreps; ++h)
      for(int i = 0; i < sopi[h]; ++i)
        evals.push_back( make_pair(epsilon->get(h,i),h) );

    // Sort the eigenvalues by energy
    sort(evals.begin(),evals.end());

    int ndocc = min(moinfo_scf->get_nael(),moinfo_scf->get_nbel()) - (reference == tcscf ? 1 : 0);
    int nactv = abs(moinfo_scf->get_nael()-moinfo_scf->get_nbel()) + (reference == tcscf ? 2 : 0);

    vector<int> new_docc;
    vector<int> new_actv;
    for(int h = 0; h < nirreps; ++h){
      new_docc.push_back(0);
      new_actv.push_back(0);
    }
    for(int i = 0; i < ndocc; ++i){
      new_docc[evals[i].second]++;
    }
    for(int i = ndocc; i < ndocc + nactv; ++i){
      new_actv[evals[i].second]++;
    }

    if((new_docc != docc) || (new_actv != actv)){
      fprintf(outfile,"\n\n  Occupation changed");
      fprintf(outfile,"\n  docc = (");
      for(int h = 0; h < nirreps; ++h)  fprintf(outfile," %d",new_docc[h]);
      fprintf(outfile," )");
      fprintf(outfile,"\n  actv = (");
      for(int h = 0; h < nirreps; ++h)  fprintf(outfile," %d",new_actv[h]);
      fprintf(outfile," )\n");
    }
    docc = new_docc;
    actv = new_actv;
  }
}

}} // End namespace
