#include <libmoinfo/libmoinfo.h>

#include "blas.h"
#include "index.h"
#include "index_iterator.h"

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

        CCIndexIterator::CCIndexIterator(std::string str)
{
  nirreps = moinfo->get_nirreps();
  ccindex = blas->get_index(str);
  startup(0,nirreps);
}

CCIndexIterator::CCIndexIterator(CCIndex* index_)
{
  nirreps = moinfo->get_nirreps();
  ccindex = index_;
  startup(0,nirreps);
}

CCIndexIterator::CCIndexIterator(std::string str,int select_irrep)
{
  nirreps = moinfo->get_nirreps();
  ccindex = blas->get_index(str);
  startup(select_irrep,select_irrep + 1);
}

CCIndexIterator::CCIndexIterator(CCIndex* index_,int select_irrep)
{
  nirreps = moinfo->get_nirreps();
  ccindex = index_;
  startup(select_irrep,select_irrep + 1);
}

CCIndexIterator::~CCIndexIterator()
{
}

void CCIndexIterator::startup(int min_sym,int max_sym)
{
  nelements     = ccindex->get_nelements();
  tuples        = ccindex->get_tuples();
  element_irrep = ccindex->get_element_irrep();
  min_abs       = ccindex->get_first(min_sym);
  max_abs       = ccindex->get_last(max_sym-1);

  for(int h = min_sym; h < max_sym; ++h){
    if(ccindex->get_tuplespi(h) > 0){
      block_last.push_back(ccindex->get_last(h));
      block_symmetry.push_back(h);
    }
  }
  block_last.push_back(0);
  block_symmetry.push_back(0);
}

bool CCIndexIterator::first()
{
  if(min_abs != max_abs){
    current_block = 0;
    relative      = 0;
    absolute      = min_abs;
    symmetry      = block_symmetry[current_block];
    return(true);
  }else
    return(false);
}

void CCIndexIterator::next()
{
  ++relative;
  ++absolute;
  if(absolute >= block_last[current_block]){
    ++current_block;
    relative = 0;
    symmetry = block_symmetry[current_block];
  }
}

}}
