/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <algorithm>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libciomr/libciomr.h"
#include "debugging.h"
#include "psi4/libqt/qt.h"

#include "blas.h"
#include "index.h"
#include "matrix.h"


namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager *memory_manager;

using namespace std;

CCBLAS::CCBLAS(Options &options):
        options_(options),
        full_in_core(false),
        work_size(0),
        buffer_size(0)
{
  init();
}


CCBLAS::~CCBLAS()
{
  cleanup();
}

void CCBLAS::init()
{
  add_indices();
  allocate_work();
  allocate_buffer();
}

void CCBLAS::cleanup()
{
  free_sortmap();
  free_buffer();
  free_work();
  free_matrices();
  free_indices();
}

void CCBLAS::allocate_work()
{
  // Make sure work is empty
  if(!work.empty())
    for(size_t n = 0; n < work.size(); ++n)
      if(work[n]!=NULL)
        release1(work[n]);

  for(int n=0;n<options_.get_int("CC_NUM_THREADS");n++)
    work.push_back(nullptr);
  // Compute the temporary work space size
  CCIndex* oo_pair = get_index("[oo]");
  CCIndex* vv_pair = get_index("[vv]");
  CCIndex* ff_pair = get_index("[ff]");

  work_size = 0;
  for(int h=0;h<moinfo->get_nirreps();h++){
    vector<size_t> dimension;
    dimension.push_back(oo_pair->get_pairpi(h));
    dimension.push_back(vv_pair->get_pairpi(h));
    dimension.push_back(ff_pair->get_pairpi(h));
    sort(dimension.begin(),dimension.end());
    work_size += dimension[2] * dimension[1];
  }
  // Allocate the temporary work space
  for(int n=0;n<options_.get_int("CC_NUM_THREADS");n++){
    allocate1(double,work[n],work_size);
    zero_arr(work[n],work_size);
  }
  outfile->Printf("\n  Allocated work array of size %ld (%.2f MiB)",work_size * sizeof(double),type_to_MiB<double>(work_size));
}

void CCBLAS::allocate_buffer()
{
  // Make sure buffer is empty
  if(!buffer.empty())
    for(size_t n = 0; n < buffer.size(); ++n)
      if(buffer[n]!=NULL)
        release1(buffer[n]);

  for(int n=0;n<options_.get_int("CC_NUM_THREADS");n++)
    buffer.push_back(nullptr);
  // Compute the temporary buffer space size, 101% of the actual strip size
  buffer_size = static_cast<size_t>(1.01 * CCMatrix::fraction_of_memory_for_buffer *
                                    static_cast<double>(memory_manager->get_FreeMemory()) /
                                    static_cast<double>(sizeof(double)));
  // The value used here , 0.05 is also used in

  // Allocate the temporary buffer space
  for(int n=0;n<options_.get_int("CC_NUM_THREADS");n++){
    allocate1(double,buffer[n],buffer_size);
    zero_arr(buffer[n],buffer_size);
  }
  outfile->Printf("\n  Allocated buffer array of size %ld (%.2f MiB)",buffer_size * sizeof(double),type_to_MiB<double>(buffer_size));
}

void CCBLAS::free_sortmap()
{
  for(SortMap::iterator iter=sortmap.begin();iter!=sortmap.end();++iter){
    for(int irrep=0;irrep<moinfo->get_nirreps();irrep++)
      delete[] iter->second[irrep];
    delete[] iter->second;
  }
}

void CCBLAS::free_work()
{
  // Delete the temporary work space
  for(size_t n = 0; n < work.size(); ++n){
    if(work[n]!=NULL){
      release1(work[n]);
    }
  }
}

void CCBLAS::free_buffer()
{
  // Delete the temporary buffer space
  for(size_t n = 0; n < buffer.size(); ++n){
    if(buffer[n]!=NULL){
      release1(buffer[n]);
    }
  }
}

void CCBLAS::free_matrices()
{
  for(MatrixMap::iterator iter=matrices.begin();iter!=matrices.end();++iter){
    delete iter->second;
  }
}

void CCBLAS::free_indices()
{
  for(IndexMap::iterator iter=indices.begin();iter!=indices.end();++iter){
    delete iter->second;
  }
}

void CCBLAS::add_indices()
{
  add_index("[]");
  add_index("[o]");
  add_index("[v]");
  add_index("[a]");
  add_index("[f]");
  add_index("[o>o]");
  add_index("[v>v]");
  add_index("[v>=v]");
  add_index("[oo]");
  add_index("[ov]");
  add_index("[vo]");
  add_index("[vv]");
  add_index("[aa]");

  add_index("[aaa]");
  add_index("[ooo]");
  add_index("[oov]");
  add_index("[voo]");
  add_index("[ovv]");
  add_index("[vvo]");
  add_index("[ovo]");

  // MP3 PCBS
  add_index("[fo]");
  add_index("[of]");
  add_index("[ff]");
  add_index("[vf]");
  add_index("[fv]");

  add_index("[ovf]");
  add_index("[ofv]");
  add_index("[foo]");
  add_index("[off]");

  // MP2-CCSD
  if(options_.get_str("CORR_WFN")=="MP2-CCSD"){
    add_index("[oav]");
    add_index("[ova]");
    add_index("[avo]");
    add_index("[aao]");
    add_index("[aoa]");
    add_index("[oaa]");
    add_index("[vaa]");
    add_index("[aav]");
    add_index("[ava]");
  }
  if(options_.get_str("CORR_WFN")!="PT2"){
    add_index("[vvv]");
  }

  // Mk-MRPT2
  add_index("[ao]");
  add_index("[av]");

  // Not useful
  add_index("[oa]");
  add_index("[va]");
}

// void CCBLAS::allocate_matrices_in_core()
// {
//   for(MatrixMap::iterator iter=matrices.begin();iter!=matrices.end();++iter){
//     CCMatrix* Matrix = iter->second;
//     outfile->Printf("\n%s(analyzing)",Matrix->get_label().c_str());
//
//     if(Matrix->get_out_of_core()){
//       Matrix->load();
//       outfile->Printf("\n%s <- reading from disk",Matrix->get_label().c_str());
//
//     }else if(!Matrix->is_allocated())
//       Matrix->allocate_memory();
//   }
// }

void CCBLAS::print(const char* cstr)
{
  string str(cstr);
  vector<string> names = moinfo->get_matrix_names(str);
  for(size_t n = 0; n < names.size(); ++n)
    print_ref(names[n]);
}

void CCBLAS::print_ref(string& str)
{
  get_Matrix(str)->print();
}

void CCBLAS::print_memory()
{
//  size_t total_memory_required = 0;
//  outfile->Printf("\n\n\t-----------------------------------------------------------------------------");
//  outfile->Printf("\n\tMatrix ID    Memory(bytes)   Cumulative Memory(bytes)  Accessed    Label");
//  outfile->Printf("\n\t------------------------------------------------------------------------------");

//  for(MatrixMap::iterator iter=matrices.begin();iter!=matrices.end();++iter){
//    total_memory_required += iter->second->get_memory2();
//    outfile->Printf("\n  %4d",distance(matrices.begin(),iter));
//    outfile->Printf("     %14d",iter->second->get_memory2());
//    outfile->Printf("        %14d",total_memory_required);
//    outfile->Printf("             %4d",iter->second->get_naccess());
//    outfile->Printf("         %s",iter->second->get_label().c_str());
//  }
//  outfile->Printf("\n\t------------------------------------------------------------------------------");
//  outfile->Printf("\n\n\tTotal memory required for matrices = %14d (bytes)\n",total_memory_required);

//  total_memory_required = 0;

//  outfile->Printf("\n\n\t-------------------------------------------------------------");
//  outfile->Printf("\n\tIndex ID    Memory(MB)   Cumulative Memory(MB)     Label");
//  outfile->Printf("\n\t--------------------------------------------------------------");

//  for(IndexMap::iterator iter=indices.begin();iter!=indices.end();++iter){
//    total_memory_required += iter->second->get_memory();
//    outfile->Printf("\n\t%4d",distance(indices.begin(),iter));
//    outfile->Printf("     %10.2f",iter->second->get_memory());
//    outfile->Printf("        %10.2f",total_memory_required);
//    outfile->Printf("         %s",iter->second->get_label().c_str());
//  }
//  outfile->Printf("\n\t--------------------------------------------------------------");

//  outfile->Printf("\n\n\tTotal memory required for indexing = %10.2f (Mb)\n",total_memory_required);
}

/**
 * This routine computes which quantities have to be initially stored in memory and which on disk
 */
int CCBLAS::compute_storage_strategy()
{
  outfile->Printf("\n\n  Computing storage strategy:");

  // N.B. Here I am using bytes as the basic unit
  size_t available_memory     = memory_manager->get_FreeMemory();
  double fraction_for_in_core = 0.97; // Fraction of the total available memory that may be used
  size_t storage_memory       = static_cast<size_t>(static_cast<double>(available_memory) * fraction_for_in_core);
  size_t fully_in_core_memory = 0;
  size_t integrals_memory     = 0;
  size_t fock_memory          = 0;
  size_t others_memory        = 0;

  outfile->Printf("\n    Input memory                           = %14lu bytes",
                  (unsigned long)memory_manager->get_MaximumAllowedMemory());
  outfile->Printf("\n    Free memory                            = %14lu bytes",
                  (unsigned long)available_memory);
  outfile->Printf("\n    Free memory available for matrices     = %14lu bytes (%3.0f%%)",
                  (unsigned long)storage_memory,fraction_for_in_core*100.0);

  // Gather the memory requirements for all the CCMAtrix object
  // and divide the integrals from all the other matrices.
  // At the same time compute the memory requirements for
  // a fully in-core algorithm.
  vector<pair<size_t,pair<CCMatrix*,int> > > integrals;
  vector<pair<size_t,pair<CCMatrix*,int> > > fock;
  vector<pair<size_t,pair<CCMatrix*,int> > > others;
  for(MatrixMap::iterator it=matrices.begin();it!=matrices.end();++it){
    for(int h=0;h<moinfo->get_nirreps();++h){
      size_t block_memory = it->second->get_memorypi2(h);
      if(it->second->is_integral()){
        integrals.push_back(make_pair(block_memory,make_pair(it->second,h)));
        integrals_memory += block_memory;
      }else if(it->second->is_fock()){
        fock.push_back(make_pair(block_memory,make_pair(it->second,h)));
        fock_memory += block_memory;
      }else{
        others.push_back(make_pair(block_memory,make_pair(it->second,h)));
        others_memory += block_memory;
      }
      fully_in_core_memory += block_memory;
    }
  }
  outfile->Printf("\n    Memory required by fock matrices       = %14lu bytes",(unsigned long)fock_memory);
  outfile->Printf("\n    Memory required by integrals           = %14lu bytes",(unsigned long)integrals_memory);
  outfile->Printf("\n    Memory required by other matrices      = %14lu bytes",(unsigned long)others_memory);
  outfile->Printf("\n    Memory required for in-core algorithm  = %14lu bytes",(unsigned long)fully_in_core_memory);

  // Check if you may use a fully in core algorithm
  full_in_core = false;
  int strategy = 0;
  if(fully_in_core_memory < storage_memory ){
    full_in_core = true;
    outfile->Printf("\n    PSIMRCC will perform a full in-core computation");
    strategy = 0;
  }else{
    if(others_memory < storage_memory ){
      outfile->Printf("\n    PSIMRCC will store some integrals out-of-core");
      strategy = 1;
    }else{
      outfile->Printf("\n    PSIMRCC will store all integrals and some other matrices out-of-core");
      strategy = 2;
      throw PSIEXCEPTION("CCBLAS::compute_storage_strategy(): Strategy #2 is not implemented yet");
    }
  }
  sort(integrals.begin(),integrals.end());
  sort(others.begin(),others.end());
  for(size_t i = 0; i < fock.size(); ++i){
    // Store all the fock matrices in core and allocate them
    storage_memory -= fock[i].first;
    load_irrep(fock[i].second.first,fock[i].second.second);
  }
  // Let the CCBlas class worry about allocating matrices
  int number_of_others_on_disk = 0;
  for(size_t i = 0;i < others.size(); ++i){
    // Check if this matrix can be stored in core
    if(others[i].first < storage_memory){
      storage_memory -= others[i].first;
      load_irrep(others[i].second.first,others[i].second.second);
    }else{
      number_of_others_on_disk++;
    }
  }
  int number_of_integrals_on_disk = 0;
  for(size_t i = 0; i < integrals.size(); ++i){
    // Check if this matrix can be stored in core
    if(integrals[i].first < storage_memory){
      storage_memory -= integrals[i].first;
      load_irrep(integrals[i].second.first,integrals[i].second.second);
    }else{
      number_of_integrals_on_disk++;
    }
  }

//  DEBUGGING(1,
//    outfile->Printf("\n    -------------------- Fock matrices -------------------------");
//    for(int i=0;i<fock.size();i++){
//      if(fock[i].first > 1.0e-5){
//        outfile->Printf("\n    %-32s irrep %d   %6.2f Mb --> ",fock[i].second.first->get_label().c_str(),
//                                                      fock[i].second.second,
//                                                      fock[i].first);
//        outfile->Printf("%s",fock[i].second.first->is_block_allocated(fock[i].second.second) ? "in-core" : "out-of-core");
//      }
//    }
//    outfile->Printf("\n    -------------------- Other matrices ------------------------");
//    for(int i=0;i<others.size();i++){
//      if(others[i].first > 1.0e-5){
//        outfile->Printf("\n    %-32s irrep %d   %6.2f Mb --> ",others[i].second.first->get_label().c_str(),
//                                                      others[i].second.second,
//                                                      others[i].first);
//        outfile->Printf("%s",others[i].second.first->is_block_allocated(others[i].second.second) ? "in-core" : "out-of-core");
//      }
//    }
//    outfile->Printf("\n    -------------------- Integrals -----------------------------");
//    for(int i=0;i<integrals.size();i++){
//      if(integrals[i].first > 1.0e-5){
//        outfile->Printf("\n    %-32s irrep %d   %6.2f Mb --> ",integrals[i].second.first->get_label().c_str(),
//                                                      integrals[i].second.second,
//                                                      integrals[i].first);
//        outfile->Printf("%s",integrals[i].second.first->is_block_allocated(integrals[i].second.second) ? "in-core" : "out-of-core");
//      }
//    }
//    outfile->Printf("\n\n");
//  );

  if(!full_in_core){
    outfile->Printf("\n    Out-of-core algorithm will store %d other matrices on disk",number_of_others_on_disk);
    outfile->Printf("\n    Out-of-core algorithm will store %d integrals on disk",number_of_integrals_on_disk);
  }
  return(strategy);
}

/**
 * This routine computes which quantities have to be initially stored in memory and which on disk
 */
void CCBLAS::show_storage()
{
//  DEBUGGING(1,
//    outfile->Printf("\n    ----------------------------------");
//    outfile->Printf("\n               Show Storage ");
//    outfile->Printf("\n    ----------------------------------");

//    for(MatrixMap::iterator it=matrices.begin();it!=matrices.end();++it){
//      for(int h=0;h<moinfo->get_nirreps();++h){
//        double block_memory = it->second->get_memorypi(h);
//        outfile->Printf("\n    %-32s irrep %d   %6.2f Mb",it->second->get_label().c_str(),h,block_memory);
//        outfile->Printf(" is %s",it->second->is_block_allocated(h) ? "allocated" : "not allocated");
//        outfile->Printf("%s",it->second->is_out_of_core(h) ? "(out-of-core)" : "");

//      }
//    }
//  )
}

}} /* End Namespaces */
