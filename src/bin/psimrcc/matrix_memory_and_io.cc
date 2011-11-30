/***************************************************************************
 *   Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *   frank@ccc.uga.edu
 *   SR/MRCC Code
 ***************************************************************************/

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <libpsio/psio.hpp>

#include "debugging.h"
#include "matrix.h"

namespace psi{
    extern FILE *outfile;
    namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager *memory_manager;

using namespace std;


/*********************************************************
  Memory Allocation Routines
*********************************************************/

/**
 * Return true if all the blocks are written to core
 * @return
 */
bool CCMatrix::is_out_of_core()
{
  for(int h=0;h<moinfo->get_nirreps();++h)
    if(!out_of_core[h] && (block_sizepi[h]>0))
      return(false);
  return(true);
}

/**
 * Return true if all the blocks are allocated
 * @return
 */
bool CCMatrix::is_allocated()
{
  for(int h=0;h<moinfo->get_nirreps();++h)
    if(!is_block_allocated(h) && (block_sizepi[h]>0))
      return(false);
  return(true);
}

bool CCMatrix::is_block_allocated(int h)
{
  if(matrix[h]==NULL)
    return(false);
  else
    return(true);
}

void CCMatrix::allocate_memory()
{
  // Allocate a matrix if we are not checking
  for(int h=0;h<nirreps;h++)
    allocate_block(h);
}

void CCMatrix::allocate_block(int h)
{
  if(block_sizepi[h]>0){
    if(!is_block_allocated(h)){
      if(memorypi2[h] < memory_manager->get_FreeMemory()){
        allocate2(double,matrix[h],left_pairpi[h],right_pairpi[h]);
        DEBUGGING(2,
          fprintf(outfile,"\n  %s[%s] <- allocated",label.c_str(),moinfo->get_irr_labs(h));
          fflush(outfile);
        )
      }else{
        fprintf(outfile,"\n\nNot enough memory to allocate irrep %d of %s\n",h,label.c_str());
        fflush(outfile);
        exit(1);
      }
    }else{
      fprintf(outfile,"\n\nCCMatrix::allocate_block(): You are trying to allocate irrep %d of %s when is already allocated!!!\n",h,label.c_str());
      fflush(outfile);
      exit(EXIT_FAILURE);
    }
  }
}


/**
 * Free the memory used to store the matrix elements
 */
void CCMatrix::free_memory()
{
  for(int h=0;h<nirreps;h++)
    free_block(h);
}

void CCMatrix::free_block(int h)
{
  if(block_sizepi[h]>0){
    if(is_block_allocated(h)){
      release2(matrix[h]);
      DEBUGGING(2,
        fprintf(outfile,"\n  %s[%s] <- deallocated",label.c_str(),moinfo->get_irr_labs(h));
        fflush(outfile);
      )
    }
  }
}

/*********************************************************
  I/O Routines
*********************************************************/

/**
 * Write the matrix to disk and free the memory.
 */
void CCMatrix::dump_to_disk()
{
  dump_to_disk(0,moinfo->get_nirreps());
}

/**
 * Write the matrix to disk and free the memory
 */
void CCMatrix::dump_to_disk(int first_irrep,int last_irrep)
{
  for(int h=first_irrep;h<last_irrep;++h)
    dump_block_to_disk(h);
}

/**
 * Write a irrep block to disk and free the memory
 * @param h irrep to write to disk
 */
void CCMatrix::dump_block_to_disk(int h)
{
  write_block_to_disk(h);
  free_block(h);
  out_of_core[h]=true;
}

/**
 * Write a irrep block to disk without freeing the memory
 * @param h irrep to write to disk
 */
void CCMatrix::write_block_to_disk(int h)
{
  if(block_sizepi[h]>0){
    // for generic matrices store the entire symmetry block on disk
    if(!is_integral()){
      char data_label[80];
      sprintf(data_label,"%s_%d",label.c_str(),h);
      _default_psio_lib_->write_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)&(matrix[h][0][0]),block_sizepi[h]*sizeof(double));
    }else{
//       fprintf(outfile,"\n    CCMatrix::write_block_to_disk(): writing %s irrep %d to disk",label.c_str(),h);
//       fprintf(outfile,"\n    This is a %d x %d block",left_pairpi[h],right_pairpi[h]);
      // for two electron integrals store strips of the symmetry block on disk
      size_t max_strip_size = static_cast<size_t>(fraction_of_memory_for_buffer * static_cast<double>(memory_manager->get_FreeMemory()));

      int    strip          = 0;
      size_t last_row       = 0;
      // Determine the size of the strip and write the strip to disk
      while(last_row < left_pairpi[h]){
        // Determine the size of the strip
        size_t first_row      = last_row;
        size_t strip_size     = 0.0; // Mb
        size_t strip_length   = 0;
        while((strip_size < max_strip_size) && (last_row < left_pairpi[h]) ){
          last_row++;
          strip_length = last_row - first_row;
          strip_size = sizeof(double) * strip_length * right_pairpi[h];
        }
//         fprintf(outfile,"\n    Writing strip %d of lenght %d (%d -> %d)",strip,strip_length,first_row,last_row);
        // Write the size of the strip
        char size_label[80];
        sprintf(size_label ,"%s_%d_%d_size",label.c_str(),h,strip);
        _default_psio_lib_->write_entry(PSIF_PSIMRCC_INTEGRALS,size_label,(char*)&(strip_length),sizeof(size_t));

        // Write the strip
        char data_label[80];
        sprintf(data_label,"%s_%d_%d",label.c_str(),h,strip);
        _default_psio_lib_->write_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)&(matrix[h][first_row][0]),
                         strip_length*right_pairpi[h]* sizeof(double));
        strip++;
      }

//       fprintf(outfile,"\n    Written %d strip%s",strip,strip>1 ? "" : "s");
      // Write the number of strips
      char nstrips_label[80];
      sprintf(nstrips_label ,"%s_%d_nstrips",label.c_str(),h);
      _default_psio_lib_->write_entry(PSIF_PSIMRCC_INTEGRALS,nstrips_label,(char*)&(strip),sizeof(int));
    }
  }
}

/**
 * A black-box version of read_from_disk() that can be called for any matrix
 */
void CCMatrix::load()
{
  if(is_out_of_core()){
    if(!is_allocated())
      read_from_disk();
  }else
    if(!is_allocated())
      allocate_memory();
}

/**
 * A black-box version of read_from_disk() that can be called for any matrix
 * @param h irrep to read from disk
 */
void CCMatrix::load_irrep(int h)
{
  if(out_of_core[h]){
    if(!is_block_allocated(h))
      read_block_from_disk(h);
  }else{
    if(!is_block_allocated(h))
      allocate_block(h);
  }
}

/**
 * Read a matrix from disk.
 */
void CCMatrix::read_from_disk()
{
  read_from_disk(0,moinfo->get_nirreps());
}

/**
 * Read irrep blocks from disk
 * @param h irrep to write to disk
 */
void CCMatrix::read_from_disk(int first_irrep,int last_irrep)
{
  for(int h=first_irrep;h<last_irrep;++h)
    read_block_from_disk(h);
}

/**
 * Read an irrep block from disk
 * @param h irrep to write to disk
 */
void CCMatrix::read_block_from_disk(int h)
{
  if(block_sizepi[h]>0){
    if(!is_block_allocated(h))
      allocate_block(h);
    // for generic matrices read the entire symmetry block on disk
    if(!is_integral()){
      char data_label[80];
      sprintf(data_label,"%s_%d",label.c_str(),h);
      _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)&(matrix[h][0][0]),block_sizepi[h]*sizeof(double));
    }else{
      // Read the number of strips
      int nstrips = 0;
      char nstrips_label[80];
      sprintf(nstrips_label ,"%s_%d_nstrips",label.c_str(),h);
      _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,nstrips_label,(char*)&(nstrips),sizeof(int));
      size_t first_row =0;
      for(int strip = 0;strip!=nstrips;++strip){
        // Read the size of the strip
        size_t strip_length =  0;
        char size_label[80];
        sprintf(size_label ,"%s_%d_%d_size",label.c_str(),h,strip);
        _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,size_label,(char*)&(strip_length),sizeof(size_t));

        // Read the strip
        char data_label[80];
        sprintf(data_label,"%s_%d_%d",label.c_str(),h,strip);
        _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)&(matrix[h][first_row][0]),
                         strip_length*right_pairpi[h]*sizeof(double));

        first_row += strip_length;
      }
    }
  }
}

/**
 * Read an irrep strip from disk and return a boolean that is true if there is strip
 * @param h irrep to write to disk
 */
size_t CCMatrix::read_strip_from_disk(int h, int strip, double* buffer)
{
  size_t strip_length =  0;
  if(block_sizepi[h]>0){
    // for generic matrices read the entire symmetry block on disk
    if(!is_integral()){
      fprintf(outfile,"\nMatrix %s is not stored in strips!!!",label.c_str());
      fflush(outfile);
      exit(EXIT_FAILURE);
    }else{
      // Read the number of strips
      int nstrips = 0;
      char nstrips_label[80];
      sprintf(nstrips_label ,"%s_%d_nstrips",label.c_str(),h);
      _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,nstrips_label,(char*)&(nstrips),sizeof(int));
      if(strip < nstrips){
        // Read the size of the strip
        char size_label[80];
        sprintf(size_label ,"%s_%d_%d_size",label.c_str(),h,strip);
        _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,size_label,(char*)&(strip_length),sizeof(size_t));

        // Read the strip
        char data_label[80];
        sprintf(data_label,"%s_%d_%d",label.c_str(),h,strip);
        _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)buffer,
                          strip_length*right_pairpi[h]*sizeof(double));
      }
    }
  }
  return(strip_length);
}

}} /* End Namespaces */
