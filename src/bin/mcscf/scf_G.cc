#include <iostream>

#include <psifiles.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libutil/libutil.h>

#include "scf.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

void SCF::construct_G(SBlockMatrix& density,
                      SBlockMatrix& G,
                      double* integrals,
                      int batch)
{
  construct_G(density,G,integrals,batch,1.0);
}


void SCF::construct_G(SBlockMatrix& density,
                      SBlockMatrix& G,
                      double* integrals,
                      int batch,
                      double factor)
{
  double* D_vector;
  double* G_vector;

  allocate1(double,D_vector,npairs);
  allocate1(double,G_vector,npairs);

  // Convert D to a vector and double the off diagonal elements
  for(int h = 0; h < nirreps; ++h){
    for(int p = 0; p < sopi[h]; ++p){
      int p_abs = p + block_offset[h];
      for(int q = 0; q <= p; ++q){
        int q_abs = q + block_offset[h];
        D_vector[ pair[p_abs][q_abs] ] = 2.0 * density->get(h,p,q);
        G_vector[ pair[p_abs][q_abs] ] = 0.0;
      }
      D_vector[ pair[p_abs][p_abs] ]  *= 0.5;
    }
  }

  // general algorithm
  double G_pq,D_pq;
  double* D_rs;
  double* G_rs;
  int     pq,rs;
  double* PK_pqrs = integrals;

  for(pq = batch_pq_min[batch]; pq < batch_pq_max[batch]; ++pq){
    G_pq = 0.0;
    D_pq = D_vector[pq];
    D_rs = &D_vector[0];
    G_rs = &G_vector[0];
    for(rs = 0; rs <= pq; ++rs){
      G_pq += *PK_pqrs * (*D_rs);
      *G_rs += *PK_pqrs * D_pq;
      ++D_rs;
      ++G_rs;
      ++PK_pqrs;
    }
    G_vector[pq] += G_pq;
  }

  // Convert G to a matrix
  for(int h = 0; h < nirreps; ++h){
    for(int p = 0; p < sopi[h]; ++p){
      int p_abs = p + block_offset[h];
      for(int q = 0; q < sopi[h]; ++q){
        int q_abs = q + block_offset[h];
        G->set(h,p,q,2.0 * factor * G_vector[ pair[p_abs][q_abs]]);
      }
    }
  }

  release1(G_vector);
  release1(D_vector);
}



void SCF::read_Raffanetti(const char* integral_type, double* integrals, int batch)
{
  // Read the PK matrix from disk
  char data_label[80];
  sprintf(data_label,"%s_%d",integral_type,batch);
  size_t buffer_size = batch_size[batch] * sizeof(double);
  psio_->read_entry(PSIF_MCSCF,data_label,(char*)integrals,buffer_size);
}

void SCF::write_Raffanetti(const char* integral_type, double* integrals, int batch)
{
  // Write the PK matrix to disk
  char data_label[80];
  sprintf(data_label,"%s_%d",integral_type,batch);
  size_t buffer_size = batch_size[batch] * sizeof(double);
  psio_->write_entry(PSIF_MCSCF,data_label,(char*)integrals,buffer_size);
}


}} /* End Namespaces */
