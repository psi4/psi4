#include <cmath>
#include <algorithm>

#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>


#define CCTRANSFORM_USE_BLAS

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include "psifiles.h"

#include "algebra_interface.h"
#include "blas.h"
#include "index.h"
#include "matrix.h"
#include "transform.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager *memory_manager;

using namespace std;

CCTransform::CCTransform() : fraction_of_memory_for_presorting(0.75)
{
  blas->add_index("[s>=s]");
  blas->add_index("[n>=n]");
  blas->add_index("[s]");
  tei_mo_indexing = blas->get_index("[n>=n]");
  tei_so_indexing = blas->get_index("[s>=s]");
  oei_so_indexing = blas->get_index("[s]");
  ioff = moinfo->get_ioff();
  first_irrep_in_core = 0;
  last_irrep_in_core  = 0;
  s_so   = NULL;
  oei_so = NULL;
  oei_mo = NULL;
  tei_so = NULL;
  tei_mo = NULL;
  tei_half_transformed = NULL;
}

CCTransform::~CCTransform()
{
  free_memory();
}

void CCTransform::read_mo_integrals()
{
  read_oei_mo_integrals();
  read_tei_mo_integrals();
}

void CCTransform::read_so_integrals()
{
  read_tei_so_integrals();
}

/*!
    \fn CCTransform::read_tei_integrals()
 */
void CCTransform::read_tei_so_integrals()
{
  // Read all the (frozen + non-frozen) TEI in Pitzer order
  // and store them in a in-core block-matrix
  CCIndex* indexing = blas->get_index("[s>=s]");

  // Allocate the tei_so matrix blocks
  allocate1(double*,tei_so,moinfo->get_nirreps());

  for(int h=0;h<moinfo->get_nirreps();h++){
    if(indexing->get_pairpi(h)>0){
      size_t block_size = INDEX(indexing->get_pairpi(h)-1,indexing->get_pairpi(h)-1)+1;
      allocate1(double,tei_so[h],block_size);
      for(size_t i=0;i<block_size;i++)
        tei_so[h][i]=0.0;
      fprintf(outfile,"\n\tCCTransform: allocated the %s block of size %lu",moinfo->get_irr_labs(h),block_size);
    }
  }

  double value;
  size_t p,q,r,s,pq,rs,pqrs,irrep;
  int ilsti,nbuf,fi,index,elements;
  elements = 0;
  struct iwlbuf ERIIN;
  iwl_buf_init(&ERIIN,PSIF_SO_TEI,0.0,1,1);
    do{
      ilsti = ERIIN.lastbuf;
      nbuf  = ERIIN.inbuf;
      fi=0;
      for(index=0;index<nbuf;index++){
        // Compute the [pq] index for this pqrs combination
        p = abs(ERIIN.labels[fi]);
        q = ERIIN.labels[fi+1];
        r = ERIIN.labels[fi+2];
        s = ERIIN.labels[fi+3];
        value = ERIIN.values[index];
        irrep = indexing->get_tuple_irrep(p,q);
        pq = indexing->get_tuple_rel_index(p,q);
        rs = indexing->get_tuple_rel_index(r,s);
        pqrs = INDEX(pq,rs);
        tei_so[irrep][pqrs]=value;
        fi+=4;
        elements++;
      }
      if(!ilsti)
        iwl_buf_fetch(&ERIIN);
    } while(!ilsti);

  fprintf(outfile,"\n    CCTransform: read %d non-zero integrals", elements);
  iwl_buf_close(&ERIIN,1);

//   for(int h=0;h<moinfo->get_nirreps();h++){
//     char label[80];
//     sprintf(label,"tei_so_%d",h);
//     psio_write_entry(MRCC_SO_INTS,label,(char*)&tei_so[h][0],INDEX(indexing->get_pairpi(h)-1,indexing->get_pairpi(h)-1)+1*sizeof(double));
//   }
//
//   for(int h=0;h<moinfo->get_nirreps();h++){
//     if(indexing->get_pairpi(h)>0){
//       delete[] tei_so[h];
//     }
//   }
//   delete[] tei_so;
}

/*!
    \fn CCTransform::transform_tei_integrals()
 */
void CCTransform::transform_tei_so_integrals()
{
  double alpha  = 1.0;
  double beta   = 0.0;
  int nirreps   = moinfo->get_nirreps();
  int pq,rs;
  int p_abs,q_abs,r_abs,s_abs,pqrs;
  int kl,i_abs,j_abs,k_abs,l_abs;
  double **A;
  double **B;
  double **D;
  double **C;
  CCIndex* rsindx = blas->get_index("[s>=s]");
  CCIndex* ijindx = blas->get_index("[n>=n]");
  CCIndex* elemindx = blas->get_index("[s]");

  // First-half transform
  fprintf(outfile,"\n\tCCTransform: beginning first-half integral trasform");
  fflush(outfile);
  for(int h_rs=0;h_rs<nirreps;h_rs++){
    for(int h_p=0;h_p<nirreps;h_p++){
      int h_q = h_rs^h_p;
      if(h_p<h_q) continue;

      int rows_A = elemindx->get_pairpi(h_q);
      int cols_A = elemindx->get_pairpi(h_p);
      int rows_B = moinfo->get_mopi(h_p);
      int cols_B = elemindx->get_pairpi(h_q);
      int rows_D = moinfo->get_mopi(h_q);
      int cols_D = moinfo->get_mopi(h_p);

      if(rows_A*cols_A*rows_B*cols_B*rows_D*cols_D>0){
        allocate2(double,A,rows_A,cols_A);
        allocate2(double,B,rows_B,cols_B);
        allocate2(double,D,rows_D,cols_D);
        for(int rs=0;rs<rsindx->get_pairpi(h_rs);rs++){
          zero_arr(&(A[0][0]),rows_A*cols_A);
          zero_arr(&(B[0][0]),rows_B*cols_B);
          zero_arr(&(D[0][0]),rows_D*cols_D);

          // Fill the A matrix
          for(int q=0;q<rows_A;q++){
            for(int p=0;p<cols_A;p++){
              p_abs = p + elemindx->get_first(h_p);
              q_abs = q + elemindx->get_first(h_q);
              p_abs >= q_abs ? (pq = rsindx->get_tuple_rel_index(p_abs,q_abs)) : (pq = rsindx->get_tuple_rel_index(q_abs,p_abs));
              pqrs = INDEX(pq,rs);
              A[q][p]=tei_so[h_rs][pqrs];
            }
          }

          // First transform
          C = moinfo->get_scf_mos(h_p);
          // B(i,q)=C(p,i)*A(q,p)
#ifdef CCTRANSFORM_USE_BLAS
          if(rows_B*cols_B*cols_A!=0)
            C_DGEMM_12(rows_B,cols_B,cols_A,alpha,&(C[0][0]),cols_A,&(A[0][0]),cols_A,beta,&(B[0][0]),cols_B);
#else
          for(int p=0;p<cols_A;p++)
            for(int q=0;q<rows_A;q++)
              for(int i=0;i<rows_B;i++)
                B[i][q]+=C[p][i]*A[q][p];
#endif

          // Second transform
          C = moinfo->get_scf_mos(h_q);
          // D(j,i)+=C(q,j)*B(i,q);
#ifdef CCTRANSFORM_USE_BLAS
          if(rows_D*cols_D*cols_B!=0)
            C_DGEMM_12(rows_D,cols_D,cols_B,alpha,&(C[0][0]),cols_B,&(B[0][0]),cols_B,beta,&(D[0][0]),cols_D);
#else
          for(int q=0;q<cols_B;q++)
            for(int i=0;i<rows_B;i++)
              for(int j=0;j<rows_D;j++)
                D[j][i]+=C[q][j]*B[i][q];
#endif
          // Store the half-transformed integrals
          for(int i=0;i<moinfo->get_mopi(h_p);i++){
            for(int j=0;j<moinfo->get_mopi(h_q);j++){
              i_abs = i + elemindx->get_first(h_p);
              j_abs = j + elemindx->get_first(h_q);
              if(i_abs >= j_abs){
                int ij = ijindx->get_tuple_rel_index(i_abs,j_abs);
                tei_half_transformed[h_rs][ij][rs]=D[j][i];
              }
            }
          }
        }
        release2(A);
        release2(B);
        release2(D);
      }
    }
  }

  // Second-half transform
  fprintf(outfile,"\n\tCCTransform: beginning second-half integral trasform");
  fflush(outfile);
  for(int h_ij=0;h_ij<nirreps;h_ij++){
    for(int h_r=0;h_r<nirreps;h_r++){
      int h_s = h_ij^h_r;
      if(h_r < h_s) continue;

      int rows_A = elemindx->get_pairpi(h_s);
      int cols_A = elemindx->get_pairpi(h_r);
      int rows_B = moinfo->get_mopi(h_r);
      int cols_B = elemindx->get_pairpi(h_s);
      int rows_D = moinfo->get_mopi(h_s);
      int cols_D = moinfo->get_mopi(h_r);

      if(rows_A*cols_A*rows_B*cols_B*rows_D*cols_D>0){
        allocate2(double,A,rows_A,cols_A);
        allocate2(double,B,rows_B,cols_B);
        allocate2(double,D,rows_D,cols_D);
        for(int ij=0;ij<ijindx->get_pairpi(h_ij);ij++){
          zero_arr(&(A[0][0]),rows_A*cols_A);
          zero_arr(&(B[0][0]),rows_B*cols_B);
          zero_arr(&(D[0][0]),rows_D*cols_D);
          // Fill the A matrix
          for(int r=0;r<elemindx->get_pairpi(h_r);r++){
            for(int s=0;s<elemindx->get_pairpi(h_s);s++){
              r_abs = r + elemindx->get_first(h_r);
              s_abs = s + elemindx->get_first(h_s);
              r_abs >= s_abs ? (rs = rsindx->get_tuple_rel_index(r_abs,s_abs)) : (rs = rsindx->get_tuple_rel_index(s_abs,r_abs));
              A[s][r]=tei_half_transformed[h_ij][ij][rs];
            }
          }

          // First transform
          C = moinfo->get_scf_mos(h_r);
          // B(k,s)=C(r,k)*A(s,r)
#ifdef CCTRANSFORM_USE_BLAS
          if(rows_B*cols_B*cols_A!=0)
            C_DGEMM_12(rows_B,cols_B,cols_A,alpha,&(C[0][0]),cols_A,&(A[0][0]),cols_A,beta,&(B[0][0]),cols_B);
#else
          for(int r=0;r<cols_A;r++)
            for(int s=0;s<rows_A;s++)
              for(int k=0;k<rows_B;k++)
                B[k][s]+=C[r][k]*A[s][r];
#endif

          // Second transform
          C = moinfo->get_scf_mos(h_s);
          // D(l,k)+=C(s,l)*B(k,s);
#ifdef CCTRANSFORM_USE_BLAS
          if(rows_D*cols_D*cols_B!=0)
            C_DGEMM_12(rows_D,cols_D,cols_B,alpha,&(C[0][0]),cols_B,&(B[0][0]),cols_B,beta,&(D[0][0]),cols_D);
#else
          for(int s=0;s<cols_B;s++)
            for(int k=0;k<rows_B;k++)
              for(int l=0;l<rows_D;l++)
                D[l][k]+=C[s][l]*B[k][s];
#endif

          // Store the half-transformed integrals
          for(int k=0;k<moinfo->get_mopi(h_r);k++){
            for(int l=0;l<moinfo->get_mopi(h_s);l++){
              k_abs = k + elemindx->get_first(h_r);
              l_abs = l + elemindx->get_first(h_s);
              if(k_abs >= l_abs){
                kl = ijindx->get_tuple_rel_index(k_abs,l_abs);
                tei_half_transformed[h_ij][ij][kl]=D[l][k];
              }
            }
          }
        }
        release2(A);
        release2(B);
        release2(D);
      }
    }
  }
//   for(int h_ij=0;h_ij<nirreps;h_ij++){
//     for(int ij=ijindx->get_first(h_ij);ij<ijindx->get_last(h_ij);ij++){
//       ij_tuple = ijindx->get_tuple(ij);
//       int ij_abs = ij - ijindx->get_first(h_ij);
//       for(int kl=ijindx->get_first(h_ij);kl<=ij;kl++){
//         kl_tuple = ijindx->get_tuple(kl);
//         int kl_abs = kl - ijindx->get_first(h_ij);
//         fprintf(outfile,"\n (%2d %2d|%2d %2d) = %15.10f",ij_tuple[0],ij_tuple[1],kl_tuple[0],kl_tuple[1],tei_half_transformed[h_ij][ij_abs][kl_abs]);
//       }
//     }
//   }
  fprintf(outfile,"\n\tCCTransform: end of integral transform");
  fflush(outfile);
}

/**
 * Read the two electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in the packed array tei_mo
 */
void CCTransform::read_tei_mo_integrals()
{
  // Read all the (frozen + non-frozen) TEI in Pitzer order
  // and store them in a in-core block-matrix
  CCIndex* mo_indexing = blas->get_index("[n>=n]");

  allocate_tei_mo();

  double value;
  size_t p,q,r,s,pq,rs,pqrs,irrep;
  size_t ilsti,nbuf,fi,index,elements;
  elements = 0;
  struct iwlbuf ERIIN;
  iwl_buf_init(&ERIIN,PSIF_MO_TEI,0.0,1,1);
    do{
      ilsti = ERIIN.lastbuf;
      nbuf  = ERIIN.inbuf;
      fi    = 0;
      for(index = 0; index < nbuf; index++){
        // Compute the [pq] index for this pqrs combination
        p = abs(ERIIN.labels[fi]);
        q = ERIIN.labels[fi+1];
        r = ERIIN.labels[fi+2];
        s = ERIIN.labels[fi+3];
        value = ERIIN.values[index];
        irrep = mo_indexing->get_tuple_irrep(p,q);
        pq    = mo_indexing->get_tuple_rel_index(p,q);
        rs    = mo_indexing->get_tuple_rel_index(r,s);
        pqrs  = INDEX(pq,rs);
        tei_mo[irrep][pqrs]=value;
        fi += 4;
        elements++;
      }
      if(!ilsti)
        iwl_buf_fetch(&ERIIN);
    } while(!ilsti);
  fprintf(outfile,"\n    CCTransform: read %lu non-zero integrals", elements);
  fflush(outfile);
  iwl_buf_close(&ERIIN,1);
}

/**
 * Read the one electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in oei_mo
 */
void CCTransform::read_oei_mo_integrals()
{
  // Read all the (frozen + non-frozen) OEI in Pitzer order
  allocate_oei_mo();

  int nmo = moinfo->get_nmo();

  double* H;
  allocate1(double,H,INDEX(nmo-1,nmo-1) + 1);

  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_MO_FZC),H,nmo*(nmo+1)/2,0,0,outfile);
//   else
//     iwl_rdone(PSIF_OEI,PSIF_MO_FZC,H,norbs*(norbs+1)/2,0,1,outfile); //TODO fix it!

  for(int i=0; i < nmo; i++)
    for(int j=0; j < nmo; j++)
      oei_mo[i][j] = H[INDEX(i,j)];
  release1(H);
}

/**
 * Free all the memory allocated by CCTransform
 */
void CCTransform::free_memory()
{
  free_oei_so();
  free_oei_mo();
  free_tei_so();
  free_tei_mo();
  free_tei_half_transformed();
  integral_map.clear();
}

/**
 * Allocate the oei_mo array
 */
void CCTransform::allocate_oei_mo()
{
  if(oei_mo==NULL){
    int nmo = moinfo->get_nmo();
    allocate2(double,oei_mo,nmo,nmo);
  }
}

/**
 * Free the oei_mo array
 */
void CCTransform::free_oei_mo()
{
  if(oei_mo!=NULL){
    int nmo = moinfo->get_nmo();
    release2(oei_mo);
    oei_mo = NULL;
  }
}

/**
 * Allocate the oei_so array
 */
void CCTransform::allocate_oei_so()
{
  if(oei_so==NULL){
    int nso = moinfo->get_nso();
    allocate2(double,oei_so,nso,nso);
  }
  if(s_so==NULL){
    int nso = moinfo->get_nso();
    allocate2(double,s_so,nso,nso);
  }
}

/**
 * Free the oei_so array
 */
void CCTransform::free_oei_so()
{
  if(oei_so!=NULL){
    release2(oei_so);
    oei_so = NULL;
  }
  if(s_so!=NULL){
    release2(s_so);
    s_so = NULL;
  }
}

/**
 * Allocate the tei_mo array and exit(EXIT_FAILURE) if there is not enough space
 */
void CCTransform::allocate_tei_mo()
{
  if(tei_mo==NULL){
    CCIndex* indexing = blas->get_index("[n>=n]");

    // Allocate the tei_mo matrix blocks
    allocate1(double*,tei_mo,moinfo->get_nirreps());

    bool failed = false;
    size_t required_size = 0;
    size_t matrix_size = 0;
    for(int h=0;h<moinfo->get_nirreps();h++){
      if(indexing->get_pairpi(h)>0){
        size_t block_size = INDEX(indexing->get_pairpi(h)-1,indexing->get_pairpi(h)-1)+1;
        matrix_size += block_size;
        if(sizeof(double) * block_size < memory_manager->get_FreeMemory()){
          allocate1(double,tei_mo[h],block_size);
          for(size_t i=0;i<block_size;i++)
            tei_mo[h][i]=0.0;
        }else{
          failed = true;
          required_size += sizeof(double) * block_size;
          tei_mo[h] = NULL;
        }
        fprintf(outfile,"\n\tCCTransform: allocated the %s block of size %lu bytes (free memory = %14lu bytes)",moinfo->get_irr_labs(h),block_size,memory_manager->get_FreeMemory());
      }
    }
    if(failed){
      fprintf(outfile,"\n\tCCTransform: not enough memory! %lu bytes extra required",required_size);
      fflush(outfile);
      exit(EXIT_FAILURE);
    }
  }
}

void CCTransform::allocate_tei_so()
{
  if(tei_so==NULL){
    CCIndex* indexing = blas->get_index("[s>=s]");

    // Allocate the tei_so matrix blocks
    allocate1(double*,tei_so,moinfo->get_nirreps());

    bool failed = false;
    size_t required_size = 0;
    size_t matrix_size = 0;
    for(int h=0;h<moinfo->get_nirreps();h++){
      if(indexing->get_pairpi(h)>0){
        int block_size = INDEX(indexing->get_pairpi(h)-1,indexing->get_pairpi(h)-1)+1;
        matrix_size += block_size;
        if(sizeof(double) * block_size < memory_manager->get_FreeMemory()){
          allocate1(double,tei_so[h],block_size);
          for(int i=0;i<block_size;i++)
            tei_so[h][i]=0.0;
        }else{
          failed = true;
          required_size += sizeof(double) * block_size;
          tei_so[h] = NULL;
        }
        fprintf(outfile,"\n\tCCTransform: allocated the %s block of size %d bytes (free memory = %14lu bytes)",moinfo->get_irr_labs(h),block_size,memory_manager->get_FreeMemory());
      }
    }
    if(failed){
      fprintf(outfile,"\n\tCCTransform: not enough memory!");
      fflush(outfile);
      exit(EXIT_FAILURE);
    }
  }
}

void CCTransform::allocate_tei_half_transformed()
{
  if(tei_half_transformed==NULL){
    CCIndex* so_indexing = blas->get_index("[s>=s]");
    CCIndex* mo_indexing = blas->get_index("[n>=n]");

    allocate1(double**,tei_half_transformed,moinfo->get_nirreps());

    bool failed = false;
    size_t required_size = 0;
    size_t matrix_size = 0;

    for(int h=0;h<moinfo->get_nirreps();h++){
      if(so_indexing->get_pairpi(h)*mo_indexing->get_pairpi(h)>0){
        allocate2(double,tei_half_transformed[h],mo_indexing->get_pairpi(h),so_indexing->get_pairpi(h));
        fprintf(outfile,"\n\tCCTransform: allocated the %s block of size %lu*%lu",moinfo->get_irr_labs(h),mo_indexing->get_pairpi(h),so_indexing->get_pairpi(h));
      }
    }
  }
}

/**
 * Free the tei_mo array
 */
void CCTransform::free_tei_mo()
{
  int nirreps = moinfo->get_nirreps();
  if(tei_mo!=NULL){
    size_t matrix_size = 0;
    CCIndex* indexing = blas->get_index("[n>=n]");
    for(int h=0;h<moinfo->get_nirreps();h++){
      if(indexing->get_pairpi(h)>0){
        size_t block_size = INDEX(indexing->get_pairpi(h)-1,indexing->get_pairpi(h)-1)+1;
        matrix_size += block_size;
        release1(tei_mo[h]);
        fprintf(outfile,"\n\tCCTransform: deallocated the %s block of size %lu",moinfo->get_irr_labs(h),block_size);
      }
    }
    release1(tei_mo);
    tei_mo=NULL;
  }
}

void CCTransform::free_tei_so()
{
  int nirreps = moinfo->get_nirreps();
  if(tei_so!=NULL){
    size_t matrix_size = 0;
    CCIndex* indexing = blas->get_index("[s>=s]");
    for(int h=0;h<moinfo->get_nirreps();h++){
      if(indexing->get_pairpi(h)>0){
        size_t block_size = INDEX(indexing->get_pairpi(h)-1,indexing->get_pairpi(h)-1)+1;
        matrix_size += block_size;
        release1(tei_so[h]);
        fprintf(outfile,"\n\tCCTransform: deallocated the %s block of size %lu",moinfo->get_irr_labs(h),block_size);
      }
    }
    release1(tei_so);
    tei_so = NULL;
  }
}

void CCTransform::free_tei_half_transformed()
{
  if(tei_half_transformed!=NULL){
    CCIndex* rsindx = blas->get_index("[s>=s]");
    CCIndex* ijindx = blas->get_index("[n>=n]");
    size_t matrix_size = 0;
    for(int h=0;h<moinfo->get_nirreps();h++){
      if(rsindx->get_pairpi(h)*ijindx->get_pairpi(h)>0){
        matrix_size += ijindx->get_pairpi(h)*rsindx->get_pairpi(h);
        release2(tei_half_transformed[h]);
        fprintf(outfile,"\n\tCCTransform: deallocated the %s block of size %lu*%lu",moinfo->get_irr_labs(h),ijindx->get_pairpi(h),rsindx->get_pairpi(h));
      }
    }
    release1(tei_half_transformed);
    tei_half_transformed = NULL;
  }
}

double CCTransform::oei(int p, int q)
{
  return(oei_mo[p][q]);
}

double CCTransform::tei(int p, int q, int r, int s)
{
  // Get the (pq|rs) integral
  return(tei_mo[tei_mo_indexing->get_tuple_irrep(MAX(p,q),MIN(p,q))][INDEX(tei_mo_indexing->get_tuple_rel_index(MAX(p,q),MIN(p,q)),tei_mo_indexing->get_tuple_rel_index(MAX(r,s),MIN(r,s)))]);
}



/*!
    \fn CCTransform::transform_oei_so_integrals()
 */
void CCTransform::transform_oei_so_integrals()
{
  fprintf(outfile,"\n  CCTransform: transforming one-electron integrals");
  fflush(outfile);

  allocate_oei_mo();

  int nso   = moinfo->get_nso();
  int nmo = moinfo->get_nmo();

  double** A;
  allocate2(double,A,nso,nmo);
  double** C = moinfo->get_scf_mos();

  // A(q,i) = H(q,p) * C(p,i)
/*#ifdef CCTRANSFORM_USE_BLAS
  C_DGEMM_12(
#else*/
  for(int q=0;q<nso;q++)
    for(int j=0;j<nmo;j++){
      A[q][j] = 0.0;
      for(int p=0;p<nso;p++)
        A[q][j] += oei_so[q][p] * C[p][j];
    }
  for(int i=0;i<nmo;i++)
    for(int j=0;j<nmo;j++){
      oei_mo[i][j] = 0.0;
      for(int q=0;q<nso;q++)
        oei_mo[i][j] += C[q][i] * A[q][j];
    }
// #endif

  release2(A);
}

}} /* End Namespaces */
