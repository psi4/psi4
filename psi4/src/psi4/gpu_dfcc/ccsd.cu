/*
 *@BEGIN LICENSE
 *
 * GPU-accelerated density-fitted coupled-cluster, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include"ccsd.h"
#include"blas.h"
#include<psi4/libmints/matrix.h>
#include<psi4/libmints/vector.h>
#include<psi4/libmints/molecule.h>
#include"gpuhelper.h"
#include<psi4/libmints/mintshelper.h>
#include<psi4/libciomr/libciomr.h>
#include<psi4/libqt/qt.h>
//#include<psi4/libparallel/process.h>
#include<omp.h>


#ifdef HAVE_MKL
    #include<mkl.h>
#else
    #define mkl_set_dynamic(a)
    #define mkl_set_num_threads(a)
    #define mkl_domain_set_num_threads(a,b)
#endif

#define NUMTHREADS 32
#define MAXBLOCKS 65535

__device__ int  GPUKernel_Position(int i,int j) {
  if (i<j){
    return j*(j+1)/2+i;
  }
  return i*(i+1)/2+j;
}
__global__ void GPUKernel_VpVm_tiled(int a, int bstart, int bsize,int v,double * in,double * outp,double * outm) {

    int blockid = blockIdx.x*gridDim.y + blockIdx.y;
    int id      = blockid*blockDim.x + threadIdx.x;

    int v2 = v*v;

    if ( id >= v2*bsize ) return;

    // id : b*v2+c*v+d

    int  d = id%v;
    int  c = (id-d)%(v*v)/v;

    if ( d > c ) return;

    //int  b = (id-d)%(v*bsize)/v;


    //int  c = (id-d-b*v)/(bsize*v);
    int  b = (id-d-c*v)/(v*v);

    if ( b + bstart < a ) return;

    int cd   = c*(c+1)/2 + d;
    int vtri = v*(v+1)/2;
    int bv2  = b*v2;

    //outp[b*vtri+cd] = in[bv2+d*v+c] + in[bv2+c*v+d];
    //outm[b*vtri+cd] = in[bv2+d*v+c] - in[bv2+c*v+d];
    outp[b*vtri+cd] = in[bv2+d*v+c] + in[id];
    outm[b*vtri+cd] = in[bv2+d*v+c] - in[id];
}

__global__ void GPUKernel_VpVm_v2(int a, int b,int v,double * in,double * outp,double * outm) {

    int blockid = blockIdx.x*gridDim.y + blockIdx.y;
    int id      = blockid*blockDim.x + threadIdx.x;

    int v2 = v*v;

    if ( id >= v2 ) return;

    int  d = id%v;
    int  c = (id-d)/v;

    if ( d > c ) return;

    int cd   = GPUKernel_Position(c,d);

    outp[cd] = in[d*v+c] + in[c*v+d];
    outm[cd] = in[d*v+c] - in[c*v+d];
}
__global__ void GPUKernel_VpVm(int a, int v,double * in,double * outp,double * outm) {

    int blockid = blockIdx.x*gridDim.y + blockIdx.y;
    int id      = blockid*blockDim.x + threadIdx.x;

    int v2 = v*v;

    if ( id >= v2*v ) return;

    int  d = id%v;
    int  b = (id-d)%(v2)/v;

    if ( b < a ) return;

    int bma = b - a;

    int  c = (id-d-b*v)/(v2);

    if ( d > c ) return;

    int cd   = GPUKernel_Position(c,d);
    int vtri = v*(v+1)/2;

    outp[bma*vtri+cd] = in[bma*v2+d*v+c] + in[bma*v2+c*v+d];
    outm[bma*vtri+cd] = in[bma*v2+d*v+c] - in[bma*v2+c*v+d];
}
__global__ void GPUKernel_Vm(int a, int v,double * in,double * out) {

    int blockid      = blockIdx.x*gridDim.y + blockIdx.y;
    int id      = blockid*blockDim.x + threadIdx.x;

    if ( id >= v*v*v ) return;

    int  d = id%v;
    int  b = (id-d)%(v*v)/v;
    int  c = (id-d-b*v)/(v*v);

    if ( b < a ) return;
    if ( d > c ) return;

    int cd   = GPUKernel_Position(c,d);
    int vtri = v*(v+1)/2;

    out[(b-a)*vtri+cd] = in[(b-a)*v*v+d*v+c] - in[(b-a)*v*v+c*v+d];
}
__global__ void GPUKernel_Vp(int a, int v,double * in,double * out) {

    int blockid      = blockIdx.x*gridDim.y + blockIdx.y;
    int id      = blockid*blockDim.x + threadIdx.x;

    if ( id >= v*v*v ) return;

    int  d = id%v;
    int  b = (id-d)%(v*v)/v;
    int  c = (id-d-b*v)/(v*v);

    if ( b < a ) return;
    if ( d > c ) return;

    int cd   = GPUKernel_Position(c,d);
    int vtri = v*(v+1)/2;

    out[(b-a)*vtri+cd] = in[(b-a)*v*v+d*v+c] + in[(b-a)*v*v+c*v+d];
}

using namespace psi;


namespace psi{namespace fnocc{

GPUDFCoupledCluster::GPUDFCoupledCluster(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options):
        DFCoupledCluster(reference_wavefunction,options)
{
    common_init();
}

GPUDFCoupledCluster::~GPUDFCoupledCluster()
{
}

// this is where we'll set up cuda/gpu stuff i suppose
void GPUDFCoupledCluster::common_init() {
    /**
      *  GPU helper class knows if we have gpus or not and how to use them.
      *  all gpu memory is allocated by the helper.  
      */
    helper_ = std::shared_ptr<GPUHelper>(new GPUHelper);

    // get device parameters, allocate gpu memory and pinned cpu memory
    helper_->ndoccact = ndoccact;
    helper_->nvirt    = nvirt;
    helper_->nmo      = nmo;
  
    helper_->CudaInit(options_);

    gpubuffer = helper_->gpubuffer;
    left      = helper_->gpumemory / 8.0;
    wasted    = helper_->extraroom / 8.0;
    num_gpus  = helper_->num_gpus;

    long int v = nvirt;
    ngputhreads=NUMTHREADS;
    num=1;
    if ((v*v*v)%ngputhreads==0)
       nblocks = (v*v*v)/ngputhreads;
    else
       nblocks = (v*v*v+ngputhreads-(v*v*v)%ngputhreads)/ngputhreads;
    if (nblocks>MAXBLOCKS){
       num = nblocks/MAXBLOCKS+1;
       nblocks = nblocks/num + 1;
    }

    ncputhreads = omp_get_max_threads();

    if (  options_.get_bool("DGEMM_TIMINGS")  ) {
        helper_->DGEMM_Timings();
    }
  
}

// accumulate results of contraction of (ac|bd) and t2
void GPUDFCoupledCluster::useVabcd1(){

  long int o = ndoccact;
  long int v = nvirt;
  long int oov = o*o*v;
  long int oo  = o*o;
  long int otri = o*(o+1)/2;
  long int vtri = v*(v+1)/2;

  std::shared_ptr<PSIO> psio(new PSIO());

  psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));

  // available gpu memory (in doubles)
  long int ndoubles = (left - wasted) - 2*otri*vtri;

  for (long int a = 0; a < v; a++) {

      // do we need to tile loop over b >= a?
      long int ntiles = 1;
      while ( ntiles < v-a )  {
          long int size = (v - a) / ntiles;
          if (size * ntiles < v - a) size++;
          long int max = (size*nQ*v+nQ*v > 2*size*vtri ? size*nQ*v + nQ*v : 2*size*vtri);;

          //if ( ndoubles >= max + 2*size*otri ) break;
          if ( ndoubles >= max + size*nQ*v ) break;
          ntiles++;
      }

      // tile dimensions
      long int * tilesize = (long int *)malloc(ntiles*sizeof(long int));
      for (long int tile = 0; tile < ntiles - 1; tile++) {
          tilesize[tile] = (v-a) / ntiles;
          if ( tilesize[tile] * ntiles < v - a) tilesize[tile]++;
      }
      tilesize[ntiles-1] = (v - a) - tilesize[0] * (ntiles - 1);

      //if (ntiles > 1) printf("%5i/%5i ntiles %5i\n",a,v,ntiles);fflush(stdout);

      for (long int tileb = 0; tileb < ntiles; tileb++) {

          long int bsize = tilesize[tileb];
          long int bstart = a + tileb*tilesize[0];
    
          // contribute to residual
          #pragma omp parallel for schedule (static)
          for (long int ij = 0; ij < o*o; ij++) {
              long int j = ij % o;
              long int i = ( ij - j ) / o;
              int sg     = ( i > j ) ? 1 : -1;
              for (long int b = bstart; b < bstart + bsize; b++) {
                  tempv[a*oo*v+b*oo+i*o+j]    +=    tempr[Position(i,j) * vtri + Position(a,b)]
                                               + sg*tempr[Position(i,j) * vtri + Position(a,b) + otri*vtri];
                  if (a!=b) {
                     tempv[b*oov+a*oo+i*o+j]  +=    tempr[Position(i,j) * vtri + Position(a,b)]
                                               - sg*tempr[Position(i,j) * vtri + Position(a,b) + otri*vtri];
                  }
              }
          }
//gohere
      }
      free(tilesize);

  }

  // contribute to residual
  psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_DCC_R2,1);
}

void GPUDFCoupledCluster::Vabcd1(){
  long int o = ndoccact;
  long int v = nvirt;
  long int oov = o*o*v;
  long int oo  = o*o;
  long int otri = o*(o+1)/2;
  long int vtri = v*(v+1)/2;

  std::shared_ptr<PSIO> psio(new PSIO());

  #pragma omp parallel for schedule (static) num_threads(num_gpus)
  for (long int i=0; i<o; i++){
      for (long int j=i; j<o; j++){
          long int ij = Position(i,j);
          for (long int a=0; a<v; a++){
              for (long int b=a; b<v; b++){
                  tempr[ij*vtri+Position(a,b)] =
                     (tb[a*oov+b*oo+i*o+j]+tb[b*oov+a*oo+i*o+j]);
                  tempr[ij*vtri+Position(a,b)+vtri*otri] =
                     (tb[a*oov+b*oo+i*o+j]-tb[b*oov+a*oo+i*o+j]);
              }
              tempr[ij*vtri+Position(a,a)] = tb[a*oov+a*oo+i*o+j];
          }
      }
  }

  if ( v > nQ ) {
      throw PsiException("GPU DFCC will break if Nv > Naux",__FILE__,__LINE__);
  }

  // available gpu memory (in doubles)
  long int ndoubles = (left - wasted) - 2*otri*vtri;
  long int ntiles_ij = 1;

  // do we need to tile ij?
  if ( ndoubles < 0 ) {
      while ( ntiles_ij < otri ) {
          ntiles_ij++;
          long int size  = otri / ntiles_ij;
          if ( size * ntiles_ij < otri ) size++;
          if ( left - wasted - size * 2*vtri ) {
              ndoubles = (left - wasted) - size * 2*vtri;
              break;
          }
      }
      outfile->Printf("    <<< warning >>> tiling composite ij index (%5li tiles)\n",ntiles_ij);
  }
  // sizes of ij tiles:
  long int * tilesize_ij = (long int *)malloc(ntiles_ij*sizeof(long int));
  for (long int tile = 0; tile < ntiles_ij - 1; tile++) {
      tilesize_ij[tile] = otri / ntiles_ij;
      if ( tilesize_ij[tile] * ntiles_ij < otri ) tilesize_ij[tile]++;
  }
  tilesize_ij[ntiles_ij-1] = otri - tilesize_ij[0] * (ntiles_ij - 1);

  for (long int tile_ij = 0; tile_ij < ntiles_ij; tile_ij++) {

      // copy this tile of t2 to the gpus
      #pragma omp parallel for schedule (static) num_threads(num_gpus)
      for (int i = 0; i < num_gpus; i++) {
          int thread = omp_get_thread_num();
          cudaSetDevice(thread);
          double * gput2 = gpubuffer[thread];
          cudaMemcpy(gput2,                    tempr + tile_ij * tilesize_ij[0] * vtri,              sizeof(double) * tilesize_ij[tile_ij] * vtri,cudaMemcpyHostToDevice);
          cudaMemcpy(gput2+tilesize_ij[0]*vtri,tempr + tile_ij * tilesize_ij[0] * vtri + otri * vtri,sizeof(double) * tilesize_ij[tile_ij] * vtri,cudaMemcpyHostToDevice);
      }

      last_a = v;
      // parallelize over multiple gpus
      #pragma omp parallel for schedule (dynamic) num_threads(num_gpus)
      for (long int a = 0; a < v; a++) {

          if (cpudone && last_a == v) { last_a = a; }

          if (last_a == v) {

          cudaStream_t stream;
          cudaEvent_t estart,estop;
          cudaEventCreate(&estart);
          cudaEventCreate(&estop);
          int thread = omp_get_thread_num();
          cudaSetDevice(thread);
          double * gput2 = gpubuffer[thread];

          // do we need to tile loop over b >= a?
          long int ntiles = 1;
          while ( ntiles < v-a )  {
              long int size = (v - a) / ntiles;
              if (size * ntiles < v - a) size++;
              long int max = (size*nQ*v+nQ*v > 2*size*vtri ? size*nQ*v + nQ*v : 2*size*vtri);

              //if ( ndoubles >= max + 2*size*otri ) break;
              if ( ndoubles >= max + size*nQ*v ) break;
              ntiles++;
          }

          // tile dimensions
          long int * tilesize = (long int *)malloc(ntiles*sizeof(long int));
          for (long int tile = 0; tile < ntiles - 1; tile++) {
              tilesize[tile] = (v-a) / ntiles;
              if ( tilesize[tile] * ntiles < v - a) tilesize[tile]++;
          }
          tilesize[ntiles-1] = (v - a) - tilesize[0] * (ntiles - 1);

          if (ntiles > 1) outfile->Printf("%5i/%5i ntiles %5i tilesize %5i\n",a,v,ntiles,tilesize[0]);fflush(stdout);

          for (long int tileb = 0; tileb < ntiles; tileb++) {

              long int bsize = tilesize[tileb];
              long int bstart = a + tileb*tilesize[0];

              // shift other buffers by 2 * tilesize_ij * vtri
              long int shift = 2L * tilesize_ij[0] * vtri;

              double * gpuVcdb = gpubuffer[thread] + shift + (bsize*nQ*v + nQ*v > 2*bsize*vtri ? bsize*nQ*v + nQ*v : 2*bsize*vtri);
              double * gpuVm   = gpubuffer[thread] + shift;
              double * gpuVp   = gpubuffer[thread] + shift + bsize*vtri;
              double * gpuA    = gpubuffer[thread] + shift + 2*bsize*vtri;
              double * gpuIqd  = gpubuffer[thread] + shift;
              double * gpuIqc  = gpubuffer[thread] + shift + bsize*nQ*v;

              long int num     = 1;
              long int nblocks = ( bsize*v*v )/ NUMTHREADS;
              if ( (bsize*v*v) % NUMTHREADS != 0 ) {
                 nblocks = (bsize*v*v+NUMTHREADS-(bsize*v*v)%NUMTHREADS)/NUMTHREADS;
              }
              if (nblocks > MAXBLOCKS){
                 num     = nblocks / MAXBLOCKS + 1;
                 nblocks = nblocks / num + 1;
              }

              dim3 dimgrid (nblocks,num);

              stream = NULL;
        
              double start2 = omp_get_wtime();
              //cudaThreadSynchronize();
              //helper_->Check_CUDA_Error(outfile,"before anything. ");
              cudaEventRecord(estart,stream);
        
                  cudaMemcpyAsync(gpuIqc,Qvv+a*nQ*v,sizeof(double)*nQ*v,cudaMemcpyHostToDevice,stream);
                  //cudaThreadSynchronize();
                  //helper_->Check_CUDA_Error(outfile,"memcpy 1");
                  cudaMemcpyAsync(gpuIqd,Qvv+bstart*nQ*v,sizeof(double)*bsize*nQ*v,cudaMemcpyHostToDevice,stream);
                  //cudaThreadSynchronize();
                  //helper_->Check_CUDA_Error(outfile,"memcpy 2");
                  cublasDgemm('t','n',v,bsize*v,nQ,1.0,gpuIqc,nQ,gpuIqd,nQ,0.0,gpuVcdb,v);
                  //cudaThreadSynchronize();
                  //helper_->Check_CUDA_Error(outfile,"building v");
       
                  GPUKernel_VpVm_tiled<<<dimgrid,NUMTHREADS>>>(a,bstart,bsize,v,gpuVcdb,gpuVp,gpuVm);
                  //cudaThreadSynchronize();
                  //helper_->Check_CUDA_Error(outfile,"building v+/v-");

                  cublasDgemm('t','n',tilesize_ij[tile_ij],bsize,vtri,0.5,gput2,                    vtri,gpuVp,vtri,0.0,gpuA,                           tilesize_ij[tile_ij]);
                  cublasDgemm('t','n',tilesize_ij[tile_ij],bsize,vtri,0.5,gput2+tilesize_ij[0]*vtri,vtri,gpuVm,vtri,0.0,gpuA+bsize*tilesize_ij[tile_ij],tilesize_ij[tile_ij]);

                  cudaMemcpyAsync(tempr2[thread],gpuA,sizeof(double)*2*bsize*tilesize_ij[tile_ij],cudaMemcpyDeviceToHost,stream);

              cudaEventRecord(estop,stream);
        
              while( cudaEventQuery(estop) == cudaErrorNotReady );
              double end2 = omp_get_wtime();
              for (int ij = 0; ij < tilesize_ij[tile_ij]; ij++) {
                  for (int b = bstart; b < bstart + bsize; b++) {
                      tempr[(ij+tile_ij*tilesize_ij[0])*vtri + Position(a,b)]           = tempr2[thread][(b-bstart)*tilesize_ij[tile_ij]+ij];
                      tempr[(ij+tile_ij*tilesize_ij[0])*vtri + Position(a,b)+otri*vtri] = tempr2[thread][(b-bstart)*tilesize_ij[tile_ij]+ij+bsize*tilesize_ij[tile_ij]];
                  }
              }
//gohere
          }
          free(tilesize);

          }

      }
  }
  free(tilesize_ij);

}
void GPUDFCoupledCluster::FinishVabcd1(){
  long int o = ndoccact;
  long int v = nvirt;
  long int oov = o*o*v;
  long int oo  = o*o;
  long int otri = o*(o+1)/2;
  long int vtri = v*(v+1)/2;

  std::shared_ptr<PSIO> psio(new PSIO());

  // need to build t2+/- for CPU to use
  #pragma omp parallel for schedule (static) num_threads(num_gpus)
  for (long int i=0; i<o; i++){
      for (long int j=i; j<o; j++){
          long int ij = Position(i,j);
          for (long int a=0; a<v; a++){
              for (long int b=a; b<v; b++){
                  tempt[ij*vtri+Position(a,b)] =
                     (tb[a*oov+b*oo+i*o+j]+tb[b*oov+a*oo+i*o+j]);
                  tempt[ij*vtri+Position(a,b)+vtri*otri] =
                     (tb[a*oov+b*oo+i*o+j]-tb[b*oov+a*oo+i*o+j]);
              }
              tempt[ij*vtri+Position(a,a)] = tb[a*oov+a*oo+i*o+j];
          }
      }
  }

  // available gpu memory (in doubles)
  long int ndoubles = (left - wasted) - 2*otri*vtri;
  long int ntiles_ij = 1;


  // available cpu memory (in doubles)
  long int nQmax = nQ > nQ_scf ? nQ : nQ_scf;
  long int dim = 2L*v*v*v;
  if (2*nQmax*o*v>dim)   dim = 2*nQmax*o*v;
  if (o*o*v*v>dim)       dim = o*o*v*v;
  if (nQmax*v*v>dim)     dim = nQmax*v*v;
  if (nQmax*nso*nso>dim) dim = nQmax*nso*nso;

  // do we need to tile ij?
  if ( ndoubles < 0 ) {
      while ( ntiles_ij < otri ) {
          ntiles_ij++;
          long int size  = otri / ntiles_ij;
          if ( size * ntiles_ij < otri ) size++;
          if ( left - wasted - size * 2*vtri ) {
              ndoubles = (left - wasted) - size * 2*vtri;
              break;
          }
      }
      //outfile->Printf("    <<< warning >>> tiling composite ij index (%5li tiles)\n",ntiles_ij);
      //outfile->Printf("    <<< warning >>> tiling composite ij index (%5li tiles)\n",ntiles_ij);
      throw PsiException("  <<< warning >>> tiling composite ij index ... feature temporarily disabled",__FILE__,__LINE__);
  }
  // sizes of ij tiles:
  long int * tilesize_ij = (long int *)malloc(ntiles_ij*sizeof(long int));
  for (long int tile = 0; tile < ntiles_ij - 1; tile++) {
      tilesize_ij[tile] = otri / ntiles_ij;
      if ( tilesize_ij[tile] * ntiles_ij < otri ) tilesize_ij[tile]++;
  }
  tilesize_ij[ntiles_ij-1] = otri - tilesize_ij[0] * (ntiles_ij - 1);


  omp_set_nested(1);
  omp_set_dynamic(0);
  mkl_set_dynamic(0);
  int nthreads = omp_get_max_threads();
  for (long int tile_ij = 0; tile_ij < ntiles_ij; tile_ij++) {

      // copy this tile of t2 to the gpus (already there)

      // parallelize over multiple gpus
      #pragma omp parallel for schedule (dynamic) num_threads(num_gpus + 1)
      for (long int a = last_a; a < v; a++) {

          int thread = omp_get_thread_num();

          if ( thread < num_gpus ) {
              cudaStream_t stream;
              cudaEvent_t estart,estop;
              cudaEventCreate(&estart);
              cudaEventCreate(&estop);
              cudaSetDevice(thread);
              double * gput2 = gpubuffer[thread];

              // do we need to tile loop over b >= a?
              long int ntiles = 1;
              while ( ntiles < v-a )  {
                  long int size = (v - a) / ntiles;
                  if (size * ntiles < v - a) size++;
                  long int max = (size*nQ*v+nQ*v > 2*size*vtri ? size*nQ*v + nQ*v : 2*size*vtri);

                  //if ( ndoubles >= max + 2*size*otri ) break;
                  if ( ndoubles >= max + size*nQ*v ) break;
                  ntiles++;
              }

              // tile dimensions
              long int * tilesize = (long int *)malloc(ntiles*sizeof(long int));
              for (long int tile = 0; tile < ntiles - 1; tile++) {
                  tilesize[tile] = (v-a) / ntiles;
                  if ( tilesize[tile] * ntiles < v - a) tilesize[tile]++;
              }
              tilesize[ntiles-1] = (v - a) - tilesize[0] * (ntiles - 1);


              for (long int tileb = 0; tileb < ntiles; tileb++) {

                  long int bsize = tilesize[tileb];
                  long int bstart = a + tileb*tilesize[0];

                  // shift other buffers by 2 * tilesize_ij * vtri
                  long int shift = 2L * tilesize_ij[0] * vtri;

                  double * gpuVcdb = gpubuffer[thread] + shift + (bsize*nQ*v + nQ*v > 2*bsize*vtri ? bsize*nQ*v + nQ*v : 2*bsize*vtri);
                  double * gpuVm   = gpubuffer[thread] + shift;
                  double * gpuVp   = gpubuffer[thread] + shift + bsize*vtri;
                  double * gpuA    = gpubuffer[thread] + shift + 2*bsize*vtri;
                  double * gpuIqd  = gpubuffer[thread] + shift;
                  double * gpuIqc  = gpubuffer[thread] + shift + bsize*nQ*v;

                  long int num     = 1;
                  long int nblocks = ( bsize*v*v )/ NUMTHREADS;
                  if ( (bsize*v*v) % NUMTHREADS != 0 ) {
                     nblocks = (bsize*v*v+NUMTHREADS-(bsize*v*v)%NUMTHREADS)/NUMTHREADS;
                  }
                  if (nblocks > MAXBLOCKS){
                     num     = nblocks / MAXBLOCKS + 1;
                     nblocks = nblocks / num + 1;
                  }

                  dim3 dimgrid (nblocks,num);

                  stream = NULL;
        
                  double start2 = omp_get_wtime();
                  //cudaThreadSynchronize();
                  //helper_->Check_CUDA_Error(outfile,"before anything. ");
                  cudaEventRecord(estart,stream);
        
                      cudaMemcpyAsync(gpuIqc,Qvv+a*nQ*v,sizeof(double)*nQ*v,cudaMemcpyHostToDevice,stream);
                      //cudaThreadSynchronize();
                      //helper_->Check_CUDA_Error(outfile,"memcpy 1");
                      cudaMemcpyAsync(gpuIqd,Qvv+bstart*nQ*v,sizeof(double)*bsize*nQ*v,cudaMemcpyHostToDevice,stream);
                      //cudaThreadSynchronize();
                      //helper_->Check_CUDA_Error(outfile,"memcpy 2");
                      cublasDgemm('t','n',v,bsize*v,nQ,1.0,gpuIqc,nQ,gpuIqd,nQ,0.0,gpuVcdb,v);
                      //cudaThreadSynchronize();
                      //helper_->Check_CUDA_Error(outfile,"building v");
       
                      GPUKernel_VpVm_tiled<<<dimgrid,NUMTHREADS>>>(a,bstart,bsize,v,gpuVcdb,gpuVp,gpuVm);
                      //cudaThreadSynchronize();
                      //helper_->Check_CUDA_Error(outfile,"building v+/v-");

                      cublasDgemm('t','n',tilesize_ij[tile_ij],bsize,vtri,0.5,gput2,                    vtri,gpuVp,vtri,0.0,gpuA,                           tilesize_ij[tile_ij]);
                      cublasDgemm('t','n',tilesize_ij[tile_ij],bsize,vtri,0.5,gput2+tilesize_ij[0]*vtri,vtri,gpuVm,vtri,0.0,gpuA+bsize*tilesize_ij[tile_ij],tilesize_ij[tile_ij]);

                      cudaMemcpyAsync(tempr2[thread],gpuA,sizeof(double)*2*bsize*tilesize_ij[tile_ij],cudaMemcpyDeviceToHost,stream);

                  cudaEventRecord(estop,stream);
        
                  while( cudaEventQuery(estop) == cudaErrorNotReady );
                  double end2 = omp_get_wtime();
                  for (int ij = 0; ij < tilesize_ij[tile_ij]; ij++) {
                      for (int b = bstart; b < bstart + bsize; b++) {
                          tempr[(ij+tile_ij*tilesize_ij[0])*vtri + Position(a,b)]           = tempr2[thread][(b-bstart)*tilesize_ij[tile_ij]+ij];
                          tempr[(ij+tile_ij*tilesize_ij[0])*vtri + Position(a,b)+otri*vtri] = tempr2[thread][(b-bstart)*tilesize_ij[tile_ij]+ij+bsize*tilesize_ij[tile_ij]];
                      }
                  }
              }
              free(tilesize);

          }else {

              // cpu work

              mkl_set_num_threads(nthreads - num_gpus);

              // do we need to tile loop over b >= a?
              long int ntiles = 1;
/*
              while ( ntiles < v-a )  {
                  long int size = (v - a) / ntiles;
                  if (size * ntiles < v - a) size++;
                  long int max = (size*nQ*v+nQ*v > 2*size*vtri ? size*nQ*v + nQ*v : 2*size*vtri);

                  //if ( ndoubles >= max + 2*size*otri ) break;
                  if ( ndoubles_cpu >= max + size*nQ*v ) break;
                  ntiles++;
              }
*/

              // tile dimensions
              long int * tilesize = (long int *)malloc(ntiles*sizeof(long int));
              for (long int tile = 0; tile < ntiles - 1; tile++) {
                  tilesize[tile] = (v-a) / ntiles;
                  if ( tilesize[tile] * ntiles < v - a) tilesize[tile]++;
              }
              tilesize[ntiles-1] = (v - a) - tilesize[0] * (ntiles - 1);

              if (ntiles > 1) outfile->Printf("%5i/%5i ntiles %5i tilesize %5i (cpu) \n",a,v,ntiles,tilesize[0]);fflush(stdout);

              for (long int tileb = 0; tileb < ntiles; tileb++) {

                  long int bsize = tilesize[tileb];
                  long int bstart = a + tileb*tilesize[0];

                  // shift other buffers by 2 * tilesize_ij * vtri
                  long int shift = 0;//2L * tilesize_ij[0] * vtri;

                  double * gpuVm   = integrals + shift;
                  double * gpuVp   = integrals + shift + bsize*vtri;
                  double * gpuA    = integrals + shift + 2*bsize*vtri;
                  double * gpuVcdb = integrals + shift + 3*bsize*vtri;//(bsize*nQ*v + nQ*v > 2*bsize*vtri ? bsize*nQ*v + nQ*v : 2*bsize*vtri);
                  //double * gpuIqd  = integrals + shift;
                  //double * gpuIqc  = integrals + shift + bsize*nQ*v;

                  double start2 = omp_get_wtime();
        
                  //cudaMemcpyAsync(gpuIqc,Qvv+a*nQ*v,sizeof(double)*nQ*v,cudaMemcpyHostToDevice,stream);
                  //cudaMemcpyAsync(gpuIqd,Qvv+bstart*nQ*v,sizeof(double)*bsize*nQ*v,cudaMemcpyHostToDevice,stream);
                  F_DGEMM('t','n',v,bsize*v,nQ,1.0,Qvv+a*nQ*v,nQ,Qvv+bstart*nQ*v,nQ,0.0,gpuVcdb,v);
       
                  #pragma omp parallel for schedule (dynamic) num_threads(nthreads - num_gpus)
                  for (int d = 0; d < v; d++) {
                      for (int c = d; c < v; c++) {
                          int cd   = c*(c+1)/2 + d;
                          for (int b = bstart; b < v; b++) {
                              int id   = d + c*v + (b-bstart)*v*v;
                              int bv2  = (b-bstart)*v*v;
                               gpuVp[(b-bstart)*vtri+cd] = gpuVcdb[bv2+d*v+c] + gpuVcdb[id];
                               gpuVm[(b-bstart)*vtri+cd] = gpuVcdb[bv2+d*v+c] - gpuVcdb[id];
                          }
                      }
                  }

                  F_DGEMM('t','n',tilesize_ij[tile_ij],bsize,vtri,0.5,tempt,                    vtri,gpuVp,vtri,0.0,gpuA,                           tilesize_ij[tile_ij]);
                  F_DGEMM('t','n',tilesize_ij[tile_ij],bsize,vtri,0.5,tempt+tilesize_ij[0]*vtri,vtri,gpuVm,vtri,0.0,gpuA+bsize*tilesize_ij[tile_ij],tilesize_ij[tile_ij]);

                  //cudaMemcpyAsync(tempr2[thread],gpuA,sizeof(double)*2*bsize*tilesize_ij[tile_ij],cudaMemcpyDeviceToHost,stream);
        
                  #pragma omp parallel for schedule (dynamic) num_threads(nthreads - num_gpus)
                  for (int ij = 0; ij < tilesize_ij[tile_ij]; ij++) {
                      for (int b = bstart; b < bstart + bsize; b++) {
                          tempr[(ij+tile_ij*tilesize_ij[0])*vtri + Position(a,b)]           = gpuA[(b-bstart)*tilesize_ij[tile_ij]+ij];
                          tempr[(ij+tile_ij*tilesize_ij[0])*vtri + Position(a,b)+otri*vtri] = gpuA[(b-bstart)*tilesize_ij[tile_ij]+ij+bsize*tilesize_ij[tile_ij]];
                      }
                  }
              }
              free(tilesize);

          }
      }
  }
  free(tilesize_ij);

  omp_set_nested(0);
  omp_set_dynamic(1);
  mkl_set_dynamic(1);
  mkl_set_num_threads(nthreads);

}


void GPUDFCoupledCluster::CudaInit(){
  num_gpus = 0;

  cublasInit();
  helper_->Check_CUDA_Error(stdout,"cudaInit");
  struct cudaDeviceProp cudaProp;
  int gpu_id;

  // how many GPUs do we have?
  cudaGetDeviceCount(&num_gpus);
  helper_->Check_CUDA_Error(stdout,"cudaGetDeviceCount");

  if ( num_gpus == 0 ) { 
      throw PsiException("    Error: no cuda capable device detected.",__FILE__,__LINE__);
  }

  if (options_["NUM_GPUS"].has_changed()) {
      num_gpus = options_.get_int("NUM_GPUS");
  }

  cudaGetDevice(&gpu_id);
  helper_->Check_CUDA_Error(stdout,"cudaGetDevice");
  cudaGetDeviceProperties( &cudaProp,gpu_id );
  helper_->Check_CUDA_Error(stdout,"cudaGetDeviceProperties");
  outfile->Printf("\n");
  outfile->Printf("  _________________________________________________________\n");
  outfile->Printf("  CUDA device properties:\n");
  outfile->Printf("  name:                 %20s\n",cudaProp.name);
  outfile->Printf("  major version:        %20d\n",cudaProp.major);
  outfile->Printf("  minor version:        %20d\n",cudaProp.minor);
  outfile->Printf("  canMapHostMemory:     %20d\n",cudaProp.canMapHostMemory);
  outfile->Printf("  totalGlobalMem:       %20lu mb\n",cudaProp.totalGlobalMem/(1024*1024));
  outfile->Printf("  sharedMemPerBlock:    %20lu\n",cudaProp.sharedMemPerBlock);
  outfile->Printf("  clockRate:            %20.3f ghz\n",cudaProp.clockRate/1.0e6);
  outfile->Printf("  regsPerBlock:         %20d\n",cudaProp.regsPerBlock);
  outfile->Printf("  warpSize:             %20d\n",cudaProp.warpSize);
  outfile->Printf("  maxThreadsPerBlock:   %20d\n",cudaProp.maxThreadsPerBlock);
  outfile->Printf("  _________________________________________________________\n");
  outfile->Printf("\n");
  //fflush(outfile);

  // device memory left after some arrays (no, now total memory)
  int v = nvirt;
  left = cudaProp.totalGlobalMem/8.;// - 3*o*o*v*v - o*v-nmo*nmo;
  wasted = 200*1024*1024/8.; // leave an extra 200 mb on there.

  ngputhreads=NUMTHREADS;
  num=1;
  if ((v*v*v)%ngputhreads==0)
     nblocks = (v*v*v)/ngputhreads;
  else
     nblocks = (v*v*v+ngputhreads-(v*v*v)%ngputhreads)/ngputhreads;
  if (nblocks>MAXBLOCKS){
     num = nblocks/MAXBLOCKS+1;
     nblocks = nblocks/num + 1;
  }

  cudaDeviceReset();
  helper_->Check_CUDA_Error(stdout,"cudaDeviceReset");
}

void GPUDFCoupledCluster::CudaFinalize(){

  #pragma omp parallel for schedule (static) num_threads(num_gpus)
  for (int i=0; i<num_gpus; i++){ 
      int thread = omp_get_thread_num();
      cudaSetDevice(thread);
      cudaFree(gpubuffer[thread]);
  }
  cudaDeviceReset();
}
void GPUDFCoupledCluster::AllocateGPUMemory(){

  gpubuffer = (double**)malloc(num_gpus*sizeof(double*));
  #pragma omp parallel for schedule (static) num_threads(num_gpus)
  for (int i=0; i<num_gpus; i++){ 
      int thread = omp_get_thread_num();

      cudaSetDevice(thread);
      helper_->Check_CUDA_Error(stdout,"cudaSetDevice");

      cudaMalloc((void**)&gpubuffer[thread],sizeof(double)*(left-wasted));
      helper_->Check_CUDA_Error(stdout,"gpu memory");

  }

}

void GPUDFCoupledCluster::AllocateMemory() {

  if (nirrep_>1){
     throw PsiException("df_ccsd requires symmetry c1",__FILE__,__LINE__);
  }

  ischolesky_ = ( options_.get_str("DF_BASIS_CC") == "CHOLESKY" );
  nQ          = (int)Process::environment.globals["NAUX (CC)"];
  nQ_scf      = (int)Process::environment.globals["NAUX (SCF)"];

  int count=0;
  eps = (double*)malloc((ndoccact+nvirt)*sizeof(double));
  std::shared_ptr<Vector> eps_test = reference_wavefunction_->epsilon_a();
  for (int h=0; h<nirrep_; h++){
      for (int norb = frzcpi_[h]; norb<doccpi_[h]; norb++){
          eps[count++] = eps_test->get(h,norb);
      }
  }
  for (int h=0; h<nirrep_; h++){
      for (int norb = doccpi_[h]; norb<nmopi_[h]-frzvpi_[h]; norb++){
          eps[count++] = eps_test->get(h,norb);
      }
  }

  long int o = ndoccact;
  long int v = nvirt;

  /*========================================================
     ccsd memory requirements:
    
     tb:     o^2v^2
     tempt:  o^2v^2+ov ( actually o(o+1)v(v+1) + ov )
     tempv:  max (o^2v^2+ov , o*v*nQ)
     integrals: max(2v^3,nQ*nso^2, o^2v^2, 2v^3, 2nQ*o*v) (this is a minimum)
     Abij (SJS v^4 result): o(o+1)v/2
     Sbij (SJS v^4 result): o(o+1)v/2
     other stuff: 2ov+2v^2+(o+v)
    
     total: 3o^2v^2 + 2v^3  + o(o+1)v + 4ov  + 2v^2 + (o+v)  or 
            4o^2v^2         + o(o+1)v + 4ov  + 2v^2 + (o+v)  or
            3o^2v^2 + 2ovnQ + o(o+1)v + 4ov  + 2v^2 + (o+v)

     compare to the requirements for the (T) part:

            2o^2v^2 + 3v^3*nthreads + o^3v + ov
    
  ========================================================*/

  // reduce available memory by the amount required by the helper class
  memory -= helper_->max_mapped_memory;

  long int nQmax = nQ > nQ_scf ? nQ : nQ_scf;

  // for the df version, the dimension of the large buffer:
  long int dim = 2L*v*v*v;
  if (2*nQmax*o*v>dim)   dim = 2*nQmax*o*v;
  if (o*o*v*v>dim)    dim = o*o*v*v;
  if (nQmax*v*v>dim)     dim = nQmax*v*v;
  if (nQmax*nso*nso>dim) dim = nQmax*nso*nso;

  double total_memory = dim+(o*o*v*v+o*v)+(o*(o+1)*v*(v+1)+o*v)+o*o*v*v+2.*o*v+2.*v*v;
  long int max = nvirt*nvirt*nQmax > (nfzv+ndocc+nvirt)*ndocc*nQmax ? nvirt*nvirt*nQmax : (nfzv+ndocc+nvirt)*ndocc*nQmax;
  double df_memory    = nQmax*(o*o+o*v)+max;

  total_memory       *= 8./1024./1024.;
  df_memory          *= 8./1024./1024.;

  outfile->Printf("  Total memory requirements:       %9.2lf mb\n",df_memory+total_memory);
  outfile->Printf("  3-index integrals:               %9.2lf mb\n",df_memory);
  outfile->Printf("  CCSD intermediates:              %9.2lf mb\n",total_memory);
  outfile->Printf("\n");

  if (1.0 * memory / 1024. / 1024. < total_memory + df_memory) {
     outfile->Printf("\n");
     outfile->Printf("  error: not enough memory for ccsd.  increase available memory by %7.2lf mb\n",total_memory+df_memory-1.0*memory/1024./1024.);
     outfile->Printf("\n");
     //fflush(outfile);
     throw PsiException("not enough memory (ccsd).",__FILE__,__LINE__);
  }
  if (options_.get_bool("COMPUTE_TRIPLES")) {
      long int nthreads = omp_get_max_threads();
      double tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+3L*v*v*v*nthreads);
      if (tempmem > memory) {
          outfile->Printf("\n  <<< warning! >>> switched to low-memory (t) algorithm\n\n");
      }
      if (tempmem > memory || options_.get_bool("TRIPLES_LOW_MEMORY")){
         throw PsiException("low-memory triples option not yet implemented",__FILE__,__LINE__);
         //DPG commented out to remove unreachable warning
	 //isLowMemory = true;
         //tempmem = 8.*(2L*o*o*v*v+o*o*o*v+o*v+5L*o*o*o*nthreads);
      }
      outfile->Printf("  memory requirements for CCSD(T): %9.2lf mb\n\n",tempmem/1024./1024.);
  }
  cudaMallocHost((void**)&Qvv,nvirt*nvirt*nQ*sizeof(double));
  cudaThreadSynchronize();
  helper_->Check_CUDA_Error(stdout,"allocate host Qvv");

  tempr = (double*)malloc(o*(o+1)*v*(v+1)/2*sizeof(double));

  cudaThreadSynchronize();
  helper_->Check_CUDA_Error(stdout,"allocate host tempr");

  // o*(o+1)*v mapped memory for each gpu:
  // for now, give the choice of using helper's or allocating more.  TODO:
  // need to figure out a cleaner way to choose the memory we want to pin
  // and Qvv needs to be considerred as well.
  if ( o*(o+1)/v*sizeof(double) < helper_->max_mapped_memory_per_thread ) {
      tempr2 = helper_->tmp;
  }else {
      tempr2 = (double**)malloc(num_gpus*sizeof(double*));
      #pragma omp parallel for schedule (static) num_threads(num_gpus)
      for (long int i=0; i<num_gpus; i++){
          long int thread = 0;
          #ifdef _OPENMP
            thread = omp_get_thread_num();
          #endif
          cudaSetDevice(thread);
          helper_->Check_CUDA_Error(stdout,"cudaSetDevice");
          cudaMallocHost((void**)&tempr2[thread],o*(o+1)*v*sizeof(double));
          helper_->Check_CUDA_Error(stdout,"cpu tempr2");
      }
  }


  // allocate some memory for 3-index tensors
  Qoo = (double*)malloc(ndoccact*ndoccact*nQmax*sizeof(double));
  Qov = (double*)malloc(ndoccact*nvirt*nQmax*sizeof(double));

  long int                       tempvdim = o*o*v*v+o*v;
  if ( nQmax * o * v > tempvdim) tempvdim = nQmax * o * v;

  integrals = (double*)malloc(dim*sizeof(double));
  tempt     = (double*)malloc((o*(o+1)*v*(v+1)+o*v)*sizeof(double));
  tempv     = (double*)malloc(tempvdim*sizeof(double));
  Abij      = (double*)malloc(o*(o+1)/2*v*sizeof(double));
  Sbij      = (double*)malloc(o*(o+1)/2*v*sizeof(double));
  tb        = (double*)malloc(o*o*v*v*sizeof(double));
  w1        = (double*)malloc(o*v*sizeof(double));
  t1        = (double*)malloc(o*v*sizeof(double));
  I1        = (double*)malloc(v*v*sizeof(double));
  I1p       = (double*)malloc(v*v*sizeof(double));

  memset((void*)integrals,'\0',dim*sizeof(double));
  memset((void*)tempv,'\0',tempvdim*sizeof(double));
  memset((void*)tempt,'\0',(o*(o+1)*v*(v+1)+o*v)*sizeof(double));
  memset((void*)tempr,'\0',(o*(o+1)*v*(v+1)/2)*sizeof(double));
  memset((void*)tb,'\0',o*o*v*v*sizeof(double));
  memset((void*)w1,'\0',o*v*sizeof(double));
  memset((void*)t1,'\0',o*v*sizeof(double));
  memset((void*)I1,'\0',v*v*sizeof(double));
  memset((void*)I1p,'\0',v*v*sizeof(double));
  memset((void*)Abij,'\0',o*(o+1)/2*v*sizeof(double));
  memset((void*)Sbij,'\0',o*(o+1)/2*v*sizeof(double));

  // DIIS:
  diisvec    = (double*)malloc(sizeof(double)*(maxdiis+1));
  memset((void*)diisvec,'\0',(maxdiis+1)*sizeof(double));

  // new 3-index stuff for t1-transformed integrals:
  Fij   = (double*)malloc(o*o*sizeof(double));
  Fia   = (double*)malloc(o*v*sizeof(double));
  Fai   = (double*)malloc(o*v*sizeof(double));
  Fab   = (double*)malloc(v*v*sizeof(double));
  Ca_R  = (double*)malloc(nso*(nmo+nfzc+nfzv)*sizeof(double));
  Ca_L  = (double*)malloc(nso*(nmo+nfzc+nfzv)*sizeof(double));

  Ca = reference_wavefunction_->Ca()->pointer();

  // one-electron integrals
  std::shared_ptr<BasisSet> basis = reference_wavefunction_->basisset();
  std::shared_ptr<MintsHelper> mints(new MintsHelper(basis,options_));
  H = mints->so_kinetic();
  H->add(mints->so_potential());

}

// GPU kernels!
__global__ void GPUKernel_Iqdb(int a,int v,int nQ,double * in,double * out) {

    int blockid = blockIdx.x*gridDim.y + blockIdx.y;
    int id      = blockid*blockDim.x + threadIdx.x;

    if ( id >= v*v*nQ ) return;

    int  q = id%nQ;
    int  d = (id-q)%(nQ*v)/nQ;
    int  b = (id-q-d*nQ)/(nQ*v);

    if ( b < a ) return;

    int id2 = (b-a)*nQ*v+d*nQ+q;
    out[id2] = in[id];

}
typedef struct {
    int id;
    GPUDFCoupledCluster * cc;
} mega;

void *doit(void*x) {

  mega * m = (mega * )x;
  m->cc->pthreadCCResidual(m->id);

  return (NULL);
}

void GPUDFCoupledCluster::pthreadCCResidual(int id) {
    bool timer = options_.get_bool("CC_TIMINGS");
    long int o = ndoccact;
    long int v = nvirt;

        //////// start gpu section! ////////
        if (id==0)
        {
// AED
            omp_set_num_threads(num_gpus);
    
            int nthreads = omp_get_max_threads();

            #pragma omp parallel for schedule (static) num_threads(num_gpus)
            for (int i = 0 ; i < num_gpus; i++) {
                int mythread = omp_get_thread_num();
                cudaSetDevice(mythread);
            }
            double start = omp_get_wtime();

            Vabcd1();

            if (last_a == v) {
                gpudone = true;
                if (timer) {
                    outfile->Printf("        A2 =      t(c,d,i,j) (ac|bd)                                    %6.2lf\n",omp_get_wtime()-start);
                }
            }
            else 
                gpudone = false;
// AED
//gpudone = false;
    
        //////// end gpu section! ////////
        } 
        //////// start cpu section! ////////
        else
        {

// AED
///*
            int mythread = omp_get_thread_num();

            // pthread has NO idea what the right number of threads is ...
            int nthreads = ncputhreads;//omp_get_max_threads();

            if (nthreads > 1 + num_gpus) nthreads -= num_gpus;

            omp_set_num_threads(nthreads);
            mkl_set_num_threads(nthreads);
            mkl_domain_set_num_threads(nthreads, MKL_DOMAIN_BLAS);

            double start;

            // C2 = -1/2 t(bc,kj) [ (ki|ac) - 1/2 t(ad,li) (kd|lc) ] 
            //      +    t(bc,ki) [ (kj|ac) - 1/2 t(ad,lj) (kd|lc) ] 
            if (timer) start = omp_get_wtime();

            if (gpudone) helper_->GPUTiledDGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
            else F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            //printf("position 1 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int i = 0; i < o; i++) {
                    for (int l = 0; l < o; l++) {
                        for (int d = 0; d < v; d++) {
                            tempt[a*o*o*v+i*o*v+l*v+d] = tb[a*o*o*v+d*o*o+l*o+i];
                        }
                    }
                }
            }
            //printf("position 2 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int l = 0; l < o; l++) {
                for (int d = 0; d < v; d++) {
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            tempv[k*o*v*v+c*o*v+l*v+d] = integrals[k*o*v*v+d*o*v+l*v+c];
                        }
                    }
                }
            }
            // hang out until the gpu finishes ...
//            double wait = omp_get_wtime();
//            double accum = 0.0;
//            do {
//                if ( omp_get_wtime() - wait > 5.0 ) {
//                    accum += omp_get_wtime() - wait;
//                    wait = omp_get_wtime();
//                    outfile->Printf("gpu has taken an extra %6.2lf s\n",accum);
//                }
//            }while(!gpudone);

            //printf("position 3 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
//            if (gpudone) helper_->GPUTiledDGEMM('t','n',o*v,o*v,o*v,-0.5,tempv,o*v,tempt,o*v,0.0,integrals,o*v);
//            else         F_DGEMM('t','n',o*v,o*v,o*v,-0.5,tempv,o*v,tempt,o*v,0.0,integrals,o*v);
            long int gpuchunk = 0;
            long int odone    = 0;
            for (int i = 0; i < o; i++) {
                if (!gpudone) {
                    F_DGEMM('t','n',o*v,v,o*v,-0.5,tempv,o*v,tempt+i*o*v*v,o*v,0.0,integrals+i*o*v*v,o*v);
                }else {
                    gpuchunk = o - i;
                    odone    = i;
                    break;
                }
            }
            if (gpudone && gpuchunk > 0) {
                helper_->GPUTiledDGEMM('t','n',o*v,gpuchunk*v,o*v,-0.5,tempv,o*v,tempt+odone*o*v*v,o*v,0.0,integrals+odone*o*v*v,o*v);
            }

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            //printf("position 4 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
            if (gpudone) helper_->GPUTiledDGEMM('t','t',v*v,o*o,nQ,1.0,Qvv,nQ,Qoo,o*o,0.0,tempv,v*v);
            else         F_DGEMM('t','t',v*v,o*o,nQ,1.0,Qvv,nQ,Qoo,o*o,0.0,tempv,v*v);

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            //printf("position 5 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
            //F_DGEMM('n','t',v*v,o*o,nQ,1.0,Qvv,v*v,Qoo,o*o,0.0,tempv,v*v);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int i = 0; i < o; i++) {
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            integrals[a*o*o*v+i*o*v+k*v+c] += tempv[k*o*v*v+i*v*v+a*v+c];
                        }
                    }
                }
            }
            //printf("position 6 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int b = 0; b < v; b++) {
                for (int j = 0; j < o; j++) {
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            tempt[b*o*o*v+j*o*v+k*v+c] = tb[b*o*o*v+c*o*o+k*o+j];
                        }
                    }
                }
            }
            //printf("position 7 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
//            if (gpudone) helper_->GPUTiledDGEMM('t','n',o*v,o*v,o*v,-1.0,integrals,o*v,tempt,o*v,0.0,tempv,o*v);
//            else         F_DGEMM('t','n',o*v,o*v,o*v,-1.0,integrals,o*v,tempt,o*v,0.0,tempv,o*v);

            gpuchunk = 0;
            odone    = 0;
            for (int i = 0; i < o; i++) {
                if (!gpudone) {
                    F_DGEMM('t','n',o*v,v,o*v,-1.0,integrals,o*v,tempt+i*o*v*v,o*v,0.0,tempv+i*o*v*v,o*v);
                }else {
                    gpuchunk = o - i;
                    odone    = i;
                    break;
                }
            }
            if (gpudone && gpuchunk > 0) {
                helper_->GPUTiledDGEMM('t','n',o*v,gpuchunk*v,o*v,-1.0,integrals,o*v,tempt+odone*o*v*v,o*v,0.0,tempv+odone*o*v*v,o*v);
            }

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            //printf("position 8 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    for (int i = 0; i < o; i++) {
                        for (int j = 0; j < o; j++) {
                            tempt[a*o*o*v+b*o*o+i*o+j] = 0.5 * tempv[b*o*o*v+j*o*v+a*o+i] + tempv[b*o*o*v+i*o*v+a*o+j];
                        }
                    }
                }
            }
            //printf("position 9 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);

            // first contribution to residual 
            std::shared_ptr<PSIO> psio(new PSIO());
            psio->open(PSIF_DCC_R2,PSIO_OPEN_NEW);
            psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
            psio->close(PSIF_DCC_R2,1);
            //printf("position 10 %20.12lf\n",omp_get_wtime()-start);fflush(stdout);
            if (timer) {
                outfile->Printf("\n");
                outfile->Printf("        C2 = -1/2 t(b,c,k,j) [ (ki|ac) - 1/2 t(a,d,l,i) (kd|lc) ]\n");
                outfile->Printf("                + t(b,c,k,i) [ (kj|ac) - 1/2 t(a,d,l,j) (kd|lc) ]       %6.2lf\n",omp_get_wtime()-start);
                start = omp_get_wtime();
            }

            // now singles residual:

            // D1: F(ai)
            C_DCOPY(o*v,Fai,1,w1,1);

            // A1 (G):  U(c,d,k,l) (ad|kc)
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int d = 0; d < v; d++) {
                for (int i = 0; i < o; i++) {
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            tempt[d*o*o*v+i*o*v+k*v+c] = (2.0*tb[c*o*o*v+d*o*o+k*o+i] - tb[c*o*o*v+d*o*o+i*o+k]);
                        }
                    }
                }
            }
            if (gpudone) helper_->GPUTiledDGEMM('t','n',o*v,nQ,o*v,1.0,tempt,o*v,Qov,o*v,0.0,tempv,o*v);
            else F_DGEMM('t','n',o*v,nQ,o*v,1.0,tempt,o*v,Qov,o*v,0.0,tempv,o*v);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int q = 0; q < nQ; q++) {
                for (int a = 0; a < v; a++) {
                    for (int b = 0; b < v; b++) {
                        integrals[q*v*v+b*v+a] = Qvv[a*v*nQ+b*nQ+q];
                    }
                }
            }
            //#pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            //for (int q = 0; q < nQ; q++) {
            //    for (int a = 0; a < v; a++) {
            //        C_DCOPY(v,Qvv+q*v*v+a*v,1,integrals+q*v*v+a,v);
            //    }
            //}
            //if (gpudone) helper_->GPUTiledDGEMM('n','t',o,v,v*nQ,1.0,tempv,o,integrals,v,1.0,w1,o);
            //else F_DGEMM('n','t',o,v,v*nQ,1.0,tempv,o,integrals,v,1.0,w1,o);
            F_DGEMM('n','t',o,v,v*nQ,1.0,tempv,o,integrals,v,1.0,w1,o);

            if (timer) {
                outfile->Printf("        A1 =      U(c,d,k,l) (ad|kc)                                    %6.2lf\n",omp_get_wtime()-start);
                start = omp_get_wtime();
            }

            // B1 (H): -U(a,c,k,l) (ki|lc)
            if (gpudone) helper_->GPUTiledDGEMM('n','t',o*v,o*o,nQ,1.0,Qov,o*v,Qoo,o*o,0.0,integrals,o*v);
            else F_DGEMM('n','t',o*v,o*o,nQ,1.0,Qov,o*v,Qoo,o*o,0.0,integrals,o*v);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int i = 0; i < o; i++) {
                for (int c = 0; c < v; c++) {
                    for (int k = 0; k < o; k++) {
                        for (int l = 0; l < o; l++) {
                            tempv[i*o*o*v+c*o*o+k*o+l] = integrals[k*o*o*v+i*o*v+l*v+c];
                        }
                    }
                }
            }
            C_DCOPY(o*o*v*v,tb,1,tempt,1);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int c = 0; c < v; c++) {
                    for (int k = 0; k < o; k++) {
                        C_DAXPY(o,-0.5,tb+a*o*o*v+c*o*o+k,o,tempt+a*o*o*v+c*o*o+k*o,1);
                    }
                }
            }
            if (gpudone) helper_->GPUTiledDGEMM('t','n',o,v,o*o*v,-2.0,tempv,o*o*v,tempt,o*o*v,1.0,w1,o);
            else F_DGEMM('t','n',o,v,o*o*v,-2.0,tempv,o*o*v,tempt,o*o*v,1.0,w1,o);

            if (timer) {
                outfile->Printf("        B1 =    - U(a,c,k,l) (ki|lc)                                    %6.2lf\n",omp_get_wtime()-start);
                start = omp_get_wtime();
            }

            // C1
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int i = 0; i < o; i++) {
                    double dum = 0.0;
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            dum += Fia[k*v+c] * (2.0*tb[a*o*o*v+c*o*o+i*o+k] - tb[a*o*o*v+c*o*o+k*o+i]);
                        }
                    }
                    w1[a*o+i] += dum;
                }
            }

            if (timer) {
                outfile->Printf("        C1 =      F(k,c) U(a,c,i,k)                                     %6.2lf\n",omp_get_wtime()-start);
                start = omp_get_wtime();
            }

            // D2: 1/2 U(b,c,j,k) [ L(a,i,k,c) + 1/2 U(a,d,i,l) L(l,d,k,c) ] 
            if (gpudone) helper_->GPUTiledDGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
            else F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            C_DCOPY(o*o*v*v,integrals,1,tempv,1);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int l = 0; l < o; l++) {
                for (int d = 0; d < v; d++) {
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            tempv[l*o*v*v+d*o*v+k*v+c] -= 0.5 * integrals[l*o*v*v+c*o*v+k*v+d];
                        }
                    }
                }
            }
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int l = 0; l < o; l++) {
                for (int d = 0; d < v; d++) {
                    for (int a = 0; a < v; a++) {
                        for (int i = 0; i < o; i++) {
                            tempt[a*o*o*v+i*o*v+l*v+d] = 2.0 * tb[a*o*o*v+d*o*o+i*o+l]-tb[a*o*o*v+d*o*o+l*o+i];
                            //tempt[l*o*v*v+d*o*v+a*o+i] = 2.0 * tb[a*o*o*v+d*o*o+i*o+l]-tb[a*o*o*v+d*o*o+l*o+i];
                        }
                    }
                }
            }


            // hang out until the gpu finishes ...
//            double wait = omp_get_wtime();
//            double accum = 0.0;
//            do {
//                if ( omp_get_wtime() - wait > 5.0 ) {
//                    accum += omp_get_wtime() - wait;
//                    wait = omp_get_wtime();
//                    outfile->Printf("gpu has taken an extra %6.2lf s\n",accum);
//                }
//            }while(!gpudone);

            //if (gpudone) helper_->GPUTiledDGEMM('n','n',o*v,o*v,o*v,1.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);
            //else F_DGEMM('n','n',o*v,o*v,o*v,1.0,tempv,o*v,tempt,o*v,0.0,integrals,o*v);
            gpuchunk = 0;
            odone    = 0;
            for (int i = 0; i < o; i++) {
                if (!gpudone) {
                    F_DGEMM('n','n',o*v,v,o*v,1.0,tempv,o*v,tempt+i*o*v*v,o*v,0.0,integrals+i*o*v*v,o*v);
                }else {
                    gpuchunk = o - i;
                    odone    = i;
                    break;
                }
            }
            if (gpudone && gpuchunk > 0) {
                helper_->GPUTiledDGEMM('n','n',o*v,gpuchunk*v,o*v,1.0,tempv,o*v,tempt+odone*o*v*v,o*v,0.0,integrals+odone*o*v*v,o*v);
            }

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
            psio->read_entry(PSIF_DCC_QSO,"qvo",(char*)&tempv[0],nQ*o*v*sizeof(double));
            psio->close(PSIF_DCC_QSO,1);
            if (gpudone) helper_->GPUTiledDGEMM('n','t',o*v,o*v,nQ,2.0,Qov,o*v,tempv,o*v,1.0,integrals,o*v);
            else F_DGEMM('n','t',o*v,o*v,nQ,2.0,Qov,o*v,tempv,o*v,1.0,integrals,o*v);

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            if (gpudone) helper_->GPUTiledDGEMM('n','n',o*o,v*v,nQ,-1.0,Qoo,o*o,Qvv,nQ,0.0,tempv,o*o);
            else F_DGEMM('n','n',o*o,v*v,nQ,-1.0,Qoo,o*o,Qvv,nQ,0.0,tempv,o*o);

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            //F_DGEMM('n','t',o*o,v*v,nQ,-1.0,Qoo,o*o,Qvv,v*v,0.0,tempv,o*o);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int i = 0; i < o; i++) {
                    for (int k = 0; k < o; k++) {
                        for (int c = 0; c < v; c++) {
                            integrals[a*o*o*v+i*o*v+k*v+c] += tempv[a*o*o*v+c*o*o+k*o+i];
                        }
                    }
                }
            }
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int k = 0; k < o; k++) {
                for (int c = 0; c < v; c++) {
                    for (int b = 0; b < v; b++) {
                        for (int j = 0; j < o; j++) {
                            tempt[k*o*v*v+c*o*v+b*o+j] = 2.0 * tb[b*o*o*v+c*o*o+j*o+k] - tb[b*o*o*v+c*o*o+k*o+j];
                        }
                    }
                }
            }
            //if (gpudone) helper_->GPUTiledDGEMM('n','n',o*v,o*v,o*v,0.5,tempt,o*v,integrals,o*v,0.0,tempv,o*v);
            //else F_DGEMM('n','n',o*v,o*v,o*v,0.5,tempt,o*v,integrals,o*v,0.0,tempv,o*v);
            gpuchunk = 0;
            odone    = 0;
            for (int i = 0; i < o; i++) {
                if (!gpudone) {
                    F_DGEMM('n','n',o*v,v,o*v,0.5,tempt,o*v,integrals+i*o*v*v,o*v,0.0,tempv+i*o*v*v,o*v);
                }else {
                    gpuchunk = o - i;
                    odone    = i;
                    break;
                }
            }
            if (gpudone && gpuchunk > 0) {
                helper_->GPUTiledDGEMM('n','n',o*v,gpuchunk*v,o*v,0.5,tempt,o*v,integrals+odone*o*v*v,o*v,0.0,tempv+odone*o*v*v,o*v);
            }

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }

            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    for (int i = 0; i < o; i++) {
                        for (int j = 0; j < o; j++) {
                            tempt[a*o*o*v+b*o*o+i*o+j] = tempv[a*o*o*v+i*o*v+b*o+j];
                        }
                    }
                }
            }
            psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
            psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
            C_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
            psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
            psio->close(PSIF_DCC_R2,1);
            if (timer) {
                outfile->Printf("        D2 =  1/2 U(b,c,j,k) [ L(a,i,k,c) + 1/2 U(a,d,i,l) L(l,d,k,c) ] %6.2lf\n",omp_get_wtime()-start);
                start = omp_get_wtime();
            }

            if (gpudone) {
                omp_set_num_threads(ncputhreads);
                mkl_set_num_threads(ncputhreads);
                mkl_domain_set_num_threads(ncputhreads, MKL_DOMAIN_BLAS);
            }


            // E2 a: t(ac,ij) [ F(bc) - U(bd,kl) (ld|kc) ]
            C_DCOPY(o*o*v*v,tb,1,tempt,1);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int b = 0; b < v; b++) {
                for (int d = 0; d < v; d++) {
                    for (int k = 0; k < o; k++) {
                        C_DAXPY(o,-0.5,tb+b*o*o*v+d*o*o+k,o,tempt+b*o*o*v+d*o*o+k*o,1);
                    }
                }
            }
            if (gpudone) helper_->GPUTiledDGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
            else F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int c = 0; c < v; c++) {
                for (int d = 0; d < v; d++) {
                    for (int k = 0; k < o; k++) {
                        for (int l = 0; l < o; l++) {
                            tempv[c*o*o*v+d*o*o+k*o+l] = integrals[l*o*v*v+d*o*v+k*v+c];
                        }
                    }
                }
            }
            // overwriting Fab here, but it gets rebuilt every iteration anyway.
            if (gpudone) helper_->GPUTiledDGEMM('t','n',v,v,o*o*v,-2.0,tempv,o*o*v,tempt,o*o*v,1.0,Fab,v);
            else F_DGEMM('t','n',v,v,o*o*v,-2.0,tempv,o*o*v,tempt,o*o*v,1.0,Fab,v);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int c = 0; c < v; c++) {
                for (int a = 0; a < v; a++) {
                    for (int i = 0; i < o; i++) {
                        for (int j = 0; j < o; j++) {
                            tempt[c*o*o*v+a*o*o+i*o+j] = tb[a*o*o*v+c*o*o+i*o+j];
                        }
                    }
                }
            }
            if (gpudone) helper_->GPUTiledDGEMM('n','n',o*o*v,v,v,1.0,tempt,o*o*v,Fab,v,0.0,tempv,o*o*v);
            else F_DGEMM('n','n',o*o*v,v,v,1.0,tempt,o*o*v,Fab,v,0.0,tempv,o*o*v);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    for (int i = 0; i < o; i++) {
                        for (int j = 0; j < o; j++) {
                            tempt[a*o*o*v+b*o*o+i*o+j] = tempv[b*o*o*v+a*o*o+i*o+j];
                        }
                    }
                }
            }
            psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
            psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
            C_DAXPY(o*o*v*v,1.0,tempv,1,tempt,1);
            psio->write_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
            psio->close(PSIF_DCC_R2,1);
            if (timer) {
                outfile->Printf("        E2 =      t(a,c,i,j) [ F(b,c) - U(b,d,k,l) (ld|kc) ]            %6.2lf\n",omp_get_wtime()-start);
                start = omp_get_wtime();
            }

            // E2 b: -t(a,b,i,k) [ F(kj) - U(c,d,l,j) (kd|lc) ]
            // note that (kd|lc) should still be in integrals buffer
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int j = 0; j < o; j++) {
                for (int d = 0; d < v; d++) {
                    for (int l = 0; l < o; l++) {
                        for (int c = 0; c < v; c++) {
                            tempt[j*o*v*v+d*o*v+l*v+c] = (2.0 * tb[c*o*o*v+d*o*o+l*o+j] - tb[c*o*o*v+d*o*o+j*o+l] );
                        }
                    }
                }
            }
            // overwriting Fij here, but it gets rebuilt every iteration anyway.
            if (gpudone) helper_->GPUTiledDGEMM('t','n',o,o,o*v*v,1.0,tempt,o*v*v,integrals,o*v*v,1.0,Fij,o);
            else F_DGEMM('t','n',o,o,o*v*v,1.0,tempt,o*v*v,integrals,o*v*v,1.0,Fij,o);

            psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
            psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
            //if (gpudone) helper_->GPUTiledDGEMM('n','n',o,o*v*v,o,-1.0,Fij,o,tb,o,1.0,tempt,o);
            //else F_DGEMM('n','n',o,o*v*v,o,-1.0,Fij,o,tb,o,1.0,tempt,o);
            F_DGEMM('n','n',o,o*v*v,o,-1.0,Fij,o,tb,o,1.0,tempt,o);

            // R2 = R2 + P(ia,jb) R2
            C_DCOPY(o*o*v*v,tempt,1,integrals,1);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    for (int i = 0; i < o; i++) {
                        for (int j = 0; j < o; j++) {
                            integrals[a*o*o*v+b*o*o+i*o+j] += tempt[b*o*o*v+a*o*o+j*o+i];
                        }
                    }
                }
            }
            psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
            psio->close(PSIF_DCC_R2,1);
            if (timer) {
                outfile->Printf("                - t(a,b,i,k) [ F(k,j) - U(c,d,l,j) (kd|lc) ]            %6.2lf\n",omp_get_wtime()-start);
                start = omp_get_wtime();
            }

            // B2 = t(ab,kl) [ (ki|lj) + t(cd,ij) (kc|ld) ]
            if (gpudone) helper_->GPUTiledDGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
            else F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int k = 0; k < o; k++) {
                for (int l = 0; l < o; l++) {
                    for (int c = 0; c < v; c++) {
                        for (int d = 0; d < v; d++) {
                            tempv[k*o*v*v+l*v*v+c*v+d] = integrals[k*v*v*o+c*o*v+l*v+d];
                        }
                    }
                }
            }
            if (gpudone) helper_->GPUTiledDGEMM('n','t',o*o,o*o,nQ,1.0,Qoo,o*o,Qoo,o*o,0.0,integrals,o*o);
            else F_DGEMM('n','t',o*o,o*o,nQ,1.0,Qoo,o*o,Qoo,o*o,0.0,integrals,o*o);
            #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
            for (int k = 0; k < o; k++) {
                for (int i = 0; i < o; i++) {
                    for (int l = 0; l < o; l++) {
                        for (int j = 0; j < o; j++) {
                            tempt[k*o*o*o+l*o*o+i*o+j] = integrals[k*o*o*o+i*o*o+l*o+j];
                        }
                    }
                }
            }
            if (gpudone) helper_->GPUTiledDGEMM('n','n',o*o,o*o,v*v,1.0,tb,o*o,tempv,v*v,1.0,tempt,o*o);
            else F_DGEMM('n','n',o*o,o*o,v*v,1.0,tb,o*o,tempv,v*v,1.0,tempt,o*o);
            if (gpudone) helper_->GPUTiledDGEMM('n','n',o*o,v*v,o*o,1.0,tempt,o*o,tb,o*o,0.0,integrals,o*o);
            else F_DGEMM('n','n',o*o,v*v,o*o,1.0,tempt,o*o,tb,o*o,0.0,integrals,o*o);

            psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
            psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempt[0],o*o*v*v*sizeof(double));
            C_DAXPY(o*o*v*v,1.0,tempt,1,integrals,1);
            psio->write_entry(PSIF_DCC_R2,"residual",(char*)&integrals[0],o*o*v*v*sizeof(double));
            psio->close(PSIF_DCC_R2,1);

            if (timer) {
                outfile->Printf("        B2 =      t(a,b,k,l) [ (ki|lj) + t(c,d,i,j) (kc|ld) ]           %6.2lf\n",omp_get_wtime()-start);
                start = omp_get_wtime();
            }

            cpudone = true;
//*/
//AED
        }
        //////// end cpu section! ////////

}

void GPUDFCoupledCluster::CCResidual(){
    bool timer = options_.get_bool("CC_TIMINGS");
    long int v = nvirt;

    // test new transposed storage of Qvv

    // qvv transpose
    #pragma omp parallel for schedule (static)
    for (int q = 0; q < nQ; q++) {
        C_DCOPY(v*v,Qvv+q*v*v,1,integrals+q,nQ);
    }
    C_DCOPY(nQ*v*v,integrals,1,Qvv,1);


    int n = 2 < omp_get_max_threads() ? 2 : omp_get_max_threads();


    pthread_attr_t pthread_custom_attr;
    pthread_t * threads = (pthread_t *)malloc(n*sizeof(*threads));
    pthread_attr_init(&pthread_custom_attr);
    
    mega * m = (mega *)malloc(sizeof(mega)*n);

    gpudone     = false;
    cpudone     = false;

    double start = omp_get_wtime();
    for (int i = 0; i < n; i++) {
        m[i].id = i;
        m[i].cc = (&(*this));
        pthread_create(&threads[i], &pthread_custom_attr, doit, (void*)(m+i));
    }

    // Synchronize the completion of each thread.
    for (int i = 0; i < n; i++) {
    	pthread_join(threads[i],NULL);
    }
    free(threads);

    // it is possible the gpu didn't finish its work.  check and finish with 
    // the CPU and the GPU together
    if (cpudone && !gpudone) {
        //outfile->Printf("cpu finished first! %5i %5i\n",v,last_a);fflush(stdout);//exit(0);
        FinishVabcd1();
        if (timer) {
            outfile->Printf("        A2 =      t(c,d,i,j) (ac|bd)                                    %6.2lf\n",omp_get_wtime()-start);
        }
    }

    // use results of contraction of (ac|bd) and t2
    useVabcd1();

}

// t1-transformed 3-index fock matrix (using 3-index integrals from SCF)
void GPUDFCoupledCluster::T1Fock(){
    long int o = ndoccact;
    long int v = nvirt;
    long int full = o+v+nfzc+nfzv;

    // Ca_L = C(1-t1^T)
    // Ca_R = C(1+t1)
    double * Catemp = (double*)malloc(nso*full*sizeof(double));
    C_DCOPY(nso*full,&Ca[0][0],1,Ca_L,1);
    C_DCOPY(nso*full,&Ca[0][0],1,Ca_R,1);
    C_DCOPY(nso*full,&Ca[0][0],1,Catemp,1);
   

    #pragma omp parallel for schedule (static)
    for (int mu = 0; mu < nso; mu++) {
        for (int a = 0; a < v; a++) {
            double dum = 0.0;
            for (int i = 0; i < o; i++) {
                dum += Catemp[mu*full+i+nfzc] * t1[a*o+i];
            }
            Ca_L[mu*full + a + ndocc] -= dum;
        }
    }
    #pragma omp parallel for schedule (static)
    for (int mu = 0; mu < nso; mu++) {
        for (int i = 0; i < o; i++) {
            double dum = 0.0;
            for (int a = 0; a < v; a++) {
                dum += Catemp[mu*full+a+ndocc] * t1[a*o+i];
            }
            Ca_R[mu*full + i + nfzc] += dum;
        }
    }
    free(Catemp);

    // (Q|rs)
    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio_address addr1  = PSIO_ZERO;
    psio_address addr2  = PSIO_ZERO;

    long int nrows = 1;
    long int rowsize = nQ_scf;
    while ( rowsize*nso*nso > o*o*v*v ) {
        nrows++;
        rowsize = nQ_scf / nrows;
        if (nrows * rowsize < nQ_scf) rowsize++;
        if (rowsize == 1) break;
    }
    long int lastrowsize = nQ_scf - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso SCF",(char*)&integrals[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,Ca_L,full,integrals,nso,0.0,tempv,full);
        for (int q = 0; q < rowdims[row]; q++) {
            for (int mu = 0; mu < nso; mu++) {
                C_DCOPY(full,tempv+q*nso*full+mu*full,1,integrals+q*nso*full+mu,nso);
            }
        }
        F_DGEMM('n','n',full,full*rowdims[row],nso,1.0,Ca_R,full,integrals,nso,0.0,tempv,full);
        // full Qmo
        psio->write(PSIF_DCC_QSO,"Qmo SCF",(char*)&tempv[0],rowdims[row]*full*full*sizeof(double),addr2,&addr2);
    }
    delete rowdims;

    // build Fock matrix

    memset((void*)Fij,'\0',o*o*sizeof(double));
    memset((void*)Fia,'\0',o*v*sizeof(double));
    memset((void*)Fai,'\0',o*v*sizeof(double));
    memset((void*)Fab,'\0',v*v*sizeof(double));

    // transform H
    double ** hp = H->pointer();
    double * h = (double*)malloc(nmo*nmo*sizeof(double));
    for (int mu = 0; mu < nso; mu++) {
        for (int p = 0; p < nmo; p++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso; nu++) {
                dum += Ca_L[nu*full + p + nfzc] * hp[nu][mu];
            }
            integrals[p*nso+mu] = dum;
        }
    }
    for (int p = 0; p < nmo; p++) {
        for (int q = 0; q < nmo; q++) {
            double dum = 0.0;
            for (int nu = 0; nu < nso; nu++) {
                dum += Ca_R[nu*full+q+nfzc] * integrals[p*nso+nu];
            }
            h[p*nmo+q] = dum;
        }
    }

    double * temp3 = (double*)malloc(full*full*sizeof(double));

    memset((void*)temp3,'\0',full*full*sizeof(double));
    psio_address addr = PSIO_ZERO;

    nrows = 1;
    rowsize = nQ_scf;
    while ( rowsize*full*full > o*o*v*v ) {
        nrows++;
        rowsize = nQ_scf / nrows;
        if (nrows * rowsize < nQ_scf) rowsize++;
        if (rowsize == 1) break;
    }
    lastrowsize = nQ_scf - (nrows - 1L) * rowsize;
    rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qmo SCF",(char*)&integrals[0],rowdims[row]*full*full*sizeof(double),addr,&addr);
        for (int q = 0; q < rowdims[row]; q++) {
            // sum k (q|rk) (q|ks)
            F_DGEMM('n','n',full,full,ndocc,-1.0,integrals+q*full*full,full,integrals+q*full*full,full,1.0,temp3,full);

            // sum k (q|kk) (q|rs)
            double dum = 0.0;
            for (int k = 0; k < ndocc; k++) {
                dum += integrals[q*full*full+k*full + k];
            }
            C_DAXPY(full*full,2.0 * dum,integrals+q*full*full,1,temp3,1);
        }
    }
    delete rowdims;
    psio->close(PSIF_DCC_QSO,1);

    // Fij
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            Fij[i*o+j] = h[i*nmo+j] + temp3[(i+nfzc)*full+(j+nfzc)];
        }
    }

    // Fia
    for (int i = 0; i < o; i++) {
        for (int a = 0; a < v; a++) {
            Fia[i*v+a] = h[i*nmo+a+o] + temp3[(i+nfzc)*full+(a+ndocc)];
        }
    }

    // Fai
    for (int a = 0; a < v; a++) {
        for (int i = 0; i < o; i++) {
            Fai[a*o+i] = h[(a+o)*nmo+i] + temp3[(a+ndocc)*full+(i+nfzc)];
        }
    }

    // Fab
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            Fab[a*v+b] = h[(a+o)*nmo+b+o] + temp3[(a+ndocc)*full+(b+ndocc)];
        }
    }

    // replace eps
    for (int i = 0; i < o; i++) {
        eps[i] = Fij[i*o+i];
    }
    for (int a = 0; a < v; a++) {
        eps[a+o] = Fab[a*v+a];
    }
    free(h);
    free(temp3);

}

// t1-transformed 3-index integrals
void GPUDFCoupledCluster::T1Integrals(){
    long int o = ndoccact;
    long int v = nvirt;
    long int full = o+v+nfzc+nfzv;

    // Ca_L = C(1-t1^T)
    // Ca_R = C(1+t1)
    double * Catemp = (double*)malloc(nso*full*sizeof(double));
    C_DCOPY(nso*full,&Ca[0][0],1,Ca_L,1);
    C_DCOPY(nso*full,&Ca[0][0],1,Ca_R,1);
    C_DCOPY(nso*full,&Ca[0][0],1,Catemp,1);

    #pragma omp parallel for schedule (static)
    for (int mu = 0; mu < nso; mu++) {
        for (int a = 0; a < v; a++) {
            double dum = 0.0;
            for (int i = 0; i < o; i++) {
                dum += Catemp[mu*full+i+nfzc] * t1[a*o+i];
            }
            Ca_L[mu*full + a + ndocc] -= dum;
        }
    }
    #pragma omp parallel for schedule (static)
    for (int mu = 0; mu < nso; mu++) {
        for (int i = 0; i < o; i++) {
            double dum = 0.0;
            for (int a = 0; a < v; a++) {
                dum += Catemp[mu*full+a+ndocc] * t1[a*o+i];
            }
            Ca_R[mu*full + i + nfzc] += dum;
        }
    }
    free(Catemp);

    // (Q|rs)
    std::shared_ptr<PSIO> psio(new PSIO());
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio_address addr1  = PSIO_ZERO;
    psio_address addrvo = PSIO_ZERO;
    long int nrows = 1;
    long int rowsize = nQ;
    while ( rowsize*nso*nso > o*o*v*v ) {
        nrows++;
        rowsize = nQ / nrows;
        if (nrows * rowsize < nQ) rowsize++;
        if ( rowsize == 1 ) break;
    }
    long int lastrowsize = nQ - (nrows - 1L) * rowsize;
    long int * rowdims = new long int [nrows];
    for (int i = 0; i < nrows-1; i++) rowdims[i] = rowsize;
    rowdims[nrows-1] = lastrowsize;
    for (int row = 0; row < nrows; row++) {
        psio->read(PSIF_DCC_QSO,"Qso CC",(char*)&integrals[0],rowdims[row]*nso*nso*sizeof(double),addr1,&addr1);
        //helper_->GPUTiledDGEMM('n','n',full,nso*rowdims[row],nso,1.0,Ca_L,full,integrals,nso,0.0,tempv,full);
        F_DGEMM('n','n',full,nso*rowdims[row],nso,1.0,Ca_L,full,integrals,nso,0.0,tempv,full);
        for (int q = 0; q < rowdims[row]; q++) {
            for (int mu = 0; mu < nso; mu++) {
                C_DCOPY(full,tempv+q*nso*full+mu*full,1,integrals+q*nso*full+mu,full);
            }
        }
        //helper_->GPUTiledDGEMM('n','n',full,full*rowdims[row],nso,1.0,Ca_R,full,integrals,nso,0.0,tempv,full);
        F_DGEMM('n','n',full,full*rowdims[row],nso,1.0,Ca_R,full,integrals,nso,0.0,tempv,full);

        // Qoo
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    Qoo[(q+rowdims[0]*row)*o*o+i*o+j] = tempv[q*full*full+(i+nfzc)*full+(j+nfzc)];
                }
            }
        }
        // Qov
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int i = 0; i < o; i++) {
                for (int a = 0; a < v; a++) {
                    Qov[(q+rowdims[0]*row)*o*v+i*v+a] = tempv[q*full*full+(i+nfzc)*full+(a+ndocc)];
                }
            }
        }
        // Qvo
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int a = 0; a < v; a++) {
                for (int i = 0; i < o; i++) {
                    integrals[q*o*v+a*o+i] = tempv[q*full*full+(a+ndocc)*full+(i+nfzc)];
                }
            }
        }
        psio->write(PSIF_DCC_QSO,"qvo",(char*)&integrals[0],rowdims[row]*o*v*sizeof(double),addrvo,&addrvo);
        // Qvv
        #pragma omp parallel for schedule (static)
        for (int q = 0; q < rowdims[row]; q++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    Qvv[(q+rowdims[0]*row)*v*v+a*v+b] = tempv[q*full*full+(a+ndocc)*full+(b+ndocc)];
                }
            }
        }
    }
    delete rowdims;
    psio->close(PSIF_DCC_QSO,1);


    // check mp2 energy
    /*double * tints = (double*)malloc(o*o*v*v*sizeof(double));
    memset((void*)tints,'\0',o*o*v*v*sizeof(double));
    F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,tints,o*v);
    double e2 = 0.0;
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    double dijab = (eps[i] + eps[j] - eps[a+o] - eps[b+o]);
                    long int iajb = i*v*v*o + a*o*v + j*v + b;
                    long int jaib = j*v*v*o + a*v*o + i*v + b;
                    e2 += (2.0 * tints[iajb] - tints[jaib]) * tints[iajb] / dijab;
                }
            }
        }
    }
    printf("mp2 energy %20.12lf\n",e2);
    //exit(0);*/



    // check ccsd energy
  /*  F_DGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,integrals,o*v);

    if (t2_on_disk){
        std::shared_ptr<PSIO> psio (new PSIO());
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&tempv[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
        tb = tempv;
    }
    double energy = 0.0;
    double mp2energy = 0.0;
    for (long int a = 0; a < v; a++){
        for (long int b = 0; b < v; b++){
            for (long int i = 0; i < o; i++){
                for (long int j = 0; j < o; j++){
                    double dijab = (eps[i] + eps[j] - eps[a+o] - eps[b+o]);
                    long int ijab = a*v*o*o + b*o*o + i*o + j;
                    long int iajb = i*v*v*o + a*v*o + j*v + b;
                    long int jaib = j*v*v*o + a*v*o + i*v + b;
                    energy += (2.0*integrals[iajb]-integrals[jaib])*tb[ijab];
                    mp2energy += (2.0*integrals[iajb]-integrals[jaib])*integrals[iajb]/dijab;
                }
            }
        }
    }
    printf("ccsd energy %20.12lf\n",energy);
    printf("mp2 energy  %20.12lf\n",mp2energy);
*/
}

double GPUDFCoupledCluster::compute_energy() {
  PsiReturnType status = Success;

  //WriteBanner();
  AllocateMemory();
  status = CCSDIterations();

  // free some memory!
  free(Fij);
  free(Fab);
  free(Abij);
  free(Sbij);
  free(integrals);
  free(w1);
  free(I1);
  free(I1p);
  free(diisvec);
  free(tempt);
  free(tempv);

  // tstart in fnocc
  tstop();

  // mp2 energy
  Process::environment.globals["MP2 CORRELATION ENERGY"] = emp2;
  Process::environment.globals["MP2 TOTAL ENERGY"] = emp2 + escf;
  Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = emp2_os;
  Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = emp2_ss;

  // ccsd energy
  Process::environment.globals["CCSD CORRELATION ENERGY"] = eccsd;
  Process::environment.globals["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = eccsd_os;
  Process::environment.globals["CCSD SAME-SPIN CORRELATION ENERGY"] = eccsd_ss;
  Process::environment.globals["CCSD TOTAL ENERGY"] = eccsd + escf;
  Process::environment.globals["CURRENT ENERGY"] = eccsd + escf;

  if (options_.get_bool("COMPUTE_TRIPLES")){
      long int o = ndoccact;
      long int v = nvirt;

      if (!isLowMemory ) {
          // write (ov|vv) integrals, formerly E2abci, for (t)
          double *tempq = (double*)malloc(v*nQ*sizeof(double));
          // the buffer integrals was at least 2v^3, so these should definitely fit.
          double *Z     = (double*)malloc(v*v*v*sizeof(double));
          double *Z2    = (double*)malloc(v*v*v*sizeof(double));
          std::shared_ptr<PSIO> psio(new PSIO());
          psio->open(PSIF_DCC_ABCI,PSIO_OPEN_NEW);
          psio_address addr2 = PSIO_ZERO;
          for (long int i=0; i<o; i++){
              #pragma omp parallel for schedule (static)
              for (long int q=0; q<nQ; q++){
                  for (long int b=0; b<v; b++){
                      tempq[q*v+b] = Qov[q*o*v+i*v+b];
                  }
              }     
              helper_->GPUTiledDGEMM('n','t',v,v*v,nQ,1.0,tempq,v,Qvv,v*v,0.0,&Z[0],v);
              #pragma omp parallel for schedule (static)
              for (long int a=0; a<v; a++){
                  for (long int b=0; b<v; b++){
                      for (long int c=0; c<v; c++){
                          Z2[a*v*v+b*v+c] = Z[a*v*v+c*v+b];
                      }
                  }
              }
              psio->write(PSIF_DCC_ABCI,"E2abci",(char*)&Z2[0],v*v*v*sizeof(double),addr2,&addr2);
          }
          psio->close(PSIF_DCC_ABCI,1);
          free(tempq);
          free(Z);
          free(Z2);
      } else {
          psio_address addr = PSIO_ZERO;
          double * temp1 = (double*)malloc(( nQ*v > o*v*v ? nQ*v : o*v*v)*sizeof(double));
          double * temp2 = (double*)malloc(o*v*v*sizeof(double));
          std::shared_ptr<PSIO> psio(new PSIO());
          psio->open(PSIF_DCC_ABCI4,PSIO_OPEN_NEW);
          for (long int a = 0; a < v; a++) {
              #pragma omp parallel for schedule (static)
              for (long int q = 0; q < nQ; q++) {
                  for (long int c = 0; c < v; c++) {
                      temp1[q*v+c] = Qvv[q*v*v+a*v+c];
                  }
              }
              helper_->GPUTiledDGEMM('n','t',o*v,v,nQ,1.0,Qov,o*v,temp1,v,0.0,temp2,o*v);
              #pragma omp parallel for schedule (static)
              for (long int b = 0; b < v; b++) {
                  for (long int i = 0; i < o; i++) {
                      for (long int c = 0; c < v; c++) {
                          temp1[b*o*v+i*v+c] = temp2[c*o*v+i*v+b];
                      }
                  }
              }
              psio->write(PSIF_DCC_ABCI4,"E2abci4",(char*)&temp1[0],o*v*v*sizeof(double),addr,&addr);
          }
          psio->close(PSIF_DCC_ABCI4,1);
          free(temp1);
          free(temp2);
      }
      cudaFreeHost(Qvv);//free(Qvv);
      double * temp1 = (double*)malloc(o*o*v*v*sizeof(double));
      double * temp2 = (double*)malloc(o*o*v*v*sizeof(double));

      // write (oo|ov) integrals, formerly E2ijak, for (t)
      helper_->GPUTiledDGEMM('n','t',o*o,o*v,nQ,1.0,Qoo,o*o,Qov,o*v,0.0,temp1,o*o);
      for (int i=0; i<o; i++){
          for (int j=0; j<o; j++){
              for (int k=0; k<o; k++){
                  for (int a=0; a<v; a++){
                      temp2[j*o*o*v+i*o*v+k*v+a] = temp1[i*o*o*v+a*o*o+j*o+k];
                  }
              }
          }
      }
      std::shared_ptr<PSIO> psio(new PSIO());
      psio->open(PSIF_DCC_IJAK,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_IJAK,"E2ijak",(char*)&temp2[0],o*o*o*v*sizeof(double));
      psio->close(PSIF_DCC_IJAK,1);

      // df (ov|ov) integrals, formerly E2klcd
      helper_->GPUTiledDGEMM('n','t',o*v,o*v,nQ,1.0,Qov,o*v,Qov,o*v,0.0,temp1,o*v);
      psio->open(PSIF_DCC_IAJB,PSIO_OPEN_NEW);
      psio->write_entry(PSIF_DCC_IAJB,"E2iajb",(char*)&temp1[0],o*o*v*v*sizeof(double));
      psio->close(PSIF_DCC_IAJB,1);

      free(Qov);
      free(Qoo);
      free(temp1);
      free(temp2);

      // triples
      tstart();

      ccmethod = 0;
      if (isLowMemory)                           status = lowmemory_triples();
      else                                       status = triples();

      if (status == Failure){
         throw PsiException(
            "Whoops, the (T) correction died.",__FILE__,__LINE__);
      }
      tstop();

      // ccsd(t) energy
      Process::environment.globals["(T) CORRECTION ENERGY"] = et;
      Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = eccsd + et;
      Process::environment.globals["CCSD(T) TOTAL ENERGY"] = eccsd + et + escf;
      Process::environment.globals["CURRENT ENERGY"] = eccsd + et + escf;
  }else {
      free(Qoo);
      free(Qov);
      free(Qvv);
  }

  // free remaining memory
  free(Fia);
  free(Fai);
  free(t1);
  free(tb);

  return Process::environment.globals["CURRENT ENERGY"];
}

void GPUDFCoupledCluster::UpdateT2(){
    long int v = nvirt;
    long int o = ndoccact;
    long int rs = nmo;

    std::shared_ptr<PSIO> psio(new PSIO());

    // df (ai|bj)
    psio->open(PSIF_DCC_QSO,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_QSO,"qvo",(char*)&tempv[0],nQ*o*v*sizeof(double));
    psio->close(PSIF_DCC_QSO,1);
    helper_->GPUTiledDGEMM('n','t',o*v,o*v,nQ,1.0,tempv,o*v,tempv,o*v,0.0,integrals,o*v);

    // residual
    psio->open(PSIF_DCC_R2,PSIO_OPEN_OLD);
    psio->read_entry(PSIF_DCC_R2,"residual",(char*)&tempv[0],o*o*v*v*sizeof(double));
    psio->close(PSIF_DCC_R2,1);

    #pragma omp parallel for schedule (static)
    for (long int a=o; a<rs; a++){
        double da = eps[a];
        for (long int b=o; b<rs; b++){
            double dab = da + eps[b];
            for (long int i=0; i<o; i++){
                double dabi = dab - eps[i];
                for (long int j=0; j<o; j++){

                    long int iajb = (a-o)*v*o*o+i*v*o+(b-o)*o+j;
                    long int ijab = (a-o)*v*o*o+(b-o)*o*o+i*o+j;

                    double dijab = dabi-eps[j];
                    double tnew  = - (integrals[iajb] + tempv[ijab])/dijab;
                    //double tnew  = - (integrals[iajb])/dijab;
                    //tempt[ijab]  = tnew;
                    tempv[ijab]  = tnew;
                }
            }
        }
    }
    // error vector is just dt
    //C_DCOPY(o*o*v*v,tempt,1,tempv,1);

    if (t2_on_disk){
        psio->open(PSIF_DCC_T2,PSIO_OPEN_OLD);
        psio->read_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        C_DAXPY(o*o*v*v,1.0,tempv,1,integrals,1);
        psio->write_entry(PSIF_DCC_T2,"t2",(char*)&integrals[0],o*o*v*v*sizeof(double));
        psio->close(PSIF_DCC_T2,1);
    }else {
        C_DAXPY(o*o*v*v,1.0,tempv,1,tb,1);
    }
}



}}
