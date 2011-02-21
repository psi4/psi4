#include "3index.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;
using namespace psi;

namespace psi { 

DFTensor::DFTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> aux) :
    psio_(psio), primary_(primary), auxiliary_(aux), metric_(new FittingMetric(aux))
{
    common_init(); 
}
DFTensor::~DFTensor()
{
    psio_->close(PSIF_DF_TENSOR, 0);
}
void DFTensor::common_init()
{
    nao_ = primary_->nbf();
    naux_ = auxiliary_->nbf();
    nocc_ = 0;
    nvir_ = 0;
    psio_->open(PSIF_DF_TENSOR, PSIO_OPEN_NEW);
}
void DFTensor::form_Aia(bool do_all)
{
    double** Cop = Co_->pointer();
    double** Cvp = Cv_->pointer();

    //Storage required for each A
    unsigned long int storage_per_row = nao_*nao_+nocc_*nao_+nocc_*nvir_ + (do_all ? nocc_*nocc_ + nvir_*nvir_ : 0L);

    int maxPshell = 0;
    for (int P = 0; P<auxiliary_->nshell(); P++)
      if (maxPshell < auxiliary_->shell(P)->nfunction())
        maxPshell = auxiliary_->shell(P)->nfunction();

    //Max rows of A to pursue
    int max_rows = memory_/storage_per_row;
    if (max_rows > naux_)
        max_rows = naux_;
    if (max_rows < maxPshell )
        max_rows = maxPshell;

    //Block sizes
    int nblocks = 2*naux_/max_rows;
    int* p_starts = init_int_array(nblocks);
    int* p_sizes = init_int_array(nblocks);
    int* block_starts = init_int_array(nblocks);
    int* block_stops = init_int_array(nblocks);
    int* block_sizes = init_int_array(nblocks);

    //Determine block sizes
    nblocks = 0;
    int fun_counter = 0;
    for (int A=0; A<auxiliary_->nshell(); A++) {
      if (A == auxiliary_->nshell() - 1) {
        block_sizes[nblocks]++;
        p_sizes[nblocks] += auxiliary_->shell(A)->nfunction();
        block_stops[nblocks] = auxiliary_->nshell();
        nblocks++;
        break;
      }
      if (fun_counter+auxiliary_->shell(A)->nfunction() > max_rows) {
        block_stops[nblocks] = A;
        block_starts[nblocks+1] = A;
        block_sizes[nblocks+1] = 1;
        p_sizes[nblocks+1] = auxiliary_->shell(A)->nfunction();
        p_starts[nblocks+1] = auxiliary_->shell(A)->function_index();
        nblocks++;
        fun_counter = auxiliary_->shell(A)->nfunction();
        continue;
      }
      p_sizes[nblocks] += auxiliary_->shell(A)->nfunction();
      block_sizes[nblocks]++;
      fun_counter += auxiliary_->shell(A)->nfunction();
    }

    //if (print_>1) {

    //  fprintf(outfile,"\n  Disk Formation of Aia, using %d blocks.\n",nblocks);
    //  for (int b = 0; b<nblocks; b++)
    //    fprintf(outfile,"  block_starts = %d, block_stops  %d, block_sizes = %d, p_start = %d, p_sizes = %d.\n",
    //      block_starts[b], block_stops[b], block_sizes[b], p_starts[b], p_sizes[b]);
    //}

    //Buffer tensors
    double** Amn = block_matrix(max_rows,nao_*nao_);
    double** Ami = block_matrix(max_rows,nao_*nocc_);
    double** Aia = block_matrix(max_rows,nocc_*nvir_);
    double** Aii;
    double** Aaa;
    if (do_all) {
        Aii = block_matrix(max_rows,nocc_*nocc_);
        Aaa = block_matrix(max_rows,nvir_*nvir_);
    }

    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif
    int rank = 0;

    //Get an ERI object for the AO three-index integrals
    shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(), primary_, primary_));
    //Get a TEI for each thread
    const double **buffer = new const double*[nthread];
    shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
      eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory->eri());
      buffer[Q] = eri[Q]->buffer();
    }
    rifactory.reset();

    psio_address next_PSIF_DFMP2_AIA = PSIO_ZERO;
    psio_address next_PSIF_DFMP2_AII = PSIO_ZERO;
    psio_address next_PSIF_DFMP2_AAA = PSIO_ZERO;

    if (do_all) {
        // Zero everything out to prevent collision 
        double *temp = init_array(nao_*nao_);
        for (int P = 0; P < naux_; P++) {
            psio_->write(PSIF_DFMP2_AIA, "OO Integrals", (char*) temp, nocc_*nocc_*sizeof(double), next_PSIF_DFMP2_AII, &next_PSIF_DFMP2_AII);    
        }
        next_PSIF_DFMP2_AII = PSIO_ZERO;
        for (int P = 0; P < naux_; P++) {
            psio_->write(PSIF_DFMP2_AIA, "OV Integrals", (char*) temp, nocc_*nvir_*sizeof(double), next_PSIF_DFMP2_AII, &next_PSIF_DFMP2_AII);    
        }
        next_PSIF_DFMP2_AII = PSIO_ZERO;
        for (int P = 0; P < naux_; P++) {
            psio_->write(PSIF_DFMP2_AIA, "VV Integrals", (char*) temp, nvir_*nvir_*sizeof(double), next_PSIF_DFMP2_AII, &next_PSIF_DFMP2_AII);    
        }
        next_PSIF_DFMP2_AII = PSIO_ZERO;
        free(temp);
    }

    //indices for three-index integrals
    int P, MU, NU, nump, nummu, numnu, p, mu, nu, op, omu, onu, index;

    shared_ptr<SchwarzSieve> schwarz(new SchwarzSieve(primary_, schwarz_cutoff_));
    long int* schwarz_shell_pairs = schwarz->get_schwarz_shells_reverse();

    //Loop over blocks of A
    for (int block = 0; block<nblocks; block++) {
      //Zero that guy out! (Schwarzing messes with it)
      memset((void*)&Amn[0][0],'\0',p_sizes[block]*nao_*(ULI)nao_);

      //Form Amn ints
      timer_on("(A|mn)");
      #pragma omp parallel for private (P, MU, NU, p, mu, nu, nump, nummu, numnu, op, omu, onu, index, rank) schedule (dynamic) num_threads(nthread)
      for (P=block_starts[block]; P < block_stops[block]; ++P) {
        #ifdef _OPENMP
           rank = omp_get_thread_num();
        #endif
        nump = auxiliary_->shell(P)->nfunction();
        for (MU=0; MU < primary_->nshell(); ++MU) {
          nummu = primary_->shell(MU)->nfunction();
          for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU)->nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
              eri[rank]->compute_shell(P, 0, MU, NU);
              for (p=0, index=0; p < nump; ++p) {
                op = auxiliary_->shell(P)->function_index() + p;
                for (mu=0; mu < nummu; ++mu) {
                  omu = primary_->shell(MU)->function_index() + mu;
                  for (nu=0; nu < numnu; ++nu, ++index) {
                    onu = primary_->shell(NU)->function_index() + nu;
                    Amn[op-p_starts[block]][omu*nao_+onu] = buffer[rank][index]; // (op | omu onu) integral
                    Amn[op-p_starts[block]][onu*nao_+omu] = buffer[rank][index]; // (op | onu omu) integral
                  }
                }
              }
            }
          }
        }
      }
      timer_off("(A|mn)");

      //fprintf(outfile, "  Amn\n");
      //print_mat(Amn,max_rows,nao_*nao_, outfile);

      //Transform to Ami
      // (A|mi) = (Amn)C_ni
      timer_on("(A|mi)");
      C_DGEMM('N', 'N', p_sizes[block]*nao_, nocc_, nao_, 1.0, &(Amn[0][0]),
          nao_, &(Cop[0][0]), nocc_, 0.0, &(Ami[0][0]), nocc_);
      timer_off("(A|mi)");


      //fprintf(outfile, "  Ami\n");
      //print_mat(Ami,max_rows,nocc_*nao_, outfile);

      #ifdef HAVE_MKL
         int mkl_nthreads = mkl_get_max_threads();
         mkl_set_num_threads(1);
      #endif

      //Transform to Aia
      // (A|ia) = C_ma(A|mi)
      timer_on("(A|ia)");
      #pragma omp parallel for
      for (int A = 0; A<p_sizes[block]; A++) {
        C_DGEMM('T', 'N', nocc_, nvir_, nao_, 1.0, &(Ami[A][0]),
          nocc_, &(Cvp[0][0]), nvir_, 0.0, &(Aia[A][0]), nvir_);
      }
      timer_off("(A|ia)");
         
      //fprintf(outfile, "  Aia\n");
      //print_mat(Aia,max_rows,nocc_*nvir_, outfile);

      //Stripe to disk
      timer_on("(A|ia) Write");
      psio_->write(PSIF_DFMP2_AIA,"OV Integrals",(char *)(&Aia[0][0]),sizeof(double)*p_sizes[block]*nocc_*(ULI)nvir_,next_PSIF_DFMP2_AIA,&next_PSIF_DFMP2_AIA);
      timer_off("(A|ia) Write");
      
      if (do_all) { 
          timer_on("(A|ii)");
          #pragma omp parallel for
          for (int A = 0; A<p_sizes[block]; A++) {
            C_DGEMM('T', 'N', nocc_, nocc_, nao_, 1.0, &(Ami[A][0]),
              nocc_, &(Cop[0][0]), nocc_, 0.0, &(Aii[A][0]), nocc_);
          }
          timer_off("(A|ii)");
         
          //fprintf(outfile, "  Aii\n");
          //print_mat(Aii,max_rows,nocc_*nocc_, outfile);

          //Stripe to disk
          timer_on("(A|ii) Write");
          psio_->write(PSIF_DFMP2_AIA,"OO Integrals",(char *)(&Aii[0][0]),sizeof(double)*p_sizes[block]*nocc_*(ULI)nocc_,next_PSIF_DFMP2_AII,&next_PSIF_DFMP2_AII);
          timer_off("(A|ii) Write");
          
          timer_on("(A|aa)");
          #pragma omp parallel for
          for (int A = 0; A<p_sizes[block]; A++) {
            C_DGEMM('T', 'N', nvir_, nvir_, nao_, 1.0, &(Ami[A][0]),
              nocc_, &(Cvp[0][0]), nvir_, 0.0, &(Aaa[A][0]), nvir_);
          }
          timer_off("(A|aa)");
         
          //fprintf(outfile, "  Aaa\n");
          //print_mat(Aaa,max_rows,nvir_*nvir_, outfile);

          //Stripe to disk
          timer_on("(A|aa) Write");
          psio_->write(PSIF_DFMP2_AIA,"VV Integrals",(char *)(&Aaa[0][0]),sizeof(double)*p_sizes[block]*nvir_*(ULI)nvir_,next_PSIF_DFMP2_AAA,&next_PSIF_DFMP2_AAA);
          timer_off("(A|aa) Write");

        }
        #ifdef HAVE_MKL
           mkl_set_num_threads(mkl_nthreads);
        #endif

    }

    delete[] buffer;
    delete[] eri;
    free(p_starts);
    free(p_sizes);
    free(block_starts);
    free(block_stops);
    free(block_sizes);
    free_block(Amn);
    free_block(Ami);
    free_block(Aia);
    if (do_all) {
        free_block(Aii);
        free_block(Aaa);
    }
}
void DFTensor::apply_fitting(const std::string& entry)
{
  //Available memory is lower due to fitting metric 
  unsigned long int available_memory = memory_-naux_*naux_;
  unsigned long int entry_size;
  if (entry == "OO Integrals")
    entry_size = nocc_*(unsigned long int) nocc_;
  else if (entry == "OV Integrals")
    entry_size = nocc_*(unsigned long int) nvir_;
  else if (entry == "VV Integrals")
    entry_size = nvir_*(unsigned long int) nvir_;

  //max_cols to attempt
  int max_cols = available_memory/((long) naux_+naux_);
  if (max_cols > entry_size)
      max_cols = entry_size;
  if (max_cols < 1)
      max_cols = 1;

  //How many blocks are there?
  int nblocks = entry_size/max_cols;
  if (max_cols*nblocks < entry_size)
    nblocks++;

  //File handlers
  psio_address next_PSIF_DFMP2_AIA = PSIO_ZERO;
  psio_address next_PSIF_DF_TENSOR = PSIO_ZERO;

  //Prestripe
  timer_on("(Q|ia) Prestripe");
  double *Prestripe = init_array(entry_size);
  for (int Q = 0; Q < naux_; Q++) {
       psio_->write(PSIF_DF_TENSOR,entry.c_str(),(char *) &(Prestripe[0]),sizeof(double)*entry_size,next_PSIF_DF_TENSOR,&next_PSIF_DF_TENSOR);
  }
  timer_off("(Q|ia) Prestripe");
  free(Prestripe);
  next_PSIF_DF_TENSOR = PSIO_ZERO;

  //Allocate buffers
  double** Aia = block_matrix(naux_, max_cols);
  double** Qia = block_matrix(max_cols, naux_);

  int current_column = 0;
  int current_columns;
  for (int block = 0; block<nblocks; block++) {

    //How many columns in this pass?
    current_columns = max_cols;
    if (current_column+current_columns>=entry_size)
        current_columns = entry_size-current_column;

    //Read Aia in
    timer_on("(A|ia) Read");
    for (int A = 0; A<naux_; A++) {
        next_PSIF_DFMP2_AIA = psio_get_address(PSIO_ZERO,(ULI)(A*(ULI)entry_size*sizeof(double)+current_column*sizeof(double)));
        psio_->read(PSIF_DFMP2_AIA,entry.c_str(),(char *)(&Aia[A][0]),sizeof(double)*(ULI)current_columns,next_PSIF_DFMP2_AIA,&next_PSIF_DFMP2_AIA);
    }
    timer_off("(A|ia) Read");

    //Embed fitting
    timer_on("(Q|ia)");
    if (fitting_algorithm_ == "CHOLESKY") {

        #pragma omp parallel for schedule(static)
        for (int A = 0; A< naux_; A++) {
            C_DCOPY(current_columns, &Aia[A][0], 1, &Qia[0][A], naux_);
        }

        int info = C_DPOTRS('L',naux_,current_columns,metric_->get_metric()->pointer()[0],naux_,&Qia[0][0],naux_);

    } else {
        C_DGEMM('T','N',current_columns,naux_,naux_,1.0, Aia[0], max_cols, metric_->get_metric()->pointer()[0], naux_,0.0, Qia[0],naux_);
    }
    timer_off("(Q|ia)");

    //fprintf(outfile,"  Nblocks = %d, Max cols = %d, current_columns = %d, current_column = %d\n",nblocks, max_cols, current_columns, current_column);

    //Write Qia out
    timer_on("(Q|ia) Write");
    if (!Qia_striping_) {
        psio_->write(PSIF_DF_TENSOR,entry.c_str(),(char *)(&Qia[0][0]),sizeof(double)*naux_*(ULI)current_columns,next_PSIF_DF_TENSOR,&next_PSIF_DF_TENSOR);
    } else {
        for (int A = 0; A<naux_; A++) {
            next_PSIF_DF_TENSOR = psio_get_address(PSIO_ZERO,(ULI)(A*(ULI)entry_size*sizeof(double)+current_column*sizeof(double)));
            psio_->read(PSIF_DF_TENSOR,entry.c_str(),(char *)(&Qia[A][0]),sizeof(double)*(ULI)current_columns,next_PSIF_DF_TENSOR,&next_PSIF_DF_TENSOR);
        }
    }
    timer_off("(Q|ia) Write");

    // Increment column offset
    current_column += current_columns;

  }

  free_block(Aia);
  free_block(Qia);
}
void DFTensor::form_MO_integrals(unsigned long int memory_doubles, shared_ptr<Matrix> Co, shared_ptr<Matrix> Cv, bool Qia_striping, const std::string& fitting_algorithm,
        double condition, double schwarz)
{
    memory_ = memory_doubles;
    Co_ = Co;
    Cv_ = Cv;
    nocc_ = Co->colspi()[0];
    nvir_ = Cv->colspi()[0];
    Qia_striping_ = Qia_striping;
    fitting_algorithm_ = fitting_algorithm;
    fitting_condition_ = condition;
    schwarz_cutoff_ = schwarz;   
 
    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_NEW);

    form_Aia(true);
    if (fitting_algorithm_ == "CHOLESKY")
        metric_->form_cholesky_inverse();
    else if (fitting_algorithm_ == "QR")
        metric_->form_QR_inverse(condition);
    else if (fitting_algorithm_ == "EIG")
        metric_->form_eig_inverse(condition);

    apply_fitting("OO Integrals"); 
    apply_fitting("OV Integrals"); 
    apply_fitting("VV Integrals"); 
    
    metric_.reset();
    psio_->close(PSIF_DFMP2_AIA, 0);
}
void DFTensor::form_OV_integrals(unsigned long int memory_doubles, shared_ptr<Matrix> Co, shared_ptr<Matrix> Cv, bool Qia_striping, const std::string& fitting_algorithm,
        double condition, double schwarz)
{
    memory_ = memory_doubles;
    Co_ = Co;
    Cv_ = Cv;
    nocc_ = Co->colspi()[0];
    nvir_ = Cv->colspi()[0];
    Qia_striping_ = Qia_striping;
    fitting_algorithm_ = fitting_algorithm;
    fitting_condition_ = condition;
    schwarz_cutoff_ = schwarz;   
    
    psio_->open(PSIF_DFMP2_AIA, PSIO_OPEN_NEW);

    form_Aia(false);
    if (fitting_algorithm_ == "CHOLESKY")
        metric_->form_cholesky_inverse();
    else if (fitting_algorithm_ == "QR")
        metric_->form_QR_inverse(condition);
    else if (fitting_algorithm_ == "EIG")
        metric_->form_eig_inverse(condition);
    apply_fitting("OV Integrals"); 
    
    metric_.reset();
    psio_->close(PSIF_DFMP2_AIA, 0);
}
shared_ptr<TensorChunk> DFTensor::get_oo_iterator(unsigned long int memory)
{
    if (Qia_striping_)
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_DF_TENSOR,
            "OO Integrals",
            naux_,
            nocc_*nocc_,
            memory
        ));
    else
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_DF_TENSOR,
            "OO Integrals",
            nocc_*nocc_,
            naux_,
            memory
        ));
}
shared_ptr<TensorChunk> DFTensor::get_ov_iterator(unsigned long int memory)
{
    if (Qia_striping_)
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_DF_TENSOR,
            "OV Integrals",
            naux_,
            nocc_*nvir_,
            memory
        ));
    else
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_DF_TENSOR,
            "OV Integrals",
            nocc_*nvir_,
            naux_,
            memory
        ));
}
shared_ptr<TensorChunk> DFTensor::get_vv_iterator(unsigned long int memory)
{
    if (Qia_striping_)
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_DF_TENSOR,
            "VV Integrals",
            naux_,
            nvir_*nvir_,
            memory
        ));
    else
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_DF_TENSOR,
            "VV Integrals",
            nvir_*nvir_,
            naux_,
            memory
        ));
}


}

