#include "3index.h"
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <utility>
#include <ctype.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {

PSTensor::PSTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, 
                   shared_ptr<BasisSet> dealias, shared_ptr<PseudoGrid> grid) :
                   psio_(psio), primary_(primary), dealias_(dealias), grid_(grid)
{
    common_init(); 
}
PSTensor::~PSTensor()
{
    psio_->close(PSIF_PS_TENSOR, 0);
}
void PSTensor::common_init()
{
    nao_ = primary_->nbf();
    naux_ = grid_->getBlock()->getTruePoints();
    nocc_ = 0;
    nvir_ = 0;
    psio_->open(PSIF_PS_TENSOR, PSIO_OPEN_NEW);
}
shared_ptr<Matrix> PSTensor::form_Q(shared_ptr<Matrix> C)
{
    int nocc = C->colspi()[0];
    shared_ptr<Matrix> Q = form_Q_AO();
    shared_ptr<Matrix> Qnew( new Matrix("Q_i^P (Fitted Collocation Matrix)", nocc, naux_));

    C_DGEMM('T', 'N', nocc, naux_, nao_, 1.0, C->pointer()[0], nocc, Q->pointer()[0], naux_, 0.0, Qnew->pointer(0)[0], naux_); 

    return Qnew;
}
shared_ptr<Matrix> PSTensor::form_X(shared_ptr<Matrix> C)
{
    int nocc = C->colspi()[0];
    shared_ptr<Matrix> Q = form_X_AO();
    shared_ptr<Matrix> Qnew( new Matrix("X_a^P (Collocation Matrix)", nocc, naux_));

    C_DGEMM('T', 'N', nocc, naux_, nao_, 1.0, C->pointer()[0], nocc, Q->pointer()[0], naux_, 0.0, Qnew->pointer(0)[0], naux_); 

    return Qnew;
}
shared_ptr<Matrix> PSTensor::form_X_AO()
{
    shared_ptr<Matrix> X(new Matrix("Primary Collocation Metric (Dirac)", primary_->nbf(), naux_));
    double** Xp = X->pointer();
    
    shared_ptr<BasisPoints> points(new BasisPoints(primary_, naux_));    
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->getBlock());
 
    // Copy the points in
    for (int i = 0; i < naux_; i++) {
        for (int Q = 0; Q < primary_->nbf(); Q++)
            Xp[Q][i] = bpoints[i][Q]; 
    }

    return X;
}
shared_ptr<Matrix> PSTensor::form_X_dealias()
{
    shared_ptr<Matrix> X(new Matrix("Dealias Collocation Metric (Dirac)", dealias_->nbf(), naux_));
    double** Xp = X->pointer();
    
    shared_ptr<BasisPoints> points(new BasisPoints(dealias_, naux_));    
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->getBlock());
 
    // Copy the points in
    for (int i = 0; i < naux_; i++) {
        for (int Q = 0; Q < dealias_->nbf(); Q++)
            Xp[Q][i] = bpoints[i][Q]; 
    }

    return X;
}
shared_ptr<Matrix> PSTensor::form_S()
{
    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_, primary_, primary_, primary_));
    shared_ptr<OneBodyAOInt> o2(fact->ao_overlap());
    shared_ptr<Matrix> S(new Matrix("Primary Overlap Matrix", primary_->nbf(), primary_->nbf()));
    o2->compute(S);

    return S;
}
shared_ptr<Matrix> PSTensor::form_S_dealias()
{
    shared_ptr<IntegralFactory> fact(new IntegralFactory(dealias_, primary_, primary_, primary_));
    shared_ptr<OneBodyAOInt> o2(fact->ao_overlap());
    shared_ptr<Matrix> S(new Matrix("Dealias Overlap Matrix", dealias_->nbf(), primary_->nbf()));
    o2->compute(S);

    return S;
}
shared_ptr<Matrix> PSTensor::form_Q_AO()
{
    int nbf = primary_->nbf();
    int ndf = dealias_->nbf();
    int ntotal = nbf + ndf;

    shared_ptr<Matrix> Q(new Matrix("Q (Least-Squares)", nbf, naux_));
    double** Qp = Q->pointer();     
    
    shared_ptr<Matrix> X = form_X_AO();
    shared_ptr<Matrix> Xd = form_X_dealias();
    shared_ptr<Matrix> S = form_S();
    shared_ptr<Matrix> Sd = form_S_dealias();

    //X->print();
    //Xd->print();
    //S->print();
    //Sd->print();

    double** Xp = X->pointer();
    double** Xdp = Xd->pointer();
    double** Sp = S->pointer();
    double** Sdp = Sd->pointer();

    // Copy S for later use (final projection)
    shared_ptr<Matrix> S2(new Matrix("S (Copy)", nbf, nbf));
    double** S2p = S->pointer();
    memcpy(static_cast<void*> (S2p[0]), static_cast<void*> (Sp[0]), nbf*nbf*sizeof(double)); 

    // Orthogonalize the primary and dealias bases
    C_DPOTRF('L', nbf, Sp[0], nbf);
    //S->print();
    C_DPOTRS('L', nbf, ndf, Sp[0], nbf, Sdp[0], nbf);
    //Sd->print();    

    // Build the Collocation matrix
    shared_ptr<Matrix> R(new Matrix("R (Full Collocation Matrix)", ntotal, naux_));
    double** Rp = R->pointer();
    double** Rdp = &Rp[nbf];
    memcpy(static_cast<void*> (Rp[0]), static_cast<void*> (Xp[0]), nbf*naux_*sizeof(double)); 
    memcpy(static_cast<void*> (Rp[nbf]), static_cast<void*> (Xdp[0]), ndf*naux_*sizeof(double)); 
    //R->print();

    // Remove the overlap from the primary basis on the dealias collocation partition  
    C_DGEMM('N','N', ndf, naux_, nbf, -1.0, Sdp[0], nbf, Xp[0], naux_, 1.0, Rdp[0], naux_); 
    //R->print();

    // Roll the square root of the weights into R (one square root w will chill until the end)
    double* w = new double[naux_];
    double* wp = grid_->getBlock()->getWeights();
    for (int P = 0; P < naux_; P++)
        w[P] = sqrt(wp[P]);

    for (int P = 0; P < naux_; P++) 
        C_DSCAL(ntotal, w[P], &Rp[0][P], naux_);
    //R->print();

    // Form (X'X)X'\sqrt(w) via QR Decomposition
    double* tau = new double[ntotal];

    // First, find out how much workspace to provide
    double work_size;
    C_DGEQRF(naux_,ntotal,Rp[0],naux_,tau,&work_size, -1);  

    // Now, do the QR decomposition
    int lwork = (int)work_size;
    double *work = new double[lwork];
    C_DGEQRF(naux_,ntotal,Rp[0],naux_,tau,work, lwork);  
    R->set_name("Q (of QR Decomposition");
    delete[] work;

    //R->print();

    // Put R in the upper triangle where it belongs
    shared_ptr<Matrix> r(new Matrix("R (of QR decomposition)", ntotal, ntotal));
    double** rp = r->pointer();
    for (int i = 0; i < ntotal; i++)
        for (int j = i; j < ntotal; j++) {
            rp[i][j] = Rp[j][i]; 
        }
 
    //r->print();
    
    // First, find out how much workspace to provide
    C_DORGQR(naux_,ntotal,ntotal,Rp[0],naux_,tau,&work_size,-1); 

    // Now, form Q
    lwork = (int)work_size;
    work = new double[lwork];
    C_DORGQR(naux_,ntotal,ntotal,Rp[0],naux_,tau,work,lwork); 
    delete[] work;
    delete[] tau;

    //R->print();
    
    // Scale Q' by sqrt w
    for (int P = 0; P < naux_; P++) 
        C_DSCAL(ntotal, w[P], &Rp[0][P], naux_);
    delete[] w;

    //R->print();
    
    // Backsolve R^-1 Q' sqrt w
    C_DTRSM('L','U','N','N', ntotal, naux_, 1.0, rp[0], ntotal, Rp[0], naux_);

    //R->print();

    // Do P R^-1 Q' sqrt w, using S (this bit is confusing, Friesner never gets it straight)
    C_DGEMM('N','N', nbf, naux_, nbf, 1.0, S2p[0], nbf, Rp[0], naux_, 0.0, Qp[0], naux_);

    //Q->print();

    return Q;
}
void PSTensor::form_Aia(bool do_all)
{
    psio_address next_PSIF_PS_AIA = PSIO_ZERO;
    psio_address next_PSIF_PS_AII = PSIO_ZERO;
    psio_address next_PSIF_PS_AAA = PSIO_ZERO;

    if (do_all) {
        // Zero everything out to prevent collision 
        double *temp = init_array(nao_*nao_);
        for (int P = 0; P < naux_; P++) {
            psio_->write(PSIF_DFMP2_AIA, "OO Integrals", (char*) temp, nocc_*nocc_*sizeof(double), next_PSIF_PS_AII, &next_PSIF_PS_AII);    
        }
        next_PSIF_PS_AII = PSIO_ZERO;
        for (int P = 0; P < naux_; P++) {
            psio_->write(PSIF_DFMP2_AIA, "OV Integrals", (char*) temp, nocc_*nvir_*sizeof(double), next_PSIF_PS_AII, &next_PSIF_PS_AII);    
        }
        next_PSIF_PS_AII = PSIO_ZERO;
        for (int P = 0; P < naux_; P++) {
            psio_->write(PSIF_DFMP2_AIA, "VV Integrals", (char*) temp, nvir_*nvir_*sizeof(double), next_PSIF_PS_AII, &next_PSIF_PS_AII);    
        }
        next_PSIF_PS_AII = PSIO_ZERO;
        free(temp);
    }

    std::vector<shared_ptr<PseudospectralInt> > ints;
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif
    ints.resize(nthread); 
   
    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_, primary_, primary_, primary_));
    for (int thread = 0; thread < nthread; thread++)
        ints[thread] = shared_ptr<PseudospectralInt>(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));
        
    shared_ptr<Matrix> Amn(new Matrix("(A|mn) Pseudospectral Integrals", nthread, nao_*nao_));
    shared_ptr<Matrix> Ami(new Matrix("(A|mi) Pseudospectral Integrals", nthread, nao_*nocc_)); 

    // TODO provide VV/OO integrals
    unsigned long int scratch = (ULI)nthread*(nao_*nao_ + nao_*nocc_);
    int max_rows = (memory_ - scratch)/ (nocc_ *(ULI) nvir_);
    if (max_rows > naux_)
        max_rows = naux_;
    if (max_rows < 1L)
        max_rows = 1L;

    shared_ptr<Matrix> Aia(new Matrix("(A|ia) Pseudospectral Integrals", max_rows, nvir_*nocc_)); 

    int nblocks = naux_ / max_rows;
    if (nblocks * max_rows != naux_)
        nblocks++;
    
    std::vector<int> block_starts;
    std::vector<int> block_sizes;
    block_starts.resize(nblocks);
    block_sizes.resize(nblocks);
    
    // Naive distribution
    block_starts[0] = 0;
    block_sizes[0] = max_rows;
    int gimp_delta = 0;
    for (int Q = 1; Q < nblocks; Q++) {
        block_starts[Q] = block_starts[Q - 1] + max_rows;
        if (block_starts[Q] + max_rows >= naux_) {
            block_sizes[Q] = naux_ - block_starts[Q];
            gimp_delta = (block_sizes[Q - 1] - block_sizes[Q]); 
        } else {
            block_sizes[Q] = max_rows;
        }
    }
   
    // Now Level the gimp out
    for (int Q = 0; Q < gimp_delta - 1; Q++) {
        block_sizes[nblocks - 2- Q] -= 1;
        block_starts[nblocks - 2 - Q] -= (gimp_delta - 1 - Q);
    }

    #ifdef _MKL
       int mkl_nthreads = mkl_get_max_threads();
       mkl_set_num_threads(1);
    #endif

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Aiap = Aia->pointer();

    double** Cop = Co_->pointer();
    double** Cvp = Cv_->pointer();

    shared_ptr<GridBlock> block = grid_->getBlock();
    double* x = block->getX();
    double* y = block->getY();
    double* z = block->getZ();

    shared_ptr<SchwarzSieve> schwarz(new SchwarzSieve(primary_,schwarz_cutoff_));
    long int* schwarz_shell_pairs = schwarz->get_schwarz_shells_reverse();

    for (int block_i = 0; block_i < nblocks; block_i++) {

        int start = block_starts[block_i];
        int size = block_sizes[block_i]; 
        
        #pragma omp parallel for schedule(static,1)
        for (int Q = 0; Q < size; Q++) {
            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif  
            const double* buffer = ints[rank]->buffer();
        
            int Q_global = Q + start;
         
            // Generation of integrals
            memset((void*)&Amnp[rank][0],'\0',nao_*(ULI)nao_*sizeof(double));
            ints[rank]->set_point(x[Q_global], y[Q_global], z[Q_global]);
            for (int MU=0; MU < primary_->nshell(); ++MU) {
              int nummu = primary_->shell(MU)->nfunction();
              for (int NU=0; NU <= MU; ++NU) {
                int numnu = primary_->shell(NU)->nfunction();
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                  ints[rank]->compute_shell(MU, NU);
                  for (int mu=0, index = 0; mu < nummu; ++mu) {
                    int omu = primary_->shell(MU)->function_index() + mu;
                    for (int nu=0; nu < numnu; ++nu, ++index) {
                      int onu = primary_->shell(NU)->function_index() + nu;
                      Amnp[rank][omu*nao_+onu] = buffer[index]; // (op | omu onu) integral
                      Amnp[rank][onu*nao_+omu] = buffer[index]; // (op | onu omu) integral
                    }
                  }
                }
              }
            } 
            
            // First half transform 
            C_DGEMM('N', 'N', nao_, nocc_, nao_, 1.0, &(Amnp[rank][0]),        
                nao_, &(Cop[0][0]), nocc_, 0.0, &(Amip[rank][0]), nocc_);
 
            // Second half transform 
            C_DGEMM('T', 'N', nocc_, nvir_, nao_, 1.0, &(Amip[rank][0]),
                nocc_, &(Cvp[0][0]), nvir_, 0.0, &(Aiap[Q][0]), nvir_);
        
            //TODO: Add VV/OO integrals    
        }
        
        psio_->write(PSIF_DFMP2_AIA, "OV Integrals", (char*) Aiap[0], (ULI)size*nvir_*nocc_*sizeof(double), next_PSIF_PS_AIA, &next_PSIF_PS_AIA);
    }
    
    #ifdef _MKL
       mkl_set_num_threads(mkl_nthreads);
    #endif
}
void PSTensor::restripe(const std::string& entry)
{
  unsigned long int available_memory = memory_;
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
  psio_address next_PSIF_PS_TENSOR = PSIO_ZERO;

  //Prestripe
  //timer_on("(Q|ia) Prestripe");
  double *Prestripe = init_array(entry_size);
  for (int Q = 0; Q < naux_; Q++) {
       psio_->write(PSIF_PS_TENSOR,entry.c_str(),(char *) &(Prestripe[0]),sizeof(double)*entry_size,next_PSIF_PS_TENSOR,&next_PSIF_PS_TENSOR);
  }
  //timer_off("(Q|ia) Prestripe");
  free(Prestripe);
  next_PSIF_PS_TENSOR = PSIO_ZERO;

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


    // Transpose
    for (int A = 0; A < naux_; A++){
        C_DCOPY(current_columns, &Aia[A][0], 1, &Qia[0][A], naux_);
    } 

    //fprintf(outfile,"  Nblocks = %d, Max cols = %d, current_columns = %d, current_column = %d\n",nblocks, max_cols, current_columns, current_column);
    //fprintf(outfile, "  Aia\n");
    //print_mat(Qia, max_cols, naux_, outfile);

    //Write Qia out
    timer_on("(Q|ia) Write");
    if (!Qia_striping_) {
        psio_->write(PSIF_PS_TENSOR,entry.c_str(),(char *)(&Qia[0][0]),sizeof(double)*naux_*(ULI)current_columns,next_PSIF_PS_TENSOR,&next_PSIF_PS_TENSOR);
    } else {
        for (int A = 0; A<naux_; A++) {
            next_PSIF_PS_TENSOR = psio_get_address(PSIO_ZERO,(ULI)(A*(ULI)entry_size*sizeof(double)+current_column*sizeof(double)));
            psio_->read(PSIF_PS_TENSOR,entry.c_str(),(char *)(&Qia[A][0]),sizeof(double)*(ULI)current_columns,next_PSIF_PS_TENSOR,&next_PSIF_PS_TENSOR);
        }
    }
    timer_off("(Q|ia) Write");

    // Increment column offset
    current_column += current_columns;

  }

  free_block(Aia);
  free_block(Qia);
}
void PSTensor::form_MO_integrals(unsigned long int memory_doubles, shared_ptr<Matrix> Co, shared_ptr<Matrix> Cv, bool Qia_striping, double schwarz)
{
    fprintf(outfile, "\n  ==> PS Tensor OO/OV/VV Integrals <==\n\n");
    fprintf(outfile, "  %s striping will be used for the integrals.\n", (Qia_striping ? "(Q|ia)" : "(ia|Q)"));
    fprintf(outfile, "  A Cauchy-Schwarz sieve with cutoff of %7.3E will be applied to AO integrals.\n", schwarz);
    
    memory_ = memory_doubles;
    Co_ = Co;
    Cv_ = Cv;
    nocc_ = Co->colspi()[0];
    nvir_ = Cv->colspi()[0];
    Qia_striping_ = Qia_striping;
    schwarz_cutoff_ = schwarz;   
 
    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);
    form_Aia(true);
    restripe("OO Integrals");
    restripe("OV Integrals");
    restripe("VV Integrals");
    psio_->close(PSIF_DFMP2_AIA,0);
}
void PSTensor::form_OV_integrals(unsigned long int memory_doubles, shared_ptr<Matrix> Co, shared_ptr<Matrix> Cv, bool Qia_striping, double schwarz)
{
    fprintf(outfile, "\n  ==> PS Tensor OV Integrals <==\n\n");
    fprintf(outfile, "  %s striping will be used for the integrals.\n", (Qia_striping ? "(Q|ia)" : "(ia|Q)"));
    fprintf(outfile, "  A Cauchy-Schwarz sieve with cutoff of %7.3E will be applied to AO integrals.\n", schwarz);

    memory_ = memory_doubles;
    Co_ = Co;
    Cv_ = Cv;
    nocc_ = Co->colspi()[0];
    nvir_ = Cv->colspi()[0];
    Qia_striping_ = Qia_striping;
    schwarz_cutoff_ = schwarz;   
    

    psio_->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);
    form_Aia(false);
    restripe("OV Integrals");
    psio_->close(PSIF_DFMP2_AIA,0);
}
shared_ptr<TensorChunk> PSTensor::get_oo_iterator(unsigned long int memory)
{
    if (Qia_striping_)
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_PS_TENSOR,
            "OO Integrals",
            naux_,
            nocc_*nocc_,
            memory
        ));
    else
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_PS_TENSOR,
            "OO Integrals",
            nocc_*nocc_,
            naux_,
            memory,
            nocc_
        ));
}
shared_ptr<TensorChunk> PSTensor::get_ov_iterator(unsigned long int memory)
{
    if (Qia_striping_)
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_PS_TENSOR,
            "OV Integrals",
            naux_,
            nocc_*nvir_,
            memory
        ));
    else
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_PS_TENSOR,
            "OV Integrals",
            nocc_*nvir_,
            naux_,
            memory,
            nvir_
        ));
}
shared_ptr<TensorChunk> PSTensor::get_vv_iterator(unsigned long int memory)
{
    if (Qia_striping_)
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_PS_TENSOR,
            "VV Integrals",
            naux_,
            nvir_*nvir_,
            memory
        ));
    else
        return shared_ptr<TensorChunk>(new TensorChunk(
            psio_,
            PSIF_PS_TENSOR,
            "VV Integrals",
            nvir_*nvir_,
            naux_,
            memory,
            nvir_
        ));
}


DirectPSTensor::DirectPSTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, 
                   shared_ptr<BasisSet> dealias, shared_ptr<PseudoGrid> grid) :
                   PSTensor(psio,primary,dealias,grid)
{
}
DirectPSTensor::~DirectPSTensor()
{
}
void DirectPSTensor::initialize_OV_integrals(shared_ptr<Matrix> Cocc, shared_ptr<Matrix> Cvir, unsigned long int memory, double schwarz_cutoff)
{
    memory_ = memory;
    Co_ = Cocc;
    Cv_ = Cvir;
    nocc_ = Cocc->colspi()[0];
    nvir_ = Cvir->colspi()[0];
    schwarz_ = shared_ptr<SchwarzSieve>(new SchwarzSieve(primary_, schwarz_cutoff));

    // Initialize indexing
    block_ = -1;

    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif
    ints_.resize(nthread); 
   
    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_, primary_, primary_, primary_));
    for (int thread = 0; thread < nthread; thread++)
        ints_[thread] = shared_ptr<PseudospectralInt>(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));
        
    Amn_ = shared_ptr<Matrix>(new Matrix("(A|mn) Pseudospectral Integrals", nthread, nao_*nao_));
    Ami_ = shared_ptr<Matrix>(new Matrix("(A|mi) Pseudospectral Integrals", nthread, nao_*nocc_)); 

    max_rows_ = memory_ / (nocc_ *(ULI) nvir_);
    if (max_rows_ > naux_)
        max_rows_ = naux_;
    if (max_rows_ < 1L)
        max_rows_ = 1L;

    Aia_ = shared_ptr<Matrix>(new Matrix("(A|ia) Pseudospectral Integrals", max_rows_, nvir_*nocc_)); 

    int nblocks = naux_ / max_rows_;
    if (nblocks * max_rows_ != naux_)
        nblocks++;
    
    block_starts_.resize(nblocks);
    block_sizes_.resize(nblocks);
    
    // Naive distribution
    block_starts_[0] = 0;
    block_sizes_[0] = max_rows_;
    int gimp_delta = 0;
    for (int Q = 1; Q < nblocks; Q++) {
        block_starts_[Q] = block_starts_[Q - 1] + max_rows_;
        if (block_starts_[Q] + max_rows_ >= naux_) {
            block_sizes_[Q] = naux_ - block_starts_[Q];
            gimp_delta = (block_sizes_[Q - 1] - block_sizes_[Q]); 
        } else {
            block_sizes_[Q] = max_rows_;
        }
    }
   
    // Now Level the gimp out
    for (int Q = 0; Q < gimp_delta - 1; Q++) {
        block_sizes_[nblocks - 2- Q] -= 1;
        block_starts_[nblocks - 2 - Q] -= (gimp_delta - 1 - Q);
    }
    
}
void DirectPSTensor::compute_block(int block_number)
{
    if (block_ = block_number) 
        return;

    block_ = block_number;

    int start = block_starts_[block_];
    int size = block_sizes_[block_]; 

    #ifdef _MKL
       int mkl_nthreads = mkl_get_max_threads();
       mkl_set_num_threads(1);
    #endif

    double** Amnp = Amn_->pointer();
    double** Amip = Ami_->pointer();
    double** Aiap = Aia_->pointer();

    double** Cop = Co_->pointer();
    double** Cvp = Cv_->pointer();

    shared_ptr<GridBlock> block = grid_->getBlock();
    double* x = block->getX();
    double* y = block->getY();
    double* z = block->getZ();

    long int* schwarz_shell_pairs = schwarz_->get_schwarz_shells_reverse();

    #pragma omp parallel for schedule(static,1)
    for (int Q = 0; Q < size; Q++) {
        int rank = 0;
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif  
        const double* buffer = ints_[rank]->buffer();
    
        int Q_global = Q + start;
     
        // Generation of integrals
        memset((void*)&Amnp[rank][0],'\0',nao_*(ULI)nao_*sizeof(double));
        ints_[rank]->set_point(x[Q_global], y[Q_global], z[Q_global]);
        for (int MU=0; MU < primary_->nshell(); ++MU) {
          int nummu = primary_->shell(MU)->nfunction();
          for (int NU=0; NU <= MU; ++NU) {
            int numnu = primary_->shell(NU)->nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
              ints_[rank]->compute_shell(MU, NU);
              for (int mu=0, index = 0; mu < nummu; ++mu) {
                int omu = primary_->shell(MU)->function_index() + mu;
                for (int nu=0; nu < numnu; ++nu, ++index) {
                  int onu = primary_->shell(NU)->function_index() + nu;
                  Amnp[rank][omu*nao_+onu] = buffer[index]; // (op | omu onu) integral
                  Amnp[rank][onu*nao_+omu] = buffer[index]; // (op | onu omu) integral
                }
              }
            }
          }
        } 
        
        // First half transform 
        C_DGEMM('N', 'N', nao_, nocc_, nao_, 1.0, &(Amnp[rank][0]),        
            nao_, &(Cop[0][0]), nocc_, 0.0, &(Amip[rank][0]), nocc_);
 
        // Second half transform 
        C_DGEMM('T', 'N', nocc_, nvir_, nao_, 1.0, &(Amip[rank][0]),
            nocc_, &(Cvp[0][0]), nvir_, 0.0, &(Aiap[Q][0]), nvir_);
        
    }

    #ifdef _MKL
       mkl_set_num_threads(mkl_nthreads);
    #endif

}

}
