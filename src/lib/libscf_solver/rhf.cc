/*
 *  rhf.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/10/08.
 *
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <libmints/basisset_parser.h>
#include "integralfunctors.h"

#include "rhf.h"
#include "pseudospectral.h"
#include "df.h"

#include "rhf_functor.h"

#define TIME_SCF
#define _DEBUG

using namespace psi;
using namespace std;

namespace psi { namespace scf {

static void sort_rows_based_on_energies(SimpleMatrix* C, double *energies, int *order_mapping)
{
    unsigned int i, j;
    int itemp;
    double dtemp;
    int length = C->nrow();

    // Populate order_mapping with original ordering.
    for (i=0; i< length; ++i)
        order_mapping[i] = i;

    // Sort using Quicksort algorithm
    for (i=0; i<length; ++i) {
        for (j = i+1; j<length; ++j) {
            if (energies[i] > energies[j]) {
                C->swap_rows(i, j);

                dtemp = energies[i];
                energies[i] = energies[j];
                energies[j] = dtemp;

                itemp = order_mapping[i];
                order_mapping[i] = order_mapping[j];
                order_mapping[j] = itemp;
            }
        }
    }
}

RHF::RHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : HF(options, psio, chkpt)
{
    common_init();
}

RHF::RHF(Options& options, shared_ptr<PSIO> psio)
    : HF(options, psio)
{
    common_init();
}

RHF::~RHF()
{
}

void RHF::common_init()
{
    Drms_ = 0.0;

    // Allocate matrix memory
    Fa_        = SharedMatrix(factory_->create_matrix("F"));
    Fb_        = Fa_;
    Ca_        = SharedMatrix(factory_->create_matrix("C"));
    Cb_        = Ca_;
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = epsilon_a_;
    Da_        = SharedMatrix(factory_->create_matrix("D"));
    Db_        = Da_;
    D_         = Da_;
    Dold_      = SharedMatrix(factory_->create_matrix("D old"));
    G_         = SharedMatrix(factory_->create_matrix("G"));
    J_         = SharedMatrix(factory_->create_matrix("J"));
    K_         = SharedMatrix(factory_->create_matrix("K"));

    //What are we using?
    fprintf(outfile, "  SCF Algorithm Type is %s.\n", scf_type_.c_str());
    // Print DIIS status
    fprintf(outfile, "  DIIS %s.\n", diis_enabled_ ? "enabled" : "disabled");

    fflush(outfile);
    // Allocate memory for PK matrix

#if CUSTOM_PK_CODE
    // PK super matrix for fast G
    pk_ = NULL;
    G_vector_ = NULL;
    if (scf_type_ == "PK") {
        pk_size_ = 0;
        pk_pairs_ = 0;
        pk_symoffset_ = new int[nirrep_];
        for (int h=0; h<nirrep_; ++h) {
            pk_symoffset_[h] = pk_pairs_;
            // Add up possible pair combinations that yield A1 symmetry
            pk_pairs_ += nsopi_[h]*(nsopi_[h] + 1)/2;
        }
        // Compute the number of pairs in PK
        pk_size_ = INDEX2(pk_pairs_-1, pk_pairs_-1) + 1;
        allocate_PK();

        // Allocate memory for threading the PK
        int nthread = 1;
#ifdef _OPENMP
        nthread = omp_get_max_threads();
#endif
        G_vector_ = block_matrix(nthread, pk_pairs_);
    }
#endif
}
void RHF::finalize()
{
#if CUSTOM_PK_CODE
    if (pk_)
        delete[] pk_;
    if (G_vector_)
        free_block(G_vector_);
#endif

    Dold_.reset();
    G_.reset();
    J_.reset();
    K_.reset();

    HF::finalize();
}
SharedMatrix RHF::Da() const
{
    return D_;
}

void RHF::save_density_and_energy()
{
    Dold_->copy(D_);  // Save previous density
    Eold_ = E_;       // Save previous energy
}

void RHF::form_G()
{
#if CUSTOM_PK_CODE
    if (scf_type_ == "PK"){
        form_G_from_PK();
        return;
    }
#endif
    J_K_Functor jk_builder(G_, K_, D_, Ca_, nalphapi_);
    process_tei<J_K_Functor>(jk_builder);
    G_->subtract(K_);
}

void RHF::save_dual_basis_projection()
{
    if (print_)
        fprintf(outfile,"\n  Computing dual basis set projection from %s to %s.\n  Results will be stored in File 100.\n", \
            options_.get_str("BASIS").c_str(),options_.get_str("DUAL_BASIS_SCF").c_str());

    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    shared_ptr<BasisSet> dual_basis = BasisSet::construct(parser, molecule_, "DUAL_BASIS_SCF");

    SharedMatrix C_2 = dualBasisProjection(Ca_,doccpi_,basisset_,dual_basis);
    C_2->set_name("DUAL BASIS MOS");
    if(print_>3)
        C_2->print(outfile);

    psio_->open(PSIF_SCF_DB_MOS,PSIO_OPEN_NEW);
    C_2->save(psio_, PSIF_SCF_DB_MOS, Matrix::SubBlocks);
    psio_->write_entry(PSIF_SCF_DB_MOS,"DB SCF Energy",(char *) &(E_),sizeof(double));
    psio_->write_entry(PSIF_SCF_DB_MOS,"DB NIRREPS",(char *) &(nirrep_),sizeof(int));
    psio_->write_entry(PSIF_SCF_DB_MOS,"DB DOCCPI",(char *) (doccpi_),8*sizeof(int));
    psio_->write_entry(PSIF_SCF_DB_MOS,"DB SOCCPI",(char *) (soccpi_),8*sizeof(int));
    psio_->write_entry(PSIF_SCF_DB_MOS,"DB NALPHAPI",(char *) (nalphapi_),8*sizeof(int));
    psio_->write_entry(PSIF_SCF_DB_MOS,"DB NBETAPI",(char *) (nbetapi_),8*sizeof(int));
    psio_->close(PSIF_SCF_DB_MOS,1);
}

void RHF::save_information()
{
}

void RHF::save_fock()
{
    if (initialized_diis_manager_ == false) {
        diis_manager_ = shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk, psio_));
        diis_manager_->set_error_vector_size(1, DIISEntry::Matrix, Fa_.get());
        diis_manager_->set_vector_size(1, DIISEntry::Matrix, Fa_.get());
        initialized_diis_manager_ = true;
    }

    // Determine error matrix for this Fock
    SharedMatrix FDS(factory_->create_matrix()), DS(factory_->create_matrix());
    SharedMatrix SDF(factory_->create_matrix()), DF(factory_->create_matrix());

    // FDS = Fa_ * D_ * S_;
    DS->gemm(false, false, 1.0, D_, S_, 0.0);
    FDS->gemm(false, false, 1.0, Fa_, DS, 0.0);
    // SDF = S_ * D_ * Fa_;
    DF->gemm(false, false, 1.0, D_, Fa_, 0.0);
    SDF->gemm(false, false, 1.0, S_, DF, 0.0);

    Matrix FDSmSDF;
    FDSmSDF.copy(FDS);
    FDSmSDF.subtract(SDF);

    // Orthonormalize the error matrix
    if (!canonical_X_)
        FDSmSDF.transform(Shalf_);
    else
        FDSmSDF.transform(X_);

    //FDSmSDF.print(outfile);

    diis_manager_->add_entry(2, &FDSmSDF, Fa_.get());
}

bool RHF::diis()
{
    return diis_manager_->extrapolate(1, Fa_.get());
}

bool RHF::test_convergency()
{
    // energy difference
    double ediff = E_ - Eold_;

    // RMS of the density
    Matrix D_rms;
    D_rms.copy(D_);
    D_rms.subtract(Dold_);
    Drms_ = D_rms.rms();

    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}


void RHF::form_F()
{
    Fa_->copy(H_);
    Fa_->add(G_);

#ifdef _DEBUG
    if (debug_) {
        Fa_->print(outfile);
    }
#endif
}

void RHF::form_C()
{
    if (!canonical_X_) {
        Matrix eigvec;
        factory_->create_matrix(eigvec);

        Fa_->transform(Shalf_);
        Fa_->diagonalize(eigvec, *epsilon_a_);

        Ca_->gemm(false, false, 1.0, Shalf_, eigvec, 0.0);

        // Save C to checkpoint file.
        //double **vectors = C_->to_block_matrix();
        //chkpt_->wt_scf(vectors);
        //free_block(vectors);

#ifdef _DEBUG
        if (debug_) {
            Ca_->eivprint(epsilon_a_);
        }
#endif
    } else {

        Ca_->zero();

        for (int h = 0; h<Ca_->nirrep(); h++) {

            int norbs = nsopi_[h];
            int nmos = nmopi_[h];

            //fprintf(outfile,"  Norbs = %d, Nmos = %d\n",norbs,nmos); fflush(outfile);

            double **X = block_matrix(norbs,nmos);
            for (int m = 0 ; m<norbs; m++)
                for (int i = 0; i<nmos; i++)
                    X[m][i] = X_->get(h,m,i);

            double **F = block_matrix(norbs,norbs);
            for (int m = 0 ; m<norbs; m++)
                for (int i = 0; i<norbs; i++)
                    F[m][i] = Fa_->get(h,m,i);

            double **C = block_matrix(norbs,nmos);
            double **Temp = block_matrix(nmos,norbs);
            double **Fp = block_matrix(nmos,nmos);
            double **Cp = block_matrix(nmos,nmos);

            //Form F' = X'FX for canonical orthogonalization
            C_DGEMM('T','N',nmos,norbs,norbs,1.0,X[0],nmos,F[0],norbs,0.0,Temp[0],norbs);
            C_DGEMM('N','N',nmos,nmos,norbs,1.0,Temp[0],norbs,X[0],nmos,0.0,Fp[0],nmos);

            //Form C' = eig(F')
            double *eigvals = init_array(nmos);
            sq_rsp(nmos, nmos, Fp,  eigvals, 1, Cp, 1.0e-14);

            for (int i = 0; i<nmos; i++)
                epsilon_a_->set(h,i,eigvals[i]);

            //fprintf(outfile,"  Canonical orbital eigenvalues");
            //for (int i = 0; i<nmos; i++)
            //    fprintf(outfile,"  %d: %10.6f\n",i+1,eigvals[i]);

            free(eigvals);

            //Form C = XC'
            C_DGEMM('N','N',norbs,nmos,nmos,1.0,X[0],nmos,Cp[0],nmos,0.0,C[0],nmos);

            for (int m = 0 ; m<norbs; m++)
                for (int i = 0; i<nmos; i++)
                    Ca_->set(h,m,i,C[m][i]);

            free_block(X);
            free_block(F);
            free_block(C);
            free_block(Temp);
            free_block(Cp);
            free_block(Fp);
        }

    }
    find_occupation();
}

void RHF::form_D()
{
    int h, i, j;
    int *opi = D_->rowspi();
    int nirreps = D_->nirrep();
    int norbs = basisset_->nbf();

    double** C = Ca_->to_block_matrix();
    double** D = block_matrix(norbs,norbs);

    int offset_R = 0;
    int offset_C = 0;
    for (h = 0; h<nirreps; ++h) {
        C_DGEMM('n','t',opi[h],opi[h],doccpi_[h],1.0,&C[offset_R][offset_C],nmo_,&C[offset_R][offset_C],nmo_,0.0,&D[offset_R][offset_R],norbs);

        for (i = 0; i<opi[h]; i++)
            for (j = 0; j<opi[h]; j++)
                D_->set(h,i,j,D[offset_R+i][offset_R+j]);

        offset_R += opi[h];
        offset_C += nmopi_[h];
    }

    free_block(C);
    free_block(D);


#ifdef _DEBUG
    if (debug_) {
        D_->print(outfile);
    }
#endif
}

double RHF::compute_initial_E()
{
    double Etotal = nuclearrep_ + D_->vector_dot(H_);
    return Etotal;
}

double RHF::compute_E()
{
    Matrix HplusF;
    HplusF.copy(H_);
    HplusF.add(Fa_);
    double Etotal = nuclearrep_ + D_->vector_dot(HplusF);
    return Etotal;
}

#if CUSTOM_PK_CODE
void RHF::allocate_PK()
{
    // The size of the pk matrix is determined in HF::form_indexing
    // Allocate memory for the PK matrix (using a vector)
    if (pk_size_ < (memory_ / sizeof(double))) {
        pk_ = new double[pk_size_];

        if (pk_ == NULL) {
            fprintf(outfile, "  Insufficient free system memory for in-core PK implementation.\n");
            fprintf(outfile, "  Switching to out-of-core algorithm.\n");
            scf_type_ = "OUT_OF_CORE";
        } else {
             // Zero out PK
            memset(pk_, 0, pk_size_*sizeof(double));
             fprintf(outfile, "  Allocated %lu elements (%lu pairs) for PK. (%5.2f MiB)\n\n", (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
        }
    }
    else {
        fprintf(outfile, "  Insufficient memory for in-core PK implementation.\n");
        fprintf(outfile, "  Would need %lu elements of double memory. (%5f MiB)\n", (unsigned long)pk_size_, pk_size_ * sizeof(double) / 1048576.0);
        fprintf(outfile, "  Switching to out-of-core algorithm.\n");
        scf_type_ = "OUT_OF_CORE";
    }
}

void RHF::form_PK()
{
    // struct iwlbuf ERIIN;
    int ilsti, nbuf;
    int i, j, k, l;
    int ii, jj, kk, ll;
    int is, js, ks, ls;
    int fi;
    size_t bra, ket, braket=0;
    int idx;
    int counter = 0;
    int pk_counter = 0;
    bool pk_flag = false;
    double value;

    // PK zeroed out during allocation
    fprintf(outfile, "  Forming PK matrix.\n");
    fflush(outfile);

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);

    do {
        ilsti = ERIIN.last_buffer();
        nbuf  = ERIIN.buffer_count();

        fi = 0;
        for (idx=0; idx<nbuf; ++idx) {
            if (ERIIN.labels()[fi] >= 0) {
                i = ERIIN.labels()[fi];
                pk_flag = false;
            }
            else {
                i = -ERIIN.labels()[fi];
                pk_flag = true;
            }
            i = ERIIN.labels()[fi] >= 0 ? ERIIN.labels()[fi] : -ERIIN.labels()[fi];
            j = ERIIN.labels()[fi+1];
            k = ERIIN.labels()[fi+2];
            l = ERIIN.labels()[fi+3];
            value = ERIIN.values()[idx];
            fi += 4;

            // Get the symmetries
            is = so2symblk_[i];
            js = so2symblk_[j];
            ks = so2symblk_[k];
            ls = so2symblk_[l];

            // Get the offset of the SO index in its symblock
            ii = so2index_[i];
            jj = so2index_[j];
            kk = so2index_[k];
            ll = so2index_[l];

            // J
            if ((is == js) && (ks == ls)) {
                bra = INDEX2(ii, jj) + pk_symoffset_[is];
                ket = INDEX2(kk, ll) + pk_symoffset_[ks];
                // _pk_symoffset corrects for the symmetry offset in the _pk vector
                braket = INDEX2(bra, ket);
                pk_[braket] += value;
                // K/2 (2nd sort)
                if ((ii != jj) && (kk != ll)) {
                    if ((is == ls) && (js == ks)) {
                        bra = INDEX2(ii, ll) + pk_symoffset_[is];
                        ket = INDEX2(jj, kk) + pk_symoffset_[js];
                        braket = INDEX2(bra, ket);
                        if ((ii == ll) || (jj == kk))
                            pk_[braket] -= 0.5 * value;
                        else
                            pk_[braket] -= 0.25 * value;
                    }
                }
            }

            // K/2 (1st sort)
            if ((is == ks) && (js == ls)) {
                bra = INDEX2(ii, kk) + pk_symoffset_[is];
                ket = INDEX2(jj, ll) + pk_symoffset_[js];
                braket = INDEX2(bra, ket);
                if ((ii == kk) || (jj == ll))
                    pk_[braket] -= 0.5 * value;
                else
                    pk_[braket] -= 0.25 * value;
            }
            pk_counter++;
            counter++;

            if (pk_flag) {
                pk_counter = 0;
            }
        }

        if (!ilsti)
            ERIIN.fetch();
    } while (!ilsti);

    // Going out of scope will close the buffer
    // iwl_buf_close(&ERIIN, 1);

    // After stage two is complete, the elements of P must be halved for the case IJ=KL.
    for (size_t ij=0; ij < pk_pairs_; ++ij)
        pk_[INDEX2(ij,ij)] *= 0.5;

    fprintf(outfile, "  Processed %d two-electron integrals.\n\n", counter);
#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "_pk:\n");
        print_array(pk_, pk_pairs_, outfile);
    }
#endif

    fflush(outfile);

    ERIIN.close();

}

void RHF::form_G_from_PK()
{
    int nirreps = factory_->nirrep();
    int *opi = factory_->rowspi();
    size_t ij;
    double *D_vector = new double[pk_pairs_];
    int nthread = 1;
    #ifdef _OPENMP
    nthread = omp_get_max_threads();
    #endif
//    double **G_vector = block_matrix(nthread, pk_pairs_);

    G_->zero();
    memset(D_vector, 0, sizeof(double) * pk_pairs_);
    memset(&(G_vector_[0][0]), 0, sizeof(double) * nthread * pk_pairs_);

    ij=0;
    for (int h=0; h<nirreps; ++h) {
        for (int p=0; p<opi[h]; ++p) {
            for (int q=0; q<=p; ++q) {
                if (p != q)
                    D_vector[ij] = 2.0 * D_->get(h, p, q);
                else
                    D_vector[ij] = D_->get(h, p, q);
                ij++;
            }
        }
    }

    double G_pq,D_pq;
    double* D_rs;
    double* G_rs;
    double* PK_block = pk_;
    #pragma omp parallel for private (G_pq, D_pq, D_rs, G_rs) schedule (dynamic)
    for(int pq = 0; pq < pk_pairs_; ++pq){
        int threadid = 0;
        double *local_PK_block = PK_block;
        #ifdef _OPENMP
        threadid = omp_get_thread_num();
        local_PK_block += ioff[pq];
        #endif
        G_pq = 0.0;
        D_pq = D_vector[pq];
        D_rs = &D_vector[0];
        G_rs = &(G_vector_[threadid][0]);
        for(int rs = 0; rs <= pq; ++rs){
            G_pq += *local_PK_block * (*D_rs);
            *G_rs += *local_PK_block * D_pq;

            ++D_rs;
            ++G_rs;
            ++local_PK_block;
        }
        G_vector_[threadid][pq] += G_pq;
    }

    for (int i = 1; i < nthread; ++i) {
        for (int j=0; j<pk_pairs_; ++j)
            G_vector_[0][j] += G_vector_[i][j];
    }

    // Convert G to a matrix
    ij = 0;
    for(int h = 0; h < nirreps; ++h){
        for(int p = 0; p < opi[h]; ++p){
            for(int q = 0; q <= p; ++q){
                G_->set(h,p,q, 2.0 * G_vector_[0][ij]);
                G_->set(h,q,p, 2.0 * G_vector_[0][ij]);
                ij++;
            }
        }
    }

    delete[](D_vector);
}
#endif

void RHF::save_sapt_info()
{
    if (factory_->nirrep() != 1)
    {
        fprintf(outfile,"Must run in C1. Period.\n"); fflush(outfile);
        abort();
    }
    if (soccpi_[0] != 0)
    {
        fprintf(outfile,"Aren't we in RHF Here? Pair those electrons up cracker!\n"); fflush(outfile);
        abort();
    }

    fprintf(outfile,"\n  Saving SAPT %s file.\n",options_.get_str("SAPT").c_str());

    int fileno;
    char* body_type = (char*)malloc(400*sizeof(char));
    char* key_buffer = (char*)malloc(4000*sizeof(char));

    string type = options_.get_str("SAPT");
    if (type == "2-DIMER") {
        fileno = PSIF_SAPT_DIMER;
        sprintf(body_type,"Dimer");
    } else if (type == "2-MONOMER_A") {
        fileno = PSIF_SAPT_MONOMERA;
        sprintf(body_type,"Monomer");
    } else if (type == "2-MONOMER_B") {
        fileno = PSIF_SAPT_MONOMERB;
        sprintf(body_type,"Monomer");
    } else if (type == "3-TRIMER") {
        fileno = PSIF_3B_SAPT_TRIMER;
        sprintf(body_type,"Trimer");
    } else if (type == "3-DIMER_AB") {
        fileno = PSIF_3B_SAPT_DIMER_AB;
        sprintf(body_type,"Dimer");
    } else if (type == "3-DIMER_BC") {
        fileno = PSIF_3B_SAPT_DIMER_BC;
        sprintf(body_type,"Dimer");
    } else if (type == "3-DIMER_AC") {
        fileno = PSIF_3B_SAPT_DIMER_AC;
        sprintf(body_type,"Dimer");
    } else if (type == "3-MONOMER_A") {
        fileno = PSIF_3B_SAPT_MONOMER_A;
        sprintf(body_type,"Monomer");
    } else if (type == "3-MONOMER_B") {
        fileno = PSIF_3B_SAPT_MONOMER_B;
        sprintf(body_type,"Monomer");
    } else if (type == "3-MONOMER_C") {
        fileno = PSIF_3B_SAPT_MONOMER_C;
        sprintf(body_type,"Monomer");
    } else {
        throw std::domain_error("SAPT Output option invalid");
    }

    psio_->open(fileno,0);

    int sapt_nso = basisset_->nbf();
    int sapt_nmo = nmo_;
    int sapt_nocc = doccpi_[0];
    int sapt_nvir = sapt_nmo-sapt_nocc;
    int sapt_ne = 2*sapt_nocc;
    double sapt_E_HF = E_;
    double sapt_E_nuc = nuclearrep_;
    double *sapt_evals = epsilon_a_->to_block_vector();
    double **sapt_C = Ca_->to_block_matrix();
    //print_mat(sapt_C,sapt_nso,sapt_nso,outfile);
    SharedMatrix potential(factory_->create_matrix("Potential Integrals"));
    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
    shared_ptr<OneBodySOInt> V(integral.so_potential());
    V->compute(potential);
    double *sapt_V_ints = potential->to_lower_triangle();
    double *sapt_S_ints = S_->to_lower_triangle();

    int errcod;

    sprintf(key_buffer,"%s NSO",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_nso, sizeof(int));
    sprintf(key_buffer,"%s NMO",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_nmo, sizeof(int));
    sprintf(key_buffer,"%s NOCC",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_nocc, sizeof(int));
    sprintf(key_buffer,"%s NVIR",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_nvir, sizeof(int));
    sprintf(key_buffer,"%s Number of Electrons",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_ne,sizeof(int));
    sprintf(key_buffer,"%s HF Energy",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_E_HF,sizeof(double));
    sprintf(key_buffer,"%s Nuclear Repulsion Energy",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_E_nuc, sizeof(double));
    sprintf(key_buffer,"%s HF Eigenvalues",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &(sapt_evals[0]),sizeof(double)*sapt_nmo);
    sprintf(key_buffer,"%s Nuclear Attraction Integrals",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &(sapt_V_ints[0]), sizeof(double)*sapt_nso*(sapt_nso+1)/2);
    sprintf(key_buffer,"%s Overlap Integrals",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &(sapt_S_ints[0]), sizeof(double)*sapt_nso*(sapt_nso+1)/2);
    sprintf(key_buffer,"%s HF Coefficients",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &(sapt_C[0][0]),sizeof(double)*sapt_nmo*sapt_nso);

    psio_->close(fileno,1);

    delete[] sapt_evals;
    delete[] sapt_V_ints;
    delete[] sapt_S_ints;
    free_block(sapt_C);

    free(body_type);
    free(key_buffer);
}
}}
