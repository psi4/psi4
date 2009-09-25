#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/matrix.h>
#include "uhf.h"

extern FILE *outfile;

using namespace std;
using namespace psi;

UHF::UHF(PSIO *psio, Chkpt *chkpt) : HF(psio, chkpt)
{
    common_init();
}

UHF::UHF(Ref<PSIO> &psio, Ref<Chkpt> &chkpt) : HF(psio, chkpt)
{
    common_init();
}

UHF::~UHF()
{
	if (p_j_)
		delete[](p_j_);
	if (p_k_)
		delete[](p_k_);
}

void UHF::common_init()
{
    use_out_of_core_ = false;
	
	Fa_ = factory_.create_matrix("F alpha");
	Fb_ = factory_.create_matrix("F beta");
	Da_ = factory_.create_matrix("D alpha");
	Db_ = factory_.create_matrix("D beta");
	Dt_ = factory_.create_matrix("D total");
	Ca_ = factory_.create_matrix("C alpha");
	Cb_ = factory_.create_matrix("C beta");
	Ga_ = factory_.create_matrix("G alpha");
	Gb_ = factory_.create_matrix("G beta");
	
	p_j_ = NULL;
	p_k_ = NULL;
	
	ip_boolean(const_cast<char*>("OUT_OF_CORE"), &(use_out_of_core_), 0);
	
	fprintf(outfile, "  DIIS not implemented for UHF, yet.\n\n");
	
	allocate_PK();
}

double UHF::compute_energy()
{
	bool converged = false;
	int iter = 0;
	
	// Do the initial work to give the iterations a starting point
	form_H();
	// find_occupation(_H, _H);
	
	if (use_out_of_core_ == false)
		form_PK();
	
	form_Shalf();
	form_initialF();
	form_C();
	form_D();
	
	// Compute an initial energy using H and D
	E_ = compute_initial_E();
	
	do {
		iter++;
		
		Eold_ = E_;
		
		if (use_out_of_core_ == false)
			form_G_from_PK();
		else
			form_G();
		
		form_F();
		
		E_ = compute_E();
//		fprintf(outfile, "  @UHF iteration %3d energy: %20.14f    %20.14f %s\n", iter, _E, _E - _Eold, diis_iter == false ? " " : "DIIS");
		fprintf(outfile, "  @UHF iteration %3d energy: %20.14f    %20.14f\n", iter, E_, E_ - Eold_);
		fflush(outfile);
		
		form_C();
		//find_occupation(Fa_, Fb_);
		form_D();
		
		converged = test_convergency();
	} while (!converged && iter < maxiter_);
	
    // Return the final RHF energy
    if (converged) {
        fprintf(outfile, "\n  Energy converged.\n");
        save_information();
        return E_;
    }
    else {
        fprintf(outfile, "\n  Failed to converge.\n");
        return 0.0;
    }
}

void UHF::save_information()
{
    // Print the final docc vector
    char **temp2 = chkpt_->rd_irr_labs();
    int nso = chkpt_->rd_nso();
    
    fprintf(outfile, "\n  Final occupation vector = (");
    for (int h=0; h<factory_.nirreps(); ++h) {
        fprintf(outfile, "%2d %3s ", doccpi_[h], temp2[h]);
    }
    fprintf(outfile, ")\n");
    
    // Needed for a couple of places.
    RefMatrix eigvectora = factory_.create_matrix();
    RefVector eigvaluesa = factory_.create_vector();
    RefMatrix eigvectorb = factory_.create_matrix();
    RefVector eigvaluesb = factory_.create_vector();
    Fa_.diagonalize(eigvectora, eigvaluesa);
    Fb_.diagonalize(eigvectorb, eigvaluesb);
    
    int print_mos = false;
    ip_boolean(const_cast<char*>("PRINT_MOS"), &(print_mos), 0);
    if (print_mos) {
        fprintf(outfile, "\n  Alpha Molecular orbitals:\n");
        Ca_.eivprint(eigvaluesa);
        
        fprintf(outfile, "\n  Beta Molecular orbitals:\n");
        Cb_.eivprint(eigvaluesb);
    }
    
    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairsa, pairsb;
    for (int h=0; h<eigvaluesa.nirreps(); ++h) {
        for (int i=0; i<eigvaluesa.dimpi()[h]; ++i) {
            pairsa.push_back(make_pair(eigvaluesa.get(h, i), h));
            pairsb.push_back(make_pair(eigvaluesb.get(h, i), h));
        }
    }
    sort(pairsa.begin(),pairsa.end());
    sort(pairsb.begin(),pairsb.end());
        
    fprintf(outfile, "\n  Orbital energies (a.u.):\n    Alpha occupied\n      ");
    for (int i=1; i<=nalpha_; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairsa[i-1].first, temp2[pairsa[i-1].second]);
        if (i % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "\n    Alpha unoccupied\n      ");
    for (int i=nalpha_+1; i<=nso; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairsa[i-1].first, temp2[pairsa[i-1].second]);
        if ((i-nalpha_) % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");
    
    fprintf(outfile, "\n    Beta occupied\n      ");
    for (int i=1; i<=nbeta_; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairsb[i-1].first, temp2[pairsb[i-1].second]);
        if (i % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "\n    Beta unoccupied\n      ");
    for (int i=nalpha_+1; i<=nso; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairsb[i-1].first, temp2[pairsb[i-1].second]);
        if ((i-nbeta_) % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");
    for (int i=0; i<eigvaluesa.nirreps(); ++i)
        free(temp2[i]);
    free(temp2);
    
    int *vec = new int[eigvaluesa.nirreps()];
    for (int i=0; i<eigvaluesa.nirreps(); ++i)
        vec[i] = 0;
        
    chkpt_->wt_nmo(nso);
    chkpt_->wt_ref(1);        // UHF
    chkpt_->wt_etot(E_);
    chkpt_->wt_escf(E_);
    chkpt_->wt_eref(E_);
    chkpt_->wt_clsdpi(doccpi_);
    chkpt_->wt_orbspi(eigvaluesa.dimpi());
    chkpt_->wt_openpi(vec);
    chkpt_->wt_phase_check(0);
    
    // Figure out frozen core orbitals
    int nfzc = chkpt_->rd_nfzc();
    int nfzv = chkpt_->rd_nfzv();
    int *frzcpi = compute_fcpi(nfzc, eigvaluesa);
    int *frzvpi = compute_fvpi(nfzv, eigvaluesa);
    chkpt_->wt_frzcpi(frzcpi);
    chkpt_->wt_frzvpi(frzvpi);
    delete[](frzcpi);
    delete[](frzvpi);
    
    // TODO: Figure out what chkpt_wt_iopen means for UHF
    chkpt_->wt_iopen(0);
    
    // Write eigenvectors and eigenvalue to checkpoint 
    double *values = eigvaluesa.to_block_vector();
    chkpt_->wt_alpha_evals(values);
    free(values);
    double **vectors = Ca_.to_block_matrix();
    chkpt_->wt_alpha_scf(vectors);
    free_block(vectors);
    values = eigvaluesb.to_block_vector();
    chkpt_->wt_beta_evals(values);
    free(values);
    vectors = Cb_.to_block_matrix();
    chkpt_->wt_beta_scf(vectors);
    free_block(vectors);
}

bool UHF::test_convergency()
{
	double ediff = E_ - Eold_;
	if (fabs(ediff) < energy_threshold_)
		return true;
	else
		return false;
}

void UHF::allocate_PK() {
	// Figure out how many pair combinations yield A1 symmetry (done in above loop)
	//   num_pair_combinations_of_A1 = ioff[_opi[0]] + ioff[_opi[1]] + ioff[_opi[2]] + ...
	// Allocate memory for the PK matrix (using a vector)
	if (pk_size_ < (memory_ / sizeof(double) / 2)) {
		p_j_ = new double[pk_size_];
		p_k_ = new double[pk_size_];

		if (p_j_ == NULL || p_k_ == NULL) {
			fprintf(outfile, "  Insufficient free system memory for in-core PK implementation.\n");
			fprintf(outfile, "  Switching to out-of-core algorithm.\n");
			use_out_of_core_ = true;
		} else {
			// Zero out PK and K
			memset(p_j_, 0, pk_size_*sizeof(double));
			memset(p_k_, 0, pk_size_*sizeof(double));

			fprintf(outfile,
				"  Allocated %lu elements (%lu pairs) for PJ. (%5f MiB)\n",
				(unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
			fprintf(outfile,
				"  Allocated %lu elements (%lu pairs) for PK.  (%5f MiB)\n\n",
				(unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
		}
	} else {
		fprintf(outfile,
				"  Insufficient memory for in-core PK implementation.\n");
		fprintf(outfile,
				"  Would need %lu elements of double memory. (%5f MiB)\n",
				(unsigned long)pk_size_*2, pk_size_ * 8.0 / 1048576.0 * 2.0);
		fprintf(outfile, "  Switching to out-of-core algorithm.\n");
		use_out_of_core_ = true;
	}
}

void UHF::form_initialF()
{
	Fa_.copy(H_);
	Fb_.copy(H_);
	
	// Transform the Focks
	Fa_.transform(Shalf_);
	Fb_.transform(Shalf_);
	
#ifdef _DEBUG
	if (debug_) {
		fprintf(outfile, "Initial Fock alpha matrix:\n");
		Fa_.print(outfile);
		fprintf(outfile, "Initial Fock beta matrix:\n");
		Fb_.print(outfile);
	}
#endif
}

void UHF::form_F() {
	Fa_.copy(H_ + Ga_);
	Fb_.copy(H_ + Gb_);
	
#ifdef _DEBUG
	if (debug_) {
		Fa_.print(outfile);
		Fb_.print(outfile);
	}
#endif
}

void UHF::form_C()
{
	RefMatrix eigvec = factory_.create_matrix();
	RefVector eigval = factory_.create_vector();
	
	Fa_.transform(Shalf_);
	Fa_.diagonalize(eigvec, eigval);
	Ca_.gemm(false, false, 1.0, Shalf_, eigvec, 0.0);
	
	Fb_.transform(Shalf_);
	Fb_.diagonalize(eigvec, eigval);
	Cb_.gemm(false, false, 1.0, Shalf_, eigvec, 0.0);
	
#ifdef _DEBUG
	if (debug_) {
		Ca_.print(outfile);
		Cb_.print(outfile);
	}
#endif
}

void UHF::form_D()
{
	int h, i, j, m;
	int *opi = Da_.rowspi();
	int nirreps = Da_.nirreps();
	double val;
	for (h=0; h<nirreps; ++h) {
		for (i=0; i<opi[h]; ++i) {
			for (j=0; j<opi[h]; ++j) {
				val = 0.0;
				for (m=0; m<nalphapi_[h]; ++m)
					val += Ca_.get(h, i, m) * Ca_.get(h, j, m);
				Da_.set(h, i, j, val);
				
				val = 0.0;
				for (m=0; m<nbetapi_[h]; ++m)
					val += Cb_.get(h, i, m) * Cb_.get(h, j, m);
				Db_.set(h, i, j, val);
			}
		}
	}

	// Form total density
	Dt_.copy(Da_ + Db_);
	
#ifdef _DEBUG
	if (debug_) {
		Da_.print(outfile);
		Db_.print(outfile);
	}
#endif
}

double UHF::compute_initial_E()
{
	double Etotal = nuclearrep_ + 0.5 * (Dt_.vector_dot(H_));
    fprintf(outfile, "\n  Initial UHF energy: %20.14f\n\n", Etotal);
    fflush(outfile);
    return Etotal;
}

double UHF::compute_E()
{
	double DH  = Dt_.vector_dot(H_);
	double DFa = Da_.vector_dot(Fa_);
	double DFb = Db_.vector_dot(Fb_);
	double Eelec = 0.5 * (DH + DFa + DFb);
	fprintf(outfile, "electronic energy = %20.14f\n", Eelec);
	double Etotal = nuclearrep_ + Eelec;
	return Etotal;
}

void UHF::form_PK()
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
    double value;
    
    // PK zeroed out during allocation
    fprintf(outfile, "  Forming PJ and PK matrices.\n");
    fflush(outfile);
    
    IWL ERIIN(psio_.pointer(), PSIF_SO_TEI, 0.0, 1, 1);
    
    do {
        ilsti = ERIIN.last_buffer();
        nbuf  = ERIIN.buffer_count();
        
        fi = 0;
        for (idx=0; idx<nbuf; ++idx) {
        	if (ERIIN.labels()[fi] >= 0) {
        		i = ERIIN.labels()[fi];
        	}
        	else {
        		i = -ERIIN.labels()[fi];
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
                // pk_symoffset_ corrects for the symmetry offset in the _pk vector
                braket = INDEX2(bra, ket);
                p_j_[braket] += value;
                // K/2 (2nd sort)
                if ((ii != jj) && (kk != ll)) {
                    if ((is == ls) && (js == ks)) {
                        bra = INDEX2(ii, ll) + pk_symoffset_[is];
                        ket = INDEX2(jj, kk) + pk_symoffset_[js];
                        braket = INDEX2(bra, ket);
                        if ((ii == ll) || (jj == kk))
                            p_k_[braket] -= 0.5 * value;
                        else
                            p_k_[braket] -= 0.25 * value;
                    }
                }
            }
            
            // K/2 (1st sort)
            if ((is == ks) && (js == ls)) {
                bra = INDEX2(ii, kk) + pk_symoffset_[is];
                ket = INDEX2(jj, ll) + pk_symoffset_[js];
                braket = INDEX2(bra, ket);
                if ((ii == kk) || (jj == ll))
                    p_k_[braket] -= 0.5 * value;
                else
                    p_k_[braket] -= 0.25 * value;
            }
            counter++;
        }

        if (!ilsti)
            ERIIN.fetch();
    } while (!ilsti);

    // Going out of scope will close the buffer
    // iwl_buf_close(&ERIIN, 1);
    
    // After stage two is complete, the elements of P must be halved for the case IJ=KL.
    for (size_t ij=0; ij < pk_pairs_; ++ij) {
        p_j_[INDEX2(ij,ij)] *= 0.5;
        p_k_[INDEX2(ij,ij)] *= 0.5;
    }
    
    fprintf(outfile, "  Processed %d two-electron integrals.\n", counter);
    #ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "p_j_:\n");
        print_array(p_j_, pk_pairs_, outfile);
        fprintf(outfile, "p_k_:\n");
        print_array(p_k_, pk_pairs_, outfile);
    }
    #endif
}

void UHF::form_G_from_PK()
{
	int nirreps = factory_.nirreps();
	int *opi = factory_.rowspi();
	size_t ij;
	double *Da_vector = new double[pk_pairs_];
	double *Db_vector = new double[pk_pairs_];
	double *Dt_vector = new double[pk_pairs_];
	double *Ga_vector = new double[pk_pairs_];
	double *Gb_vector = new double[pk_pairs_];

	Ga_.zero();
	Gb_.zero();

	memset(Da_vector, 0, sizeof(double) * pk_pairs_);
	memset(Db_vector, 0, sizeof(double) * pk_pairs_);
	memset(Dt_vector, 0, sizeof(double) * pk_pairs_);
	memset(Ga_vector, 0, sizeof(double) * pk_pairs_);
	memset(Gb_vector, 0, sizeof(double) * pk_pairs_);

	ij=0;
	for (int h=0; h<nirreps; ++h) {
		for (int p=0; p<opi[h]; ++p) {
			for (int q=0; q<=p; ++q) {
				if (p != q) {
					Da_vector[ij] = 2.0 * Da_.get(h, p, q);
					Db_vector[ij] = 2.0 * Db_.get(h, p, q);
					Dt_vector[ij] = 2.0 * Dt_.get(h, p, q);
				} else {
					Da_vector[ij] = Da_.get(h, p, q);
					Db_vector[ij] = Db_.get(h, p, q);
					Dt_vector[ij] = Dt_.get(h, p, q);
				}
				ij++;
			}
		}
	}

#ifdef _DEBUG
	if (debug_) {
		fprintf(outfile, "PK: ij = %lu\n", (unsigned long)ij);
		fflush(outfile);
		fprintf(outfile, "PK: Da matrix:\n");
		Da_.print(outfile);
		fprintf(outfile, "PK: Da vector (appears to be OK):\n");
		for (ij=0; ij<pk_pairs_; ++ij)
			fprintf(outfile, "PK: Da vector [%lu] = %20.14f\n", (unsigned long)ij, Da_vector[ij]);
		fprintf(outfile, "PK: Db matrix:\n");
		Db_.print(outfile);
		fprintf(outfile, "PK: Db vector (appears to be OK):\n");
		for (ij=0; ij<pk_pairs_; ++ij)
			fprintf(outfile, "PK: Db vector [%lu] = %20.14f\n", (unsigned long)ij, Db_vector[ij]);
	}
#endif

	/* 
	 * This code goes through the densities (Da_ and Db_), J, and K to form
	 * two G matrices. One G matrix is for Fa_ and the other for Fb_.
	 * See derivation notebook for equations.
	 */
	double Ga_pq, Da_pq;
	double Gb_pq, Db_pq;
	double Dt_pq;
	double* Da_rs;
	double* Ga_rs;
	double* Db_rs;
	double* Gb_rs;
	double* Dt_rs;
	int pq, rs;
	double* J_block = p_j_;
	double* K_block = p_k_;
	int ts_pairs = pk_pairs_;
	for (pq = 0; pq < ts_pairs; ++pq) {
		Ga_pq = 0.0;
		Da_pq = Da_vector[pq];
		Da_rs = &Da_vector[0];
		Ga_rs = &Ga_vector[0];
		Gb_pq = 0.0;
		Db_pq = Db_vector[pq];
		Db_rs = &Db_vector[0];
		Gb_rs = &Gb_vector[0];
		Dt_pq = Dt_vector[pq];
		Dt_rs = &Dt_vector[0];
		for (rs = 0; rs <= pq; ++rs) {
			// D_{rs}^{c} * PK_{pqrs}         Also found in RHF
			// Doing F_mn_a about to add the K term
			Ga_pq  += *J_block * (*Dt_rs) + *K_block * (*Da_rs);
			*Ga_rs += *J_block * Dt_pq + *K_block * Da_pq;
			
			Gb_pq  += *J_block * (*Dt_rs) + *K_block * (*Db_rs);
			*Gb_rs += *J_block * Dt_pq + *K_block * Db_pq;
						
			++Da_rs;
			++Ga_rs;
			++Db_rs;
			++Gb_rs;
			++Dt_rs;
			++J_block;
			++K_block;
		}
		Ga_vector[pq] += Ga_pq;
		Gb_vector[pq] += Gb_pq;
	}

	// Convert G to a matrix
	ij = 0;
	for (int h = 0; h < nirreps; ++h) {
		for (int p = 0; p < opi[h]; ++p) {
			for (int q = 0; q <= p; ++q) {
				Ga_.set(h, p, q, 2.0 * Ga_vector[ij]);
				Ga_.set(h, q, p, 2.0 * Ga_vector[ij]);
				Gb_.set(h, p, q, 2.0 * Gb_vector[ij]);
				Gb_.set(h, q, p, 2.0 * Gb_vector[ij]);
				ij++;
			}
		}
	}

#ifdef _DEBUG
	if (debug_) {
		Ga_.print(outfile);
		Gb_.print(outfile);
	}
#endif

	delete[](Da_vector);
	delete[](Db_vector);
	delete[](Dt_vector);
	delete[](Ga_vector);
	delete[](Gb_vector);
}

void UHF::form_G()
{
	fprintf(stderr, "UHF out-of-core algorithm is not implemented yet!\n");
	abort();
}