/*
 *@BEGIN LICENSE
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

#include "infsapt.h"
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <lib3index/3index.h>
#include <libfock/jk.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include <physconst.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {
namespace dftsapt {

INFSAPT::INFSAPT()
{
    common_init();
}
INFSAPT::~INFSAPT()
{
}
void INFSAPT::common_init()
{
    print_ = 1;
    debug_ = 0;
    bench_ = 0;
}
boost::shared_ptr<INFSAPT> INFSAPT::build(boost::shared_ptr<Wavefunction> d,
                                          std::vector<boost::shared_ptr<Wavefunction> > m)
{
    INFSAPT* sapt = new INFSAPT();

    Options& options = Process::environment.options; 

    sapt->print_ = options.get_int("PRINT");
    sapt->debug_ = options.get_int("DEBUG");
    sapt->bench_ = options.get_int("BENCH");

    sapt->memory_ = (unsigned long int)(Process::environment.get_memory() * options.get_double("SAPT_MEM_FACTOR") * 0.125);

    sapt->schwarz_ = options.get_double("INTS_TOLERANCE");

    sapt->cluster_ = d;
    sapt->monomers_ = m;

    sapt->primary_  = d->basisset();
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    sapt->mp2fit_ = BasisSet::construct(parser, sapt->cluster_->molecule(), "DF_BASIS_SAPT");

    return boost::shared_ptr<INFSAPT>(sapt);
}
double INFSAPT::compute_energy()
{
    energies_["Total"] = 0.0; 
    energies_["Elst"]  = 0.0;
    energies_["Exch"]  = 0.0;
    energies_["Ind"]   = 0.0;
    energies_["Disp"]  = 0.0;

    print_header();

    scf_terms();
    
    pt2_terms();

    print_trailer();
    
    return 0.0;
}
void INFSAPT::print_header() const 
{
    outfile->Printf( "\t --------------------------------------------------------\n");
    outfile->Printf( "\t                         INF-SAPT0                       \n");
    outfile->Printf( "\t               Rob Parrish and Ed Hohenstein             \n");
    outfile->Printf( "\t --------------------------------------------------------\n");
    outfile->Printf( "\n");

    outfile->Printf( "  ==> Sizes <==\n");
    outfile->Printf( "\n");

    outfile->Printf( "   => Parameters <=\n\n");

    int threads = 1;
    #ifdef _OPENMP
        threads = omp_get_max_threads();
    #endif
    
    outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
    outfile->Printf( "    Threads:           %11d\n", threads);
    outfile->Printf( "    Schwarz Cutoff:    %11.3E\n", schwarz_);
    outfile->Printf( "\n");

    outfile->Printf("   => Molecular Cluster <=\n\n");
    cluster_->molecule()->print_cluster();

    outfile->Printf( "   => Orbital Ranges (Pauli Blockade Space) <=\n\n");

    int nso = cluster_->nso();
    int nmo = cluster_->nmo();
    int nocc = cluster_->doccpi().sum();
    int nvir = nmo - nocc;
    int nfocc = cluster_->frzcpi().sum();
    int nfvir = cluster_->frzvpi().sum();
    int naocc = nocc - nfocc;
    int navir = nvir - nfvir;

    outfile->Printf("    ---------------------------------------------------------------\n");
    outfile->Printf("    %-7s %6s %6s %6s %6s %6s %6s %6s %6s\n", 
        "Unit", "Nso", "Nmo", "Nocc", "Nvir", "Nfocc", "Naocc", "Navir", "Nfvir");

    outfile->Printf("    ---------------------------------------------------------------\n");
    outfile->Printf("    %-7s %6d %6d %6d %6d %6d %6d %6d %6d\n", 
        "Cluster", nso, nmo, nocc, nvir, nfocc, naocc, navir, nfvir);
    
    for (int A = 0; A < monomers_.size(); A++) {
        boost::shared_ptr<Wavefunction> m = monomers_[A];
        nocc = m->doccpi().sum();
        nfocc = m->frzcpi().sum();
        naocc = nocc - nfocc;
        
        outfile->Printf("    %-7d %6d %6d %6d %6d %6d %6d %6d %6d\n", 
            A+1, nso, nmo, nocc, nvir, nfocc, naocc, navir, nfvir);
    }
    outfile->Printf("    ---------------------------------------------------------------\n");
    outfile->Printf( "\n");

    outfile->Printf( "   => Primary Cluster Basis Set <=\n\n");
    primary_->print_by_level("outfile", print_);

    
}
void INFSAPT::scf_terms()
{
    outfile->Printf( "  SCF TERMS:\n\n");

    // ==> Setup <== //

    // Number of bodies in this job
    int N = monomers_.size();
    // Number of AO basis functions
    int nso = primary_->nbf();

    // => PB Object <= //
    
    boost::shared_ptr<PB> pb = PB::build(monomers_);
    pb->initialize();
    pb->print_header();
    
    boost::shared_ptr<JK> jk = pb->jk();

    // D matrices
    const std::vector<SharedMatrix>& D = jk->D();
    // J matrices
    const std::vector<SharedMatrix>& J = jk->J();
    // K matrices
    const std::vector<SharedMatrix>& K = jk->K();
    // V matrices
    const std::vector<SharedMatrix>& V = pb->V();

    // ==> (HF) Terms <== //

    // E_HF^(0)
    double E_HF = 0.0;
    for (int A = 0; A < monomers_.size(); A++) {
        E_HF += pb->E_HF_0()->get(A);
    }
    
    // ==> (0) Terms <== //

    // Total Elst^(10)
    double E_elst = 0.0;
    // Elst^(10)
    boost::shared_ptr<Matrix> EJ_0(new Matrix("EJ_0", N, N)); 
    // Tunneling term (Heitler-London)
    boost::shared_ptr<Matrix> EK_0(new Matrix("EK_0", N, N)); 
    // Total Heitler-London term
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double EJ = 0.0;
            double EK = 0.0;
            EJ += 2.0 * D[A]->vector_dot(V[B]);
            EJ += 2.0 * D[B]->vector_dot(V[A]);
            EJ += 4.0 * D[A]->vector_dot(J[B]);
            EJ += 1.0 * monomers_[A]->molecule()->pairwise_nuclear_repulsion_energy(monomers_[B]->molecule());  
            EK -= 2.0 * D[A]->vector_dot(K[B]); 
            EJ_0->set(A,B,EJ);
            EK_0->set(A,B,EK);
            E_elst += EJ;
        }
    } 

    //EJ_0->print();
    //EK_0->print(); 
    
    // ==> (A) Terms <== //

    // Run the Pauli Blockade
    pb->compute();

    // Tunneling term
    double E_tunnel = 0.0;
    // Pauli term
    double E_Pauli = 0.0;
    // Cross term
    double E_cross = 0.0;
    // Total Exch
    double E_exch = 0.0;
    // Pauli Repulsion term (diagonal) and rearrangement term (off-diagonal)
    boost::shared_ptr<Matrix> EK_A(new Matrix("EK_A", N, N)); 
    // Pauli Repulsion term
    for (int A = 0; A < monomers_.size(); A++) {
        double EK = pb->E_HF_A()->get(A) - pb->E_HF_0()->get(A);
        EK_A->set(A,A,EK);
        E_exch += EK;
        E_Pauli += EK;
    }
    // Rearrangement term
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double EK = 0.0;
            EK += 2.0 * D[A]->vector_dot(V[B]);
            EK += 2.0 * D[B]->vector_dot(V[A]);
            EK += 4.0 * D[A]->vector_dot(J[B]);
            EK += 1.0 * monomers_[A]->molecule()->pairwise_nuclear_repulsion_energy(monomers_[B]->molecule());  
            EK -= 2.0 * D[A]->vector_dot(K[B]); 
            EK -= EK_0->get(A,B);
            EK -= EJ_0->get(A,B);
            EK_A->set(A,B,EK);
            E_exch += EK;
            E_exch += EK_0->get(A,B);
            E_cross  += EK; 
            E_tunnel += EK_0->get(A,B);
        }
    } 

    //EK_A->print(); 

    // ==> (I) Terms <== //

    double E_ind = cluster_->reference_energy() - E_HF - E_elst - E_exch; 

    // ==> Printing <== //

    outfile->Printf( "  ==> SCF Analysis <==\n\n");
    
    outfile->Printf( "   => Electrostatics <=\n\n");

    outfile->Printf("    Total:\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EJ_0->get(A,B);
            outfile->Printf( "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_elst, pc_hartree2kcalmol*E_elst);
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");

    outfile->Printf("   => Exchange <=\n\n");

    outfile->Printf("    Pauli Term:\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
            double Eval = EK_A->get(A,A);
            outfile->Printf( "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, A+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
    }    
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf( "    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_Pauli, pc_hartree2kcalmol*E_Pauli);
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");

    outfile->Printf("    Tunneling Term:\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EK_0->get(A,B);
            outfile->Printf( "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf( "    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_tunnel, pc_hartree2kcalmol*E_tunnel);
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");

    outfile->Printf("    Rearrangement Term:\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EK_A->get(A,B);
            outfile->Printf( "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_cross, pc_hartree2kcalmol*E_cross);
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");

    outfile->Printf("    Total:\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A; B < monomers_.size(); B++) {
            double Eval = EK_0->get(A,B) + EK_A->get(A,B);
            outfile->Printf( "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_exch, pc_hartree2kcalmol*E_exch);
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");

    outfile->Printf( "   => Induction <=\n\n");

    outfile->Printf("    Total:\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_ind, pc_hartree2kcalmol*E_ind);
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");

    energies_["Elst"] = E_elst;
    energies_["Exch"] = E_exch;
    energies_["Ind"] = E_ind;

    

    // => Setup PT2 Terms <= //

    std::vector<boost::shared_ptr<Matrix> > Co = pb->Cocc(); 
    std::vector<boost::shared_ptr<Matrix> > Cv = pb->Cvir(); 
    std::vector<boost::shared_ptr<Vector> > eo = pb->eps_occ(); 
    std::vector<boost::shared_ptr<Vector> > ev = pb->eps_vir(); 

    pb.reset();

    Caocc_.clear();
    Cavir_.clear();
    eps_aocc_.clear();
    eps_avir_.clear();
    for (int A = 0; A < monomers_.size(); A++) {
        int nvir = cluster_->nmo() - cluster_->doccpi().sum();
        int nfvir = cluster_->frzvpi().sum();
        int navir = nvir - nfvir;
        
        boost::shared_ptr<Wavefunction> m = monomers_[A];
        int nocc = m->doccpi().sum();
        int nfocc = m->frzcpi().sum();
        int naocc = nocc - nfocc;

        Caocc_.push_back(boost::shared_ptr<Matrix>(new Matrix("Caocc",nso,naocc)));
        for (int a = 0; a < naocc; a++) {
            C_DCOPY(nso,&Co[A]->pointer()[0][a + nfocc],nocc,&Caocc_[A]->pointer()[0][a],naocc);
        }

        Cavir_.push_back(boost::shared_ptr<Matrix>(new Matrix("Cavir",nso,navir)));
        for (int a = 0; a < navir; a++) {
            C_DCOPY(nso,&Cv[A]->pointer()[0][a],nvir,&Cavir_[A]->pointer()[0][a],navir);
        }

        eps_aocc_.push_back(boost::shared_ptr<Vector>(new Vector("eps_aocc",naocc)));
        C_DCOPY(naocc,&eo[A]->pointer()[nfocc],1,eps_aocc_[A]->pointer(),1);

        eps_avir_.push_back(boost::shared_ptr<Vector>(new Vector("eps_avir",navir)));
        C_DCOPY(navir,&ev[A]->pointer()[0],1,eps_avir_[A]->pointer(),1);
    }    
}
void INFSAPT::pt2_terms()
{
    outfile->Printf( "  PT2 TERMS:\n\n");

    outfile->Printf( "   => MP2FIT Auxiliary Basis Set <=\n\n");
    mp2fit_->print_by_level("outfile", print_);

    

    // ==> Sizing <== //

    int N = monomers_.size();
    int nso   = cluster_->nso();
    int naux  = mp2fit_->nbf();
    int naocc = cluster_->doccpi().sum() - cluster_->frzcpi().sum();
    int navir = Cavir_[0]->colspi()[0];
    std::vector<int> aocc_starts;
    aocc_starts.push_back(0);
    for (int A = 0; A < N; A++) {
        aocc_starts.push_back(Caocc_[A]->colspi()[0] + aocc_starts[A]);
    }

    // ==> Integral Generation <== //
    
    build_iaQ(primary_,mp2fit_,Caocc_,Cavir_);

    // ==> PT2 computation <== //

    // => Targets <= //

    // Uncoupled Casimir-Polder Dispersion
    double EJ = 0.0;
    boost::shared_ptr<Matrix> EJ_PT2(new Matrix("EJ_PT2", N, N)); 
    // Tunneling Exchange-Dispersion
    double EK = 0.0;
    boost::shared_ptr<Matrix> EK_PT2(new Matrix("EK_PT2", N, N)); 

    // => Memory sizing <= //

    // Threads need some overhead
    int threads = 1;
    #ifdef _OPENMP
        threads = omp_get_max_threads();
    #endif

    // Largest block of occupied orbitals at a time
    size_t doubles = memory_;
    doubles -= threads * (size_t) navir * navir;
    size_t size_per_row = 2L * navir * naux;
    size_t max_val = (doubles / size_per_row);

    int max_naocc = 0;
    for (int A = 0; A < N; A++) {
        int n = Caocc_[A]->colspi()[0];
        max_naocc = (n > max_naocc ? n : max_naocc);
    }

    int max_block = (max_val > max_naocc ? max_naocc : max_val);  
    max_block = (max_block < 1 ? 1 : max_block);

    // => Tensor Blocks <= //
    
    // Disk blocks of (ar|Q) integrals
    boost::shared_ptr<Matrix> arQ(new Matrix("(ar|Q)", max_block, navir*naux));
    boost::shared_ptr<Matrix> bsQ(new Matrix("(bs|Q)", max_block, navir*naux));
    double** arQp = arQ->pointer();
    double** bsQp = bsQ->pointer();

    // Thread-level blocks for GEMM to integrals
    std::vector<boost::shared_ptr<Matrix> > I;
    std::vector<double**> Ipointers;
    for (int thread = 0; thread < threads; thread++) {
        I.push_back(boost::shared_ptr<Matrix>(new Matrix("I", navir, navir)));
        Ipointers.push_back(I[thread]->pointer());
    } 

    boost::shared_ptr<PSIO> psio(new PSIO()); 
    psio_address next_AIA = PSIO_ZERO;

    psio->open(PSIF_DFMP2_AIA,PSIO_OPEN_OLD);

    for (int A = 0; A < N-1; A++) {

        int astart = aocc_starts[A];
        int astop  = aocc_starts[A+1];
        int adelta = astop - astart;

        double* eap = eps_aocc_[A]->pointer();
        double* erp = eps_avir_[A]->pointer();

        for (int a = 0; a < adelta; a += max_block) {

            int na = (a + max_block < adelta ? max_block : adelta - a); 
            next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(a + astart)*navir*naux);
            psio->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)arQp[0],sizeof(double)*na*navir*naux,next_AIA,&next_AIA);

            for (int B = A + 1; B < N; B++) {

                int bstart = aocc_starts[B];
                int bstop  = aocc_starts[B+1];
                int bdelta = bstop - bstart;
                
                double* ebp = eps_aocc_[B]->pointer();
                double* esp = eps_avir_[B]->pointer();

                for (int b = 0; b < bdelta; b += max_block) {

                    int nb = (b + max_block < bdelta ? max_block : bdelta - b); 
                    next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(b + bstart)*navir*naux);
                    psio->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)bsQp[0],sizeof(double)*nb*navir*naux,next_AIA,&next_AIA);

                    int nab = na * nb;

                    double EJval = 0.0;
                    double EKval = 0.0;

                    #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: EJval, EKval)
                    for (int ab = 0; ab < nab; ab++) {

                        int aval = ab / nb;
                        int bval = ab % nb; 

                        int thread = 0; 
                        #ifdef _OPENMP
                            thread = omp_get_thread_num();
                        #endif
                        double** Ip = Ipointers[thread];

                        C_DGEMM('N','T',navir,navir,naux,1.0,arQp[aval],naux,bsQp[bval],naux,0.0,Ip[0],navir);

                        double ea = eap[aval];
                        double eb = ebp[bval];

                        for (int r = 0; r < navir; r++) {
                            double er = erp[r];
                            for (int s = 0; s < navir; s++) {
                                double es = esp[s];
                                double Irs = Ip[r][s];    
                                double Isr = Ip[s][r];    
                                double denom = 1.0 / (er + es - ea - eb);
                                EJval += -4.0 * Irs * Irs * denom; 
                                EKval +=  2.0 * Irs * Isr * denom; 
                            }
                        }
                        
                    }

                    EJ_PT2->set(A,B,EJval + EJ_PT2->get(A,B));
                    EK_PT2->set(A,B,EKval + EK_PT2->get(A,B));

                    EJ += EJval;
                    EK += EKval;

                }
            }
        }
    }

    psio->close(PSIF_DFMP2_AIA,0);

    // ==> Printing <== //

    outfile->Printf( "  ==> PT2 Analysis <==\n\n");
    
    outfile->Printf("    Dispersion (Uncoupled, Casimir-Polder):\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EJ_PT2->get(A,B);
            outfile->Printf( "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*EJ, pc_hartree2kcalmol*EJ);
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");
    
    outfile->Printf("    Exchange-Dispersion (Tunneling):\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EK_PT2->get(A,B);
            outfile->Printf( "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*EK, pc_hartree2kcalmol*EK);
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");
    
    outfile->Printf("    Total:\n\n");
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    outfile->Printf("    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EJ_PT2->get(A,B) + EK_PT2->get(A,B);
            outfile->Printf( "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*(EJ+EK), pc_hartree2kcalmol*(EJ+EK));
    outfile->Printf("    -------------------------------------------------\n");
    outfile->Printf("\n");

    energies_["Disp"] = EJ + EK;

    
}
void INFSAPT::print_trailer()
{
    outfile->Printf( "  ALL TERMS:\n\n");

    energies_["Total"] = energies_["Elst"] + energies_["Exch"] + energies_["Ind"] + energies_["Disp"];
    
    outfile->Printf("    ----------------------------------------------------\n");
    outfile->Printf("    %-14s %18s %18s\n",
        "Term", "mH", "kcal mol^-1");
    outfile->Printf("    ----------------------------------------------------\n");
    outfile->Printf("    %-14s %18.12f %18.12f\n",
        "Electrostatics", 1.0E3*energies_["Elst"], pc_hartree2kcalmol*energies_["Elst"]);
    outfile->Printf("    %-14s %18.12f %18.12f\n",
        "Exchange", 1.0E3*energies_["Exch"], pc_hartree2kcalmol*energies_["Exch"]);
    outfile->Printf("    %-14s %18.12f %18.12f\n",
        "Induction", 1.0E3*energies_["Ind"], pc_hartree2kcalmol*energies_["Ind"]);
    outfile->Printf("    %-14s %18.12f %18.12f\n",
        "Dispersion", 1.0E3*energies_["Disp"], pc_hartree2kcalmol*energies_["Disp"]);
    outfile->Printf("    ----------------------------------------------------\n");
    outfile->Printf("    %-14s %18.12f %18.12f\n",
        "Total", 1.0E3*energies_["Total"], pc_hartree2kcalmol*energies_["Total"]);
    outfile->Printf("    ----------------------------------------------------\n");
    outfile->Printf("\n");

    outfile->Printf( "    \"You can observe a lot by just watching.\"\n");
    outfile->Printf( "                      --Yogi Berra\n");
    outfile->Printf( "\n");

    
}

void INFSAPT::build_iaQ(
        boost::shared_ptr<BasisSet> primary,
        boost::shared_ptr<BasisSet> auxiliary,
        std::vector<boost::shared_ptr<Matrix> >& Caocc, 
        std::vector<boost::shared_ptr<Matrix> >& Cavir)
{
    // PSIO object
    boost::shared_ptr<PSIO> psio(new PSIO()); 

    boost::shared_ptr<Matrix> Caocc2 = Matrix::horzcat(Caocc);

    // ==> Target <== //

    psio->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);

    // => Global Sizing <= //

    int nso = primary->nbf();
    int naux = auxiliary->nbf();
    int naocc = Caocc2->colspi()[0];
    int navir = Cavir[0]->colspi()[0];

    // ==> (A|ia) <== //

    { // Begin Aia

    // Schwarz Sieve 
    boost::shared_ptr<ERISieve> sieve(new ERISieve(primary,schwarz_));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // ERI objects
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary,BasisSet::zero_ao_basis_set(),
        primary,primary));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int thread = 0; thread < nthread; thread++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri()));
        buffer.push_back(eri[thread]->buffer());
    }

    // Sizing
    int maxQ = auxiliary->max_function_per_shell();

    // Max block size in naux
    ULI Amn_cost_per_row = nso * (ULI) nso;
    ULI Ami_cost_per_row = nso * (ULI) naocc;
    ULI Aia_cost_per_row = naocc * (ULI) navir;
    ULI total_cost_per_row = Amn_cost_per_row + Ami_cost_per_row + Aia_cost_per_row;
    ULI doubles = memory_;
    ULI max_temp = doubles / (total_cost_per_row);
    int max_naux = (max_temp > (ULI) naux ? naux : max_temp);
    max_naux = (max_naux < maxQ ? maxQ : max_naux);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < auxiliary->nshell(); Q++) {
        int nQ = auxiliary->shell(Q).nfunction();
        if (counter + nQ > max_naux) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(auxiliary->nshell());

    // Tensor blocks
    SharedMatrix Amn(new Matrix("(A|mn) Block", max_naux, nso * (ULI) nso));
    SharedMatrix Ami(new Matrix("(A|mi) Block", max_naux, nso * (ULI) naocc));
    SharedMatrix Aia(new Matrix("(A|ia) Block", max_naux, naocc * (ULI) navir));
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Aiap = Aia->pointer();

    // C Matrices
    double** Caoccp = Caocc2->pointer();

    psio_address next_AIA = PSIO_ZERO;

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++) {

        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop  = block_Q_starts[block+1];
        int qoff   = auxiliary->shell(Qstart).function_index();
        int nrows  = (Qstop == auxiliary->nshell() ?
                     auxiliary->nbf() -
                     auxiliary->shell(Qstart).function_index() :
                     auxiliary->shell(Qstop).function_index() -
                     auxiliary->shell(Qstart).function_index());

        // Clear Amn for Schwarz sieve
        ::memset((void*) Amnp[0], '\0', sizeof(double) * nrows * nso * nso);

        // Compute TEI tensor block (A|mn)
        timer_on("INFSAPT (A|mn)");
        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int QMN = 0L; QMN < (Qstop - Qstart) * (ULI) npairs; QMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int Q =  QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = auxiliary->shell(Q).nfunction();
            int nm = primary->shell(M).nfunction();
            int nn = primary->shell(N).nfunction();

            int sq =  auxiliary->shell(Q).function_index();
            int sm =  primary->shell(M).function_index();
            int sn =  primary->shell(N).function_index();

            eri[thread]->compute_shell(Q,0,M,N);

            for (int oq = 0; oq < nq; oq++) {
                for (int om = 0; om < nm; om++) {
                    for (int on = 0; on < nn; on++) {
                        Amnp[sq + oq - qoff][(om + sm) * nso + (on + sn)] =
                        Amnp[sq + oq - qoff][(on + sn) * nso + (om + sm)] =
                        buffer[thread][oq * nm * nn + om * nn + on];
                    }
                }
            }
        }
        timer_off("INFSAPT (A|mn)");

        // Compute (A|mi) tensor block (A|mn) C_ni
        timer_on("INFSAPT (A|mn)C_mi");
        C_DGEMM('N','N',nrows*(ULI)nso,naocc,nso,1.0,Amnp[0],nso,Caoccp[0],naocc,0.0,Amip[0],naocc);
        timer_off("INFSAPT (A|mn)C_mi");

        // Compute (A|ia) tensor block (A|ia) = (A|mi) C_ma
        timer_on("INFSAPT (A|mi)C_na");
        #pragma omp parallel for
        for (int row = 0; row < nrows; row++) {
            double* Aia2p = Aiap[row];
            int aoff = 0;
            for (int A = 0; A < Caocc_.size(); A++) {
                double** Cap = Cavir_[A]->pointer();
                int na = Caocc_[A]->colspi()[0];
                C_DGEMM('T','N',na,navir,nso,1.0,&Amip[row][aoff],naocc,Cap[0],navir,0.0,Aia2p,navir);
                Aia2p += na * navir;
                aoff += na;
            }
        }
        timer_off("INFSAPT (A|mi)C_na");

        // Stripe (A|ia) out to disk
        timer_on("INFSAPT Aia Write");
        psio->write(PSIF_DFMP2_AIA,"(A|ia)",(char*)Aiap[0],sizeof(double)*nrows*naocc*navir,next_AIA,&next_AIA);
        timer_off("INFSAPT Aia Write");
    }

    } // End Aia

    // ==> Fitting Metric <== //

    // Fitting metric is needed below
    SharedMatrix Jm12;
    { // Begin Metric

    // Form the inverse metric manually
    timer_on("INFSAPT Metric");
    boost::shared_ptr<FittingMetric> metric(new FittingMetric(auxiliary, true));
    metric->form_eig_inverse(1.0E-10);
    Jm12 = metric->get_metric();
    timer_off("INFSAPT Metric");

    } // End Metric

    // ==> (ia|Q) <== //

    { // Begin iaQ

    // Number of columns in the (A|ia) disk tensor
    ULI nia = naocc * navir; 

    // Memory constraints
    ULI Jmem = naux * naux;
    ULI doubles = memory_; 
    if (doubles < 2L * Jmem) {
        throw PSIEXCEPTION("INFSAPT: More memory required for tractable disk transpose");
    }
    ULI rem = (doubles - Jmem) / 2L;
    ULI max_nia = (rem / naux);
    max_nia = (max_nia > nia ? nia : max_nia);
    max_nia = (max_nia < 1L ? 1L : max_nia);

    // Block sizing
    std::vector<ULI> ia_starts;
    ia_starts.push_back(0);
    for (ULI ia = 0L; ia < nia; ia+=max_nia) {
        if (ia + max_nia >= nia) {
            ia_starts.push_back(nia);
        } else {
            ia_starts.push_back(ia + max_nia);
        }
    }

    // Tensor blocks
    SharedMatrix Aia(new Matrix("Aia", naux, max_nia));
    SharedMatrix Qia(new Matrix("Qia", max_nia, naux));
    double** Aiap = Aia->pointer();
    double** Qiap = Qia->pointer();
    double** Jp   = Jm12->pointer();

    // Loop through blocks
    psio_address next_AIA = PSIO_ZERO;
    psio_address next_QIA = PSIO_ZERO;
    for (int block = 0; block < ia_starts.size() - 1; block++) {

        // Sizing
        ULI ia_start = ia_starts[block];
        ULI ia_stop  = ia_starts[block+1];
        ULI ncols = ia_stop - ia_start;

        // Read Aia
        timer_on("INFSAPT Aia Read");
        for (ULI Q = 0; Q < naux; Q++) {
            next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(Q*nia+ia_start));
            psio->read(PSIF_DFMP2_AIA,"(A|ia)",(char*)Aiap[Q],sizeof(double)*ncols,next_AIA,&next_AIA);
        }
        timer_off("INFSAPT Aia Read");

        // Apply Fitting
        timer_on("INFSAPT (Q|A)(A|ia)");
        C_DGEMM('T','N',ncols,naux,naux,1.0,Aiap[0],max_nia,Jp[0],naux,0.0,Qiap[0],naux);
        timer_off("INFSAPT (Q|A)(A|ia)");

        // Write Qia
        timer_on("INFSAPT Qia Write");
        psio->write(PSIF_DFMP2_AIA,"(Q|ia)",(char*)Qiap[0],sizeof(double)*ncols*naux,next_QIA,&next_QIA);
        timer_off("INFSAPT Qia Write");

    }

    } // End iaQ

    // ==> Cleanup <== //

    psio->close(PSIF_DFMP2_AIA,1);
}

PB::PB()
{
    common_init();
}
PB::~PB()
{
}
void PB::common_init()
{
    print_ = 1;
    debug_ = 0;
    bench_ = 0;
}
boost::shared_ptr<PB> PB::build(std::vector<boost::shared_ptr<Wavefunction> > m)
{
    PB* pb = new PB();

    Options& options = Process::environment.options; 

    pb->print_ = options.get_int("PRINT");
    pb->debug_ = options.get_int("DEBUG");
    pb->bench_ = options.get_int("BENCH");

    pb->lambda_ = options.get_double("PB_LAMBDA");

    pb->S_cutoff_ = options.get_double("S_TOLERANCE");

    pb->maxiter_ = options.get_int("MAXITER");
    pb->E_convergence_ = options.get_double("E_CONVERGENCE");
    pb->D_convergence_ = options.get_double("D_CONVERGENCE");

    pb->diis_ = options.get_bool("DIIS");
    pb->diis_min_ = options.get_int("DIIS_MIN_VECS");
    pb->diis_max_ = options.get_int("DIIS_MAX_VECS");
    pb->diis_start_ = options.get_int("DIIS_START");

    pb->memory_ = (unsigned long int)(Process::environment.get_memory() * options.get_double("SAPT_MEM_FACTOR") * 0.125);

    pb->monomers_ = m;

    return boost::shared_ptr<PB>(pb);
}
void PB::initialize()
{
    if (monomers_.size() == 0) {
        throw PSIEXCEPTION("What is going on here, no monomers?");
    }

    // Sizing
    int N = monomers_.size();
    int nso = monomers_[0]->basisset()->nbf();
    occ_starts_.clear();
    occ_starts_.push_back(0);
    for (int A = 0; A < N; A++) {
        occ_starts_.push_back(monomers_[A]->doccpi().sum() + occ_starts_[A]);
    }

    // S/T matrix
    S_ = boost::shared_ptr<Matrix>(new Matrix("S",nso,nso));
    T_ = boost::shared_ptr<Matrix>(new Matrix("T",nso,nso));
    boost::shared_ptr<IntegralFactory> Tfact(new IntegralFactory(monomers_[0]->basisset()));   
    boost::shared_ptr<OneBodyAOInt> Tint(Tfact->ao_kinetic());
    boost::shared_ptr<OneBodyAOInt> Sint(Tfact->ao_overlap());
    Tint->compute(T_);
    Sint->compute(S_);
    Tint.reset();
    Sint.reset();
    Tfact.reset();

    // X matrix
    X_ = S_->canonical_orthogonalization(S_cutoff_);
    if (X_->colspi()[0] != monomers_[0]->nmo()) {
        throw PSIEXCEPTION("Orthogonalization techniques are different, call Rob");
    }
    int nmo = X_->colspi()[0];
    int nvir = nmo - occ_starts_[N];

    // V matrices
    V_.clear();
    for (int A = 0; A < monomers_.size(); A++) {
        V_.push_back(boost::shared_ptr<Matrix>(new Matrix("V",nso,nso)));
        boost::shared_ptr<IntegralFactory> Vfact(new IntegralFactory(monomers_[A]->basisset()));   
        boost::shared_ptr<OneBodyAOInt> Vint(Vfact->ao_potential());
        Vint->compute(V_[A]);
    } 

    // JK object (hopefully DFJK already on disk)
    jk_ = JK::build_JK();
    long int jk_memory = (long int)memory_; 
    jk_memory -= 2L * N * nso * nso;
    jk_memory -= 1L * N * nmo * nmo;
    jk_memory -= 7L * nso * nso;
    if (jk_memory < 0L) {
        throw PSIEXCEPTION("Too little static memory for INFSAPT::scf_terms");
    }
    jk_->set_memory((unsigned long int )jk_memory);

    jk_->set_do_J(true);
    jk_->set_do_K(true);
    jk_->initialize();

    // Initial Cocc matrices and other allocation
    Cocc_.clear();
    Cvir_.clear();
    eps_occ_.clear();
    eps_vir_.clear();
    for (int A = 0; A < monomers_.size(); A++) {
        Cocc_.push_back(monomers_[A]->Ca_subset("AO", "OCC"));
        Cvir_.push_back(boost::shared_ptr<Matrix>(new Matrix("Cvir",nso,nvir)));
        eps_occ_.push_back(boost::shared_ptr<Vector>(new Vector("eps_occ",Cocc_[A]->colspi()[0])));
        eps_vir_.push_back(boost::shared_ptr<Vector>(new Vector("eps_vir",nvir)));
    }

    // => Initial J/K matrices [for E_HF^(0)] <= //

    // Cl matrices
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    // Cr matrices
    std::vector<SharedMatrix>& Cr = jk_->C_right();

    Cl.clear();
    Cr.clear();

    for (int A = 0; A < monomers_.size(); A++) {
        Cl.push_back(Cocc_[A]);
        Cr.push_back(Cocc_[A]);
    }

    jk_->compute();
    
    E_HF_0_ = HF_energy();
    E_HF_0_->set_name("E_HF^(0)");
    E_HF_A_ = E_HF_0_;
    E_0_ = E_;

    // F/Q matrices
    F_.clear();
    Q_.clear();
    for (int A = 0; A < monomers_.size(); A++) {
        F_.push_back(boost::shared_ptr<Matrix>(new Matrix("F",nso,nso)));
        Q_.push_back(boost::shared_ptr<Matrix>(new Matrix("Q",nmo,nmo)));

    }

    S_0_ = compute_Soo();
}
void PB::print_header() const 
{
    outfile->Printf( "  ==> Pauli Blockade <==\n\n");

    outfile->Printf( "    Lambda:               %11.3E\n", lambda_);
    outfile->Printf( "    S cutoff:             %11.3E\n", S_cutoff_);
    outfile->Printf( "    Maximum iterations:   %11d\n", maxiter_);
    outfile->Printf( "    Energy threshold:     %11.3E\n", E_convergence_);
    outfile->Printf( "    Commutator threshold: %11.3E\n", D_convergence_);
    outfile->Printf( "    DIIS:                 %11s\n", diis_ ? "Yes" : "No");
    outfile->Printf( "\n");

    jk_->print_header();

    
}
void PB::compute()
{
    outfile->Printf("    Unconstrained Total Energy: %24.14f\n\n", E_);
    
    outfile->Printf("    Occupied guess via symmetric orthogonalization.\n\n");

    // Do a Pauli explosion to start from a good spot
    symmetric_orthogonalize();

    // ==> Master PBD Loop <== //    

    bool converged = false; 
    bool extrapolated = false;
    E_HF_L_ = boost::shared_ptr<Vector>(E_HF_A_->clone());
    E_HF_L_->zero();

    outfile->Printf("   => Pauli Blockade Iterations <=\n\n");

    outfile->Printf( "    %-14s %24s %11s %11s %s\n",
        "", "Total Energy", "DeltaE", "Q RMS", "");
    outfile->Printf( "\n");
    

    for (int iter = 0; iter < maxiter_ || converged; iter++) {

        // Build the new J/K matrices
        jk_->compute();

        // Compute the new HF energy
        E_HF_A_ = HF_energy();

        // Build the full Fock matrices with penalties
        build_fock();

        // Build the commutators
        build_commutators();

        // Check that the energy and commutators are converged
        converged = check_convergence();
        
        // Message to the user
        outfile->Printf( "    @PB Iter %4d: %24.14f %11.3E %11.3E %s\n",
            iter + 1, E_, dE_, dQ_, (extrapolated ? "DIIS" : ""));
        

        // New Fock matrix components are built, we are self-consistent
        if (converged) break; 

        // DIIS (modifies the Fock matrices)
        extrapolated = false;
        if (iter + 1 >= diis_start_) {
            extrapolated = diis();
        }

        // Diagonalize the Fock matrices
        diagonalize(); 

    }

    if (!converged) {
        throw PSIEXCEPTION("PB did not converge");
    }

    outfile->Printf("\n    Pauli Blockade Converged.\n\n");

    outfile->Printf("    Constrained Total Energy: %24.14f\n", E_);
    outfile->Printf("    Pauli Total Energy:       %24.14f\n", E_ - E_0_);
    outfile->Printf("\n");

    // Remove the penalty function from the orbitals
    purify_eigenvalues();

    // Give the user an indication of what happened
    print_summary();
    print_errors();
}
void PB::symmetric_orthogonalize()
{
    // Metrix matrix
    boost::shared_ptr<Matrix> Soo = compute_Soo();

    if (debug_) {
        outfile->Printf( "    > Symmetric Orthogonalization <\n\n");
        Soo->print();
    }

    // Symmetric orthogonalization transformation
    boost::shared_ptr<Matrix> X(Soo->clone());
    X->copy(Soo);
    X->set_name("X"); 
    X->power(-1.0/2.0); 

    // Orbital transformation
    boost::shared_ptr<Matrix> Cocc1 = Matrix::horzcat(Cocc_);
    boost::shared_ptr<Matrix> Cocc2(Cocc1->clone());

    int nso =  Cocc1->rowspi()[0];
    int nocc = Cocc1->colspi()[0];

    double** Xp =  X->pointer();
    double** C1p = Cocc1->pointer();
    double** C2p = Cocc2->pointer();

    C_DGEMM('N','N',nso,nocc,nocc,1.0,C1p[0],nocc,Xp[0],nocc,0.0,C2p[0],nocc); 

    // Orbital selection (Based on <p|a><a|p>)

    // The selections should remain unique (pronounced "assignable")
    std::set<int> taken;

    for (int A = 0; A < Cocc_.size(); A++) {
        int na = Cocc_[A]->colspi()[0];
        boost::shared_ptr<Matrix> T(new Matrix("T",nso,na));
        boost::shared_ptr<Matrix> R(new Matrix("R",nocc,na));
        boost::shared_ptr<Vector> P(new Vector("P",nocc));
        double** Cp = Cocc_[A]->pointer();
        double** Sp = S_->pointer();
        double** Tp = T->pointer();
        double** Rp = R->pointer();    
        double* Pp = P->pointer();

        // Compute <p|a><a|p> for this a space
        C_DGEMM('N','N',nso,na,nso,1.0,Sp[0],nso,Cp[0],na,0.0,Tp[0],na);
        C_DGEMM('T','N',nocc,na,nso,1.0,C2p[0],nocc,Tp[0],na,0.0,Rp[0],na);
        std::vector<std::pair<double,int> > order;
        for (int a = 0; a < nocc; a++) {
            Pp[a] = C_DDOT(na,Rp[a],1,Rp[a],1);
            order.push_back(std::pair<double,int>(Pp[a],a));
        } 
        
        if (debug_) {
            P->print();
        }

        // Sort the projections by p descending
        std::sort(order.begin(),order.end(),std::greater<std::pair<double, int> >());

        // Assign the best na p to be the chosen a
        for (int a = 0; a < na; a++) {
            int ind = order[a].second;
            if (taken.count(ind)) {
                throw PSIEXCEPTION("Symmetric orthogonalization of monomers is not assignable.");
            } else {
                taken.insert(ind);
                C_DCOPY(nso,&C2p[0][ind],nocc,&Cp[0][a],na);
            }
        }
    }
}
void PB::build_fock()
{
    // D matrices
    const std::vector<SharedMatrix>& D = jk_->D();
    // J matrices
    const std::vector<SharedMatrix>& J = jk_->J();
    // K matrices
    const std::vector<SharedMatrix>& K = jk_->K();

    int nso = Cocc_[0]->rowspi()[0];

    // Standard Fock Matrix
    for (int A = 0; A < Cocc_.size(); A++) {
        double** Fp =  F_[A]->pointer();
        C_DAXPY(nso * nso, 1.0, T_->pointer()[0], 1, Fp[0], 1);
        C_DAXPY(nso * nso, 1.0, V_[A]->pointer()[0], 1, Fp[0], 1);
        C_DAXPY(nso * nso, 2.0, J[A]->pointer()[0], 1, Fp[0], 1);
        C_DAXPY(nso * nso,-1.0, K[A]->pointer()[0], 1, Fp[0], 1);
    }

    // Penalty operator
    boost::shared_ptr<Matrix> R(new Matrix("R", nso, nso));
    double** Rp = R->pointer();
    for (int B = 0; B < Cocc_.size(); B++) {
        int na = Cocc_[B]->colspi()[0];
        boost::shared_ptr<Matrix> T(new Matrix("T", na, nso));
        double** Tp = T->pointer();
        double** Sp = S_->pointer();
        double** Cp = Cocc_[B]->pointer(); 
        C_DGEMM('T','N',na,nso,nso,1.0,Cp[0],na,Sp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('T','N',nso,nso,na,1.0,Tp[0],nso,Tp[0],nso,0.0,Rp[0],nso);
        for (int A = 0; A < Cocc_.size(); A++) {
            if (A == B) continue;
            double** Fp =  F_[A]->pointer();
            C_DAXPY(nso * nso, lambda_, Rp[0], 1, Fp[0], 1);
        }
    } 
}
void PB::build_commutators()
{
    // D matrices
    const std::vector<SharedMatrix>& D = jk_->D();

    int nso = X_->rowspi()[0];
    int nmo = X_->colspi()[0];

    double** Xp = X_->pointer();
    double** Sp = S_->pointer();

    // X[FDS-SDF]X
    boost::shared_ptr<Matrix> FD(new Matrix("FD",nso,nso));
    boost::shared_ptr<Matrix> FDS(new Matrix("FDS",nso,nso));
    boost::shared_ptr<Matrix> XFDS(new Matrix("XFDS",nmo,nso));

    double** FDp = FD->pointer();
    double** FDSp = FDS->pointer();
    double** XFDSp = XFDS->pointer();

    for (int A = 0; A < Cocc_.size(); A++) {
         
        double** Fp = F_[A]->pointer();
        double** Dp = D[A]->pointer();
        double** Qp = Q_[A]->pointer();

        C_DGEMM('N','N',nso,nso,nso,1.0,Fp[0],nso,Dp[0],nso,0.0,FDp[0],nso);
        C_DGEMM('N','N',nso,nso,nso,1.0,FDp[0],nso,Sp[0],nso,0.0,FDSp[0],nso);
        C_DGEMM('T','N',nmo,nso,nso,1.0,Xp[0],nmo,FDSp[0],nso,0.0,XFDSp[0],nso);
        C_DGEMM('N','N',nmo,nmo,nso,1.0,XFDSp[0],nso,Xp[0],nmo,0.0,Qp[0],nmo);

        for (int p = 0; p < nmo; p++) {
            for (int q = p; q < nmo; q++) {
                double val = Qp[p][q] - Qp[q][p];
                Qp[p][q] = val;
                Qp[q][p] = -val;    
            }
        }
    } 
}
bool PB::check_convergence()
{
    dE_ = 0.0;
    dQ_ = 0.0;
    for (int A = 0; A < Q_.size(); A++) {
        dE_ += (E_HF_A_->get(A) - E_HF_L_->get(A));
        double val = Q_[A]->rms();
        dQ_ = (val > dQ_ ? val : dQ_);
    }

    bool converged = false;
    //if (abs(dE_) < E_convergence_ && abs(dQ_) < D_convergence_) {
    if (abs(dE_) < E_convergence_) {
        converged = true;
    }

    E_HF_L_ = E_HF_A_;
    
    return converged;
}
bool PB::diis()
{
    if (!diis_manager_) {   
        diis_manager_ = boost::shared_ptr<DIISN>(new DIISN("PB DIIS",diis_min_,diis_max_,F_,Q_));
    }

    return diis_manager_->diis(F_,Q_);
}
void PB::diagonalize()
{
    int nso = X_->rowspi()[0];
    int nmo = X_->colspi()[0];

    boost::shared_ptr<Matrix> T(new Matrix("T",nso,nmo));
    boost::shared_ptr<Matrix> F2(new Matrix("F2",nmo,nmo));
    boost::shared_ptr<Matrix> C(new Matrix("C",nmo,nmo));
    boost::shared_ptr<Vector> e(new Vector("e",nmo));

    double** Tp = T->pointer(); 
    double** F2p = F2->pointer(); 
    double** Xp = X_->pointer();
    double** Cp = C->pointer();
    double*  ep = e->pointer();

    for (int A = 0; A < Cocc_.size(); A++) {
        double** Fp = F_[A]->pointer();
        
        // Transform to orthonormal basis
        C_DGEMM('N','N',nso,nmo,nso,1.0,Fp[0],nso,Xp[0],nmo,0.0,Tp[0],nmo);
        C_DGEMM('T','N',nmo,nmo,nso,1.0,Xp[0],nmo,Tp[0],nmo,0.0,F2p[0],nmo);
    
        // Diagonalize
        F2->diagonalize(C,e);

        // Select low-energy occ/vir orbitals
        int na = Cocc_[A]->colspi()[0];
        int nv = Cvir_[A]->colspi()[0];

        double** Cap = Cocc_[A]->pointer();
        double** Cvp = Cvir_[A]->pointer();
        double*  eap = eps_occ_[A]->pointer();
        double*  evp = eps_vir_[A]->pointer();

        C_DGEMM('N','N',nso,na,nmo,1.0,Xp[0],nmo,&Cp[0][0],nmo,0.0,Cap[0],na);
        C_DGEMM('N','N',nso,nv,nmo,1.0,Xp[0],nmo,&Cp[0][na],nmo,0.0,Cvp[0],nv);

        C_DCOPY(na,&ep[0],1,eap,1);
        C_DCOPY(nv,&ep[na],1,evp,1);
    }
}
void PB::purify_eigenvalues()
{
    // J matrices
    const std::vector<SharedMatrix>& J = jk_->J();
    // K matrices
    const std::vector<SharedMatrix>& K = jk_->K();

    int nso = Cocc_[0]->rowspi()[0];
    int nvir = Cvir_[0]->colspi()[0];

    boost::shared_ptr<Matrix> F(new Matrix("F",nso,nso));
    boost::shared_ptr<Matrix> T2(new Matrix("T",nvir,nso)); 

    for (int A = 0; A < Cocc_.size(); A++) {
        int na = Cocc_[0]->colspi()[0];
        boost::shared_ptr<Matrix> T1(new Matrix("T",na,nso)); 

        double** Fp =  F->pointer();
        double** T1p = T1->pointer();
        double** T2p = T2->pointer();
        double** C1p = Cocc_[A]->pointer();
        double** C2p = Cvir_[A]->pointer();
        double*  e1p = eps_occ_[A]->pointer();
        double*  e2p = eps_vir_[A]->pointer();
 
        // Build the true Fock matrix
        F->zero();
        C_DAXPY(nso * nso, 1.0, T_->pointer()[0], 1, Fp[0], 1);
        C_DAXPY(nso * nso, 1.0, V_[A]->pointer()[0], 1, Fp[0], 1);
        C_DAXPY(nso * nso, 2.0, J[A]->pointer()[0], 1, Fp[0], 1);
        C_DAXPY(nso * nso,-1.0, K[A]->pointer()[0], 1, Fp[0], 1);

        // <a|F_A|a>
        C_DGEMM('T','N',na,nso,nso,1.0,C1p[0],na,Fp[0],nso,0.0,T1p[0],nso);
        for (int a = 0; a < na; a++) {
            e1p[a] = C_DDOT(nso,T1p[a],1,&C1p[0][a],na);
        }

        // <v|F_A|v>
        C_DGEMM('T','N',nvir,nso,nso,1.0,C2p[0],nvir,Fp[0],nso,0.0,T2p[0],nso);
        for (int a = 0; a < nvir; a++) {
            e2p[a] = C_DDOT(nso,T2p[a],1,&C2p[0][a],nvir);
        }
    }
}
boost::shared_ptr<Matrix> PB::compute_Soo()
{
    boost::shared_ptr<Matrix> Cocc = Matrix::horzcat(Cocc_);

    int nso = Cocc->rowspi()[0];
    int nocc = Cocc->colspi()[0];

    boost::shared_ptr<Matrix> T(new Matrix("T",nso,nocc));
    boost::shared_ptr<Matrix> Soo(new Matrix("Soo",nocc,nocc));

    double** Cp = Cocc->pointer();
    double** Sp = S_->pointer();
    double** Tp = T->pointer();
    double** Soop = Soo->pointer();

    C_DGEMM('N','N',nso,nocc,nso,1.0,Sp[0],nso,Cp[0],nocc,0.0,Tp[0],nocc);
    C_DGEMM('T','N',nocc,nocc,nso,1.0,Cp[0],nocc,Tp[0],nocc,0.0,Soop[0],nocc);
    
    return Soo;
}
boost::shared_ptr<Vector> PB::HF_energy()
{
    // D matrices
    const std::vector<SharedMatrix>& D = jk_->D();
    // J matrices
    const std::vector<SharedMatrix>& J = jk_->J();
    // K matrices
    const std::vector<SharedMatrix>& K = jk_->K();

    boost::shared_ptr<Vector> E(new Vector("E_HF", monomers_.size())); 
    
    E_ = 0.0;
    for (int A = 0; A < monomers_.size(); A++) {
        double E_HF = 0.0;
        E_HF += 2.0 * D[A]->vector_dot(T_);
        E_HF += 2.0 * D[A]->vector_dot(V_[A]);
        E_HF += 2.0 * D[A]->vector_dot(J[A]);
        E_HF -= 1.0 * D[A]->vector_dot(K[A]);
        E_HF += monomers_[A]->molecule()->nuclear_repulsion_energy(); 
        E->set(A,E_HF);
        E_ += E_HF;
    } 

    return E;
}
void PB::print_summary()
{
    outfile->Printf("   => Pauli Blockade Results <=\n\n");
    
    for (int A = 0; A < Cocc_.size(); A++) {
        outfile->Printf("    Monomer %d:\n\n", A+1);

        outfile->Printf("    Unconstrained Energy: %24.14f\n", E_HF_0_->get(A));
        outfile->Printf("    Constrained Energy:   %24.14f\n", E_HF_A_->get(A));
        outfile->Printf("    Pauli Energy:         %24.14f\n", E_HF_A_->get(A) - E_HF_0_->get(A));
        outfile->Printf("\n");

        int count;
        int n;  
        int offset;
        double* ep;

        outfile->Printf( "\t%-70s\n\n\t", "Purified Occupied Orbital Eigenvalues:");
        count = 0;
        n = eps_occ_[A]->dimpi()[0];
        offset = 0;
        ep = eps_occ_[A]->pointer();
        for (int i = 0; i < n; i++) {
            outfile->Printf( "%4d%-4s%11.6f  ", i + offset + 1, "A", ep[i]);
            if (count++ % 3 == 2 && count != n)
                outfile->Printf( "\n\t");
        }
        outfile->Printf( "\n\n");

        outfile->Printf( "\t%-70s\n\n\t", "Purified Virtual Orbital Eigenvalues:");
        count = 0;
        n = eps_vir_[A]->dimpi()[0];
        offset = eps_occ_[A]->dimpi()[0];
        ep = eps_vir_[A]->pointer();
        for (int i = 0; i < n; i++) {
            outfile->Printf( "%4d%-4s%11.6f  ", i + offset + 1, "A", ep[i]);
            if (count++ % 3 == 2 && count != n)
                outfile->Printf( "\n\t");
        }
        outfile->Printf( "\n\n");
        
    }    

    
}
void PB::print_errors()
{
    outfile->Printf("   => Pauli Blockade Overlap Errors <=\n\n");
    
    boost::shared_ptr<Matrix> Soo = compute_Soo();

    double** SAp = Soo->pointer();
    double** S0p = S_0_->pointer();

    double SAG = 0.0;
    double S0G = 0.0;
    
    outfile->Printf( "    -----------------------------------------------\n");
    outfile->Printf( "    %-11s %11s %11s %11s\n",
        "Interaction", "HF^(A)", "HF^(0)", "Ratio");
    outfile->Printf( "    -----------------------------------------------\n");

    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            
            double SA = 0.0;
            double S0 = 0.0;
            for (int a = occ_starts_[A]; a < occ_starts_[A+1]; a++) {
                for (int b = occ_starts_[B]; b < occ_starts_[B+1]; b++) {
                    SA = (SA > SAp[a][b] ? SA : SAp[a][b]);
                    S0 = (S0 > S0p[a][b] ? S0 : S0p[a][b]);
                }
            }
            
            SAG = (SAG > SA ? SAG : SA);
            S0G = (S0G > S0 ? S0G : S0);

            outfile->Printf( "    %-3d <-> %3d %11.3E %11.3E %11.3E\n",
                A+1, B+1, SA, S0, SA / S0);
        }
    }

    outfile->Printf( "    -----------------------------------------------\n");
    outfile->Printf( "    %-11s %11.3E %11.3E %11.3E\n",
        "Overall", SAG, S0G, SAG / S0G);
    outfile->Printf( "    -----------------------------------------------\n");
    outfile->Printf( "\n");

    
}

DIISN::DIISN(const std::string& name,
    int min_vecs,
    int max_vecs,
    const std::vector<boost::shared_ptr<Matrix> >& x,
    const std::vector<boost::shared_ptr<Matrix> >& f) :
    name_(name), min_vecs_(min_vecs), max_vecs_(max_vecs)
{
    N_ = x.size();

    if (N_ != f.size()) {
        throw PSIEXCEPTION("DIISN: x.size() != f.size()");
    }
    if (N_ <= 0) {
        throw PSIEXCEPTION("DIISN: x.size() <= 0?");
    }
    if (x[0]->nirrep() != 1) {
        throw PSIEXCEPTION("DIISN: x is not C1");
    }
    if (f[0]->nirrep() != 1) {
        throw PSIEXCEPTION("DIISN: f is not C1");
    }

    x_numel_ = x[0]->rowspi()[0] * (size_t) x[0]->colspi()[0];
    f_numel_ = f[0]->rowspi()[0] * (size_t) f[0]->colspi()[0];

    for (int n = 0; n < N_; n++) {
        if (x[n]->rowspi()[0] * (size_t) x[n]->colspi()[0] != x_numel_) {
            throw PSIEXCEPTION("DIISN: x must all be the same size");
        }
        if (f[n]->rowspi()[0] * (size_t) f[n]->colspi()[0] != f_numel_) {
            throw PSIEXCEPTION("DIISN: f must all be the same size");
        }
    }

    metric_ = boost::shared_ptr<Matrix>(new Matrix("B", max_vecs, max_vecs)); 

    fh_ = fopen(filename().c_str(), "wb+");

    iter_ = 0;
} 
DIISN::~DIISN()
{
    fclose(fh_);
    remove(filename().c_str());
}
std::string DIISN::filename() const
{
    std::stringstream ss;
    ss <<  PSIOManager::shared_object()->get_default_path();
    ss <<  "/";
    ss << psi_file_prefix;
    ss << ".";
    ss << getpid(); 
    ss << ".";
    ss << PSIO::get_default_namespace();
    ss << ".";
    ss << name_;
    ss << ".dat";
    return ss.str();
}
bool DIISN::diis(std::vector<boost::shared_ptr<Matrix> >& x,
                std::vector<boost::shared_ptr<Matrix> >& f)
{
    // ==> Add Phase <== //

    // => Index determination <= //
    int nvec;
    int ind;

    if (iter_ < max_vecs_) {
        // Pattern is not full, add to the next open space
        nvec = iter_ + 1;
        ind = iter_;
    } else {
        // Pattern is full, overwrite the worst residual
        double** Bp = metric_->pointer(); 
        double maxB = 0.0;
        ind = 0;
        for (int v = 0; v < max_vecs_; v++) {
            if (Bp[v][v] >= maxB) {
                maxB = Bp[v][v];
                ind = v;    
            }
        }
        nvec = max_vecs_;
    }

    // => Entry export <= //

    fseek(fh_, ind * (x_numel_ + f_numel_) * N_ * sizeof(double), SEEK_SET);
    for (int n = 0; n < N_; n++) {
        fwrite((void*) x[n]->pointer()[0],sizeof(double),x_numel_,fh_);
    } 
    for (int n = 0; n < N_; n++) {
        fwrite((void*) f[n]->pointer()[0],sizeof(double),f_numel_,fh_);
    } 

    // => B matrix update <= //

    double** Bp = metric_->pointer();
    boost::shared_ptr<Matrix> f2(f[0]->clone()); 

    for (int v = 0; v < nvec; v++) {
        fseek(fh_, (v * (x_numel_ + f_numel_) * N_ + x_numel_ * N_) * sizeof(double), SEEK_SET);
        double Bval = 0.0;
        for (int n = 0; n < N_; n++) {
            fread((void*) f2->pointer()[0],sizeof(double),f_numel_,fh_);
            Bval += f2->vector_dot(f[n]); 
        } 
        Bp[v][ind] = Bp[ind][v] = Bval;
    }

    iter_++;
    
    // ==> Extrapolate Phase <== //

    if (nvec < min_vecs_) {
        return false;
    }

    // => Extrapolation Coefficients <= //

    std::vector<std::pair<double,int> > order;
    for (int v = 0; v < nvec; v++) {
        order.push_back(std::pair<double, int>(Bp[v][v],v));
    }
    std::sort(order.begin(), order.end(), std::greater<std::pair<double, int> >());
    std::vector<int> sorted;
    for (int v = 0; v < nvec; v++) {
        sorted.push_back(order[v].second);
    }

    // > Sorted Basis C matrix < //

    boost::shared_ptr<Matrix> C(new Matrix("C", nvec + 1, nvec + 1));
    boost::shared_ptr<Matrix> D(new Matrix("D", nvec + 1, nvec + 1));
    boost::shared_ptr<Vector> r(new Vector("r", nvec + 1));
    boost::shared_ptr<Vector> l(new Vector("l", nvec + 1));

    double** Cp = C->pointer();
    double** Dp = D->pointer();
    double*  lp = l->pointer();
    double*  rp = r->pointer();

    rp[nvec] = -1.0; 

    for (int v = 0; v < nvec; v++) {
        for (int w = 0; w < nvec; w++) {
            Cp[v][w] = Bp[sorted[v]][sorted[w]];
        }
        Cp[v][nvec] = Cp[nvec][v] = -1.0;
    }

    // > Conditioned Inversion < //

    int* pivots = new int[nvec+1];
    int offset = 0;
    while (true) {
        // Destructive LU 
        D->copy(C);
        // Do the LU
        int info = C_DGETRF(nvec+1 - offset, nvec+1 - offset, &Dp[offset][offset], nvec+1, pivots);
        // LU failing is usually what you call a clue
        if (info) {
            offset++;
            continue;
        }
        // You are guaranteed to get here at some point (WCS is [X 1; 1 0], with the smallest X in the set)
        C_DGETRS('N',nvec+1 - offset, 1, &Dp[offset][offset], nvec+1, pivots, &rp[offset], nvec+1);
        for (int v = offset; v < nvec; v++) {
            lp[sorted[v]] = rp[v];
        }
        break;
    }
    delete[] pivots;

    // => Extrapolation <= //

    boost::shared_ptr<Matrix> x2(x[0]->clone()); 

    for (int n = 0; n < N_; n++) {
        x[n]->zero();
    }

    for (int v = 0; v < nvec; v++) {
        fseek(fh_, (v * (x_numel_ + f_numel_) * N_) * sizeof(double), SEEK_SET);
        double lval = lp[v];
        for (int n = 0; n < N_; n++) {
            fread((void*) x2->pointer()[0],sizeof(double),x_numel_,fh_);
            C_DAXPY(x_numel_,lval,x2->pointer()[0],1,x[n]->pointer()[0],1);
        } 
    }    
    
    return true;
}

}}
