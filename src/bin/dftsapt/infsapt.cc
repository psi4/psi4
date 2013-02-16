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
    print_header();

    scf_terms();
    
    pt2_terms();

    print_trailer();
    
    return 0.0;
}
void INFSAPT::print_header() const 
{
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t                         INF-SAPT0                       \n");
    fprintf(outfile, "\t               Rob Parrish and Ed Hohenstein             \n");
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  ==> Sizes <==\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "   => Parameters <=\n\n");

    int threads = 1;
    #ifdef _OPENMP
        threads = omp_get_max_threads();
    #endif
    
    fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
    fprintf(outfile, "    Threads:           %11d\n", threads);
    fprintf(outfile, "    Schwarz Cutoff:    %11.3E\n", schwarz_);
    fprintf(outfile, "\n");

    fprintf(outfile,"   => Molecular Cluster <=\n\n");
    cluster_->molecule()->print_cluster();

    fprintf(outfile, "   => Orbital Ranges (Pauli Blockade Space) <=\n\n");

    int nso = cluster_->nso();
    int nmo = cluster_->nmo();
    int nocc = cluster_->doccpi().sum();
    int nvir = nmo - nocc;
    int nfocc = cluster_->frzcpi().sum();
    int nfvir = cluster_->frzvpi().sum();
    int naocc = nocc - nfocc;
    int navir = nvir - nfvir;

    fprintf(outfile,"    ---------------------------------------------------------------\n");
    fprintf(outfile,"    %-7s %6s %6s %6s %6s %6s %6s %6s %6s\n", 
        "Unit", "Nso", "Nmo", "Nocc", "Nvir", "Nfocc", "Naocc", "Navir", "Nfvir");

    fprintf(outfile,"    ---------------------------------------------------------------\n");
    fprintf(outfile,"    %-7s %6d %6d %6d %6d %6d %6d %6d %6d\n", 
        "Cluster", nso, nmo, nocc, nvir, nfocc, naocc, navir, nfvir);
    
    for (int A = 0; A < monomers_.size(); A++) {
        boost::shared_ptr<Wavefunction> m = monomers_[A];
        nocc = m->doccpi().sum();
        nfocc = m->frzcpi().sum();
        naocc = nocc - nfocc;
        
        fprintf(outfile,"    %-7d %6d %6d %6d %6d %6d %6d %6d %6d\n", 
            A+1, nso, nmo, nocc, nvir, nfocc, naocc, navir, nfvir);
    }
    fprintf(outfile,"    ---------------------------------------------------------------\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "   => Primary Cluster Basis Set <=\n\n");
    primary_->print_by_level(outfile, print_);

    fflush(outfile);
}
void INFSAPT::scf_terms()
{
    fprintf(outfile, "  SCF TERMS:\n\n");

    // ==> Setup <== //

    // Number of bodies in this job
    int N = monomers_.size();
    // Number of AO basis functions
    int nso = primary_->nbf();
    
    // JK object (hopefully DFJK already on disk)
    boost::shared_ptr<JK> jk = JK::build_JK();
    long int jk_memory = (long int)memory_; 
    if (jk_memory < 0L) {
        throw PSIEXCEPTION("Too little static memory for INFSAPT::scf_terms");
    }
    jk->set_memory((unsigned long int )jk_memory);

    jk->set_do_J(true);
    jk->set_do_K(true);
    jk->initialize();
    jk->print_header();

    fflush(outfile);

    // T matrix
    boost::shared_ptr<Matrix> T(new Matrix("T",nso,nso));
    boost::shared_ptr<IntegralFactory> Tfact(new IntegralFactory(primary_));   
    boost::shared_ptr<OneBodyAOInt> Tint(Tfact->ao_kinetic());
    Tint->compute(T);
    Tint.reset();
    Tfact.reset();

    // V matrices
    std::vector<boost::shared_ptr<Matrix> > V;
    for (int A = 0; A < monomers_.size(); A++) {
        V.push_back(boost::shared_ptr<Matrix>(new Matrix("V",nso,nso)));
        boost::shared_ptr<IntegralFactory> Vfact(new IntegralFactory(monomers_[A]->basisset()));   
        boost::shared_ptr<OneBodyAOInt> Vint(Vfact->ao_potential());
        Vint->compute(V[A]);
    } 

    // C matrices
    std::vector<boost::shared_ptr<Matrix> > C;
    for (int A = 0; A < monomers_.size(); A++) {
        C.push_back(monomers_[A]->Ca_subset("AO", "OCC"));
    }

    // Cl matrices
    std::vector<SharedMatrix>& Cl = jk->C_left();
    // Cr matrices
    std::vector<SharedMatrix>& Cr = jk->C_right();

    // D matrices
    const std::vector<SharedMatrix>& D = jk->D();
    // J matrices
    const std::vector<SharedMatrix>& J = jk->J();
    // K matrices
    const std::vector<SharedMatrix>& K = jk->K();

    // F matrices
    std::vector<boost::shared_ptr<Matrix> > F;
    for (int A = 0; A < monomers_.size(); A++) {
        F.push_back(boost::shared_ptr<Matrix>(new Matrix("F",nso,nso)));
    }

    // ==> (HF) Terms <== //

    // E_HF^(0)
    double E_HF = 0.0;
    for (int A = 0; A < monomers_.size(); A++) {
        E_HF += monomers_[A]->reference_energy();
    }
    
    // ==> (0) Terms <== //

    Cl.clear();
    Cr.clear();

    for (int A = 0; A < monomers_.size(); A++) {
        Cl.push_back(C[A]);
        Cr.push_back(C[A]);
    }

    jk->compute();

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

    // TODO: PBD => (F, C, and \epsilon)

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
        double EK = 0.0;
        EK += 2.0 * D[A]->vector_dot(T);
        EK += 2.0 * D[A]->vector_dot(V[A]);
        EK += 2.0 * D[A]->vector_dot(F[A]);
        EK -= monomers_[A]->reference_energy(); 
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

    fprintf(outfile, "  ==> SCF Analysis <==\n\n");
    
    fprintf(outfile, "   => Electrostatics <=\n\n");

    fprintf(outfile,"    Total:\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EJ_0->get(A,B);
            fprintf(outfile, "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_elst, pc_hartree2kcalmol*E_elst);
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");

    fprintf(outfile,"   => Exchange <=\n\n");

    fprintf(outfile,"    Pauli Term:\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
            double Eval = EK_A->get(A,A);
            fprintf(outfile, "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, A+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
    }    
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile, "    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_Pauli, pc_hartree2kcalmol*E_Pauli);
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");

    fprintf(outfile,"    Tunneling Term:\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EK_0->get(A,B);
            fprintf(outfile, "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile, "    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_tunnel, pc_hartree2kcalmol*E_tunnel);
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");

    fprintf(outfile,"    Rearrangement Term:\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EK_A->get(A,B);
            fprintf(outfile, "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_cross, pc_hartree2kcalmol*E_cross);
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");

    fprintf(outfile,"    Total:\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A; B < monomers_.size(); B++) {
            double Eval = EK_0->get(A,B) + EK_A->get(A,B);
            fprintf(outfile, "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_exch, pc_hartree2kcalmol*E_exch);
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");

    fprintf(outfile, "   => Induction <=\n\n");

    fprintf(outfile,"    Total:\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*E_ind, pc_hartree2kcalmol*E_ind);
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");

    energies_["Elst"] = E_elst;
    energies_["Exch"] = E_exch;
    energies_["Ind"] = E_ind;

    fflush(outfile);
}
void INFSAPT::pt2_terms()
{
    fprintf(outfile, "  PT2 TERMS:\n\n");

    fprintf(outfile, "   => MP2FIT Auxiliary Basis Set <=\n\n");
    mp2fit_->print_by_level(outfile, print_);

    fflush(outfile);

    // ==> Sizing <== //

    int N = monomers_.size();
    int nso   = cluster_->nso();
    int naux  = mp2fit_->nbf();
    int naocc = cluster_->doccpi().sum() - cluster_->frzcpi().sum();
    int navir = Cavir_->colspi()[0];
    std::vector<int> aocc_starts;
    aocc_starts.push_back(0);
    for (int A = 0; A < N; A++) {
        aocc_starts.push_back(Caocc_[A]->colspi()[0] + aocc_starts[A]);
    }

    // ==> Integral Generation <== //
    
    boost::shared_ptr<Matrix> Caocc = Matrix::horzcat(Caocc_);
    build_iaQ(primary_,mp2fit_,Caocc,Cavir_);
    Caocc.reset();

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
    std::vector<boost::shared_ptr<Matrix> > I(threads); 
    std::vector<double**> Ipointers;
    for (int thread = 0; thread <= threads; thread++) {
        I.push_back(boost::shared_ptr<Matrix>(new Matrix("I", navir, navir)));
        Ipointers[threads] = I[thread]->pointer();
    } 

    // Pointer to shared virtual orbitals
    double* erp = eps_avir_->pointer();

    boost::shared_ptr<PSIO> psio(new PSIO()); 
    psio_address next_AIA = PSIO_ZERO;

    psio->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);

    for (int A = 0; A < N-1; A++) {

        int astart = aocc_starts[A];
        int astop  = aocc_starts[A+1];
        int adelta = astop - astart;

        double* eap = eps_aocc_[A]->pointer();

        for (int a = 0; a < adelta; a += max_block) {

            int na = (a + max_block < adelta ? max_block : adelta - a); 
            next_AIA = psio_get_address(PSIO_ZERO,sizeof(double)*(a + astart)*navir*naux);
            psio->read(PSIF_DFMP2_AIA,"(Q|ia)",(char*)arQp[0],sizeof(double)*na*navir*naux,next_AIA,&next_AIA);

            for (int B = A + 1; B < N; B++) {

                int bstart = aocc_starts[B];
                int bstop  = aocc_starts[B+1];
                int bdelta = bstop - bstart;
                
                double* ebp = eps_aocc_[B]->pointer();

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

                        C_DGEMM('N','T',navir,navir,naux,1.0,arQp[aval + astart],naux,bsQp[bval + bstart],naux,0.0,Ip[0],navir);

                        double ea = eap[aval + astart];
                        double eb = ebp[bval + bstart];

                        for (int r = 0; r < navir; r++) {
                            double er = erp[r];
                            for (int s = 0; s < navir; s++) {
                                double es = erp[s];
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

    fprintf(outfile, "  ==> PT2 Analysis <==\n\n");
    
    fprintf(outfile,"    Dispersion (Uncoupled, Casimir-Polder):\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EJ_PT2->get(A,B);
            fprintf(outfile, "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*EJ, pc_hartree2kcalmol*EJ);
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");
    
    fprintf(outfile,"    Exchange-Dispersion (Tunneling):\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EK_PT2->get(A,B);
            fprintf(outfile, "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*EK, pc_hartree2kcalmol*EK);
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");
    
    fprintf(outfile,"    Total:\n\n");
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18s %18s\n",
        "Interaction", "mH", "kcal mol^-1");
    fprintf(outfile,"    -------------------------------------------------\n");
    for (int A = 0; A < monomers_.size(); A++) {
        for (int B = A+1; B < monomers_.size(); B++) {
            double Eval = EJ_PT2->get(A,B) + EK_PT2->get(A,B);
            fprintf(outfile, "    %-3d <-> %3d %18.12f %18.12f\n",
                A+1, B+1, 1.0E3*Eval, pc_hartree2kcalmol*Eval);
        }
    }    
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"    %-11s %18.12f %18.12f\n",
        "Total", 1.0E3*(EJ+EK), pc_hartree2kcalmol*(EJ+EK));
    fprintf(outfile,"    -------------------------------------------------\n");
    fprintf(outfile,"\n");

    energies_["Disp"] = EJ + EK;

    fflush(outfile);
}
void INFSAPT::print_trailer()
{
    fprintf(outfile, "  ALL TERMS:\n\n");

    energies_["Total"] = energies_["Elst"] + energies_["Exch"] + energies_["Ind"] + energies_["Disp"];
    
    fprintf(outfile,"    ----------------------------------------------------\n");
    fprintf(outfile,"    %-14s %18s %18s\n",
        "Term", "mH", "kcal mol^-1");
    fprintf(outfile,"    ----------------------------------------------------\n");
    fprintf(outfile,"    %-14s %18.12f %18.12f\n",
        "Electrostatics", 1.0E3*energies_["Elst"], pc_hartree2kcalmol*energies_["Elst"]);
    fprintf(outfile,"    %-14s %18.12f %18.12f\n",
        "Exchange", 1.0E3*energies_["Exch"], pc_hartree2kcalmol*energies_["Exch"]);
    fprintf(outfile,"    %-14s %18.12f %18.12f\n",
        "Induction", 1.0E3*energies_["Ind"], pc_hartree2kcalmol*energies_["Ind"]);
    fprintf(outfile,"    %-14s %18.12f %18.12f\n",
        "Dispersion", 1.0E3*energies_["Disp"], pc_hartree2kcalmol*energies_["Disp"]);
    fprintf(outfile,"    ----------------------------------------------------\n");
    fprintf(outfile,"    %-14s %18.12f %18.12f\n",
        "Total", 1.0E3*energies_["Total"], pc_hartree2kcalmol*energies_["Total"]);
    fprintf(outfile,"    ----------------------------------------------------\n");
    fprintf(outfile,"\n");

    fprintf(outfile, "    \"You can observe a lot by just watching.\"\n");
    fprintf(outfile, "                      --Yogi Berra\n");
    fprintf(outfile, "\n");

    fflush(outfile);
}

void INFSAPT::build_iaQ(
        boost::shared_ptr<BasisSet> primary,
        boost::shared_ptr<BasisSet> auxiliary,
        boost::shared_ptr<Matrix> Caocc, 
        boost::shared_ptr<Matrix> Cavir)
{
    // PSIO object
    boost::shared_ptr<PSIO> psio(new PSIO()); 

    // ==> Target <== //

    psio->open(PSIF_DFMP2_AIA,PSIO_OPEN_NEW);

    // => Global Sizing <= //

    int nso = primary->nbf();
    int naux = auxiliary->nbf();
    int naocc = Caocc->colspi()[0];
    int navir = Cavir->colspi()[0];

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
    double** Caoccp = Caocc->pointer();
    double** Cavirp = Cavir->pointer();

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
            C_DGEMM('T','N',naocc,navir,nso,1.0,Amip[row],naocc,Cavirp[0],navir,0.0,Aiap[row],navir);
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

    } // End Metrix

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
    psio->open(PSIF_DFMP2_AIA, PSIO_OPEN_OLD);
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

}}
