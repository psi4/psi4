#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cstring>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <lib3index/3index.h>

#ifdef HAVE_MADNESS

#include "mad_mp2.h"

#ifdef _OPENMP
#include <omp.h>
#endif


#define ORDER_PRINT_START for(int dummyproc = 0; dummyproc < nproc_; ++dummyproc){\
                              Communicator::world->sync();\
                              if(dummyproc == rank_){
//                              MPI_Barrier(MPI_COMM_WORLD);\

#define ORDER_PRINT_END  }Communicator::world->sync();}
//#define ORDER_PRINT_START ;
//#define ORDER_PRINT_END ;

using namespace std;
using namespace psi;
using namespace boost;

namespace psi{ namespace mad_mp2 {

MAD_MP2::MAD_MP2(Options& options, boost::shared_ptr<PSIO> psio) :
    Wavefunction(options, psio), madness::WorldObject<MAD_MP2>(*Communicator::world->get_madworld())
{
    nproc_   = Communicator::world->nproc();
    mad_nthread_ = Communicator::world->nthread();
    rank_    = Communicator::world->me();
    comm_    = Communicator::world->communicator();
#ifdef HAVE_MADNESS
    madworld_ = Communicator::world->get_madworld();
    mutex_ = Communicator::world->get_mutex();
#endif

    sleep(options.get_int("MADMP2_SLEEP"));
    common_init();
    parallel_init();

    if (debug_ > 2) {
        ORDER_PRINT_START
        printf("Printing C from node %d\n", rank_);
        Ca_->print();
        AO2USO_->print();
        epsilon_a_->print();
        fflush(outfile);
        printf("Start of MadMP2\n");
        ORDER_PRINT_END
    }

    process_pending();

}
MAD_MP2::~MAD_MP2()
{
}
void MAD_MP2::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("MADMP2_DEBUG");

    omp_nthread_ = 1;
    #ifdef _OPENMP
        omp_nthread_ = omp_get_max_threads();
    #endif

    scale_ss_ = options_.get_double("SCALE_SS");
    scale_os_ = options_.get_double("SCALE_OS");

    reference_ = Process::environment.reference_wavefunction();

    if (!reference_.get()) {

        Eref_ = chkpt_->rd_escf();
        nirrep_ = chkpt_->rd_nirreps();
        nso_ = chkpt_->rd_nso();
        nmo_ = chkpt_->rd_nmo();

        int* doccpi = chkpt_->rd_clsdpi();
        int* nmopi  = chkpt_->rd_orbspi();
        int* nsopi  = chkpt_->rd_sopi();
        int* frzcpi = chkpt_->rd_frzcpi();
        int* frzvpi = chkpt_->rd_frzvpi();

        for (int h = 0; h < nirrep_; h++) {
            doccpi_[h] = doccpi[h];
            soccpi_[h] = 0;
            nalphapi_[h] = doccpi[h];
            nbetapi_[h] = doccpi[h];
            nsopi_[h] = nsopi[h];
            nmopi_[h] = nmopi[h];
            frzcpi_[h] = frzcpi[h];
            frzvpi_[h] = frzvpi[h];
        }

        free(doccpi);
        free(nmopi);
        free(nsopi);
        free(frzcpi);
        free(frzvpi);

        Ca_ = boost::shared_ptr<Matrix>(new Matrix("Ca", nirrep_, nsopi_, nmopi_));
        epsilon_a_ = boost::shared_ptr<Vector>(new Vector("Evals a", nirrep_, nmopi_));

        double **tempmat = chkpt_->rd_scf();
        Ca_->set(tempmat);
        free_block(tempmat);
        double *temparr = chkpt_->rd_evals();
        epsilon_a_->set(temparr);
        free(temparr);
        if (chkpt_->rd_iopen()) {
            throw PSIEXCEPTION("DFMP2 MADNESS is only closed-shell for now");
        }

    }else{
        // Might be good place for a copy constructor
        Ca_ = SharedMatrix(reference_->Ca()->clone());
        Cb_ = SharedMatrix(reference_->Cb()->clone());
        epsilon_a_ = reference_->epsilon_a();
        epsilon_b_ = reference_->epsilon_b();

        Ca_->bcast(0);
        Cb_->bcast(0);
        epsilon_a_->bcast(0);
        epsilon_b_->bcast(0);

        nirrep_ = reference_->nirrep();
        nso_    = reference_->nso();
        nmo_    = reference_->nmo();
        Eref_   = reference_->reference_energy();

        ::memcpy(static_cast<void*>(doccpi_), static_cast<void*>(reference_->doccpi()), 8 * sizeof(int));
        ::memcpy(static_cast<void*>(soccpi_), static_cast<void*>(reference_->soccpi()), 8 * sizeof(int));
        ::memcpy(static_cast<void*>(nsopi_), static_cast<void*>(reference_->nsopi()), 8 * sizeof(int));
        ::memcpy(static_cast<void*>(nmopi_), static_cast<void*>(reference_->nmopi()), 8 * sizeof(int));
        ::memcpy(static_cast<void*>(nalphapi_), static_cast<void*>(reference_->nalphapi()), 8 * sizeof(int));
        ::memcpy(static_cast<void*>(nbetapi_), static_cast<void*>(reference_->nbetapi()), 8 * sizeof(int));
        ::memcpy(static_cast<void*>(frzcpi_), static_cast<void*>(reference_->frzcpi()), 8 * sizeof(int));
        ::memcpy(static_cast<void*>(frzvpi_), static_cast<void*>(reference_->frzvpi()), 8 * sizeof(int));
        // End copy constructor opportunity
        if (!reference_->restricted()) {
            throw PSIEXCEPTION("DFMP2 MADNESS is only closed-shell for now");
        }
    }


    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral));
    AO2USO_ = boost::shared_ptr<Matrix>(pet->aotoso());



    nfocc_  = 0;
    nfvir_  = 0;
    naocc_  = 0;
    navir_  = 0;
    nalpha_ = 0;
    nbeta_  = 0;
    for (int h = 0; h < nirrep_; h++) {
        naocc_  += doccpi_[h] - frzcpi_[h];
        navir_  += nmopi_[h] - doccpi_[h] - frzvpi_[h];
        nfocc_  += frzcpi_[h];
        nfvir_  += frzvpi_[h];
        nalpha_ += nalphapi_[h];
        nbeta_  += nbetapi_[h];
    }

    ::memset(static_cast<void*> (naoccpi_), '\0', 8 * sizeof(int));
    ::memset(static_cast<void*> (navirpi_), '\0', 8 * sizeof(int));
    ::memset(static_cast<void*> (offset_aocc_), '\0', 8 * sizeof(int));
    ::memset(static_cast<void*> (offset_avir_), '\0', 8 * sizeof(int));

    for (int h = 0; h < nirrep_; h++) {
        naoccpi_[h] = doccpi_[h] - frzcpi_[h];
        navirpi_[h] = nmopi_[h] - doccpi_[h] - frzvpi_[h];
    }

    for (int h = 1; h < nirrep_; h++) {
        offset_aocc_[h] = naoccpi_[h - 1] + offset_aocc_[h - 1];
        offset_avir_[h] = navirpi_[h - 1] + offset_avir_[h - 1];
    }

    if (debug_) {
        for (int h = 0; h < nirrep_; h++) {
            fprintf(outfile, "  h = %1d: naocc = %4d aocc_off = %4d navir = %4d avir_off = %4d\n",
                h, naoccpi_[h], offset_aocc_[h], navirpi_[h], offset_avir_[h]);
        }
    }

    Caocc_ = boost::shared_ptr<Matrix>(new Matrix("Caocc", nso_, naocc_));
    Cavir_ = boost::shared_ptr<Matrix>(new Matrix("Cavir", nso_, navir_));
    eps_aocc_ = boost::shared_ptr<Vector>(new Vector("eps_aocc", naocc_));
    eps_avir_ = boost::shared_ptr<Vector>(new Vector("eps_avir", navir_));
    irrep_aocc_ = boost::shared_ptr<IntVector>(new IntVector("irrep_aocc", naocc_));
    irrep_avir_ = boost::shared_ptr<IntVector>(new IntVector("irrep_avir", navir_));

    int occ_counter = 0;
    int vir_counter = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int i = frzcpi_[h]; i < doccpi_[h]; i++) {
            eps_aocc_->set(0, occ_counter, epsilon_a_->get(h, i));
            C_DGEMV('N',nso_,nsopi_[h],1.0,AO2USO_->pointer(h)[0],nsopi_[h],&Ca_->pointer(h)[0][i],nmopi_[h],0.0,&Caocc_->pointer(0)[0][occ_counter],naocc_);
            irrep_aocc_->set(0, occ_counter, h);
            occ_counter++;
        }
        for (int a = doccpi_[h]; a < nmopi_[h] - frzvpi_[h]; a++) {
            eps_avir_->set(0, vir_counter, epsilon_a_->get(h, a));
            C_DGEMV('N',nso_,nsopi_[h],1.0,AO2USO_->pointer(h)[0],nsopi_[h],&Ca_->pointer(h)[0][a],nmopi_[h],0.0,&Cavir_->pointer(0)[0][vir_counter],navir_);
            irrep_avir_->set(0, vir_counter, h);
            vir_counter++;
        }
    }
    AO2USO_.reset();

    if (debug_ > 2) {
        Caocc_->print();
        Cavir_->print();
        eps_aocc_->print();
        eps_avir_->print();
        irrep_aocc_->print(outfile);
        irrep_avir_->print(outfile);
    }


    // Auxiliary basis information
    auxiliary_automatic_ = false;
    if (options_.get_str("RI_BASIS_MP2") == "") {
        auxiliary_automatic_ = true;
        basisset_->molecule()->set_basis_all_atoms(options_.get_str("BASIS") + "-RI", "RI_BASIS_MP2");
    }

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    auxiliary_ = BasisSet::construct(parser, basisset_->molecule(), "RI_BASIS_MP2");
    parser.reset();

    naux_ = auxiliary_->nbf();

    boost::shared_ptr<IntegralFactory> integral2(new IntegralFactory(auxiliary_,auxiliary_,auxiliary_,auxiliary_));
    boost::shared_ptr<PetiteList> pet2(new PetiteList(auxiliary_, integral2));
    AO2USO_aux_ = boost::shared_ptr<Matrix>(pet2->aotoso());

    if (debug_ > 2) {
        AO2USO_aux_->print();
        fflush(outfile);
    }

    Dimension dim = pet2->SO_basisdim();
    ::memset(static_cast<void*>(nauxpi_), '\0', 8 * sizeof(int));
    for (int h = 0; h < nirrep_; h++) {
        nauxpi_[h] = dim[h];
    }

    zero_ = BasisSet::zero_ao_basis_set();

    max_nauxpi_ = 0;
    max_naoccpi_ = 0;
    max_navirpi_ = 0;

    for (int h = 0; h < nirrep_; h++) {
        max_nauxpi_ += nauxpi_[h];
        max_naoccpi_ += naoccpi_[h];
        max_navirpi_ += navirpi_[h];
    }
}
void MAD_MP2::parallel_init()
{
    nia_ = naocc_ * (ULI) navir_;

    int na_delta;
    int na_delta_extra;
    int na_deltac;
    int na_deltac_extra;
    std::vector<int> na_delta_per;
    std::vector<int> na_deltac_per;
    std::vector<int> a_delta_owner;
    std::vector<int> a_deltac_owner;

    // Check for hipsters
    if (nproc_ > naocc_ * (ULI) navir_)
        throw PSIEXCEPTION("What are you thinking?");

    // Maximum blocks
    if (nproc_ <= naocc_) {
        // nproc_ <= naocc_ Case
        int ni_per_proc = naocc_ / nproc_;
        int ni_extra = naocc_ % nproc_;

        std::vector<int> ni_per;
        for (int ind = 0; ind < nproc_; ind++) {
            if (ind < ni_extra)
                ni_per.push_back(ni_per_proc + 1);
            else
                ni_per.push_back(ni_per_proc);
        }

        std::vector<int> i_owner;
        for (int ind = 0; ind < nproc_; ind++) {
            for (int ind2 = 0; ind2 < ni_per[ind]; ind2++) {
                i_owner.push_back(ind);
            }
        }

        ULI counter = 0;
        for (ULI ia = 0; ia < nia_; ia++) {
            int i = ia / navir_;
            int a = ia % navir_;

            ia_owner_.push_back(i_owner[i]);
            if (i_owner[i] == rank_) {
                ia_local_to_global_.push_back(ia);
                ia_global_to_local_[ia]= counter++;
            }
        }

        int counter3 = 0;
        for (int i = 0; i < naocc_; i++) {
            if (i_owner[i] == rank_)
                i_global_to_local_[i] = counter3++;

            std::vector<int> blank1;
            blank1.push_back(i_owner[i]);
            ablock_owner_.push_back(blank1);

            std::vector<int> blank2;
            blank2.push_back(0);
            ablock_start_.push_back(blank2);

            std::vector<int> blank3;
            blank3.push_back(navir_);
            ablock_size_.push_back(blank3);
        }

    } else {
        // nproc_ > naocc_ Case
        int na_per_i = nproc_ / naocc_ + (nproc_ % naocc_ == 0 ? 0 : 1);
        int ni_delta = na_per_i * naocc_ - nproc_;

        na_delta = navir_ / (na_per_i - 1);
        na_delta_extra = navir_ % (na_per_i - 1);
        na_deltac = navir_ / (na_per_i);
        na_deltac_extra = navir_ % (na_per_i);

        for (int ind = 0; ind < na_per_i - 1; ind++) {
            if (ind < na_delta_extra)
                na_delta_per.push_back(na_delta + 1);
            else
                na_delta_per.push_back(na_delta);
        }
        for (int ind = 0; ind < na_per_i; ind++) {
            if (ind < na_deltac_extra)
                na_deltac_per.push_back(na_deltac + 1);
            else
                na_deltac_per.push_back(na_deltac);
        }

        for (int ind = 0; ind < na_per_i - 1; ind++) {
            for (int ind2 = 0; ind2 < na_delta_per[ind]; ind2++) {
                a_delta_owner.push_back(ind);
            }
        }

        for (int ind = 0; ind < na_per_i; ind++) {
            for (int ind2 = 0; ind2 < na_deltac_per[ind]; ind2++) {
                a_deltac_owner.push_back(ind);
            }
        }

        ULI counter = 0L;
        int start_proc = 0;
        for (int i = 0 ; i < naocc_; i++) {
            for (int a = 0; a < navir_; a++) {
                int owner_proc;
                if (i < ni_delta) {
                    owner_proc = start_proc + a_delta_owner[a];
                    // In the N - 1 block region
                } else {
                    // In the N block region
                    owner_proc = start_proc + a_deltac_owner[a];
                }

                if (owner_proc == rank_)
                    i_global_to_local_[i] = 0;

                ia_owner_.push_back(owner_proc);
                if (rank_ == owner_proc) {
                    ia_local_to_global_.push_back(i * navir_ + a);
                    ia_global_to_local_[i * navir_ + a] = counter++;
                }
            }

            std::vector<int> blank1;
            std::vector<int> blank2;
            std::vector<int> blank3;

            if (i < ni_delta) {
                int counter2 = 0;
                for (int ind = 0; ind < na_delta_per.size(); ind++) {
                    blank1.push_back(start_proc + ind);
                    blank2.push_back(counter2 + na_delta_per[ind]);
                    blank3.push_back(na_delta_per[ind]);
                    counter += na_delta_per[ind];
                }
            } else {
                int counter2 = 0;
                for (int ind = 0; ind < na_deltac_per.size(); ind++) {
                    blank1.push_back(start_proc + ind);
                    blank2.push_back(counter2 + na_deltac_per[ind]);
                    blank3.push_back(na_deltac_per[ind]);
                    counter += na_deltac_per[ind];
                }
            }

            ablock_owner_.push_back(blank1);
            ablock_start_.push_back(blank2);
            ablock_size_.push_back(blank3);

            if (i < ni_delta) start_proc += na_per_i - 1;
            else start_proc += na_per_i;
        }
    }
    nia_local_ = ia_local_to_global_.size();

    std::set<int> unique_i;
    std::set<int> unique_a;
    for (int ia_local = 0 ; ia_local < nia_local_; ia_local++) {
        int ia_global = ia_local_to_global_[ia_local];
        int i = ia_global / navir_;
        int a = ia_global % navir_;

        unique_i.insert(i);
        unique_a.insert(a);
    }
    naocc_local_ = unique_i.size();
    navir_local_ = unique_a.size();

    for (std::set<int>::iterator it = unique_i.begin(); it != unique_i.end(); it++)
        aocc_local_.push_back(*it);
    for (std::set<int>::iterator it = unique_a.begin(); it != unique_a.end(); it++)
        avir_local_.push_back(*it);

    std::sort(aocc_local_.begin(), aocc_local_.end());
    std::sort(avir_local_.begin(), avir_local_.end());
    if (debug_) {
        ORDER_PRINT_START

        printf("Number of processors: %d\n", nproc_);
        printf("Number of threads:    %d\n", omp_nthread_);
        printf("Current processors:   %d\n", rank_);
        printf("Communicator Type:    %s\n", comm_.c_str());
        printf("na_delta:             %d\n", na_delta);
        printf("na_delta_extra:       %d\n", na_delta_extra);
        printf("na_deltac:            %d\n", na_deltac);
        printf("na_deltac_extr        %d\n", na_deltac_extra);

        std::vector<int>::const_iterator iter = na_delta_per.begin();
        printf("na_delta_per: ");
        for(; iter != na_delta_per.end(); ++iter)
            printf("%d ", *iter);
        printf("\n\n");

        iter = na_deltac_per.begin();
        printf("na_deltac_per: ");
        for(; iter != na_deltac_per.end(); ++iter)
            printf("%d ", *iter);
        printf("\n\n");

        iter = a_delta_owner.begin();
        printf("a_delta_owner: ");
        for(; iter != a_delta_owner.end(); ++iter)
            printf("%d ", *iter);
        printf("\n\n");

        iter = a_deltac_owner.begin();
        printf("a_deltac_owner: ");
        for(; iter != a_deltac_owner.end(); ++iter)
            printf("%d ", *iter);
        printf("\n\n");

        printf("IA Owner array: ");
        for(int ia = 0; ia < nia_; ia++)
            printf("%d ", ia_owner_[ia]);
        printf("\n\n");

        std::map<ULI, int>::const_iterator global_it = ia_global_to_local_.begin();
        printf("Global to local array: ");
        for(; global_it != ia_global_to_local_.end(); ++global_it)
            printf("%ld -> %d ", global_it->first, global_it->second);
        printf("\n\n");

        std::vector<ULI>::const_iterator local_it = ia_local_to_global_.begin();
        printf("Local to global array: ");
        for(; local_it != ia_local_to_global_.end(); ++local_it)
            printf("%ld ", *local_it);
        printf("\n\n");

        printf("There are %d local occ and %d local vir orbitals, giving %ld local pairs\n",
              naocc_local_, navir_local_, nia_local_);

        std::vector<int>::const_iterator aocc_it = aocc_local_.begin();
        printf("aocc_local values: ");
        for(; aocc_it != aocc_local_.end(); ++aocc_it)
            printf("%d ", *aocc_it);
        printf("\n\n");

        std::vector<int>::const_iterator avir_it = avir_local_.begin();
        printf("avir_local values: ");
        for(; avir_it != avir_local_.end(); ++avir_it)
            printf("%d ", *avir_it);
        printf("\n\n");

        std::map<int, int>::const_iterator i_global_it = i_global_to_local_.begin();
        printf("i Global to local array: ");
        for(; i_global_it != i_global_to_local_.end(); ++i_global_it)
            printf("%d -> %d ", i_global_it->first, i_global_it->second);
        printf("\n\n");

        ORDER_PRINT_END
    }
}
void MAD_MP2::print_header()
{
    fprintf(outfile, "\n");
    fprintf(outfile, "         ------------------------------------------------------------\n");
    fprintf(outfile, "                                DFMP2 MADNESS\n");
    fprintf(outfile, "                          %8s Implementation\n", (options_.get_bool("PARALLEL") ? "MADNESS" : "SERIAL"));
    fprintf(outfile, "                              %6s Algorithm\n", options_.get_str("MP2_ALGORITHM").c_str());
    fprintf(outfile, "                        %3d Threads, %6ld MiB Core\n", omp_nthread_, memory_ / 1000000L);
    fprintf(outfile, "          Ben Mintz, Rob Parrish, Andy Simmonnett, and Justin Turney\n");
    fprintf(outfile, "         ------------------------------------------------------------\n\n");

    fprintf(outfile, " ==> Geometry <==\n\n");
    molecule_->print();
    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    fprintf(outfile, "  Reference energy  = %20.15f\n\n", Eref_);

    fprintf(outfile, "  ==> Primary Basis <==\n\n");
    basisset_->print_by_level(outfile, print_);

    fprintf(outfile, "  ==> Auxiliary Basis <==\n\n");
    if (auxiliary_automatic_) {
        fprintf(outfile, "  No auxiliary basis selected, defaulting to %s-RI\n\n", options_.get_str("BASIS").c_str());
    }
    auxiliary_->print_by_level(outfile, print_);

    fprintf(outfile, "  ==> Orbital Dimensions <==\n\n");
    CharacterTable ct = molecule_->point_group()->char_table();
    fprintf(outfile, "   ------------------------------------------------------------------\n");
    fprintf(outfile, "    Irrep    %6s %6s %6s %6s %6s %6s %6s %6s\n",
            "Nso", "Nmo", "Nfocc", "Nocc", "Naocc", "Navir", "Nvir", "Nfvir");
    fprintf(outfile, "   ------------------------------------------------------------------\n");
    for (int h= 0; h < nirrep_; h++) {
        fprintf(outfile, "    %-3s      %6d %6d %6d %6d %6d %6d %6d %6d\n",
            ct.gamma(h).symbol(), nsopi_[h], nmopi_[h], frzcpi_[h], doccpi_[h],
            doccpi_[h] - frzcpi_[h], nmopi_[h] - doccpi_[h] - frzvpi_[h], nmopi_[h] - doccpi_[h],
            frzvpi_[h]);
    }
    fprintf(outfile, "   ------------------------------------------------------------------\n");
    fprintf(outfile, "    Total    %6d %6d %6d %6d %6d %6d %6d %6d\n",
        nso_, nmo_, nfocc_, nfocc_ + naocc_, naocc_, navir_, navir_ + nfvir_, nfvir_);
    fprintf(outfile, "   ------------------------------------------------------------------\n\n");
    fflush(outfile);
}
void MAD_MP2::check_memory()
{
    #define MEM_SAFETY 0.8

    fprintf(outfile, "  ==> Memory Checking <==\n\n");

    int max_pshell = auxiliary_->max_function_per_shell();

    // Required memory, doubles
    ULI required = 0L;
    ULI overhead = 0L;

    // Tensor Sizes
    ULI Qia_mem_AO = naux_ * (ULI) naocc_ * navir_;
    ULI J_mem_AO   = naux_ * (ULI) naux_;

    ULI J_mem_SO   = 0L;
    for (int h = 0; h < nirrep_; h++) {
        J_mem_SO += nauxpi_[h] * (ULI) nauxpi_[h];
    }
    ULI Qia_mem_SO = 0L;
    for (int hQ = 0; hQ < nirrep_; hQ++) {
        for (int hi = 0; hi < nirrep_; hi++) {
            int ha = hi ^ hQ;
            Qia_mem_SO += nauxpi_[hQ] * naoccpi_[hi] * navirpi_[ha];
        }
    }
    ULI I_SO = 0L;
    for (int hi = 0; hi < nirrep_; hi++) {
        for (int hj = 0; hj < nirrep_; hj++) {
            ULI trial = 0L;
            for (int hQ = 0; hQ < nirrep_; hQ++) {
                int ha = hi ^ hQ;
                int hb = hj ^ hQ;
                trial += navirpi_[ha] * (ULI) navirpi_[hb];
            }
            I_SO = (I_SO > trial ? I_SO : trial);
        }
    }

    // Rate Limiting Operations
    ULI Amn_mem = naux_ * (ULI) (nso_ * (ULI) nso_ +
                                 nso_ * (ULI) naocc_) +
        Qia_mem_AO; // (A|ia)
    ULI Aia_USO_mem = Qia_mem_SO + Qia_mem_AO +
        max_naoccpi_ * (ULI) max_navirpi_ * naux_; // (A'|ia)
    ULI Aia_mem = Qia_mem_SO + J_mem_SO +
        max_nauxpi_ * (ULI) max_nauxpi_; // (Q'|ia)
    ULI J_USO_mem = J_mem_AO + J_mem_SO +
        naux_ * (ULI) max_nauxpi_ + Qia_mem_SO; // (A'|B')
    ULI I_mem = omp_nthread_ * I_SO + Qia_mem_SO; // Conventional energy
    if (options_.get_str("MP2_ALGORITHM") == "DFMP2J")
        I_mem = 2L * Qia_mem_SO + J_mem_SO; // Rough estimate for now

    required = (required > Amn_mem ? required : Amn_mem);
    if (nirrep_ > 1) {
        required = (required > Aia_USO_mem ? required : Aia_USO_mem);
        required = (required > J_USO_mem ? required : J_USO_mem);
    }
    required = (required > Aia_mem ? required : Aia_mem);
    required = (required > I_mem ? required : I_mem);

    //Overhead
    overhead += J_mem_AO; // AO2USO_aux_
    required += overhead;

    fprintf(outfile, "    ------------------------------------------------------------\n");
    fprintf(outfile, "     %-20s %18s %18s\n", "Tensor/Operation", "Doubles", "MiB");
    fprintf(outfile, "    ------------------------------------------------------------\n");
    fprintf(outfile, "     %-20s %18ld %18ld\n", "(A|B) Tensor AO", J_mem_AO, J_mem_AO * 8L / 1000000L);
    fprintf(outfile, "     %-20s %18ld %18ld\n", "(A|B) Tensor USO", J_mem_SO, J_mem_SO * 8L / 1000000L);
    fprintf(outfile, "     %-20s %18ld %18ld\n", "(Q|ia) Tensor AO", Qia_mem_AO, Qia_mem_AO * 8L / 1000000L);
    fprintf(outfile, "     %-20s %18ld %18ld\n", "(Q|ia) Tensor USO", Qia_mem_SO, Qia_mem_SO * 8L / 1000000L);
    fprintf(outfile, "     %-20s %18ld %18ld\n", "I_ab Tensor USO", I_SO, I_SO * 8L / 1000000L);
    fprintf(outfile, "     %-20s %18ld %18ld\n", "MO Transform (A|ia)", Amn_mem, Amn_mem * 8L / 1000000L);
    if (nirrep_ > 1) {
        fprintf(outfile, "     %-20s %18ld %18ld\n", "USO Transform (A|ia)", Aia_USO_mem, Amn_mem * 8L / 1000000L);
        fprintf(outfile, "     %-20s %18ld %18ld\n", "USO Transform (A|B)", J_USO_mem, J_USO_mem * 8L / 1000000L);
    }
    fprintf(outfile, "     %-20s %18ld %18ld\n", "Fitting Transform", Aia_mem, Aia_mem * 8L / 1000000L);
    fprintf(outfile, "     %-20s %18ld %18ld\n", "Energy Computation", I_mem, I_mem * 8L / 1000000L);
    fprintf(outfile, "     %-20s %18ld %18ld\n", "Overhead", overhead, overhead * 8L / 1000000L);
    fprintf(outfile, "    ------------------------------------------------------------\n");
    fprintf(outfile, "     %-20s %18ld %18ld\n", "Required Memory", required, required * 8L / 1000000L);
    fprintf(outfile, "     %-20s %18ld %18ld\n", "Available Memory", ((ULI)( MEM_SAFETY * memory_)) / 8L, ((ULI) (MEM_SAFETY * memory_)) / 1000000L);
    fprintf(outfile, "    ------------------------------------------------------------\n\n");
    fflush(outfile);

    if ((double)required > MEM_SAFETY * ((double) memory_ / 8L)) {
        throw PSIEXCEPTION("MAD_MP2 does not have enough memory");
    }
}
double MAD_MP2::compute_energy()
{
    print_header();

    if (options_.get_str("MP2_ALGORITHM") == "DFMP2J") {
        denominator();
    }

    check_memory();
    J();
    Jm12();

    Aia();

    Communicator::world->sync();

    if (options_.get_str("MP2_ALGORITHM") == "DFMP2")
        I();
    else if (options_.get_str("MP2_ALGORITHM") == "DFMP2J")
        IJ();

    print_energy();

    return energy_;
}
void MAD_MP2::J()
{
    // => AO Basis J <= //

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(auxiliary_,zero_,auxiliary_,zero_));
    boost::shared_ptr<TwoBodyAOInt> eri(fact->eri());
    const double* buffer = eri->buffer();

    boost::shared_ptr<Matrix> J(new Matrix("J^-1/2", naux_, naux_));
    double** Jp = J->pointer();

    timer_on("MP2 J AO");
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        for (int Q = 0; Q <= P ; Q++) {
            eri->compute_shell(P,0,Q,0);
            int nP     = auxiliary_->shell(P)->nfunction();
            int Pstart = auxiliary_->shell(P)->function_index();
            int nQ     = auxiliary_->shell(Q)->nfunction();
            int Qstart = auxiliary_->shell(Q)->function_index();
            for (int op = 0, index = 0; op < nP; op++) {
                for (int oq = 0; oq < nQ; oq++, index++) {
                    Jp[Pstart + op][Qstart + oq] = buffer[index];
                    Jp[Qstart + oq][Pstart + op] = buffer[index];
                }
            }
        }
    }
    timer_off("MP2 J AO");

    if (debug_ > 2) {
        fprintf(outfile, "  ==> After J Generation <==\n\n");
        J->print();
    }

    Jm12_ = J;
}
void MAD_MP2::Jm12()
{
    timer_on("MP2 J^-1/2");
    Jm12_->power(-1.0/2.0);
    timer_off("MP2 J^-1/2");

    if (debug_ > 2) {
        fprintf(outfile, "  ==> After J^-1/2 <==\n\n");
        Jm12_->print();
    }
}
void MAD_MP2::Aia()
{
    // => AO Basis A: (A|ia) <= //

    boost::shared_ptr<Matrix> Aia(new Matrix("(A|ia)", naux_, nia_local_));
    double** Aiap = Aia->pointer();

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(auxiliary_,zero_,basisset_,basisset_));
    boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[omp_nthread_];
    const double** buffer = new const double*[omp_nthread_];
    for (int thread = 0; thread < omp_nthread_; thread++) {
        eri[thread] = boost::shared_ptr<TwoBodyAOInt>(fact->eri());
        buffer[thread] = eri[thread]->buffer();
    }

    int max_pshell = auxiliary_->max_function_per_shell();

    boost::shared_ptr<Matrix> Amn(new Matrix("(A|mn) Chunk", max_pshell, nso_ * (ULI) nso_));
    double** Amnp = Amn->pointer();
    boost::shared_ptr<Matrix> Ami(new Matrix("(A|mi) Chunk", max_pshell, nso_ * (ULI) naocc_local_));
    double** Amip = Ami->pointer();

    boost::shared_ptr<Matrix> Ci(new Matrix("C_mi local", nso_, naocc_local_));
    double** Cip = Ci->pointer();
    for (int ind = 0; ind < naocc_local_; ind++) {
        int i = aocc_local_[ind];
        C_DCOPY(nso_, &Caocc_->pointer(0)[0][i], naocc_, &Cip[0][ind], naocc_local_);
    }

    boost::shared_ptr<Matrix> Ca(new Matrix("C_na local", nso_, navir_local_));
    double** Cap = Ca->pointer();
    for (int ind = 0; ind < navir_local_; ind++) {
        int a = avir_local_[ind];
        C_DCOPY(nso_, &Cavir_->pointer(0)[0][a], navir_, &Cap[0][ind], navir_local_);
    }

    if(debug_ > 1){
        ORDER_PRINT_START
        printf("Ci for node %d\n", rank_);
        Ci->print();
        printf("Ca for node %d\n", rank_);
        Ca->print();
        ORDER_PRINT_END
    }
    // Schwarz Sieve object
    boost::shared_ptr<SchwarzSieve> schwarz(new SchwarzSieve(basisset_, options_.get_double("SCHWARZ_CUTOFF")));
    long int* schwarz_shells = schwarz->get_schwarz_shells_reverse();

    timer_on("MP2 AO -> MO");
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P)->nfunction();
        int pstart = auxiliary_->shell(P)->function_index();

        ::memset(static_cast<void*>(Amnp[0]), '\0', nP * (ULI) nso_ * nso_ * sizeof(double));
        ::memset(static_cast<void*>(Amip[0]), '\0', nP * (ULI) nso_ * naocc_local_ * sizeof(double));

        // (A|mn) block
        timer_on("MP2 (A|mn)");
        #pragma omp parallel for schedule(guided)
        for (int M = 0; M < basisset_->nshell(); M++) {

            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int nM = basisset_->shell(M)->nfunction();
            int mstart = basisset_->shell(M)->function_index();

            for (int N = 0; N <= M; N++) {
                int nN = basisset_->shell(N)->nfunction();
                int nstart = basisset_->shell(N)->function_index();

                // Schwarz
                if (schwarz_shells[M * (M + 1) / 2 + N] == -1) continue;

                eri[rank]->compute_shell(P,0,M,N);

                for (int op = 0, index = 0; op < nP; op++) {
                    for (int om = 0; om < nM; om++) {
                        for (int on = 0; on < nN; on++, index++) {
                            Amnp[op][(om + mstart) * nso_ + (on + nstart)] = buffer[rank][index];
                            Amnp[op][(on + nstart) * nso_ + (om + mstart)] = buffer[rank][index];
                        }
                    }
                }
            }
        }
        timer_off("MP2 (A|mn)");

        // (A|mi) block
        timer_on("MP2 (A|mi)");
        C_DGEMM('N','N',nP * (ULI) nso_, naocc_local_, nso_, 1.0, Amnp[0], nso_, Cip[0], naocc_local_, 0.0, Amip[0], naocc_local_);
        timer_off("MP2 (A|mi)");


        // (A|ia) block
        timer_on("MP2 (A|ia)");
        #pragma omp parallel for
        for (int op = 0; op < nP; op++) {
            C_DGEMM('T','N',naocc_local_,navir_local_,nso_,1.0,Amip[op],naocc_local_,Cap[0],navir_local_,0.0,Aiap[op + pstart], navir_local_);
        }


//        #pragma omp parallel for
//        for (int index = 0; index < nia_local_ * nP; index++) {
//            int ia_local = index / nP;
//            int op = index % nP;
//            int ia_global = ia_local_to_global_[ia_local];
//            int i = ia_global / navir_;
//            int a = ia_global % navir_;

//            Aiap[op + pstart][ia_local] = C_DDOT(nso_, &Amip[op][i], naocc_local_, &Cavir_->pointer(0)[0][a], navir_);
//        }
        timer_off("MP2 (A|ia)");
    }
    timer_off("MP2 AO -> MO");
    Amn.reset();
    Ami.reset();
    Ci.reset();
    Ca.reset();

    if (debug_ > 2) {
        ORDER_PRINT_START
        printf("  ==> After MO Transform, node %d <==\n\n", rank_);
        Aia->print();
        ORDER_PRINT_END
    }

    // => J^-1/2 (A|ia) Fitting <= //

    int n2 = (naux_ < nia_local_ ? naux_: nia_local_);
    boost::shared_ptr<Matrix> T(new Matrix("Temp", naux_, n2));
    double** Tp = T->pointer();
    double** Jp = Jm12_->pointer();

    for (int ia = 0; ia < nia_local_; ia+= n2)
    {
        int nmult = (ia + n2 > nia_local_ ? nia_local_ - ia : n2);
        for (int Q = 0; Q < naux_; Q++) {
            ::memcpy(static_cast<void*>(Tp[Q]), static_cast<void*>(&Aiap[Q][ia]), nmult * sizeof(double));
        }
        C_DGEMM('N','N',naux_,nmult,naux_,1.0,Jp[0],naux_,Tp[0],n2,0.0,&Aiap[0][ia],nia_local_);
    }

    if (debug_ > 2) {
        ORDER_PRINT_START
        printf("  ==> After Fitting <==\n\n");
        printf("A matrix from node %d\n", rank_);
        Aia->print();
        ORDER_PRINT_END
    }
    Aia_ = Aia;
}

madness::Future<std::vector<double> > MAD_MP2::fetch_Qia_block(const int& i, const int& ablock)
{
    mutex_->lock();

    int i_local = i_global_to_local_[i];

    if(debug_){
        printf("Copy Start, Rank %d: i = %d, i_local = %d from node %d\n", rank_, i, i_local, rank_);
        fflush(stdout);
    }

    int na = ablock_size_[i][ablock];
    std::vector<double> block(na * naux_);

    double** Qp = Aia_->pointer();

    ULI index = 0L;
    for (int Q = 0; Q < naux_; Q++) {
        for (int a = 0; a < na; a++, index++) {
            block[index] = Qp[Q][i_local * na + a];
        }
    }

    if(debug_){
        printf("Copy Done, Rank %d: i = %d, i_local = %d from node %d\n", rank_, i, i_local, rank_);
        fflush(stdout);
    }

    mutex_->unlock();

    return madness::Future<std::vector<double> >(block);
}
madness::Void MAD_MP2::unpack_Qia_block(const std::vector<double>& block, SharedMatrix Q, const int& astart, const int& asize, const int & i)
{
    mutex_->lock();


    if(debug_){
        printf("Unpacking i = %d on node %d block size %ld max index %d,"
               " astart %d asize %d\n",
               i, rank_, block.size(), naux_*asize -1, astart, asize);
        fflush(stdout);
    }

    double** Qp = Q->pointer();
    for (int P = 0; P < naux_; P++) {
        for (int a = 0; a < asize; a++) {
            Qp[P][astart + a] = block[P * asize + a];
        }
    }

    mutex_->unlock();

    return madness::None;
}

madness::Future<SharedMatrix> MAD_MP2::build_Qa(const int &i) {

    mutex_->lock();

    SharedMatrix Qa(new Matrix("Qa", naux_, navir_));
    double** Qap = Qa->pointer();
    double** Qiap = Aia_->pointer();

    std::vector<madness::Future<std::vector<double> > > fut(ablock_owner_[i].size());
    for (int ind = 0; ind < ablock_owner_[i].size(); ind++) {
        if (debug_)
            printf("rank = %3d, i = %3d: Must fetch from process %3d, a block of %4d a, starting at a = %4d\n",
                   rank_, i, ablock_owner_[i][ind], ablock_size_[i][ind], ablock_start_[i][ind]);

        int astart = ablock_start_[i][ind];
        int asize  = ablock_size_[i][ind];
        int aowner = ablock_owner_[i][ind];

        if(aowner == rank_){
            int i_local = i_global_to_local_[i];
            for (int Q = 0; Q < naux_; Q++)
                ::memcpy(&(Qap[Q][astart]), &(Qiap[Q][i_local*navir_local_]),
                         asize * sizeof(double));
        }else{
            fut[ind] = task(aowner, &MAD_MP2::fetch_Qia_block, i, ind);
            task(rank_, &MAD_MP2::unpack_Qia_block, fut[ind], Qa, astart, asize, i);
//            for (int P = 0; P < naux_; P++) {
//                for (int a = 0; a < asize; a++) {
//                    Qap[P][astart + a] = fut.get()[P * asize + a];
//                }
//            }
        }
    }

    mutex_->unlock();

    return madness::Future<SharedMatrix> (Qa);
}

madness::Future<SharedMatrix> MAD_MP2::build_Qb(const int &j) {

    mutex_->lock();

    SharedMatrix Qb(new Matrix("Qb", naux_, navir_));
    double** Qbp = Qb->pointer();
    double** Qiap = Aia_->pointer();

    std::vector<madness::Future<std::vector<double> > > fut(ablock_owner_[j].size());
    for (int ind = 0; ind < ablock_owner_[j].size(); ind++) {
        if (debug_)
            printf("rank = %3d, i = %3d: Must fetch from process %3d, a block of %4d a, starting at a = %4d\n",
                   rank_, j, ablock_owner_[j][ind], ablock_size_[j][ind], ablock_start_[j][ind]);

        int bstart = ablock_start_[j][ind];
        int bsize  = ablock_size_[j][ind];
        int bowner = ablock_owner_[j][ind];

        if(bowner == rank_){
            int j_local = i_global_to_local_[j];
            for (int Q = 0; Q < naux_; Q++)
                ::memcpy(&(Qbp[Q][bstart]), &(Qiap[Q][j_local*navir_local_]),
                         bsize * sizeof(double));
        }else{
            fut[ind] = task(bowner,&MAD_MP2::fetch_Qia_block,j,ind);
            task(rank_, &MAD_MP2::unpack_Qia_block,fut[ind], Qb, bstart, bsize, j);
//            for (int P = 0; P < naux_; P++) {
//                for (int b = 0; b < bsize; b++) {
//                    Qbp[P][bstart + b] = fut.get()[P * bsize + b];
//                }
//            }
        }
    }

    mutex_->unlock();

    return madness::Future<SharedMatrix> (Qb);
}
madness::Future<SharedMatrix> MAD_MP2::build_I(SharedMatrix Qa, SharedMatrix Qb)
{
    mutex_->lock();

    boost::shared_ptr<Matrix> I(new Matrix("I", navir_, navir_));

    C_DGEMM('T','N',navir_,navir_,naux_, 1.0, Qa->pointer()[0], navir_,
            Qb->pointer()[0], navir_, 0.0, I->pointer()[0], navir_);

    mutex_->unlock();

    return madness::Future<SharedMatrix>(I);
}

void MAD_MP2::I()
{

    timer_on("Parallel_I");
    E_MP2J_ = 0.0;
    E_MP2K_ = 0.0;

    boost::shared_ptr<Matrix> I(new Matrix("I", navir_, navir_));
    double** Ip = I->pointer();
    boost::shared_ptr<Matrix> Qa(new Matrix("Qa", naux_, navir_));
    double** Qap = Qa->pointer();
    boost::shared_ptr<Matrix> Qb(new Matrix("Qb", naux_, navir_));
    double** Qbp = Qb->pointer();

    double** Qiap = Aia_->pointer();

    std::map<int,int> ij_i_map, ij_j_map;
    int counter = 0L;
    for (int i = 0, ij=0 ; i < naocc_; i++) {
        for (int j = 0; j <= i; j++, ij++) {
            ij_i_map.insert(std::pair<int,int>(counter, i));
            ij_j_map.insert(std::pair<int,int>(counter, j));
            counter++;
        }
    }

    for (int ij = 0; ij < naocc_ * (naocc_ + 1) / 2; ij++) {

        if (ij % nproc_ == rank_) {

            int i = ij_i_map[ij];
            int j = ij - i * (i+1) / 2;

            madness::Future<SharedMatrix> Qa_fut = task(rank_, &MAD_MP2::build_Qa, i);
            madness::Future<SharedMatrix> Qb_fut = task(rank_, &MAD_MP2::build_Qb, j);


            boost::shared_ptr<Matrix> I_fut = task(rank_, &MAD_MP2::build_I, Qa_fut, Qb_fut);

            if (debug_ > 2) {
                Qa->print();
                Qb->print();
            }

            E_MP2J_ += task(rank_, &MAD_MP2::energy_j, I_fut, i,j);
            E_MP2K_ += task(rank_, &MAD_MP2::energy_k, I_fut, i,j);

        }
    }

    Communicator::world->sync();
    Communicator::world->sum(&E_MP2J_, 1);
    Communicator::world->sum(&E_MP2K_, 1);

    timer_off("Parallel_I");

}
madness::Future<double> MAD_MP2::energy_j(const SharedMatrix I, const int &i, const int &j) {

    mutex_->lock();

    double e_mp2j;
    for (int a = 0; a < navir_; a++) {
        for (int b = 0; b < navir_; b++) {
            double iajb = I->get(0,a,b);
            double denom = (i == j ? 1.0 : 2.0) /
                    (eps_avir_->get(a) +
                     eps_avir_->get(b) -
                     eps_aocc_->get(i) -
                     eps_aocc_->get(j));
            e_mp2j -= 2.0 * denom * iajb * iajb;
        }
    }

    mutex_->unlock();

    return madness::Future<double>(e_mp2j);
}
madness::Future<double> MAD_MP2::energy_k(const SharedMatrix I, const int &i, const int &j)
{
    mutex_->lock();

    double e_mp2k=0.0;
    for (int a = 0; a < navir_; a++) {
        for (int b = 0; b < navir_; b++) {

            double iajb = I->get(0,a,b);
            double ibja = I->get(0,b,a);
            double denom = (i == j ? 1.0 : 2.0) /
                    (eps_avir_->get(a) +
                     eps_avir_->get(b) -
                     eps_aocc_->get(i) -
                     eps_aocc_->get(j));
            e_mp2k += 1.0 * denom * iajb * ibja;
        }
    }

    mutex_->unlock();

    return madness::Future<double>(e_mp2k);
}


void MAD_MP2::denominator()
{
    denom_ = boost::shared_ptr<Denominator>(Denominator::buildDenominator(
        options_.get_str("DENOMINATOR_ALGORITHM"), eps_aocc_, eps_avir_,
        options_.get_double("DENOMINATOR_DELTA")));
}
void MAD_MP2::IJ()
{
#if 0
    E_MP2J_ = 0.0;
    E_MP2K_ = 0.0;

    int nw = denom_->nvector();
    double** tau = denom_->denominator()->pointer();

    timer_on("MP2J Energy");
    for (int hQ = 0; hQ < nirrep_; hQ++) {
        for (int hi = 0; hi < nirrep_; hi++) {
            int ha = hi ^ hQ;

            int nQ = nauxpi_[hQ];
            int ni = naoccpi_[hi];
            int na = navirpi_[ha];

            if (!nQ || !ni || !na) continue;

            double** Qiap = Aia_[std::pair<int,int>(hQ,hi)]->pointer();
            boost::shared_ptr<Matrix> Qiaw(new Matrix("(Q|ia)^w)", nQ, ni * (ULI) na));
            double** Qiawp = Qiaw->pointer();

            boost::shared_ptr<Matrix> Z(new Matrix("Z^QQ", nQ, nQ));
            double** Zp = Z->pointer();

            for (int w = 0; w < nw; w++) {
                ::memcpy(static_cast<void*>(Qiawp[0]), static_cast<void*>(Qiap[0]), nQ * (ULI) na * ni * sizeof(double));
                for (int i = 0; i < ni; i++) {
                    for (int a = 0; a < na; a++) {
                        C_DSCAL(nQ, tau[w][(offset_aocc_[hi] + i) * navir_ + (offset_avir_[ha] + a)],
                            &Qiawp[0][i * na + a], ni * (ULI) na);
                    }
                }
                C_DGEMM('N','T',nQ,nQ,ni*(ULI)na,1.0,Qiawp[0],ni*(ULI)na,Qiap[0],ni*(ULI)na,0.0,Zp[0],nQ);
                E_MP2J_ -= 2.0 * C_DDOT(nQ * (ULI) nQ, Zp[0], 1, Zp[0], 1);
            }
        }
    }
    timer_off("MP2J Energy");
#endif
}
void MAD_MP2::print_energy()
{
    fprintf(outfile, "  ==> Energies <==\n\n");
    energies_["Reference Energy"]         = Eref_;
    energies_["MP2J Energy"]              = E_MP2J_;
    energies_["MP2K Energy"]              = E_MP2K_;
    energies_["Same-Spin Energy"]         = 0.5 * E_MP2J_ + E_MP2K_;
    energies_["Opposite-Spin Energy"]     = 0.5 * E_MP2J_;
    energies_["Correlation Energy"]       = E_MP2J_ + E_MP2K_;
    energies_["Total Energy"]             = Eref_ + E_MP2J_ + E_MP2K_;
    energies_["SCS Same-Spin Energy"]     = scale_ss_ * (0.5 * E_MP2J_ + E_MP2K_);
    energies_["SCS Opposite-Spin Energy"] = scale_os_ * (0.5 * E_MP2J_);
    energies_["SCS Correlation Energy"]   = energies_["SCS Opposite-Spin Energy"] + energies_["SCS Same-Spin Energy"];
    energies_["SCS Total Energy"]         = Eref_ + energies_["SCS Correlation Energy"];

    energy_ = energies_["Total Energy"];

    Process::environment.globals["CURRENT ENERGY"] = energy_;
    Process::environment.globals["E_MP2J"] = E_MP2J_;
    Process::environment.globals["E_MP2K"] = E_MP2K_;

    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ====================> MP2 Energies <==================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Reference Energy",         energies_["Reference Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "MP2J Energy",              energies_["MP2J Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "MP2K Energy",              energies_["MP2K Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Same-Spin Energy",         energies_["Same-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Opposite-Spin Energy",     energies_["Opposite-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Correlation Energy",       energies_["Correlation Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Total Energy",             energies_["Total Energy"]);
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ==================> SCS-MP2 Energies <================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [-]\n", "SCS Same-Spin Scale",      scale_ss_);
    fprintf(outfile, "\t %-25s = %24.16f [-]\n", "SCS Opposite-Spin Scale",  scale_os_);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Same-Spin Energy",     energies_["SCS Same-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Opposite-Spin Energy", energies_["SCS Opposite-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Correlation Energy",   energies_["SCS Correlation Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Total Energy",         energies_["SCS Total Energy"]);
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\n");
    fflush(outfile);

}

}} // End Namespaces
#endif // MADNESS

