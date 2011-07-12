/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <liblmp2_solver/lmp2.h>
#include <libqt/qt.h>


namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

class BasisSet;
class Options;
class PSIO;
class Chkpt;

namespace lmp2 {


double LMP2::compute_energy()
{    


//    sleep(90);
    params();
    moinfo();
    setup_factories();
    reference();
    overlap();

    localize_pipek_mezey();
    transform_fock();
    build_domains();
    set_ij_owner_local();
    build_pairdomains();
    init_global_mat();
    projection();
    integral_direct_transformation();
    allocate_T2_memory();

    setup_ij_kj_ik_maps();


    if (me_ == 0)
        fprintf(outfile, "\n  ====> Begin LMP2 Iterations <====\n\n");
    int iter = 0;
    int conv = 0;
    while (conv != 1) {

        if (iter > 0) Elmp2_old_ = Elmp2_;
        Elmp2_ = 0.0;
        Drms_T2_ = 0.0;

        setup_diis(iter);

        if (comm_ == "MADNESS")
            conv = compute_T2_energy(iter);

        print_results(iter);
        fflush(outfile);
        iter++;

//        if (iter > 0)
//            conv = 1;
    }

    if (me_ == 0)
        fprintf(outfile, "\n  =================================\n\n");

    print_summary();

}

void LMP2::setup_ij_kj_ik_maps() {
    // Set up the ij to kj/ik maps
    for (int ij = 0; ij < ij_pairs_; ij++) {

        ij_kj_map_.push_back(std::vector<int>(ndocc_, 0));
        ij_ik_map_.push_back(std::vector<int>(ndocc_, 0));

        int i = ij_i_map_[ij];
        int j = ij_j_map_[ij];

        for (int k = 0; k < ndocc_; k++) {
            if (k > j)
                ij_kj_map_[ij][k] = (k * (k + 1)) / 2 + j;
            else
                ij_kj_map_[ij][k] = (j * (j + 1)) / 2 + k;

            if (i > k)
                ij_ik_map_[ij][k] = (i * (i + 1)) / 2 + k;
            else
                ij_ik_map_[ij][k] = (k * (k + 1)) / 2 + i;
        }
    }
}
#if HAVE_MADNESS == 1
LMP2::LMP2(Options& options, boost::shared_ptr<PSIO> psio)
    : Wavefunction(options, psio), madness::WorldObject<LMP2>(*Communicator::world->get_madworld())
#else
LMP2::LMP2(Options& options, boost::shared_ptr<PSIO> psio)
    : Wavefunction(options, psio)
#endif
{

//    options_ = options;
//    psio_ = psio;

    common_init();

#if HAVE_MADNESS == 1
    process_pending();
#endif
}

void LMP2::common_init() {

    me_ = Communicator::world->me();
    nproc_ = Communicator::world->nproc();
    nthread_ = Communicator::world->nthread();
    comm_ = Communicator::world->communicator();
//    molecule_ = Process::environment.molecule();
    wfn_ = Process::environment.reference_wavefunction();
//    basis_ = wfn_->basisset();

    // after update remove integral_: we get if from
    // inheriting from wavefunction
    integral_ = boost::shared_ptr<IntegralFactory>(
                new IntegralFactory(basisset_, basisset_,
                                    basisset_, basisset_));

    for (int i=0; i < nthread_; i++) {
        eri_.push_back( boost::shared_ptr<TwoElectronInt>(
                           static_cast<TwoElectronInt*>(integral_->eri())) );
    }

#if HAVE_MADNESS == 1
    if (comm_ == "MADNESS") {
        madworld_ = Communicator::world->get_madworld();
        print_mutex = Communicator::world->get_mutex();
        mutex_ = Communicator::world->get_mutex();
        F_mutex_ = Communicator::world->get_mutex();

    }
#endif

    print_header();

}

void LMP2::print_header() const
{
    if (me_ == 0) {
        fprintf(outfile, "\n");
        fprintf(outfile, "         ---------------------------------------------------------\n");
        fprintf(outfile, "                                 LMP2\n");
        fprintf(outfile, "                           written by Ben Mintz\n");
        fprintf(outfile, "         ---------------------------------------------------------\n");
        fprintf(outfile, "\n");
    }

    Communicator::world->print();
}



LMP2::~LMP2()
{

    std::cout << "Finished LMP2" << std::endl;
    ao_start_.clear();
    ao_start_.clear();
    eri_.clear();

    // clear domain vectors
    for (int i=0; i < domain_.size(); i++) {
        domain_[i].clear();
    }
    for (int i=0; i < pair_domain_.size(); i++) {
        pair_domain_[i].clear();
    }
    pairdom_exist_.clear();
    domain_.clear();
    domain_len_.clear();
    pair_domain_.clear();
    pair_domain_len_.clear();
    pair_domain_nr_len_.clear();

    // clear ij map vectors
    ij_owner_.clear();
    ij_local_.clear();
    ij_i_map_.clear();
    ij_j_map_.clear();
    ij_map_neglect_.clear();

    // clear global distributed matrices
    evals_.clear();
    S_virt_.clear();
    F_virt_.clear();
    W_.clear();

    // Clear MN_owner and local
    MN_Owner_.clear();
    MN_local_.clear();

    // Clear eri_2_MN and eri_ij
    eri_2_MN_.clear();
    eri_ij_.clear();

    // Clear the T2 amplitudes
    if (diis_ == 1) {
//        for (int i=0; i < max_diis_vectors_; i++) {
//            T2_ext_[i].clear();
//            error_[i].clear();
//        }
//        T2_old_.clear();
//        T2_ext_.clear();
//        error_.clear();

    }
    else {
        T2_amp_.clear();
    }
    std::cout << "finished lmp2 destructor" << std::endl;

}

void LMP2::moinfo() {
    // The integrals code below does not use symmetry, so we need to accumulate the
    // orbital info for each irrep
    nirreps_ = wfn_->nirrep();
    int *clsdpi = wfn_->doccpi();
    int *orbspi = wfn_->nmopi();
    int *frzcpi = wfn_->frzcpi();
    int *frzvpi = wfn_->frzvpi();
    ndocc_ = 0;
    nvirt_ = 0;
    nfocc_ = 0;
    nfvir_ = 0;
    nso_ = wfn_->nso();
    nact_docc_ = 0;
    nact_virt_ = 0;
    for(int h=0; h < nirreps_; ++h){
        nfocc_     += frzcpi[h];
        nfvir_     += frzvpi[h];
        ndocc_     += clsdpi[h];
        nact_docc_ += clsdpi[h] - frzcpi[h];
        nvirt_     += orbspi[h] - clsdpi[h];
        nact_virt_ += orbspi[h] - frzvpi[h] - clsdpi[h];
    }
    natom_ = molecule_->natom();
    nshell_ = basisset_->nshell();

    print_moinfo();
}

void LMP2::print_moinfo() const {

    if (nirreps_ != 1) {
        std::string symm = molecule_->sym_label();
        std::cout << "symm_from_input = " << symm << std::endl;
        throw InputException("Local MP2 is only valid in C1 symmetry", symm.c_str(), __FILE__, __LINE__);
    }

    if (me_ == 0) {
        fprintf(outfile, "\n  ====> Orbital Information <====\n\n");
        fprintf(outfile, "  Irreps\t\t= %d\n", nirreps_);
        fprintf(outfile, "  AO's\t\t\t= %d\n", nso_);
        fprintf(outfile, "  Doubly Occupied\t= %d\n", ndocc_);
        fprintf(outfile, "  Active Occupied\t= %d\n", nact_docc_);
        fprintf(outfile, "  Frozen Occupied\t= %d\n", nfocc_);
        fprintf(outfile, "  Total Virtuals\t= %d\n", nvirt_);
        fprintf(outfile, "  Active Virtuals\t= %d\n", nact_virt_);
        fprintf(outfile, "  Frozen Virtuals\t= %d\n", nfvir_);
        fprintf(outfile, "\n  ===============================\n\n");
    }

}

void LMP2::setup_factories()
{
    nso_nso_ = boost::shared_ptr<MatrixFactory> (new MatrixFactory());
    occ_occ_ = boost::shared_ptr<MatrixFactory> (new MatrixFactory());

    nso_nso_->init_with(nirreps_, &nso_, &nso_);
    occ_occ_->init_with(nirreps_, &ndocc_, &ndocc_);
}

void LMP2::reference() {
    C_ = wfn_->Ca();
    D_AO_ = wfn_->Da();
    F_AO_ = wfn_->Fa();

    escf_ = wfn_->reference_energy();
    enuc_ = molecule_->nuclear_repulsion_energy();

    C_->set_name("MO Coefficients (AO)");
    D_AO_->set_name("Density Matrix (AO)");
    F_AO_->set_name("Fock Matrix (AO)");

    if (me_ == 0) {
        fprintf(outfile, "\n  ====> Reference WFN Information <====\n\n");
        fprintf(outfile,"  Nuclear repusion\t= %5.15f\n",enuc_);
        fprintf(outfile,"  SCF energy\t\t= %5.15f\n",escf_);

        if (print_ >= 2) {
            C_->print();
            D_AO_->print();
            F_AO_->print();
        }

        fprintf(outfile, "\n  =====================================\n\n");

        fflush(outfile);
    }
}

void LMP2::overlap()
{
    boost::shared_ptr<OneBodyAOInt> S(integral_->ao_overlap());
    S_ = SharedMatrix(nso_nso_->create_matrix("Overlap Matrix (AO)"));

    S->compute(S_);

    if (print_ >= 2)
        S_->print();

}

void LMP2::transform_fock()
{
    // Transform the Fock matrix from the AO to LO basis
    F_LO_ = SharedMatrix(nso_nso_->create_matrix("Fock Matrix (LO)"));
    F_LO_->transform(F_AO_, C_);

    if (print_ >= 2) {
        F_AO_->print();
        F_LO_->print();
    }

}

void LMP2::set_ij_maps()
{

    int counter;

    counter = 0;
    for (int i = 0, ij=0 ; i < ndocc_; i++) {
        for (int j = 0; j <= i; j++, ij++) {
            ij_map_neglect_.insert(std::pair<int,int>(counter, ij));
            if (pairdom_exist_[ij]) {
                ij_i_map_.insert(std::pair<int,int>(counter, i));
                ij_j_map_.insert(std::pair<int,int>(counter, j));
                counter++;
            }
        }
    }    
}

void LMP2::set_ij_owner_local()
{
    ij_pairs_per_proc_ = 0;
    ij_owner_.clear();
    ij_local_.clear();
    int v = 0, count = 0;
    for(int ij=0; ij < ij_pairs_; ij++) {
        ij_local_.push_back(count);
        ij_owner_.push_back(v % nproc_);
        if (v%nproc_ == nproc_-1) count++;
        if (me_ == ij_owner_[ij])
            ij_pairs_per_proc_++;
        v++;
    }
}

void LMP2::build_pairdomains()
{
    pair_domain_.clear();
    pair_domain_len_.clear();
    for (int i=0; i < ij_pairs_; i++) {
        pair_domain_.push_back(std::vector<int>(natom_,0));
        pair_domain_len_.push_back(0);
    }

    for(int m=0; m < ij_pairs_; m++) {
        int i = ij_i_map_[m];
        int j = ij_j_map_[m];
        int ij = (i * (i + 1)) / 2 + j;
        if (pairdom_exist_[ij]) {
            for (int A=0; A < natom_; A++) {
                if (domain_[i][A] || domain_[j][A]) {
                    pair_domain_[m][A] = 1;
                    pair_domain_len_[m] += ao_stop_[A] - ao_start_[A];
                }
            }
        }
    }
}

void LMP2::allocate_T2_memory() {

    if (diis_ == 1) {
        // Allocate memory for the amplitudes
        for (int i=0; i < max_diis_vectors_; i++) {
            T2_amp_.push_back( std::vector< SharedMatrix > () );
            T2_ext_.push_back( std::vector< SharedMatrix > () );
            error_.push_back( std::vector< SharedMatrix > () );

            for (int ij=0; ij < ij_pairs_; ij++) {
                if (me_ == ij_owner_[ij]) {
                    T2_amp_[i].push_back( SharedMatrix(new Matrix("T2 Amplitudes", nirreps_,
                                                              &pair_domain_len_[ij],
                                                              &pair_domain_len_[ij])) );
                    T2_ext_[i].push_back( SharedMatrix(new Matrix("T2 Extrapolated Amplitudes", nirreps_,
                                                                  &pair_domain_len_[ij],
                                                                  &pair_domain_len_[ij])) );
                    error_[i].push_back( SharedMatrix(new Matrix("DIIS Error Matrix", nirreps_,
                                                                 &pair_domain_len_[ij],
                                                                 &pair_domain_len_[ij])) );
                }
            }
        }
    }
    else {
        for (int i=0; i < 2; i++) {
            T2_amp_.push_back(std::vector<SharedMatrix> (ij_pairs_per_proc_));
            for (int ij=0; ij < ij_pairs_; ij++) {
                if (me_ == ij_owner_[ij]) {
                    T2_amp_[i][ij_local_[ij]] = SharedMatrix(new Matrix("T2 Old Amplitudes", nirreps_,
                                                                        &pair_domain_len_[ij],
                                                                        &pair_domain_len_[ij]));
                }
            }
        }
    }
}

void LMP2::setup_diis(const int &iter) {

    if(diis_ == 1) {
      //  **** Set up variables for glob.diis_ extrapolation ****
      if(iter < max_diis_vectors_-1) div_ = iter;
      else div_ = iter%max_diis_vectors_;

      if(iter <= max_diis_vectors_) matsize_ = iter - 1;
      else matsize_ = max_diis_vectors_;

      if(div_ == 0) dmat1_ = max_diis_vectors_ - 1;
      else dmat1_ = div_ - 1;

      if(div_ == 0) dmat2_ = max_diis_vectors_ - 2;
      else if(div_ == 1) dmat2_ = max_diis_vectors_ - 1;
      else dmat2_ = div_ - 2;

      nmat_ = iter%2;
      if(nmat_ == 1) omat_ = 0;
      else omat_ = 1;
    }
    else {
      div_ = iter%2;
      if(div_ == 1) dmat1_ = 0;
      else dmat1_ = 1;
    }

}


void LMP2::print_summary() const {
    if (me_ == 0) {
        fprintf(outfile, "\n  ====> LMP2 Summary <====\n\n");
        fprintf(outfile, "  SCF energy\t\t = %20.12f\n",escf_);
        fprintf(outfile, "  LMP2 Corr. Energy\t = %20.12f\n", Elmp2_);
        fprintf(outfile, "  LMP2 Total Energy\t = %20.12f\n", escf_ + Elmp2_);
        fprintf(outfile, "\n  ========================\n\n");
        fflush(outfile);
    }
}

int LMP2::print_matrix(const SharedMatrix mat, const std::string &str) const {
    if (mat->name().c_str() == NULL) mat->print(outfile, str.c_str());
    else {
        mat->set_name(str);
        mat->print();
    }
    return 0;
}



}}

