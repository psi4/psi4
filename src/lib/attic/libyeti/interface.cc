#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#if HAVE_MPQC
#include "yeti.h"
#include "schwarz.h"

using namespace sc;
using namespace std;
using namespace yeti;
using namespace mpqc;

static ClassDesc MPQCMole_cd(typeid(MPQCMole), "MPQCMole", 1, "public Wavefunction",
                        0, create<MPQCMole>, create<MPQCMole>);

Ref<GaussianBasisSet> MPQC::unit_basis = 0;

MPQCMole::MPQCMole(const Ref<KeyVal> &keyval):Wavefunction(keyval) {
    ref_ << keyval->describedclassvalue("reference");
    if(ref_.null()) {
        throw InputError("require a OneBodyWavefunction object",
                         __FILE__, __LINE__, "reference", 0,
                         class_desc());
    }
    obs_ = ref_->basis();

    MPQC::unit_basis = new GaussianBasisSet(GaussianBasisSet::Unit);
    nrepeat_ = keyval->intvalue("nrepeat", KeyValValueint(10));
    nthread_ = keyval->intvalue("nthread", KeyValValueint(1));
    sleep_ = keyval->intvalue("sleep", KeyValValueint(0));
    memory_ = keyval->intvalue("memory", KeyValValueint(10000000));
}

MPQCMole::MPQCMole(StateIn &statein):Wavefunction(statein)
{
    ref_ << SavableState::restore_state(statein);
}

void
MPQCMole::save_data_state(StateOut &stateout) {
    Wavefunction::save_data_state(stateout);
    SavableState::save_state(ref_.pointer(),stateout);
}

void
MPQCMole::compute(void)
{
    if(gradient_needed()) {
        throw sc::FeatureNotImplemented("no gradients yet",
                                    __FILE__    , __LINE__, class_desc());
    }

    set_energy(0.0);

    timer::Timer::stop("all");

    timer::Timer::print();
}

void
MPQCMole::obsolete(void) {
    Wavefunction::obsolete();
    ref_->obsolete();
}

int
MPQCMole::nelectron(void) {
    return ref_->nelectron();
}

RefSymmSCMatrix
MPQCMole::density(void) {
    throw sc::FeatureNotImplemented("no density yet",
                                __FILE__, __LINE__, class_desc());
    return 0;
}

int
MPQCMole::spin_polarized(void) {
    return 0;
}

int
MPQCMole::value_implemented(void) const {
    return 1;
}

Ref<CLHF>
MPQCMole::ref()
{
    return ref_;
}

//#endif


MPQCShellComputeFunctor::MPQCShellComputeFunctor(
    Ref<TwoBodyInt> tbint
) : tbint_(tbint),
    buffer_(tbint->buffer())
{

}

void
MPQCShellComputeFunctor::operator()(
    uli i,
    uli j,
    uli a,
    uli b
) const
{
    tbint_->compute_shell(i, j, a, b);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//  TwoElectronIntegralComputer class
////////////////////////////////////////////////////////////////////////////////////////////////////

TwoElectronIntegralComputer::TwoElectronIntegralComputer(
    const Ref<TwoBodyIntDescr>& descr,
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2,
    const Ref<GaussianBasisSet>& bs3,
    const Ref<GaussianBasisSet>& bs4
 ) :
    bs1_(bs1),
    bs2_(bs2),
    bs3_(bs3),
    bs4_(bs4),
    shmap1_(0),
    shmap2_(0),
    shmap3_(0),
    shmap4_(0),
    IF_(IF),
    buffer_(0),
    tbint_(0),
    tbint_descr_(descr),
    multishell_(false),
    sort_(0),
    TwoElectronEstimableComputer()
{
    init();
}

TwoElectronIntegralComputer::TwoElectronIntegralComputer(
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2,
    const Ref<GaussianBasisSet>& bs3,
    const Ref<GaussianBasisSet>& bs4
) :
    bs1_(bs1),
    bs2_(bs2),
    bs3_(bs3),
    bs4_(bs4),
    shmap1_(0),
    shmap2_(0),
    shmap3_(0),
    shmap4_(0),
    IF_(IF),
    buffer_(0),
    tbint_(0),
    tbint_descr_(0),
    multishell_(false),
    sort_(0),
    TwoElectronEstimableComputer()
{
    init();
}

TensorElementComputer*
TwoElectronIntegralComputer::copy() const
{
    TensorElementComputer* comp = 0;
    if (tbint_descr_)
    {
        comp = new TwoElectronIntegralComputer(
            tbint_descr_, IF_, bs1_, bs2_, bs3_, bs4_
        );
    }
    else
    {
        comp = new TwoElectronIntegralComputer(
            IF_, bs1_, bs2_, bs3_, bs4_
        );
    }

    comp->set_index_descr(descr_);
    return comp;
}


TwoElectronIntegralComputer::~TwoElectronIntegralComputer()
{
    for(int i = 0; i < estimators_.size(); ++i) {
        delete estimators_[i];
    }

    if (sort_)
        delete sort_;
}

void
TwoElectronIntegralComputer::sort(Permutation* p)
{
    if (sort_)
        yeti_throw(SanityCheckError, "cannot resort ao integral recomputed tensor");

    sort_ = new Sort(p);
}

void
TwoElectronIntegralComputer::compute(
    const uli* indices,
    double* data,
    uli nelements
)
{
    double* tmp_data = data;
    if (sort_)
    {
        tmp_data = reinterpret_cast<double*>(sort_->data_buffer());
        sort_->get_permutation()->permute(indices, permuted_indices_);
        indices = permuted_indices_;
    }

    if (!multishell_)
    {
        tbint_->compute_shell(
            aobs1_->program_shell_number(indices[0]),
            aobs2_->program_shell_number(indices[1]),
            aobs3_->program_shell_number(indices[2]),
            aobs4_->program_shell_number(indices[3])
        );
        ::memcpy(tmp_data, buffer_, nelements * sizeof(double));
    }
    else
    {
        YetiRuntime::start_timer("compute shell");
        uli nsh1, nsh2, nsh3, nsh4;
        uli shstart4, shstart3, shstart2, shstart1;
        uli nfxn1, nfxn2, nfxn3, nfxn4;

        uli stride = 1;
        uli stride4 = stride;
        uli idx4 = indices[3];
        nsh4 = shmap4_->nshell(idx4);
        shstart4 = shmap4_->shell_start(idx4);
        nfxn4 = shmap4_->nfxn(idx4);
        stride *= nfxn4;

        uli stride3 = stride;
        uli idx3 = indices[2];
        nsh3 = shmap3_->nshell(idx3);
        shstart3 = shmap3_->shell_start(idx3);
        nfxn3 = shmap3_->nfxn(idx3);
        stride *= nfxn3;

        uli stride2 = stride;
        uli idx2 = indices[1];
        nsh2 = shmap2_->nshell(idx2);
        shstart2 = shmap2_->shell_start(idx2);
        nfxn2 = shmap2_->nfxn(idx2);
        stride *= nfxn2;

        uli stride1 = stride;
        uli idx1 = indices[0];
        nsh1 = shmap1_->nshell(idx1);
        shstart1 = shmap1_->shell_start(idx1);
        nfxn1 = shmap1_->nfxn(idx1);

        uli _sh1, _sh2, _sh3, _sh4;
        uli nfxn_sh1, nfxn_sh2, nfxn_sh3, nfxn_sh4;

        double *sh1_dptr, *sh2_dptr, *sh3_dptr, *sh4_dptr;
        double *idx1_dptr, *idx2_dptr, *idx3_dptr, *idx4_dptr;

        sh1_dptr = tmp_data;
        for (uli sh1 = 0; sh1 < nsh1; ++sh1)
        {
            _sh1 = sh1 + shstart1;
            sh2_dptr = sh1_dptr;
            nfxn_sh1 = shmap1_->shell_size(_sh1);
            for (uli sh2 = 0; sh2 < nsh2; ++sh2)
            {
                _sh2 = sh2 + shstart2;
                sh3_dptr = sh2_dptr;
                nfxn_sh2 = shmap2_->shell_size(_sh2);
                for (uli sh3 = 0; sh3 < nsh3; ++sh3)
                {
                    _sh3 = sh3 + shstart3;
                    sh4_dptr = sh3_dptr;
                    nfxn_sh3 = shmap3_->shell_size(_sh3);
                    for (uli sh4 = 0; sh4 < nsh4; ++sh4)
                    {
                        _sh4 = sh4 + shstart4;
                        nfxn_sh4 = shmap4_->shell_size(_sh4);
                        tbint_->compute_shell(
                            aobs1_->program_shell_number(_sh1),
                            aobs2_->program_shell_number(_sh2),
                            aobs3_->program_shell_number(_sh3),
                            aobs4_->program_shell_number(_sh4)
                        );
                        idx1_dptr = sh4_dptr;
                        const double* intptr = buffer_;

                        size_t offset = (char*) sh4_dptr - (char*) tmp_data;
                        offset /= sizeof(double);

                        for (usi idx1=0; idx1 < nfxn_sh1; ++idx1, idx1_dptr += stride1)
                        {
                            idx2_dptr = idx1_dptr;
                            for (usi idx2=0; idx2 < nfxn_sh2; ++idx2, idx2_dptr += stride2)
                            {
                                idx3_dptr = idx2_dptr;
                                for (usi idx3=0; idx3 < nfxn_sh3; ++idx3, idx3_dptr += stride3)
                                {
                                    idx4_dptr = idx3_dptr;
                                    for (usi idx4=0; idx4 < nfxn_sh4; ++idx4, ++idx4_dptr, ++intptr)
                                    {
                                        *idx4_dptr = *intptr;
                                    } //end idx4 loop
                                } //end idx3 loop
                            } //end idx2 loop
                        } //end idx1 loop
                        sh4_dptr += stride4 * nfxn_sh4;
                    } //end sh4 loop
                    sh3_dptr += stride3 * nfxn_sh3;
                } //end sh3 loop
                sh2_dptr += stride2 * nfxn_sh2;
            } //end sh2 loop
            sh1_dptr += stride1 * nfxn_sh1;
        } //end sh1 loop
        YetiRuntime::stop_timer("compute shell");
    }

    if (sort_) //go to physicist notation
    {
        descr_->get_nelements(indices, nelements_);
        sort_->configure(nelements_);
        sort_->sort<double>(tmp_data, data);
    }
}

void
TwoElectronIntegralComputer::init()
{
    if (!aobs1_) aobs1_ = YetiRuntime::get_basis(bs1_->name());
    if (!aobs2_) aobs2_ = YetiRuntime::get_basis(bs2_->name());
    if (!aobs3_) aobs3_ = YetiRuntime::get_basis(bs3_->name());
    if (!aobs4_) aobs4_ = YetiRuntime::get_basis(bs4_->name());
    shmap1_ = aobs1_->get_multishell_map();
    shmap2_ = aobs2_->get_multishell_map();
    shmap3_ = aobs3_->get_multishell_map();
    shmap4_ = aobs4_->get_multishell_map();

    IF_->set_basis(bs1_, bs2_, bs3_, bs4_);

    multishell_ = shmap1_;

    if (tbint_descr_.nonnull())
    {
        tbint_ = tbint_descr_->inteval();
        unsigned int t = 0;
        buffer_ = tbint_->buffer(tbint_descr_->intset(t));
    }
    else
    {
        tbint_ = IF_->electron_repulsion();
        buffer_ = tbint_->buffer();
    }

}



void
TwoElectronIntegralComputer::init_cauchy_schwarz()
{
    CauchySchwarzValueEstimater* last_est = new CauchySchwarzValueEstimater(new MPQCShellComputeFunctor(tbint_), aobs1_, descr_.get());
    estimators_.push_back(last_est);
    for(int i = 2; i <= descr_->depth(); ++ i) {
        last_est = new CauchySchwarzValueEstimater(last_est, descr_.get(), i);
        estimators_.push_back(last_est);
    }
}




yeti::TemplateInfo::type_t
mpqc::TwoElectronIntegralComputer::element_type(
    const uli* indices,
    usi depth
)
{
    return yeti::TemplateInfo::double_type;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
//  DijabElementComputer class
////////////////////////////////////////////////////////////////////////////////////////////////////

DijabElementComputer::DijabElementComputer(
    const RefDiagSCMatrix& ei,
    const RefDiagSCMatrix& ej,
    const RefDiagSCMatrix& ea,
    const RefDiagSCMatrix& eb
)
    :
    ei_(ei),
    ej_(ej),
    ea_(ea),
    eb_(eb)
{
}

void
DijabElementComputer::compute(
    const uli* indices,
    double* data,
    uli nelements
)
{
    ni_ = descr_->get(0)->nelements(indices[0]);
    istart_ = descr_->get(0)->index_start(indices[0]);
    nj_ = descr_->get(1)->nelements(indices[1]);
    jstart_ = descr_->get(1)->index_start(indices[1]);
    na_ = descr_->get(2)->nelements(indices[2]);
    astart_ = descr_->get(2)->index_start(indices[2]);
    nb_ = descr_->get(3)->nelements(indices[3]);
    bstart_ = descr_->get(3)->index_start(indices[3]);

    uli istop = istart_ + ni_;
    uli jstop = jstart_ + nj_;
    uli astop = astart_ + na_;
    uli bstop = bstart_ + nb_;

    double* dataptr = data;
    for (uli i=istart_; i < istop; ++i)
    {
        double ei = ei_.get_element(i);
        for (uli j=jstart_; j <  jstop; ++j)
        {
            double ej = ej_.get_element(j);
            for (uli a=astart_; a < astop; ++a)
            {
                double ea = ea_.get_element(a);
                for (uli b=bstart_; b < bstop; ++b, ++dataptr)
                {
                    double eb = eb_.get_element(b);
                    (*dataptr) = ei + ej - ea - eb;
                }
            }
        }
    }
}

TensorElementComputer*
DijabElementComputer::copy() const
{
    TensorElementComputer* comp = new DijabElementComputer(ei_, ej_, ea_, eb_);
    comp->set_index_descr(descr_);
    return comp;
}

yeti::TemplateInfo::type_t
mpqc::DijabElementComputer::element_type(
    const uli* indices,
    usi depth
)
{
    return yeti::TemplateInfo::double_type;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


double
mpqc::min_exponent(
    const GaussianShell& shell
)
{
    uli ncxn = shell.ncontraction();
    uli nprim = shell.nprimitive();
    double minexp = 999999999.0;
    for (uli cxn=0; cxn < ncxn; ++cxn)
    {
        double invexp = 0;
        for (uli prim=0; prim < nprim; ++prim)
        {
            double exp = shell.exponent(prim);
            double coef = shell.coefficient_norm(cxn,prim);
            invexp += coef * 1.0 / (exp * exp);
        }
        double exp = 1.0/sqrt(invexp);
        if (minexp > exp)
            minexp = exp;
    }
    return minexp;
}

void
mpqc::build_ao_range(
    const Ref<GaussianBasisSet>& obs,
    const std::string& descr,
    uli nfxn_min_per_tile,
    const char* id1,
    const char* id2,
    const char* id3,
    const char* id4,
    PartitioningPolicyPtr policy
)
{
    Ref<Molecule> mol = obs->molecule();
    //set up the yeti index range
    AOBasisPtr aobasis = new AOBasis(obs->name(), descr);
    uli nbasis = obs->nbasis();
    uli natoms = obs->ncenter();

    for (int atomnum=0; atomnum < natoms; ++atomnum)
    {
        AtomPtr atom(
            new Atom(
                mol->atom_symbol(atomnum),
                mol->atom_name(atomnum),
                mol->r(atomnum, 0),
                mol->r(atomnum, 1),
                mol->r(atomnum, 2),
                atomnum
            )
        );
        aobasis->add_atom(atom);

        #define AM_MAX 10
        uli counts[AM_MAX];
        for (uli i=0; i < AM_MAX; ++i)
            counts[i] = 0;

        int start = obs->shell_on_center(atomnum,0);
        int stop = start + obs->nshell_on_center(atomnum);
        for (int shellnum=start; shellnum < stop; ++shellnum)
        {
            const GaussianShell& gshell = obs->shell(shellnum);
            usi am = gshell.min_angular_momentum();
            double minexp = min_exponent(gshell);
            uli ncxn = gshell.ncontraction();
            uli nfxn = gshell.nfunction();
            uli fxn_start = obs->shell_to_function(shellnum);
            ShellPtr shell = new Shell(atom, am, ncxn, nfxn, minexp, shellnum, fxn_start, counts[am]);
            //compute the ''average'' exponent for each contraction
            aobasis->add_shell(shell);

            ++counts[am];
        }
    }

    aobasis->configure_min_nfxn_per_tile(nfxn_min_per_tile);

    if(policy == 0) {
        uli nfxn_per_atom = nbasis / natoms;
        if      (natoms > 10 && nfxn_per_atom > 20)
            aobasis->configure_many_atoms_large_basis();
        else if (natoms <= 10 && nfxn_per_atom > 20)
            aobasis->configure_few_atoms_large_basis();
        else if (natoms > 10 && nfxn_per_atom < 20)
            aobasis->configure_many_atoms_small_basis();
        else if (natoms <= 10 && nfxn_per_atom < 20)
            aobasis->configure_few_atoms_small_basis();
    }
    else
    {
        policy->set_min_fxn_per_tile(nfxn_min_per_tile);
        aobasis->configure(policy);
    }

    YetiRuntime::register_basis_set(aobasis);

    YetiRuntime::register_index_descr(
        aobasis->get_index_descr(),
        id1, id2, id3, id4
    );
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//  OneElectronIntegralComputer class
////////////////////////////////////////////////////////////////////////////////////////////////////

OneElectronIntegralComputer::OneElectronIntegralComputer(
    MPQC::obint_type_t obint_type,
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2
) :
    bs1_(bs1),
    bs2_(bs2),
    shmap1_(0),
    shmap2_(0),
    IF_(IF),
    buffer_(0),
    obint_(0),
    obint_type_(obint_type),
    multishell_(false)
{
    if (!bs2_)
        bs2_ = bs1_;

    aobs1_ = YetiRuntime::get_basis(bs1->name());
    aobs2_ = YetiRuntime::get_basis(bs2->name());
    shmap1_ = aobs1_->get_multishell_map();
    shmap2_ = aobs2_->get_multishell_map();

    //shmap1_ = YetiRuntime::get_basis(bs1->name())->get_multishell_map();
    //shmap2_ = YetiRuntime::get_basis(bs2->name())->get_multishell_map();

    multishell_ = shmap1_;

    IF->set_basis(bs1_, bs2_);

    switch(obint_type_)
    {
        case MPQC::hcore:
            obint_ = IF->hcore();
            break;
        case MPQC::overlap:
            obint_ = IF->overlap();
            break;
    }
    buffer_ = obint_->buffer();
}

TensorElementComputer*
OneElectronIntegralComputer::copy() const
{
    return new OneElectronIntegralComputer(
            obint_type_, IF_, bs1_, bs2_
    );
}

OneElectronIntegralComputer::~OneElectronIntegralComputer()
{
}

void
OneElectronIntegralComputer::compute(
    const uli* indices,
    double* data,
    uli nelements
)
{
    if (!multishell_)
    {
        obint_->compute_shell(
            aobs1_->program_shell_number(indices[0]),
            aobs2_->program_shell_number(indices[1])
        );
        ::memcpy(data, buffer_, nelements * sizeof(double));
    }
    else
    {
        uli nsh1, nsh2;
        uli shstart2, shstart1;
        uli nfxn1, nfxn2;

        uli stride = 1;

        uli stride2 = stride;
        uli idx2 = indices[1];
        nsh2 = shmap2_->nshell(idx2);
        shstart2 = shmap2_->shell_start(idx2);
        nfxn2 = shmap2_->nfxn(idx2);
        stride *= nfxn2;

        uli stride1 = stride;
        uli idx1 = indices[0];
        nsh1 = shmap1_->nshell(idx1);
        shstart1 = shmap1_->shell_start(idx1);
        nfxn1 = shmap1_->nfxn(idx1);

        uli _sh1, _sh2;
        uli nfxn_sh1, nfxn_sh2;

        double *sh1_dptr, *sh2_dptr;
        double *idx1_dptr, *idx2_dptr;

        sh1_dptr = data;
        for (uli sh1 = 0; sh1 < nsh1; ++sh1)
        {
            _sh1 = sh1 + shstart1;
            sh2_dptr = sh1_dptr;
            nfxn_sh1 = shmap1_->shell_size(_sh1);
            for (uli sh2 = 0; sh2 < nsh2; ++sh2)
            {
                _sh2 = sh2 + shstart2;
                nfxn_sh2 = shmap2_->shell_size(_sh2);
                obint_->compute_shell(
                    aobs1_->program_shell_number(_sh1),
                    aobs2_->program_shell_number(_sh2)
                );
                idx1_dptr = sh2_dptr;
                const double* intptr = buffer_;
                for (usi idx1=0; idx1 < nfxn_sh1; ++idx1, idx1_dptr += stride1)
                {
                    idx2_dptr = idx1_dptr;
                    for (usi idx2=0; idx2 < nfxn_sh2; ++idx2, ++idx2_dptr, ++intptr)
                    {
                        *idx2_dptr = *intptr;
                    } //end idx2 loop
                } //end idx1 loop
                sh2_dptr += stride2 * nfxn_sh2;
            } //end sh2 loop
            sh1_dptr += stride1 * nfxn_sh1;
        } //end sh1 loop
    }
}

yeti::TemplateInfo::type_t
mpqc::OneElectronIntegralComputer::element_type(
    const uli* indices,
    usi depth
)
{
    return yeti::TemplateInfo::double_type;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//  TwoElectronIntegralComputer2Index class
////////////////////////////////////////////////////////////////////////////////////////////////////


TwoElectronIntegralComputer2Index::TwoElectronIntegralComputer2Index(
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2
) :
    bs1_(bs1),
    bs2_(bs2),
    shmap1_(0),
    shmap2_(0),
    IF_(IF),
    buffer_(0),
    tbint_(0),
    tbint_descr_(0)
{
    init();
}

TwoElectronIntegralComputer2Index::TwoElectronIntegralComputer2Index(
    const Ref<TwoBodyIntDescr>& descr,
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2
) :
    bs1_(bs1),
    bs2_(bs2),
    shmap1_(0),
    shmap2_(0),
    IF_(IF),
    buffer_(0),
    tbint_(0),
    tbint_descr_(descr),
    multishell_(false)
{
    init();
}

TensorElementComputer*
TwoElectronIntegralComputer2Index::copy() const
{
    if (tbint_descr_)
    {
        return new TwoElectronIntegralComputer2Index(
            tbint_descr_, IF_, bs1_, bs2_
        );
    }
    else
    {
        return new TwoElectronIntegralComputer2Index(
            IF_, bs1_, bs2_
        );
    }

}

TwoElectronIntegralComputer2Index::~TwoElectronIntegralComputer2Index()
{
}

void
TwoElectronIntegralComputer2Index::compute(
    const uli* indices,
    double* data,
    uli nelements
)
{
    if (!multishell_)
    {
        tbint_->compute_shell(
            aobs1_->program_shell_number(indices[0]),
            0,
            aobs2_->program_shell_number(indices[1]),
            0
        );
        ::memcpy(data, buffer_, nelements * sizeof(double));
    }
    else
    {
        uli nsh1, nsh2;
        uli shstart2, shstart1;
        uli nfxn1, nfxn2;

        uli stride = 1;

        uli stride2 = stride;
        uli idx2 = indices[1];
        nsh2 = shmap2_->nshell(idx2);
        shstart2 = shmap2_->shell_start(idx2);
        nfxn2 = shmap2_->nfxn(idx2);
        stride *= nfxn2;

        uli stride1 = stride;
        uli idx1 = indices[0];
        nsh1 = shmap1_->nshell(idx1);
        shstart1 = shmap1_->shell_start(idx1);
        nfxn1 = shmap1_->nfxn(idx1);

        uli _sh1, _sh2;
        uli nfxn_sh1, nfxn_sh2;

        double *sh1_dptr, *sh2_dptr;
        double *idx1_dptr, *idx2_dptr;

        sh1_dptr = data;
        for (uli sh1 = 0; sh1 < nsh1; ++sh1)
        {
            _sh1 = sh1 + shstart1;
            sh2_dptr = sh1_dptr;
            nfxn_sh1 = shmap1_->shell_size(_sh1);
            for (uli sh2 = 0; sh2 < nsh2; ++sh2)
            {
                _sh2 = sh2 + shstart2;
                nfxn_sh2 = shmap2_->shell_size(_sh2);
                tbint_->compute_shell(
                    aobs1_->program_shell_number(_sh1),
                    0,
                    aobs2_->program_shell_number(_sh2),
                    0
                );
                idx1_dptr = sh2_dptr;
                const double* intptr = buffer_;
                for (usi idx1=0; idx1 < nfxn_sh1; ++idx1, idx1_dptr += stride1)
                {
                    idx2_dptr = idx1_dptr;
                    for (usi idx2=0; idx2 < nfxn_sh2; ++idx2, ++idx2_dptr, ++intptr)
                    {
                        *idx2_dptr = *intptr;
                    } //end idx2 loop
                } //end idx1 loop
                sh2_dptr += stride2 * nfxn_sh2;
            } //end sh2 loop
            sh1_dptr += stride1 * nfxn_sh1;
        } //end sh1 loop
    }
}

yeti::TemplateInfo::type_t
mpqc::TwoElectronIntegralComputer2Index::element_type(
    const uli* indices,
    usi depth
)
{
    return yeti::TemplateInfo::double_type;
}

void
TwoElectronIntegralComputer2Index::init()
{
    if (!aobs1_) aobs1_ = YetiRuntime::get_basis(bs1_->name());
    if (!aobs2_) aobs2_ = YetiRuntime::get_basis(bs2_->name());
    shmap1_ = aobs1_->get_multishell_map();
    shmap2_ = aobs2_->get_multishell_map();
    multishell_ = shmap1_;

    IF_->set_basis(bs1_, MPQC::unit_basis, bs2_, MPQC::unit_basis);

    if (tbint_descr_.nonnull())
    {
        tbint_ = tbint_descr_->inteval();
        unsigned int t = 0;
        buffer_ = tbint_->buffer(tbint_descr_->intset(t));
    }
    else
    {
        tbint_ = IF_->electron_repulsion();
        buffer_ = tbint_->buffer();
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//  TwoElectronIntegralComputer3Index class
////////////////////////////////////////////////////////////////////////////////////////////////////


TwoElectronIntegralComputer3Index::TwoElectronIntegralComputer3Index(
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2,
    const Ref<GaussianBasisSet>& bs3
 ) :
    bs1_(bs1),
    bs2_(bs2),
    bs3_(bs3),
    shmap1_(0),
    shmap2_(0),
    shmap3_(0),
    IF_(IF),
    buffer_(0),
    tbint_(0),
    tbint_descr_(0),
    multishell_(0)
{
    init();
}

TwoElectronIntegralComputer3Index::TwoElectronIntegralComputer3Index(
    const Ref<TwoBodyIntDescr>& descr,
    const Ref<Integral>& IF,
    const Ref<GaussianBasisSet>& bs1,
    const Ref<GaussianBasisSet>& bs2,
    const Ref<GaussianBasisSet>& bs3
 ) :
    bs1_(bs1),
    bs2_(bs2),
    bs3_(bs3),
    shmap1_(0),
    shmap2_(0),
    shmap3_(0),
    IF_(IF),
    buffer_(0),
    tbint_(0),
    tbint_descr_(descr),
    multishell_(0)
{
    init();
}

TensorElementComputer*
TwoElectronIntegralComputer3Index::copy() const
{
    if (tbint_descr_)
    {
        return new TwoElectronIntegralComputer3Index(
            tbint_descr_, IF_, bs1_, bs2_, bs3_
        );
    }
    else
    {
        return new TwoElectronIntegralComputer3Index(
            IF_, bs1_, bs2_, bs3_
        );
    }

}

TwoElectronIntegralComputer3Index::~TwoElectronIntegralComputer3Index()
{
}

void
TwoElectronIntegralComputer3Index::compute(
    const uli* indices,
    double* data,
    uli nelements
)
{
    if (!multishell_)
    {
        tbint_->compute_shell(
            aobs1_->program_shell_number(indices[0]),
            aobs2_->program_shell_number(indices[1]),
            aobs3_->program_shell_number(indices[2]),
            0
        );
        ::memcpy(data, buffer_, nelements * sizeof(double));
    }
    else
    {
        uli nsh1, nsh2, nsh3;
        uli shstart3, shstart2, shstart1;
        uli nfxn1, nfxn2, nfxn3;

        uli stride = 1;

        uli stride3 = stride;
        uli idx3 = indices[2];
        nsh3 = shmap3_->nshell(idx3);
        shstart3 = shmap3_->shell_start(idx3);
        nfxn3 = shmap3_->nfxn(idx3);
        stride *= nfxn3;

        uli stride2 = stride;
        uli idx2 = indices[1];
        nsh2 = shmap2_->nshell(idx2);
        shstart2 = shmap2_->shell_start(idx2);
        nfxn2 = shmap2_->nfxn(idx2);
        stride *= nfxn2;

        uli stride1 = stride;
        uli idx1 = indices[0];
        nsh1 = shmap1_->nshell(idx1);
        shstart1 = shmap1_->shell_start(idx1);
        nfxn1 = shmap1_->nfxn(idx1);

        uli _sh1, _sh2, _sh3;
        uli nfxn_sh1, nfxn_sh2, nfxn_sh3;

        double *sh1_dptr, *sh2_dptr, *sh3_dptr;
        double *idx1_dptr, *idx2_dptr, *idx3_dptr;

        sh1_dptr = data;
        for (uli sh1 = 0; sh1 < nsh1; ++sh1)
        {
            _sh1 = sh1 + shstart1;
            sh2_dptr = sh1_dptr;
            nfxn_sh1 = shmap1_->shell_size(_sh1);
            for (uli sh2 = 0; sh2 < nsh2; ++sh2)
            {
                _sh2 = sh2 + shstart2;
                sh3_dptr = sh2_dptr;
                nfxn_sh2 = shmap2_->shell_size(_sh2);
                for (uli sh3 = 0; sh3 < nsh3; ++sh3)
                {
                    _sh3 = sh3 + shstart3;
                    nfxn_sh3 = shmap3_->shell_size(_sh3);
                    tbint_->compute_shell(
                        aobs1_->program_shell_number(_sh1),
                        aobs2_->program_shell_number(_sh2),
                        aobs3_->program_shell_number(_sh3),
                        0
                    );
                    idx1_dptr = sh3_dptr;
                    const double* intptr = buffer_;
                    for (usi idx1=0; idx1 < nfxn_sh1; ++idx1, idx1_dptr += stride1)
                    {
                        idx2_dptr = idx1_dptr;
                        for (usi idx2=0; idx2 < nfxn_sh2; ++idx2, idx2_dptr += stride2)
                        {
                            idx3_dptr = idx2_dptr;
                            for (usi idx3=0; idx3 < nfxn_sh3; ++idx3, ++idx3_dptr, ++intptr)
                            {
                                *idx3_dptr = *intptr;
                            } //end idx3 loop
                        } //end idx2 loop
                    } //end idx1 loop
                    sh3_dptr += stride3 * nfxn_sh3;
                } //end sh3 loop
                sh2_dptr += stride2 * nfxn_sh2;
            } //end sh2 loop
            sh1_dptr += stride1 * nfxn_sh1;
        } //end sh1 loop
    }
}

yeti::TemplateInfo::type_t
mpqc::TwoElectronIntegralComputer3Index::element_type(
    const uli* indices,
    usi depth
)
{
    return yeti::TemplateInfo::double_type;
}

void
TwoElectronIntegralComputer3Index::init()
{
    if (!aobs1_) aobs1_ = YetiRuntime::get_basis(bs1_->name());
    if (!aobs2_) aobs2_ = YetiRuntime::get_basis(bs2_->name());
    if (!aobs3_) aobs3_ = YetiRuntime::get_basis(bs3_->name());
    shmap1_ = aobs1_->get_multishell_map();
    shmap2_ = aobs2_->get_multishell_map();
    shmap3_ = aobs3_->get_multishell_map();
    multishell_ = shmap1_;

    IF_->set_basis(bs1_, bs2_, bs3_, MPQC::unit_basis);

    if (tbint_descr_.nonnull())
    {
        tbint_ = tbint_descr_->inteval();
        unsigned int t = 0;
        buffer_ = tbint_->buffer(tbint_descr_->intset(t));
    }
    else
    {
        tbint_ = IF_->electron_repulsion();
        buffer_ = tbint_->buffer();
    }
}


#if HAVE_CHEMISTRY_QC_ZAPTR12_ROHFWFN_H


////////////////////////////////////////////////////////////////////////////////////////////////////
//  IntegralTransformElementComputer class
////////////////////////////////////////////////////////////////////////////////////////////////////


IntegralTransformElementComputer::IntegralTransformElementComputer(
    const Ref<IntegralTransform>& tform,
    uli offset_i,
    uli offset_j,
    uli offset_a,
    uli offset_b,
    uli nindex4,
    TwoBodyInt::tbint_type inttype
) :
    tform_(tform),
    nindex4_(nindex4),
    offset_i_(offset_i),
    offset_j_(offset_j),
    offset_a_(offset_a),
    offset_b_(offset_b),
    inttype_(inttype)
{
}

void
IntegralTransformElementComputer::compute(
    const uli* indices,
    double* data,
    uli nelements
)
{
    ni_ = descr_->get(0)->nelements(indices[0]);
    istart_ = descr_->get(0)->index_start(indices[0]) - offset_i_;
    nj_ = descr_->get(1)->nelements(indices[1]);
    jstart_ = descr_->get(1)->index_start(indices[1]) - offset_j_;
    na_ = descr_->get(2)->nelements(indices[2]);
    astart_ = descr_->get(2)->index_start(indices[2]) - offset_a_;
    nb_ = descr_->get(3)->nelements(indices[3]);
    bstart_ = descr_->get(3)->index_start(indices[3]) - offset_b_;

    uli istop = istart_ + ni_;
    uli jstop = jstart_ + nj_;
    uli astop = astart_ + na_;
    uli bstop = bstart_ + nb_;



    double* dataptr = data;
    uli stride = nindex4_ - nb_;
    for (uli i=istart_; i < istop; ++i)
    {
        for (uli j=jstart_; j <  jstop; ++j)
        {
            const double* ints = tform_->retrieve_pair_block(i, j, inttype_);
            const double* intptr = ints;
            intptr += astart_ * nindex4_ + bstart_;
            for (uli a=astart_; a < astop; ++a, intptr += stride)
            {
                for (uli b=bstart_; b < bstop; ++b, ++dataptr, ++intptr)
                {
                    (*dataptr) = (*intptr);
                }
            }
            tform_->release_pair_block(i,j,inttype_);
        }
    }
}


TensorElementComputer*
IntegralTransformElementComputer::copy() const
{
    return new IntegralTransformElementComputer(
        tform_,
        offset_i_,
        offset_j_,
        offset_a_,
        offset_b_,
        nindex4_
     );
}

yeti::TemplateInfo::type_t
mpqc::IntegralTransformElementComputer::element_type(
    const uli* indices,
    usi depth
)
{
    return yeti::TemplateInfo::double_type;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


Ref<MOIndexSpace>
sc::get_ao_space(
    Ref<GaussianBasisSet> obs,
    Ref<Integral> IF,
    const std::string& id
)
{
    RefSCDimension aodim = new SCDimension(obs->nbasis());
    RefSCMatrix aocoefs = make_block_matrix(MatrixKits::defaultkit->matrix(aodim, aodim));
    aocoefs.assign(0.0);
    for (int i=0; i < aocoefs.nrow(); ++i)
        aocoefs.set_element(i,i,1.0);
    Ref<MOIndexSpace> aospace(new MOIndexSpace(id, "ao", aocoefs, obs, IF));
    return aospace;
}

#endif // ZAPT_WAVEFUNCTION

#elif HAVE_PSI

#include "libmints/mints.h"
#include "libmints/coordentry.h"
#include "boost/shared_ptr.hpp"
#include "libdpd/dpd.h"
#include "yeti.h"

using namespace psi;

double
psi::min_exponent(const psi::GaussianShell& shell) {
    uli nprim = shell.nprimitive();
    double minexp = 1.0E10;
    double invexp = 0;
    for (uli prim=0; prim < nprim; ++prim)
    {
        double exp  = shell.exp(prim);
        double coef = shell.coef(prim);
        invexp += coef * 1.0 / (exp * exp);
    }
    double exp = 1.0/sqrt(invexp);
    if (minexp > exp)
        minexp = exp;
    return minexp;
}

void
psi::build_ao_range(const boost::shared_ptr<BasisSet>& obs,
               const char* id1, const char* id2, const char* id3, const char* id4)
{
    boost::shared_ptr<Molecule> mol = Process::environment.molecule();
    AOBasisPtr aobasis              = new AOBasis(obs->name(), obs->name());
    uli nbasis                      = obs->nbf();
    uli natoms                      = mol->natom();
    int shellnum                    = 0;
    #define AM_MAX 20
    uli counts[AM_MAX];
    for (uli i=0; i < AM_MAX; ++i)
        counts[i] = 0;
    for (int atomnum=0; atomnum < natoms; ++atomnum)
    {
        AtomPtr atom(
            new Atom(
                mol->atom_entry(atomnum)->symbol(),
                mol->atom_entry(atomnum)->label(),
                mol->x(atomnum),
                mol->y(atomnum),
                mol->z(atomnum),
                atomnum
            )
        );
        aobasis->add_atom(atom);

        for (int sh = 0; sh < obs->nshell_on_center(atomnum); ++sh) {
            const GaussianShell& gshell = obs->shell(shellnum);
            usi am         = gshell.am();
            double minexp  = min_exponent(gshell);
            uli ncxn       = 1;
            uli nfxn       = gshell.nfunction();
            uli fxn_start  = obs->shell_to_basis_function(shellnum);
            ShellPtr shell = new Shell(atom, am, ncxn, nfxn, minexp, shellnum, fxn_start, counts[am]);
            aobasis->add_shell(shell);
            ++counts[am];
            ++shellnum;
        }
    }

    uli nfxn_min_per_tile = 20;
    if (nbasis > 200)
        nfxn_min_per_tile = 40;
    aobasis->configure_min_nfxn_per_tile(nfxn_min_per_tile);

    uli nfxn_per_atom = nbasis / natoms;
    if      (natoms > 10 && nfxn_per_atom > 20)
        aobasis->configure_many_atoms_large_basis();
    else if (natoms <= 10 && nfxn_per_atom > 20)
        aobasis->configure_few_atoms_large_basis();
    else if (natoms > 10 && nfxn_per_atom < 20)
        aobasis->configure_many_atoms_small_basis();
    else if (natoms <= 10 && nfxn_per_atom < 20)
        aobasis->configure_few_atoms_small_basis();

    YetiRuntime::register_basis_set(aobasis);

    YetiRuntime::register_index_descr(aobasis->get_index_descr(), id1, id2, id3, id4);
}

TensorElementComputer*
TwoElectronIntegralComputer::copy() const
{
    return new TwoElectronIntegralComputer(eri_);
}

TwoElectronIntegralComputer::TwoElectronIntegralComputer(boost::shared_ptr<TwoBodyAOInt> eri):
   eri_(eri)
{
    buffer_ = eri_->buffer();
    aobs1_  = YetiRuntime::get_basis(eri_->basis1()->name());
    aobs2_  = YetiRuntime::get_basis(eri_->basis2()->name());
    aobs3_  = YetiRuntime::get_basis(eri_->basis3()->name());
    aobs4_  = YetiRuntime::get_basis(eri_->basis4()->name());
    shmap1_ = aobs1_->get_multishell_map();
    shmap2_ = aobs2_->get_multishell_map();
    shmap3_ = aobs3_->get_multishell_map();
    shmap4_ = aobs4_->get_multishell_map();

    multishell_ = shmap1_;
}

void
TwoElectronIntegralComputer::compute(const uli* indices, double* data, uli nelements) {
    if (!multishell_){
        uli sh1 = aobs1_->program_shell_number(indices[0]);
        uli sh2 = aobs2_->program_shell_number(indices[1]);
        uli sh3 = aobs3_->program_shell_number(indices[2]);
        uli sh4 = aobs4_->program_shell_number(indices[3]);
        eri_->compute_shell(sh1, sh2, sh3, sh4);
        ::memcpy(data, buffer_, nelements * sizeof(double));
    }else{
        uli nsh1, nsh2, nsh3, nsh4;
        uli shstart4, shstart3, shstart2, shstart1;
        uli nfxn1, nfxn2, nfxn3, nfxn4;

        uli stride = 1;
        uli stride4 = stride;
        uli idx4 = indices[3];
        nsh4 = shmap4_->nshell(idx4);
        shstart4 = shmap4_->shell_start(idx4);
        nfxn4 = shmap4_->nfxn(idx4);
        stride *= nfxn4;

        uli stride3 = stride;
        uli idx3 = indices[2];
        nsh3 = shmap3_->nshell(idx3);
        shstart3 = shmap3_->shell_start(idx3);
        nfxn3 = shmap3_->nfxn(idx3);
        stride *= nfxn3;

        uli stride2 = stride;
        uli idx2 = indices[1];
        nsh2 = shmap2_->nshell(idx2);
        shstart2 = shmap2_->shell_start(idx2);
        nfxn2 = shmap2_->nfxn(idx2);
        stride *= nfxn2;

        uli stride1 = stride;
        uli idx1 = indices[0];
        nsh1 = shmap1_->nshell(idx1);
        shstart1 = shmap1_->shell_start(idx1);
        nfxn1 = shmap1_->nfxn(idx1);

        uli _sh1, _sh2, _sh3, _sh4;
        uli nfxn_sh1, nfxn_sh2, nfxn_sh3, nfxn_sh4;

        double *sh1_dptr, *sh2_dptr, *sh3_dptr, *sh4_dptr;
        double *idx1_dptr, *idx2_dptr, *idx3_dptr, *idx4_dptr;

        sh1_dptr = data;
        for (uli sh1 = 0; sh1 < nsh1; ++sh1)
        {
            _sh1 = sh1 + shstart1;
            sh2_dptr = sh1_dptr;
            nfxn_sh1 = shmap1_->shell_size(_sh1);
            for (uli sh2 = 0; sh2 < nsh2; ++sh2) {
                _sh2 = sh2 + shstart2;
                sh3_dptr = sh2_dptr;
                nfxn_sh2 = shmap2_->shell_size(_sh2);
                for (uli sh3 = 0; sh3 < nsh3; ++sh3) {
                    _sh3 = sh3 + shstart3;
                    sh4_dptr = sh3_dptr;
                    nfxn_sh3 = shmap3_->shell_size(_sh3);
                    for (uli sh4 = 0; sh4 < nsh4; ++sh4) {
                        _sh4 = sh4 + shstart4;
                        nfxn_sh4 = shmap4_->shell_size(_sh4);
                        eri_->compute_shell(
                            aobs1_->program_shell_number(_sh1),
                            aobs2_->program_shell_number(_sh2),
                            aobs3_->program_shell_number(_sh3),
                            aobs4_->program_shell_number(_sh4)
                        );
                        idx1_dptr = sh4_dptr;
                        const double* intptr = buffer_;

                        size_t offset = (char*) sh4_dptr - (char*) data;
                        offset /= sizeof(double);
                        for (usi idx1=0; idx1 < nfxn_sh1; ++idx1, idx1_dptr += stride1) {
                            idx2_dptr = idx1_dptr;
                            for (usi idx2=0; idx2 < nfxn_sh2; ++idx2, idx2_dptr += stride2) {
                                idx3_dptr = idx2_dptr;
                                for (usi idx3=0; idx3 < nfxn_sh3; ++idx3, idx3_dptr += stride3) {
                                    idx4_dptr = idx3_dptr;
                                    for (usi idx4=0; idx4 < nfxn_sh4; ++idx4, ++idx4_dptr, ++intptr) {
                                        *idx4_dptr = *intptr;
                                    } //end idx4 loop
                                } //end idx3 loop
                            } //end idx2 loop
                        } //end idx1 loop
                        sh4_dptr += stride4 * nfxn_sh4;
                    } //end sh4 loop
                    sh3_dptr += stride3 * nfxn_sh3;
                } //end sh3 loop
                sh2_dptr += stride2 * nfxn_sh2;
            } //end sh2 loop
            sh1_dptr += stride1 * nfxn_sh1;
        } //end sh1 loop
    }
}



MatrixFiller::MatrixFiller(SharedMatrix matrix, yeti::uli poff, yeti::uli qoff) :
    callcount_(0),
    poff_(poff),
    qoff_(qoff),
    matrix_(matrix)
{
}

TensorElementComputer*
MatrixFiller::copy() const
{
    return new MatrixFiller(matrix_, poff_, qoff_);
}

void
MatrixFiller::compute(const uli* indices, double* data, uli n){
    double* dptr = data;
    uli pStart = descr_->get(0)->index_start(indices[0]) - poff_;
    uli qStart = descr_->get(1)->index_start(indices[1]) - qoff_;
    uli np     = descr_->get(0)->nelements(indices[0]);
    uli nq     = descr_->get(1)->nelements(indices[1]);
    uli pEnd   = pStart + np;
    uli qEnd   = qStart + nq;
    for(uli p = pStart; p < pEnd; ++p){
        for(uli q = qStart; q < qEnd; ++q){
            *dptr = matrix_->get(0, p, q);
            ++dptr;
        }
    }
}


VectorFiller::VectorFiller(SharedVector vector, yeti::uli offset) :
    callcount_(0),
    offset_(offset),
    vector_(vector)
{
}

TensorElementComputer*
VectorFiller::copy() const
{
    return new VectorFiller(vector_, offset_);
}

void
VectorFiller::compute(const uli* indices, double* data, uli n) {
    double* dptr = data;
    uli pStart = descr_->get(0)->index_start(indices[0]) - offset_;
    uli np     = descr_->get(0)->nelements(indices[0]);
    uli pEnd   = pStart + np;
    for(uli p = pStart; p < pEnd; ++p){
        *dptr = vector_->get(0, p);
        ++dptr;
    }
}

DPDFiller::DPDFiller(dpdbuf4 *buf, yeti::uli poff, yeti::uli qoff, yeti::uli roff, yeti::uli soff) :
    callcount_(0),
    poff_(poff),
    qoff_(qoff),
    roff_(roff),
    soff_(soff),
    buf_(buf)
{
}

TensorElementComputer*
DPDFiller::copy() const
{
    return new DPDFiller(buf_, poff_, qoff_, roff_, soff_);
}

void
DPDFiller::compute(const uli* indices, double* data, yeti::uli n) {
    double* dptr = data;
    uli pStart = descr_->get(0)->index_start(indices[0]) - poff_;
    uli qStart = descr_->get(1)->index_start(indices[1]) - qoff_;
    uli rStart = descr_->get(2)->index_start(indices[2]) - roff_;
    uli sStart = descr_->get(3)->index_start(indices[3]) - soff_;
    uli np     = descr_->get(0)->nelements(indices[0]);
    uli nq     = descr_->get(1)->nelements(indices[1]);
    uli nr     = descr_->get(2)->nelements(indices[2]);
    uli ns     = descr_->get(3)->nelements(indices[3]);
    uli pEnd   = pStart + np;
    uli qEnd   = qStart + nq;
    uli rEnd   = rStart + nr;
    uli sEnd   = sStart + ns;
    for(uli p = pStart; p < pEnd; ++p){
        for(uli q = qStart; q < qEnd; ++q){
            for(uli r = rStart; r < rEnd; ++r){
                for(uli s = sStart; s < sEnd; ++s){
                    size_t pq = buf_->params->rowidx[p][q];
                    size_t rs = buf_->params->colidx[r][s];
                    *dptr = buf_->matrix[0][pq][rs];
                    ++dptr;
                }
            }
        }
    }
}

#else  // have package
#error You must link to an external package and defining HAVE_PSI, or HAVE_MPQC
#endif //end have package

