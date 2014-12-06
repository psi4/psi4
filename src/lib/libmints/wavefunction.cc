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

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <liboptions/liboptions.h>
#include <libparallel/parallel.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libparallel/parallel.h>
#include <libmints/view.h>

#include "mints.h"
#include "orbitalspace.h"

#include <psi4-dec.h>

#include <boost/python.hpp>
#include <boost/python/call.hpp>

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

using namespace boost;
using namespace psi;

// Globals
size_t ioff[MAX_IOFF];
double df[MAX_DF];
double bc[MAX_BC][MAX_BC];
double fac[MAX_FAC];

Wavefunction::Wavefunction(Options & options, boost::shared_ptr<PSIO> psio) :
    options_(options), psio_(psio)
{
    chkpt_ = boost::shared_ptr<Chkpt>(new Chkpt(psio.get(), PSIO_OPEN_OLD));
    common_init();
}

Wavefunction::Wavefunction(Options & options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
    options_(options), psio_(psio), chkpt_(chkpt)
{
    common_init();
}

Wavefunction::~Wavefunction()
{
}

void Wavefunction::copy(boost::shared_ptr<Wavefunction> other)
{
    name_ = other->name_;
    basisset_ = other->basisset_;
    sobasisset_ = other->sobasisset_;
    AO2SO_ = other->AO2SO_;
    S_ = other->S_;
    H_ = other->H_;
    molecule_ = other->molecule_;
    psio_ = other->psio_;
    chkpt_ = other->chkpt_;
    integral_ = other->integral_;
    factory_ = other->factory_;
    reference_wavefunction_ = other;
    memory_ = other->memory_;
    print_ = other->print_;
    debug_ = other->debug_;
    density_fitted_ = other->density_fitted_;
    energy_threshold_ = other->energy_threshold_;
    density_threshold_ = other->density_threshold_;
    nalpha_ = other->nalpha_;
    nbeta_ = other->nbeta_;
    nfrzc_ = other->nfrzc_;

    doccpi_ = other->doccpi_;
    soccpi_ = other->soccpi_;
    frzcpi_ = other->frzcpi_;
    frzvpi_ = other->frzvpi_;
    nalphapi_ = other->nalphapi_;
    nbetapi_ = other->nbetapi_;
    nsopi_ = other->nsopi_;
    nmopi_ = other->nmopi_;

    energy_ = other->energy_;

    nso_ = other->nso_;
    nmo_ = other->nmo_;
    nirrep_ = other->nirrep_;

    Ca_ = other->Ca_;
    Cb_ = other->Cb_;
    Da_ = other->Da_;
    Db_ = other->Db_;
    Fa_ = other->Fa_;
    Fb_ = other->Fb_;
    epsilon_a_ = other->epsilon_a_;
    epsilon_b_ = other->epsilon_b_;

    gradient_ = other->gradient_;
    tpdm_gradient_contribution_ = other->tpdm_gradient_contribution_;

    QLMO_to_LMO_ = other->QLMO_to_LMO_;
    CIM_orbital_factors_ = other->CIM_orbital_factors_;
    CIM_orbital_energies_ = other->CIM_orbital_energies_;
    CIM_nactive_occupied_pointer_ = &other->CIM_nactive_occupied_;
    CIM_nactive_virtual_pointer_ = &other->CIM_nactive_virtual_;
    isCIM_ = other->isCIM_;
}

void Wavefunction::common_init()
{
    Wavefunction::initialize_singletons();

    // Take the molecule from the environment
    molecule_ = Process::environment.molecule();

    // Load in the basis set
    basisset_ = BasisSet::pyconstruct_orbital(molecule_, "BASIS", options_.get_str("BASIS"));

    // Check the point group of the molecule. If it is not set, set it.
    if (!molecule_->point_group()) {
        molecule_->set_point_group(molecule_->find_point_group());
    }

    // Create an SO basis...we need the point group for this part.
    integral_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
    sobasisset_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, integral_));

    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
    AO2SO_ = pet->aotoso();

    // Obtain the dimension object to initialize the factory.
    const Dimension dimension = sobasisset_->dimension();
    factory_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    factory_->init_with(dimension, dimension);

    nirrep_ = dimension.n();

    S_ = factory_->create_shared_matrix("S");
    boost::shared_ptr<OneBodySOInt> Sint(integral_->so_overlap());
    Sint->compute(S_);

    // Initialize array that hold dimensionality information
    nsopi_    = Dimension(nirrep_, "SOs per irrep");
    nmopi_    = Dimension(nirrep_, "MOs per irrep");
    nalphapi_ = Dimension(nirrep_, "Alpha electrons per irrep");
    nbetapi_  = Dimension(nirrep_, "Beta electrons per irrep");
    doccpi_   = Dimension(nirrep_, "Doubly occupied orbitals per irrep");
    soccpi_   = Dimension(nirrep_, "Singly occupied orbitals per irrep");
    frzcpi_   = Dimension(nirrep_, "Frozen core orbitals per irrep");
    frzvpi_   = Dimension(nirrep_, "Frozen virtual orbitals per irrep");

    // Obtain memory amount from the environment
    memory_ = Process::environment.get_memory();

    nso_ = basisset_->nbf();
    nmo_ = basisset_->nbf();
    for (int k = 0; k < nirrep_; k++) {
        nsopi_[k] = dimension[k];
        nmopi_[k] = 0;
        doccpi_[k] = 0;
        soccpi_[k] = 0;
        nalphapi_[k] = 0;
        nbetapi_[k] = 0;
    }

    energy_ = 0.0;

    // Read in the debug flag
    debug_ = options_.get_int("DEBUG");
    print_ = options_.get_int("PRINT");

    // Read in energy convergence threshold
    energy_threshold_ = options_.get_double("E_CONVERGENCE");

    // Read in density convergence threshold
    density_threshold_ = options_.get_double("D_CONVERGENCE");;

    density_fitted_ = false;

    // not a CIM computation by default
    isCIM_ = false;
}

void Wavefunction::map_irreps(std::vector<int*> &arrays)
{
    boost::shared_ptr<PointGroup> full = Process::environment.parent_symmetry();
    // If the parent symmetry hasn't been set, no displacements have been made
    if(!full) return;
    boost::shared_ptr<PointGroup> sub = molecule_->point_group();
    // Build the correlation table between full, and subgroup
    CorrelationTable corrtab(full, sub);
    int nirreps = corrtab.n();
    std::vector<int*>::iterator iter = arrays.begin();
    for(; iter != arrays.end(); ++iter){
        int *array = *iter;
        int temp[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        for(int h = 0; h < nirreps; ++h){
            int target = corrtab.gamma(h, 0);
            temp[target] += array[h];
        }
        for(int h = 0; h < nirreps; ++h)
            array[h] = temp[h];
    }
}

void Wavefunction::map_irreps(int* &array)
{
    std::vector<int*> vec;
    vec.push_back(array);
    map_irreps(vec);
}

void Wavefunction::map_irreps(Dimension &array)
{
    int *int_array = array;
    std::vector<int*> vec;
    vec.push_back(int_array);
    map_irreps(vec);
}

double Wavefunction::compute_energy()
{
    return 0.0;
}

void Wavefunction::initialize_singletons()
{
    static bool done = false;

    if (done)
        return;

    ioff[0] = 0;
    for (size_t i=1; i<MAX_IOFF; ++i)
        ioff[i] = ioff[i-1] + i;

    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for (int i=3; i<MAX_DF; ++i)
        df[i] = (i-1)*df[i-2];

    for (int i=0; i<MAX_BC; ++i)
        for (int j=0; j<=i; ++j)
            bc[i][j] = combinations(i, j);

    fac[0] = 1.0;
    for (int i=1; i<MAX_FAC; ++i)
        fac[i] = i*fac[i-1];

    done = true;
}

void Wavefunction::semicanonicalize()
{
    throw PSIEXCEPTION("This type of wavefunction cannot be semicanonicalized!");
}

boost::shared_ptr<Molecule> Wavefunction::molecule() const
{
    return molecule_;
}

boost::shared_ptr<PSIO> Wavefunction::psio() const
{
    return psio_;
}

Options& Wavefunction::options() const
{
    return options_;
}

boost::shared_ptr<IntegralFactory> Wavefunction::integral() const
{
    return integral_;
}

boost::shared_ptr<BasisSet> Wavefunction::basisset() const
{
    return basisset_;
}

boost::shared_ptr<SOBasisSet> Wavefunction::sobasisset() const
{
    return sobasisset_;
}

boost::shared_ptr<MatrixFactory> Wavefunction::matrix_factory() const
{
    return factory_;
}

boost::shared_ptr<Wavefunction> Wavefunction::reference_wavefunction() const
{
    return reference_wavefunction_;
}

void Wavefunction::set_reference_wavefunction(const boost::shared_ptr<Wavefunction> wfn)
{
    reference_wavefunction_ = wfn;
}

void Wavefunction::add_preiteration_callback(PyObject *pyobject)
{
    if (pyobject != Py_None)
        precallbacks_.push_back(pyobject);
}

void Wavefunction::add_postiteration_callback(PyObject *pyobject)
{
    if (pyobject != Py_None)
        postcallbacks_.push_back(pyobject);
}

void Wavefunction::call_preiteration_callbacks()
{
    std::vector<void*>::const_iterator iter;
    for (iter = precallbacks_.begin(); iter != precallbacks_.end(); ++iter) {
        if ((PyObject*)*iter != Py_None)
            boost::python::call<void>((PyObject*)*iter, Process::environment.wavefunction());
    }
}
void Wavefunction::call_postiteration_callbacks()
{
    std::vector<void*>::const_iterator iter;
    for (iter = postcallbacks_.begin(); iter != postcallbacks_.end(); ++iter) {
        if ((PyObject*)*iter != Py_None)
            boost::python::call<void>((PyObject*)*iter, Process::environment.wavefunction());
    }
}

SharedMatrix Wavefunction::Ca() const {
    if (!Ca_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Wavefunction::Ca: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Ca();
    }

    return Ca_;
}

SharedMatrix Wavefunction::Cb() const {
    if (!Cb_) {
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Wavefunction::Cb: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Cb();
    }

    return Cb_;
}

std::vector<std::vector<int> > Wavefunction::subset_occupation(const Dimension& noccpi, const std::string& subset)
{
    if (!(subset == "FROZEN_OCC" ||
          subset == "FROZEN_VIR" ||
          subset == "ACTIVE_OCC" ||
          subset == "ACTIVE_VIR" ||
          subset == "FROZEN" ||
          subset == "ACTIVE" ||
          subset == "OCC" ||
          subset == "VIR" ||
          subset == "ALL"))
        throw PSIEXCEPTION("Orbital subset is not defined, should be ALL, OCC, VIR, FROZEN, ACTIVE, or look at doxygen");

    // Vector of relevent positions by irrep
    std::vector<std::vector<int> > positions;
    positions.resize(nirrep_);

    for (int h = 0; h < nirrep_; h++) {
        for (int i = 0; i < frzcpi_[h]; i++) {
            if (subset == "ALL" || subset == "FROZEN_OCC" || subset == "FROZEN" || subset == "OCC")
                positions[h].push_back(i);
        }
        for (int i = frzcpi_[h]; i < noccpi[h]; i++) {
            if (subset == "ALL" || subset == "ACTIVE_OCC" || subset == "ACTIVE" || subset == "OCC")
                positions[h].push_back(i);
        }
        for (int i = noccpi[h]; i < nmopi_[h] - frzvpi_[h]; i++) {
            if (subset == "ALL" || subset == "ACTIVE_VIR" || subset == "ACTIVE" || subset == "VIR")
                positions[h].push_back(i);
        }
        for (int i = nmopi_[h] - frzvpi_[h]; i < nmopi_[h]; i++) {
            if (subset == "ALL" || subset == "FROZEN_VIR" || subset == "FROZEN" || subset == "VIR")
                positions[h].push_back(i);
        }
    }
    return positions;
}
SharedMatrix Wavefunction::C_subset_helper(SharedMatrix C, const Dimension& noccpi, SharedVector epsilon, const std::string& basis, const std::string& subset)
{
    std::vector<std::vector<int> > positions = subset_occupation(noccpi, subset);

    Dimension nmopi(nirrep_);
    for (int h = 0; h < positions.size(); h++) {
        nmopi[h] = positions[h].size();
    }
    SharedMatrix C2(new Matrix("C " + basis + " " + subset, nsopi_, nmopi));
    for (int h = 0; h < positions.size(); h++) {
        for (int i = 0; i < positions[h].size(); i++) {
            C_DCOPY(nsopi_[h], &C->pointer(h)[0][positions[h][i]], nmopi_[h], &C2->pointer(h)[0][i], nmopi[h]);
        }
    }

    if (basis == "AO") {

        SharedMatrix C3(new Matrix("C " + basis + " " + subset, nso_, nmopi.sum()));
        boost::swap(C2,C3);

        std::vector<boost::tuple<double, int, int> > order;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < positions[h].size(); i++) {
                order.push_back(boost::tuple<double,int,int>(epsilon->get(h,positions[h][i]),i,h));
            }
        }

        std::sort(order.begin(), order.end(), std::less<boost::tuple<double,int,int> >());

        for (int index = 0; index < order.size(); index++) {
            int i = boost::get<1>(order[index]);
            int h = boost::get<2>(order[index]);

            int nao = nso_;
            int nso = nsopi_[h];

            if (!nso) continue;

            C_DGEMV('N',nao,nso,1.0,AO2SO_->pointer(h)[0],nso,&C3->pointer(h)[0][i],nmopi[h],0.0,&C2->pointer()[0][index],nmopi.sum());
        }

    } else if (basis == "SO" || basis == "MO") {
        // Already done
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }

    return C2;
}
SharedVector Wavefunction::epsilon_subset_helper(SharedVector epsilon, const Dimension& noccpi, const std::string& basis, const std::string& subset)
{
    std::vector<std::vector<int> > positions = subset_occupation(noccpi, subset);

    Dimension nmopi(nirrep_);
    for (int h = 0; h < positions.size(); h++) {
        nmopi[h] = positions[h].size();
    }

    SharedVector C2;

    if (basis == "AO") {

        C2 = SharedVector(new Vector("Epsilon " + basis + " " + subset, nmopi.sum()));

        std::vector<boost::tuple<double, int, int> > order;
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < positions[h].size(); i++) {
                order.push_back(boost::tuple<double,int,int>(epsilon->get(h,positions[h][i]),i,h));
            }
        }

        std::sort(order.begin(), order.end(), std::less<boost::tuple<double,int,int> >());

        for (int index = 0; index < order.size(); index++) {
            C2->set(0,index,boost::get<0>(order[index]));
        }

    } else if (basis == "SO" || basis == "MO") {

        C2 = SharedVector(new Vector("Epsilon " + basis + " " + subset, nmopi));
        for (int h = 0; h < positions.size(); h++) {
            for (int i = 0; i < positions[h].size(); i++) {
                C2->set(h,i,epsilon->get(h,positions[h][i]));
            }
        }

    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }

    return C2;
}
SharedMatrix Wavefunction::D_subset_helper(SharedMatrix D, SharedMatrix C, const std::string& basis)
{
    if (basis == "AO") {
        double* temp = new double[AO2SO_->max_ncol() * AO2SO_->max_nrow()];
        SharedMatrix D2 = SharedMatrix(new Matrix("D (AO basis)", basisset_->nbf(), basisset_->nbf()));
        int symm = D->symmetry();
        for (int h = 0; h < AO2SO_->nirrep(); ++h) {
            int nao = AO2SO_->rowspi()[0];
            int nsol = AO2SO_->colspi()[h];
            int nsor = AO2SO_->colspi()[h^symm];
            if (!nsol || !nsor) continue;
            double** Ulp = AO2SO_->pointer(h);
            double** Urp = AO2SO_->pointer(h^symm);
            double** DSOp = D->pointer(h^symm);
            double** DAOp = D2->pointer();
            C_DGEMM('N','T',nsol,nao,nsor,1.0,DSOp[0],nsor,Urp[0],nsor,0.0,temp,nao);
            C_DGEMM('N','N',nao,nao,nsol,1.0,Ulp[0],nsol,temp,nao,1.0,DAOp[0],nao);
        }
        delete[] temp;
        return D2;
    } else if (basis == "SO") {
        return SharedMatrix(D->clone());
    } else if (basis == "MO") {
        SharedMatrix D2(new Matrix("D (MO Basis)", C->colspi(), C->colspi()));

        int symm = D->symmetry();
        int nirrep = D->nirrep();

        double* SC = new double[C->max_ncol() * C->max_nrow()];
        double* temp = new double[C->max_ncol() * C->max_nrow()];
        for (int h = 0; h < nirrep; h++) {
            int nmol = C->colspi()[h];
            int nmor = C->colspi()[h^symm];
            int nsol = C->rowspi()[h];
            int nsor = C->rowspi()[h^symm];
            if (!nmol || !nmor || !nsol || !nsor) continue;
            double** Slp = S_->pointer(h);
            double** Srp = S_->pointer(h^symm);
            double** Clp = C->pointer(h);
            double** Crp = C->pointer(h^symm);
            double** Dmop = D2->pointer(h);
            double** Dsop = D->pointer(h);

            C_DGEMM('N','N',nsor,nmor,nsor,1.0,Srp[0],nsor,Crp[0],nmor,0.0,SC,nmor);
            C_DGEMM('N','N',nsol,nmor,nsor,1.0,Dsop[0],nsor,SC,nmor,0.0,temp,nmor);
            C_DGEMM('N','N',nsol,nmol,nsol,1.0,Slp[0],nsol,Clp[0],nmol,0.0,SC,nmol);
            C_DGEMM('T','N',nmol,nmor,nsol,1.0,SC,nmol,temp,nmor,0.0,Dmop[0],nmor);
        }
        delete[] temp;
        delete[] SC;
        return D;
    } else {
        throw PSIEXCEPTION("Invalid basis requested, use AO, SO, or MO");
    }
}

OrbitalSpace Wavefunction::alpha_orbital_space(const std::string& id, const std::string &basis, const std::string &subset)
{
    return OrbitalSpace(id, subset, Ca_subset(basis, subset), epsilon_a_subset(basis, subset), basisset_, integral_);
}

OrbitalSpace Wavefunction::beta_orbital_space(const std::string& id, const std::string &basis, const std::string &subset)
{
    return OrbitalSpace(id, subset, Cb_subset(basis, subset), epsilon_b_subset(basis, subset), basisset_, integral_);
}

SharedMatrix Wavefunction::Ca_subset(const std::string& basis, const std::string& subset)
{
    return C_subset_helper(Ca_, nalphapi_, epsilon_a_, basis, subset);
}
SharedMatrix Wavefunction::Cb_subset(const std::string& basis, const std::string& subset)
{
    return C_subset_helper(Cb_, nbetapi_, epsilon_b_, basis, subset);
}
SharedMatrix Wavefunction::Da_subset(const std::string& basis)
{
    return D_subset_helper(Da_, Ca_, basis);
}
SharedMatrix Wavefunction::Db_subset(const std::string& basis)
{
    return D_subset_helper(Db_, Cb_, basis);
}
SharedVector Wavefunction::epsilon_a_subset(const std::string& basis, const std::string& subset)
{
    return epsilon_subset_helper(epsilon_a_, nalphapi_, basis, subset);
}
SharedVector Wavefunction::epsilon_b_subset(const std::string& basis, const std::string& subset)
{
    return epsilon_subset_helper(epsilon_b_, nbetapi_, basis, subset);
}

SharedMatrix Wavefunction::Fa() const
{
    return Fa_;
}

SharedMatrix Wavefunction::Fb() const
{
    return Fb_;
}

boost::shared_ptr<Matrix> Wavefunction::Lagrangian() const
{
    return Lagrangian_;
}

boost::shared_ptr<Matrix> Wavefunction::tpdm_gradient_contribution() const
{
    throw PSIEXCEPTION("This type of wavefunction has not defined a TPDM gradient contribution!");
}

boost::shared_ptr<Vector> Wavefunction::epsilon_a() const
{
    return epsilon_a_;
}

boost::shared_ptr<Vector> Wavefunction::epsilon_b() const
{
    return epsilon_b_;
}

const SharedMatrix Wavefunction::Da() const
{
    return Da_;
}

SharedMatrix Wavefunction::Db() const
{
    return Db_;
}

SharedMatrix Wavefunction::X() const
{
    return Lagrangian_;
}

SharedMatrix Wavefunction::gradient() const
{
    return gradient_;
}

void Wavefunction::set_gradient(SharedMatrix& grad)
{
    gradient_ = grad;
}

boost::shared_ptr<Vector> Wavefunction::frequencies() const
{
    return frequencies_;
}

void Wavefunction::set_frequencies(boost::shared_ptr<Vector>& freqs)
{
    frequencies_ = freqs;
}

void Wavefunction::save() const
{
}

SharedMatrix Wavefunction::CIMTransformationMatrix()
{
    if (!isCIM_){
       throw PSIEXCEPTION("This is not a CIM computation!");
    }
    return QLMO_to_LMO_;
}

SharedVector Wavefunction::CIMOrbitalEnergies()
{
    if (!isCIM_){
       throw PSIEXCEPTION("This is not a CIM computation!");
    }
    return CIM_orbital_energies_;
}

SharedVector Wavefunction::CIMOrbitalFactors()
{
    if (!isCIM_){
       throw PSIEXCEPTION("This is not a CIM computation!");
    }
    return CIM_orbital_factors_;
}

int Wavefunction::CIMActiveVirtual()
{
    if (!isCIM_){
       throw PSIEXCEPTION("This is not a CIM computation!");
    }
    return CIM_nactive_virtual_;
}

int Wavefunction::CIMActiveOccupied()
{
    if (!isCIM_){
       throw PSIEXCEPTION("This is not a CIM computation!");
    }
    return CIM_nactive_occupied_;
}

void Wavefunction::CIMSet(bool value,int nactive_occupied)
{
    isCIM_ = value;
    if (!value) return;

    QLMO_to_LMO_.reset();
    CIM_orbital_factors_.reset();
    CIM_orbital_energies_.reset();

    QLMO_to_LMO_ = 
        SharedMatrix (new Matrix("CIM Rii",nactive_occupied,nactive_occupied));
    CIM_orbital_factors_ = 
        SharedVector (new Vector("CIM Orbital Factors",nactive_occupied));
    CIM_orbital_energies_ = 
        SharedVector (new Vector("CIM Orbital Energies",nmo_));
}

bool Wavefunction::isCIM()
{
    return isCIM_;
}

void Wavefunction::load_values_from_chkpt() 
{
/*
 Below is a list of all quantities that could be stored in the chkpt file as of August 2014.
 Presumably, not all of them are still in use.  At present, I am bothering to code the
 reading of variables into the wavefunction object necessary for single-point CC
 computations.  Other variables may have to be added to support other kinds of restarts.
- RAK
*/

// Symmetry information
    nirrep_ = chkpt_->rd_nirreps();

// ** Orbital information (closed, frozen, etc., per irrep)
//    int nfzc
//    int nfzv
//    nalphapi_ = Dimension(nirrep_, "Alpha electrons per irrep");
//    nbetapi_  = Dimension(nirrep_, "Beta electrons per irrep");

    int *v; 
    if (chkpt_->exist_add_prefix("SO's per irrep")) {
      v = chkpt_->rd_sopi();
      for (int i=0; i<nirrep_; ++i) nsopi_[i] = v[i];
      free(v);
    }

    if (chkpt_->exist_add_prefix("MO's per irrep")) {
      v = chkpt_->rd_orbspi();
      for (int i=0; i<nirrep_; ++i) nmopi_[i] = v[i];
      free(v);
    }

    if (chkpt_->exist_add_prefix("Open shells per irrep")) {
      v = chkpt_->rd_openpi();
      for (int i=0; i<nirrep_; ++i) soccpi_[i] = v[i];
      free(v);
    }

    if (chkpt_->exist_add_prefix("Closed shells per irrep")) {
      v = chkpt_->rd_clsdpi();
      for (int i=0; i<nirrep_; ++i) doccpi_[i] = v[i];
      free(v);
    }

    if (chkpt_->exist_add_prefix("Frozen DOCC per irrep")) {
      v = chkpt_->rd_frzcpi();
      for (int i=0; i<nirrep_; ++i) frzcpi_[i] = v[i];
      free(v);
    }

    if (chkpt_->exist_add_prefix("Frozen UOCC per irrep")) {
      v = chkpt_->rd_frzvpi();
      for (int i=0; i<nirrep_; ++i) frzvpi_[i] = v[i];
      free(v);
    }

// Necessary for restarting complex calculations (e.g., multiple irrep ones)
/*
int ref
char *prefix             const char label               int disp
int disp_irrep           int ncalcs                     int override_occ
int *rd_statespi(void)   double *rd_vib_freqs(void)
*/

// Variables which store current energies, gradient
    energy_ = chkpt_->rd_escf();
/*
double eref	double ecorr
double efzc	double etot
double eccsd	double e_t	double emp2
double *grad	double e_labeled (varies, const char *key)
*/

// Geometry information
/*
    nuclearrep_ = chkpt_->rd_enuc();
int *atom_dummy
int *atom_position
char **felement
double **fgeom
double **geom
int iopen
int nallatom
int natom
int num_unique_atom
double **rref
struct z_entry *rd_zmat(void);
double *rd_zvals(void);
int rottype
int rd_nfragment(void);
int *rd_natom_per_fragment(void);
int *rd_nallatom_per_fragment(void);
int *rd_nref_per_fragment(void);
double *rd_rotconst(void);
double ***rd_fragment_coeff(void);
int rd_rot_symm_num(void);
*/

    nmo_ = chkpt_->rd_nmo();
    nso_ = chkpt_->rd_nso();

// basis-set/symmetry info
/*
nao  ; initialized by basis set object we hope ?
irrep_labels
int *am2canon_shell_order
double **cartrep
double **cdsalc2cd
double *contr
double **contr_full
double *exps (varies, const char *key)
int **ict
char **irr_labels
char **irr_labels_lowercase
int max_am(varies, const char *key2 = "");
bool puream(varies, const char *key2 = "");
int nprim
int nshell
int nsymhf
int num_unique_shell
int **rd_shell_transm(const char *key2 = "");
int *rd_sloc(const char *key2 = "");
int *rd_shells_per_am(const char *key2 = "");
int *rd_sloc_new(const char *key2 = "");
int *rd_snuc(const char *key2 = "");
int *rd_snumg(const char *key2 = "");
int *rd_ua2a(void);
int *rd_us2s(const char *key2 = "");
double **rd_usotao(const char *key2 = "");
double **rd_usotbf(const char *key2 = "");
int *rd_cdsalcpi();
int *rd_sprim(const char *key2 = "");
int *rd_stype(const char *key2 = "");
char *rd_sym_label(void);
int *rd_symoper(void);
*/

// scf solution
/*
double **ccvecs
double *evals
double *alpha_evals
double *beta_evals
double *fock
double **lagr
double **alpha_lagr
double **beta_lagr
double **rd_scf(void);
double **rd_alpha_scf(void);
double **rd_beta_scf(void);
double **rd_local_scf(void);
double **rd_scf_irrep(int irrep);
double **rd_alpha_scf_irrep(int irrep);
double **rd_beta_scf_irrep(int irrep);
double **set_mo_phases(double**, int, int);
int phase_check
*/
}

