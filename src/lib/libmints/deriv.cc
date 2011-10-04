/*
 * deriv.cc
 *
 *  Created on: Feb 24, 2009
 *      Author: jturney
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <boost/foreach.hpp>

#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>

#include "mints.h"

using namespace std;

namespace psi {

size_t counter;

class CorrelatedRestrictedFunctor
{
    boost::shared_ptr<Wavefunction> wavefunction_;

    dpdbuf4 G_;
public:
    int nthread;
    std::vector<SharedVector> result;

    CorrelatedRestrictedFunctor() {
        throw PSIEXCEPTION("CorrelatedRestrictedFunctor(): Default constructor called. This shouldn't happen.");
    }

    CorrelatedRestrictedFunctor(SharedVector results, boost::shared_ptr<Wavefunction> wave, IntegralTransform& ints_transform)
        : wavefunction_(wave)
    {
        _default_psio_lib_->open(PSIF_TPDM_HALFTRANS, PSIO_OPEN_OLD);
        dpd_buf4_init(&G_, PSIF_TPDM_HALFTRANS, 0,
                      ints_transform.DPD_ID("[n,n]"), ints_transform.DPD_ID("[n,n]"),
                      ints_transform.DPD_ID("[n>=n]+"), ints_transform.DPD_ID("[n>=n]+"),
                      0, "SO Basis TPDM (nn|nn)");

        for (int h=0; h<wavefunction_->nirrep(); ++h) {
            dpd_buf4_mat_irrep_init(&G_, h);
            dpd_buf4_mat_irrep_rd(&G_, h);
        }

        nthread = Communicator::world->nthread();
        result.push_back(results);
        for (int i=1; i<nthread; ++i)
            result.push_back(SharedVector(result[0]->clone()));
    }

    void finalize() {
        for (int h=0; h<wavefunction_->nirrep(); ++h)
            dpd_buf4_mat_irrep_close(&G_, h);

        dpd_buf4_close(&G_);
        _default_psio_lib_->close(PSIF_TPDM_HALFTRANS, 1);

        // Do summation over threads
        for (int i=1; i<nthread; ++i) {
            result[0]->add(result[i]);
        }
        // Do MPI global summation
        result[0]->sum();
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs,
                    int pirrep, int pso,
                    int qirrep, int qso,
                    int rirrep, int rso,
                    int sirrep, int sso,
                    double value)
    {
        int thread = Communicator::world->thread_id(pthread_self());

        double prefactor = 8.0;

        if (pirrep ^ qirrep ^ rirrep ^ sirrep)
            return;

        int h = pirrep ^ qirrep;

        if (pabs == qabs)
            prefactor *= 0.5;
        if (rabs == sabs)
            prefactor *= 0.5;
        if ((pabs == rabs && qabs == sabs) || (pabs == sabs && qabs == rabs))
            prefactor *= 0.5;

        int PQ = G_.params->colidx[pabs][qabs];   // pabs, qabs?
        int RS = G_.params->rowidx[rabs][sabs];   // pabs, qabs?

        result[thread]->add(salc, prefactor * G_.matrix[h][PQ][RS] * value);
    }
};

class DFMP2RestrictedFunctor
{
    SharedMatrix P_2_;
    SharedMatrix D_;

public:
    int nthread;
    std::vector<SharedVector> result;

    DFMP2RestrictedFunctor() {
        throw PSIEXCEPTION("DFMP2RestrictedFunctor(): Default constructor called. This shouldn't happen.");
    }

    DFMP2RestrictedFunctor(SharedVector results,
                           boost::shared_ptr<Matrix> P_2,
                           boost::shared_ptr<Matrix> D)
        : P_2_(P_2), D_(D)
    {
        counter=0;
        nthread = Communicator::world->nthread();
        result.push_back(results);

        for (int i=1; i<nthread; ++i)
            result.push_back(SharedVector(result[0]->clone()));
    }
    ~DFMP2RestrictedFunctor() {
    }

    void finalize() {
        // Do summation over threads
        for (int i=1; i<nthread; ++i) {
            result[0]->add(result[i]);
        }
        // Do MPI global summation
        result[0]->sum();
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs,
                    int pirrep, int pso,
                    int qirrep, int qso,
                    int rirrep, int rso,
                    int sirrep, int sso,
                    double value)
    {
        int thread = Communicator::world->thread_id(pthread_self());

        // Previously, we applied a factor of 4 after the fact...apply it from the beginning now.
        double prefactor = 4.0;

        if (pabs == qabs)
            prefactor *= 0.5;
        if (rabs == sabs)
            prefactor *= 0.5;
        if (INDEX2(pabs, qabs) == INDEX2(rabs, sabs))
            prefactor *= 0.5;

        double four_index_D = 0.0;

        if (pirrep == qirrep && rirrep == sirrep)
            four_index_D = P_2_->get(pirrep, pso, qso) * D_->get(rirrep, rso, sso);
        if (pirrep == rirrep && qirrep == sirrep)
            four_index_D -= P_2_->get(pirrep, pso, rso) * D_->get(qirrep, qso, sso);
        if (pirrep == sirrep && qirrep == rirrep)
            four_index_D -= P_2_->get(pirrep, pso, sso) * D_->get(qirrep, qso, rso);

        four_index_D *= prefactor;

        result[thread]->add(salc, four_index_D * value);
        counter++;
    }
};

class ScfRestrictedFunctor
{
    SharedMatrix D_;

public:
    int nthread;
    std::vector<SharedVector> result;

    ScfRestrictedFunctor() {
        throw PSIEXCEPTION("ScfRestrictedFunctor(): Default constructor called. This shouldn't happen.");
    }

    // Added for debugging MADNESS.
//    ScfRestrictedFunctor(const ScfRestrictedFunctor&) {
//        throw PSIEXCEPTION("ScfRestrictedFunctor(): Copy constructor called.\n");
//    }

//    ScfRestrictedFunctor& operator=(const ScfRestrictedFunctor&) {
//        throw PSIEXCEPTION("ScfRestrictedFunctor(): Assignment operator called. This shouldn't happen.");
//        return *this;
//    }

    ScfRestrictedFunctor(SharedVector results, boost::shared_ptr<Matrix> D)
        : D_(D)
    {
        counter=0;
        nthread = Communicator::world->nthread();
        result.push_back(results);

        for (int i=1; i<nthread; ++i)
            result.push_back(SharedVector(result[0]->clone()));
    }
    ~ScfRestrictedFunctor() {
    }

    void finalize() {
        // Do summation over threads
        for (int i=1; i<nthread; ++i) {
            result[0]->add(result[i]);
        }
        // Do MPI global summation
        result[0]->sum();
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs,
                    int pirrep, int pso,
                    int qirrep, int qso,
                    int rirrep, int rso,
                    int sirrep, int sso,
                    double value)
    {
        int thread = Communicator::world->thread_id(pthread_self());

        // Previously, we applied a factor of 4 after the fact...apply it from the beginning now.
        double prefactor = 4.0;

        if (pabs == qabs)
            prefactor *= 0.5;
        if (rabs == sabs)
            prefactor *= 0.5;
        if (INDEX2(pabs, qabs) == INDEX2(rabs, sabs))
            prefactor *= 0.5;

        double four_index_D = 0.0;

        if (pirrep == qirrep && rirrep == sirrep)
            four_index_D = 4.0 * D_->get(pirrep, pso, qso) * D_->get(rirrep, rso, sso);
        if (pirrep == rirrep && qirrep == sirrep)
            four_index_D -= D_->get(pirrep, pso, rso) * D_->get(qirrep, qso, sso);
        if (pirrep == sirrep && qirrep == rirrep)
            four_index_D -= D_->get(pirrep, pso, sso) * D_->get(qirrep, qso, rso);

        four_index_D *= prefactor;

        result[thread]->add(salc, four_index_D * value);
        counter++;
    }
};

class ScfUnrestrictedFunctor
{
    boost::shared_ptr<Wavefunction> wavefunction_;
    SharedMatrix Da_;
    SharedMatrix Db_;

public:
    int nthread;
    std::vector<SharedVector> result;

    ScfUnrestrictedFunctor() { throw PSIEXCEPTION("ScfUnrestrictedFunctor(): Oh come on!!!"); }

    ScfUnrestrictedFunctor(SharedVector results, boost::shared_ptr<Wavefunction> wave)
        : wavefunction_(wave),
          Da_(wave->Da()),
          Db_(wave->Db())
    {
        nthread = Communicator::world->nthread();
        result.push_back(results);
        for (int i=1; i<nthread; ++i) {
            result.push_back(SharedVector(result[0]->clone()));
        }
    }
    ~ScfUnrestrictedFunctor() {
    }

    void finalize() {
        // Do summation over threads
        for (int i=1; i<nthread; ++i) {
            result[0]->add(result[i]);
        }
        // Do MPI global summation
        result[0]->sum();
    }

    void operator()(int salc, int pabs, int qabs, int rabs, int sabs,
                    int pirrep, int pso,
                    int qirrep, int qso,
                    int rirrep, int rso,
                    int sirrep, int sso,
                    double value)
    {
        int thread = Communicator::world->thread_id(pthread_self());
        double prefactor = 1.0;

        if (pabs == qabs)
            prefactor *= 0.5;
        if (rabs == sabs)
            prefactor *= 0.5;
        if ((pabs == rabs && qabs == sabs) || (pabs == sabs && qabs == rabs))
            prefactor *= 0.5;

        double four_index_D = 0.0;

        if (pirrep == qirrep && rirrep == sirrep) {
            four_index_D = 4.0 * (Da_->get(pirrep, pso, qso) + Db_->get(pirrep, pso, qso)) *
                                 (Da_->get(rirrep, rso, sso) + Db_->get(rirrep, rso, sso));
        }
        if (pirrep == rirrep && qirrep == sirrep) {
            four_index_D -= 2.0 * ((Da_->get(pirrep, pso, rso) * Da_->get(qirrep, qso, sso))
                                 + (Db_->get(pirrep, pso, rso) * Db_->get(qirrep, qso, sso)));
        }
        if (pirrep == sirrep && qirrep == rirrep) {
            four_index_D -= 2.0 * ((Da_->get(pirrep, pso, sso) * Da_->get(rirrep, rso, qso))
                                 + (Db_->get(pirrep, pso, sso) * Db_->get(rirrep, rso, qso)));
        }
//        four_index_D *= prefactor;
        value *= prefactor;

        result[thread]->add(salc, four_index_D * value);
    }
};

Deriv::Deriv(const boost::shared_ptr<Wavefunction>& wave,
             char needed_irreps,
             bool project_out_translations,
             bool project_out_rotations)
    : wavefunction_(wave),
      cdsalcs_(wave->molecule(),
          wave->matrix_factory(),
          needed_irreps,
          project_out_translations,
          project_out_rotations)
{
    integral_ = wave->integral();
    basis_    = wave->basisset();
    sobasis_  = wave->sobasisset();
    factory_  = wave->matrix_factory();
    molecule_ = wave->molecule();
    natom_    = molecule_->natom();

    // Results go here.
    opdm_contr_ = factory_->create_shared_matrix("One-electron contribution to gradient", natom_, 3);
    x_contr_    = factory_->create_shared_matrix("Lagrangian contribution to gradient", natom_, 3);
    tpdm_contr_ = factory_->create_shared_matrix("Two-electron contribution to gradient", natom_, 3);
    gradient_   = factory_->create_shared_matrix("Total gradient", natom_, 3);
}

Deriv::Deriv(const boost::shared_ptr<BasisSet>& basis,
             const boost::shared_ptr<Matrix>& P_2,
             const boost::shared_ptr<Matrix>& W_2,
             const boost::shared_ptr<Matrix>& SCF_D,
             const boost::shared_ptr<MatrixFactory>& factory,
             char needed_irreps,
             bool project_out_translations,
             bool project_out_rotations)
    : basis_(basis),
      cdsalcs_(basis->molecule(),
          factory,
          needed_irreps,
          project_out_translations,
          project_out_rotations),
      P_2_(P_2),
      W_2_(W_2),
      SCF_D_(SCF_D),
      factory_(factory)
{
    integral_ = boost::shared_ptr<IntegralFactory> (new IntegralFactory(basis,basis,basis,basis));
    sobasis_  = boost::shared_ptr<SOBasisSet> (new SOBasisSet(basis,integral_));
    molecule_ = basis->molecule();
    natom_    = molecule_->natom();

    // Results go here.
    opdm_contr_ = factory_->create_shared_matrix("One-electron contribution to gradient", natom_, 3);
    x_contr_    = factory_->create_shared_matrix("Lagrangian contribution to gradient", natom_, 3);
    tpdm_contr_ = factory_->create_shared_matrix("Two-electron contribution to gradient", natom_, 3);
    gradient_   = factory_->create_shared_matrix("Total gradient", natom_, 3);
}

SharedMatrix Deriv::compute()
{
    molecule_->print_in_bohr();

    // Initialize an ERI object requesting derivatives.
    std::vector<boost::shared_ptr<TwoBodyAOInt> > ao_eri;
    for (int i=0; i<Communicator::world->nthread(); ++i)
        ao_eri.push_back(boost::shared_ptr<TwoBodyAOInt>(integral_->eri(1)));
    TwoBodySOInt so_eri(ao_eri, integral_, cdsalcs_);

    // A certain optimization can be used if we know we only need totally symmetric
    // derivatives.
    so_eri.set_only_totally_symmetric(false);

    // Compute one-electron derivatives.
    vector<SharedMatrix> s_deriv = cdsalcs_.create_matrices("S'");
    vector<SharedMatrix> h_deriv = cdsalcs_.create_matrices("H'");

    boost::shared_ptr<OneBodySOInt> s_int(integral_->so_overlap(1));
    boost::shared_ptr<OneBodySOInt> t_int(integral_->so_kinetic(1));
    boost::shared_ptr<OneBodySOInt> v_int(integral_->so_potential(1));

    s_int->compute_deriv1(s_deriv, cdsalcs_);
    t_int->compute_deriv1(h_deriv, cdsalcs_);
    v_int->compute_deriv1(h_deriv, cdsalcs_);

    double *Xcont = new double[cdsalcs_.ncd()];
    double *Dcont = new double[cdsalcs_.ncd()];
    double *TPDMcont = new double[cdsalcs_.ncd()];

    int ncd = cdsalcs_.ncd();
    SharedVector TPDMcont_vector(new Vector(1, &ncd));

    if (!wavefunction_ || (wavefunction_ && !wavefunction_->reference_wavefunction())) {
        if (P_2_ || (wavefunction_ && wavefunction_->restricted())) {
            SharedMatrix D = wavefunction_ ? wavefunction_->Da() : P_2_;
            SharedMatrix X = wavefunction_ ? wavefunction_->X() : W_2_;

            // Check the incoming matrices.
            if (!D)
                throw PSIEXCEPTION("Deriv::compute: Unable to access OPDM.");
            if (!X)
                throw PSIEXCEPTION("Deriv::compute: Unable to access Lagrangian.");

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                double temp = 2.0 * D->vector_dot(h_deriv[cd]);
                Dcont[cd] = temp;
                fprintf(outfile, "    SALC #%d One-electron contribution: %+lf\n", cd, temp);
            }

            fprintf(outfile, "\n");

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                double temp = -2.0 * X->vector_dot(s_deriv[cd]);
                Xcont[cd] = temp;
                fprintf(outfile, "    SALC #%d Lagrandian contribution:   %+lf\n", cd, temp);
            }

            fprintf(outfile, "\n");

            if (wavefunction_) {
                ScfRestrictedFunctor functor(TPDMcont_vector, D);
                so_eri.compute_integrals_deriv1(functor);
                functor.finalize();
            }
            else {
                DFMP2RestrictedFunctor functor(TPDMcont_vector, D, SCF_D_);
                so_eri.compute_integrals_deriv1(functor);
                functor.finalize();
            }

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                TPDMcont[cd] = TPDMcont_vector->get(cd);
                fprintf(outfile, "    SALC #%d TPDM contribution:         %+lf\n", cd, TPDMcont[cd]);
            }
            fflush(outfile);
        }
        else /* unrestricted */ {
            SharedMatrix Da = wavefunction_->Da();
            SharedMatrix Db = wavefunction_->Db();
            SharedMatrix X = wavefunction_->X();

            // Check the incoming matrices.
            if (!Da || !Db)
                throw PSIEXCEPTION("Deriv::compute: Unable to access OPDM.");
            if (!X)
                throw PSIEXCEPTION("Deriv::compute: Unable to access Lagrangian.");

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                double temp = Da->vector_dot(h_deriv[cd]);
                temp       += Db->vector_dot(h_deriv[cd]);
                Dcont[cd] = temp;
                fprintf(outfile, "    SALC #%d One-electron contribution: %+lf\n", cd, temp);
            }

            fprintf(outfile, "\n");

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                double temp = X->vector_dot(s_deriv[cd]);
                Xcont[cd] = -temp;
                fprintf(outfile, "    SALC #%d Lagrangian contribution:   %+lf\n", cd, temp);
            }

            fprintf(outfile, "\n");

            ScfUnrestrictedFunctor functor(TPDMcont_vector, wavefunction_);
            so_eri.compute_integrals_deriv1(functor);
            functor.finalize();

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                TPDMcont[cd] = TPDMcont_vector->get(cd);
                fprintf(outfile, "    SALC #%d TPDM contribution:         %+lf\n", cd, TPDMcont[cd]);
            }
            fflush(outfile);
        }
    }
    else { /* correlated */
        if (wavefunction_->restricted()) {
            // Define the MO orbital space we need
            vector<boost::shared_ptr<MOSpace> > spaces;
            spaces.push_back(MOSpace::all);

            // Uses a different constructor of IntegralTransform.
            IntegralTransform ints_transform(wavefunction_,
                                             spaces,
                                             IntegralTransform::Restricted, // Transformation type
                                             IntegralTransform::DPDOnly,    // Output buffer
                                             IntegralTransform::QTOrder,    // MO ordering
                                             IntegralTransform::None);      // Frozen orbitals?

            dpd_set_default(ints_transform.get_dpd_id());

            ints_transform.backtransform_density();

            SharedMatrix D = factory_->create_shared_matrix("SO-basis OPDM");
            D->load(_default_psio_lib_, PSIF_AO_OPDM);

            SharedMatrix X = factory_->create_shared_matrix("SO-basis Lagrangian");
            X->load(_default_psio_lib_, PSIF_AO_OPDM);

            // Check the incoming matrices.
            if (!D)
                throw PSIEXCEPTION("Deriv::compute: Unable to access OPDM.");
            if (!X)
                throw PSIEXCEPTION("Deriv::compute: Unable to access Lagrangian.");

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                double temp = D->vector_dot(h_deriv[cd]);
                Dcont[cd] = temp;
                fprintf(outfile, "    SALC #%d One-electron contribution: %+lf\n", cd, temp);
            }

            fprintf(outfile, "\n");

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                double temp = -0.5 * X->vector_dot(s_deriv[cd]);
                Xcont[cd] = temp;
                fprintf(outfile, "    SALC #%d Lagrandian contribution:   %+lf\n", cd, temp);
            }

            fprintf(outfile, "\n");

            CorrelatedRestrictedFunctor functor(TPDMcont_vector, wavefunction_, ints_transform);
            so_eri.compute_integrals_deriv1(functor);
            functor.finalize();

            for (int cd=0; cd < cdsalcs_.ncd(); ++cd) {
                TPDMcont[cd] = TPDMcont_vector->get(cd);
                fprintf(outfile, "    SALC #%d TPDM contribution:         %+lf\n", cd, TPDMcont[cd]);
            }
            fflush(outfile);
        }
        else {
            throw PSIEXCEPTION("Unrestricted correlated gradients not implemented.");
        }
    }

    // Transform the SALCs back to cartesian space
    SharedMatrix st = cdsalcs_.matrix();
    double **B = st->pointer(0);
    double *cart = new double[3*natom_];

    // B^t g_q^t = g_x^t -> g_q B = g_x
    C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
            1.0, Dcont, cdsalcs_.ncd(), B[0],
            3*natom_, 0.0, cart, 3*natom_);

    for (int a=0; a<natom_; ++a)
        for (int xyz=0; xyz<3; ++xyz)
            opdm_contr_->set(a, xyz, cart[3*a+xyz]);

    // B^t g_q^t = g_x^t -> g_q B = g_x
    C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
            1.0, Xcont, cdsalcs_.ncd(), B[0],
            3*natom_, 0.0, cart, 3*natom_);

    for (int a=0; a<natom_; ++a)
        for (int xyz=0; xyz<3; ++xyz)
            x_contr_->set(a, xyz, cart[3*a+xyz]);

    // B^t g_q^t = g_x^t -> g_q B = g_x
    C_DGEMM('n', 'n', 1, 3*natom_, cdsalcs_.ncd(),
            1.0, TPDMcont, cdsalcs_.ncd(), B[0],
            3*natom_, 0.0, cart, 3*natom_);

    for (int a=0; a<natom_; ++a)
        for (int xyz=0; xyz<3; ++xyz)
            tpdm_contr_->set(a, xyz, cart[3*a+xyz]);

    delete[] Dcont;
    delete[] Xcont;
    delete[] TPDMcont;

    // Obtain nuclear repulsion contribution from the wavefunction
    Matrix enuc = molecule_->nuclear_repulsion_energy_deriv1();

    // Print things out
    enuc.print_atom_vector();
    opdm_contr_->print_atom_vector();
    x_contr_->print_atom_vector();
    tpdm_contr_->print_atom_vector();

    // Add everything up into a temp.
    Matrix temp("Temp SCF gradient", molecule_->natom(), 3);
    temp.add(&enuc);
    temp.add(opdm_contr_);
    temp.add(x_contr_);
    temp.add(tpdm_contr_);

    // Symmetrize the gradients to remove any noise:
    CharacterTable ct = molecule_->point_group()->char_table();

    // Obtain atom mapping of atom * symm op to atom
    int **atom_map = compute_atom_map(molecule_);

    // Symmetrize the gradients to remove any noise
    for (int atom=0; atom<molecule_->natom(); ++atom) {
        for (int g=0; g<ct.order(); ++g) {

            int Gatom = atom_map[atom][g];

            SymmetryOperation so = ct.symm_operation(g);

            gradient_->add(atom, 0, so(0, 0) * temp(Gatom, 0) / ct.order());
            gradient_->add(atom, 1, so(1, 1) * temp(Gatom, 1) / ct.order());
            gradient_->add(atom, 2, so(2, 2) * temp(Gatom, 2) / ct.order());
        }
    }

    // Delete the atom map.
    delete_atom_map(atom_map, molecule_);

    // Print the atom vector
    gradient_->print_atom_vector();

    // Save the gradient to the wavefunction so that optking can optimize with it
    if (wavefunction_) wavefunction_->set_gradient(gradient_);

    return gradient_;
}

}

#ifdef HAVE_MADNESS
namespace madness {
namespace archive {

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::ScfRestrictedFunctor> {
    static void store(const Archive &ar, const psi::ScfRestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::ScfUnrestrictedFunctor> {
    static void store(const Archive &ar, const psi::ScfUnrestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveStoreImpl< Archive, psi::CorrelatedRestrictedFunctor> {
    static void store(const Archive &ar, const psi::CorrelatedRestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::ScfRestrictedFunctor> {
    static void load(const Archive &ar, const psi::ScfRestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::ScfUnrestrictedFunctor> {
    static void load(const Archive &ar, const psi::ScfUnrestrictedFunctor &t) {
    }
};

template <class Archive>
struct ArchiveLoadImpl< Archive, psi::CorrelatedRestrictedFunctor> {
    static void load(const Archive &ar, const psi::CorrelatedRestrictedFunctor &t) {
    }
};

}
}
#endif // HAVE_MADNESS
