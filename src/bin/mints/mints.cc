#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>

#include <psi4-dec.h>

using namespace boost;

namespace psi { namespace mints {

static int determine_unique_shell_quartets(int usii, int usjj, int uskk, int usll,
                                           int* usi_arr,
                                           int* usj_arr,
                                           int* usk_arr,
                                           int* usl_arr);

class IWLWriter {
    IWL& writeto_;
    size_t count_;

    public:

    IWLWriter(IWL& writeto) : writeto_(writeto), count_(0)
      { }

    void operator()(int i, int j, int k, int l, int , int , int , int , int , int , int , int , double value)
    {
        writeto_.write_value(i, j, k, l, value, 0, NULL, 0);
        count_++;
    }

    size_t count() const { return count_; }
};

PsiReturnType mints(Options & options)
{
    tstart();

    shared_ptr<PSIO> psio = PSIO::shared_object();

    fprintf(outfile, " MINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // We'll only be working with the active molecule.
    shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!\n   Mints wrapper is not meant to be run with IPV1 inputs.");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    // Make sure molecule is valid.
    molecule->update_geometry();

    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basisset = BasisSet::construct(parser, molecule, options.get_str("BASIS"));

    // Create integral factory.
    shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset, basisset, basisset, basisset));

    // Get the SO basis object.
    shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basisset, integral));

    // Get ERI object
    shared_ptr<TwoBodyAOInt> tb(integral->eri());
    shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, integral));

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:                %4d\n", molecule->natom());
    fprintf(outfile, "      Number of AO shells:            %4d\n", basisset->nshell());
    fprintf(outfile, "      Number of SO shells:            %4d\n", sobasis->nshell());
    fprintf(outfile, "      Number of primitives:           %4d\n", basisset->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals:      %4d\n", basisset->nao());
    fprintf(outfile, "      Number of basis functions:      %4d\n\n", basisset->nbf());
    fprintf(outfile, "      Number of irreps:               %4d\n", sobasis->nirrep());
    fprintf(outfile, "      Number of functions per irrep: [");
    for (int i=0; i<sobasis->nirrep(); ++i) {
        fprintf(outfile, "%4d ", sobasis->nfunction_in_irrep(i));
    }
    fprintf(outfile, "]\n\n");

    // Compute and dump one-electron SO integrals.
    // Obtain dimensions from the sobasis
    const Dimension dimension = sobasis->dimension();

    // Create a matrix factory and initialize it
    shared_ptr<MatrixFactory> factory(new MatrixFactory());
    factory->init_with(dimension, dimension);

    // Overlap
    shared_ptr<OneBodySOInt> overlap(integral->so_overlap());
    shared_ptr<Matrix>       overlap_mat(factory->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio, PSIF_OEI);

    // Kinetic
    shared_ptr<OneBodySOInt> kinetic(integral->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(factory->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(psio, PSIF_OEI);

    // Potential
    shared_ptr<OneBodySOInt> potential(integral->so_potential());
    shared_ptr<Matrix>       potential_mat(factory->create_matrix(PSIF_SO_V));
    potential->compute(potential_mat);
    potential_mat->save(psio, PSIF_OEI);

    // Open the IWL buffer where we will store the integrals.
    IWL ERIOUT(psio.get(), PSIF_SO_TEI, 0.0, 0, 0);
    IWLWriter writer(ERIOUT);

    // How many SO shells do we have?
    int nsoshell = sobasis->nshell();

    // Let the user know what we're doing.
    fprintf(outfile, "      Computing integrals..."); fflush(outfile);

    // Begin new iterator code...works for all centers being equal
    int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];

    for (int usii=0; usii<nsoshell; ++usii) {
        for (int usjj=0; usjj<=usii; ++usjj) {
            for (int uskk=0; uskk<=usjj; ++uskk) {
                for (int usll=0; usll<=uskk; ++usll) {
                    int num_unique_pk = determine_unique_shell_quartets(usii, usjj, uskk, usll,
                                                                        usi_arr, usj_arr, usk_arr, usl_arr);

                    // For each num_unique_pk we need to call TwoBodySOInt::compute
                    for (int upk=0; upk<num_unique_pk; ++upk) {
                        eri->compute_shell(usi_arr[upk],
                                           usj_arr[upk],
                                           usk_arr[upk],
                                           usl_arr[upk],
                                           writer);
                    }

                    // If we are making a PK-matrix the end of this loop marks
                    // the end of a PK block.
                }
            }
        }
    }

    // Flush out buffers.
    ERIOUT.flush(1);

    // We just did all this work to create the file, let's keep it around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    fprintf(outfile, "done\n\n"); fflush(outfile);

    fprintf(outfile, "      Computed %lu integrals.\n\n", writer.count());

    tstop();

    return Success;
}

static int determine_unique_shell_quartets(int usii, int usjj, int uskk, int usll,
                                           int* usi_arr,
                                           int* usj_arr,
                                           int* usk_arr,
                                           int* usl_arr)
{
    //===--- Decide what shell quarters out of (ij|kl), (ik|jl), and (il|jk) are unique ---===
    usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;

    if (usii == usjj && usii == uskk || usjj == uskk && usjj == usll)
        return 1;
    else if (usii == uskk || usjj == usll) {
        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
        return 2;
    }
    else if (usjj == uskk) {
        usi_arr[1] = usii; usj_arr[1] = usll; usk_arr[1] = usjj; usl_arr[1] = uskk;
        return 2;
    }
    else if (usii == usjj || uskk == usll) {
        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
        return 2;
    }
    else {
        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
        usi_arr[2] = usii; usj_arr[2] = usll; usk_arr[2] = usjj; usl_arr[2] = uskk;
        return 3;
    }
}

}}
