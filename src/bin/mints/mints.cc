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

PsiReturnType mints(Options & options)
{
    tstart();

    shared_ptr<PSIO> psio(new PSIO);
    psiopp_ipv1_config(psio);           // Do we really need this?

    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    fprintf(outfile, " MINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // We'll only be working with the active molecule.
    shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!\n   Mints wrapper is not meant to be run with IPV1 inputs.");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basisset = BasisSet::construct(parser, molecule, options.get_str("BASIS"));

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:           %4d\n", molecule->natom());
    fprintf(outfile, "      Number of shells:          %4d\n", basisset->nshell());
    fprintf(outfile, "      Number of primitives:      %4d\n", basisset->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals: %4d\n", basisset->nao());
    fprintf(outfile, "      Number of basis functions: %4d\n\n", basisset->nbf());

    // Create integral factory.
    IntegralFactory integral(basisset, basisset, basisset, basisset);

    // Get ERI object
    shared_ptr<TwoBodyInt> eri(integral.eri());

    // Buffer where the integrals will be
    const double *buffer = eri->buffer();

    // Open the IWL buffer where we will store the integrals.
    IWL ERIOUT(psio.get(), PSIF_SO_TEI, 0.0, 0, 0);

    fprintf(outfile, "      Computing integrals..."); fflush(outfile);

    // Handy variables
    int P, Q, R, S;
    int i, j, k, l;
    int index;
    double value;
    size_t count=0;

    ShellCombinationsIterator iter = integral.shells_iterator();
    for (iter.first(); !iter.is_done(); iter.next()) {
        P = iter.p();
        Q = iter.q();
        R = iter.r();
        S = iter.s();

        // Compute quartet
        eri->compute_shell(P, Q, R, S);

        // From the quartet get all the integrals
        IntegralsIterator int_iter = iter.integrals_iterator();
        for (int_iter.first(); !int_iter.is_done(); int_iter.next()) {
            i = int_iter.i();
            j = int_iter.j();
            k = int_iter.k();
            l = int_iter.l();
            index = int_iter.index();
            value = buffer[index];

            if (fabs(value) > 1.0e-14) {
                // Note: I do not handle the -i which we always had to test for with cints.
                // Write it to the IWL buffer.
                ERIOUT.write_value(i, j, k, l, value, 0, NULL, 0);
                count++;
            }
        }
    }

    // We just did all this work to create the file, let's keep it around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    fprintf(outfile, "done\n\n"); fflush(outfile);

    fprintf(outfile, "      Computed %lu integrals.\n\n", count);

    tstop();

    return Success;
}

}}
