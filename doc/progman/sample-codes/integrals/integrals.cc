#include "psi4-dec.h"
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>

namespace psi{ namespace integrals{

PsiReturnType
integrals(Options &options, int argc, char *argv[])
{
    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");
    // This will print out all of the user-provided options for this module
    options.print();
    // Set up some essential arrays in the Wavefunction object
    Wavefunction::initialize_singletons();
    // Make a (reference-counted) psio object for I/O operations and set the parser up with it
    shared_ptr<PSIO> psio(new PSIO);
    psiopp_ipv1_config(psio);
    // Open the existing checkpoint object by creating a Chkpt object
    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));
    // The matrix factory can create matrices of the correct dimensions...
    shared_ptr<MatrixFactory> factory(new MatrixFactory);
    // ...it needs to raid the checkpoint file to find those dimensions.
    factory->init_with_chkpt(chkpt);
    // The basis set is also created from the information stored in the checkpoint file
    shared_ptr<BasisSet> basis(new BasisSet(chkpt));
    // The integral factory oversees the creation of integral objects
    shared_ptr<IntegralFactory> integral(new IntegralFactory
            (basis, basis, basis, basis));

    // Form the one-electron integral objects from the integral factory
    shared_ptr<OneBodyAOInt> sOBI(integral->overlap());
    shared_ptr<OneBodyAOInt> tOBI(integral->kinetic());
    shared_ptr<OneBodyAOInt> vOBI(integral->potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));
    
    // Compute the one electron integrals, telling each object where to store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);

    if(print > 5){
        sMat->print();
    }
    if(print > 3){
        tMat->print();
        vMat->print();
    }
    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);
    hMat->print();

    if(doTei){
        // Here's an example of an exception - we throw it if there are more than 99 basis
        // functions, because only 2 digits are used in the print function below
        if(chkpt->rd_nso() > 99) 
            throw PsiException("This code can only handle fewer than 100 basis functions", __FILE__, __LINE__);
        // Now, the two-electron integrals
        shared_ptr<TwoBodyAOInt> eri(integral->eri());
        // The buffer will hold the integrals for each shell, as they're computed
        const double *buffer = eri->buffer();
        // The iterator conveniently lets us iterate over functions within shells
        AOShellCombinationsIterator shellIter = integral->shells_iterator();
        int count=0;
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
            // Compute quartet
            eri->compute_shell(shellIter);
            // From the quartet get all the integrals
            AOIntegralsIterator intIter = shellIter.integrals_iterator();
            for (intIter.first(); intIter.is_done() == false; intIter.next()) {
                int p = intIter.i();
                int q = intIter.j();
                int r = intIter.k();
                int s = intIter.l();
                fprintf(outfile, "\t(%2d %2d | %2d %2d) = %20.15f\n",
                    p, q, r, s, buffer[intIter.index()]);
                ++count;
            }
        }
        fprintf(outfile, "\n\tThere are %d unique integrals\n\n", count);
    }

    return Success;   
}

}} // End Namespaces
