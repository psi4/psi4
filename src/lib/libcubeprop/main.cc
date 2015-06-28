#include "psi4-dec.h"
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <libplugin/plugin.h>
#include <libqt/qt.h>
#include <boost/shared_ptr.hpp>
#include "prop.h"

INIT_PLUGIN

namespace psi{ 

extern "C" int read_options(std::string name, Options &options)
{
    if (name == "PROP" || options.read_globals()) {
        // => Printing Metadata <= //

        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- Debug level -*/ 
        options.add_int("DEBUG", 0);

        // => CubicScalarGrid options <= //

        /*- CubicScalarGrid basis cutoff. !expert -*/
        options.add_double("CUBIC_BASIS_TOLERANCE", 1.0E-12);
        /*- CubicScalarGrid maximum number of grid points per evaluation block. !expert -*/
        options.add_int("CUBIC_BLOCK_MAX_POINTS",1000);
        /*- CubicScalarGrid overages in bohr [O_X, O_Y, O_Z]. Defaults to 2.0 bohr each. -*/ 
        options.add("CUBIC_GRID_OVERAGE", new ArrayType());
        /*- CubicScalarGrid spacing in bohr [D_X, D_Y, D_Z]. Defaults to 0.2 bohr each. -*/ 
        options.add("CUBIC_GRID_SPACING", new ArrayType());
    
        // => Properties tasks <= //

        /*- Grid property data filepath -*/ 
        options.add_str_i("PROPERTY_FILEPATH", ".");
        /*-
        Properties to compute. Valid tasks include:
            DENSITY - Da, Db, Dt, Ds
            ESP - Dt, ESP 
            ORBITALS - Psi_a_N, Psi_b_N 
            BASIS_FUNCTIONS - Phi_N
            LOL - LOLa, LOLb
            ELF - ELFa, ELFb
        -*/
        options.add("PROPERTY_TASKS", new ArrayType());
        /*- List of desired orbital indices (1-based, + for alpha, - for beta). All orbitals computed if empty.-*/
        options.add("PROPERTY_ORBITALS", new ArrayType());
        /*- List of desired basis function indices (1-based). All basis functions computed if empty.-*/
        options.add("PROPERTY_BASIS_FUNCTIONS", new ArrayType());
    }

    return true;
}

extern "C" PsiReturnType plugin_prop(Options &options)
{
    tstart();

    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<Properties> prop(new Properties(wfn));
    prop->compute_properties();

    tstop();

    return Success;
}

} // end namespaces
