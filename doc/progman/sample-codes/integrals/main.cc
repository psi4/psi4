#include <libpsio/psio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

#include "psi4-def.h"

using namespace psi;

namespace psi {
    FILE *infile;

    namespace integrals{
        PsiReturnType integrals(Options &options, int argc, char **argv);
    }

    int 
    read_options(std::string name, Options &options)
    {
        ip_cwk_clear();
        ip_cwk_add(":BASIS");
        ip_cwk_add(":DEFAULT");
        ip_cwk_add(":PSI");
        ip_set_uppercase(1);
        options.clear();
        if(name == "INTEGRALS") {
            ip_cwk_add(":INTEGRALS");
            /*- The amount of information printed
                to the output file -*/
            options.add_int("PRINT", 1);
            /*- Whether to compute two-electron integrals -*/
            options.add_bool("DO_TEI", false);
        }
        options.read_ipv1();
        
        return true;
    }
}


int
main(int argc, char *argv[])
{
    int num_unparsed, i;
    char *argv_unparsed[100];

    for (i=1, num_unparsed=0; i<argc; ++i)
        argv_unparsed[num_unparsed++] = argv[i];

    psi_start(&infile,&outfile,&psi_file_prefix,num_unparsed, argv_unparsed, 0);
    psio_init();
    psio_ipv1_config();
    Options options;

    read_options("INTEGRALS", options);
    psi::integrals::integrals(options, argc, argv);

    psi_stop(infile, outfile, psi_file_prefix);
    return (EXIT_SUCCESS);
}
