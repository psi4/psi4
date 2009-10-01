#include <libpsio/psio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

#include "psi4-def.h"
#include "libmints/mints.h"

using namespace psi;

namespace psi{
    FILE *infile;
    int read_options(std::string name, Options &options)
    {
        ip_cwk_clear();
        ip_cwk_add(":BASIS");
        ip_cwk_add(":DEFAULT");
        ip_cwk_add(":PSI");
        ip_set_uppercase(1);
        options.clear();
        if(name == "DF-MP2") {
            ip_cwk_add(":DF-MP2");
            /*- The amount of information printed
                to the output file -*/
            options.add_int("PRINT", 1);
            /*- Whether to compute the SCS energy -*/
            options.add_bool("DO_SCS", true);
            /*- Whether to compute the SCS-N energy -*/
            options.add_bool("DO_SCS-N", true);
            /*- The name of the orbital basis set -*/
            options.add_str("BASIS", "");
            /*- The name of the auxilliary basis set -*/
            options.add_str("RI_BASIS", "");
            /*- The opposite-spin scale factor for the SCS energy -*/
            options.add_double("SCALE_OS", 6.0/5.0);
            /*- The same-spin scale factor for the SCS energy -*/
            options.add_double("SCALE_SS", 1.0/3.0);
        }
        options.read_ipv1();
    }
    namespace dfmp2{
        PsiReturnType dfmp2(Options &options, int argc, char **argv);
        
        void title()
        {
            fprintf(outfile, "\t\t\t*************************\n");
            fprintf(outfile, "\t\t\t*                       *\n");
            fprintf(outfile, "\t\t\t*         DF-MP2        *\n");
            fprintf(outfile, "\t\t\t*                       *\n");
            fprintf(outfile, "\t\t\t*************************\n");
            fflush(outfile);
        }
    } // End Namespace df-mp2
} // End namespace psi


int
main(int argc, char *argv[])
{
    int num_unparsed, i;
    char *argv_unparsed[100];

    for (i=1, num_unparsed=0; i<argc; ++i)
        argv_unparsed[num_unparsed++] = argv[i];

    psi_start(&infile, &outfile, &psi_file_prefix,
                num_unparsed, argv_unparsed, 0);
    psio_init();
    psio_ipv1_config();
    Options options;

    psi::dfmp2::title();
    read_options("DF-MP2", options);
    psi::dfmp2::dfmp2(options, argc, argv);

    psi_stop(infile, outfile, psi_file_prefix);
    return (EXIT_SUCCESS);
}