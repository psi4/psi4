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

/*!
** \file
** \brief Initialize input, output, file prefix, etc.
** \ingroup
*/

#include "psi4.h"
#include <psi4-dec.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <boost/filesystem.hpp>
#include <getopt.h>
#include <psifiles.h>
#include <psiconfig.h>
#include <libplugin/plugin.h>
#include <libparallel/parallel.h>
#include "libparallel/ParallelPrinter.h"
using namespace std;

namespace psi {

void create_new_plugin(std::string plugin_name, const std::string& template_name);
void print_version(std::string);
void print_usage();

/*!
** psi_start():
** This function initializes the input, output files, file prefix, etc.,
** by checking command line arguments and environmental variables. It also
** initializes the input parsing library.
**
** \param argc       = number of command-line arguments passed
** \param argv       = command-line arguments
**
** Returns: one of standard PSI error codes
** \ingroup PSI4
*/
int psi_start(int argc, char *argv[])
{
    int next_option;
    // Defaults:
    std::string ifname;
    std::string ofname;
    std::string fprefix;
    std::string templatename;
    std::string pluginname;
#if defined(PSIDEBUG)
    bool debug          = true;
#else
    bool debug          = false;
#endif
    bool append         = false;
    check_only          = false;
    clean_only          = false;
    verbose             = false;
    skip_input_preprocess = false;
    interactive_python  = false;

    // A string listing of valid short option letters
    const char* const short_options = "ahvVdcwo:p:i:l:s:n:r:mkt";
    const struct option long_options[] = {
        { "append",  0, NULL, 'a' },
        { "help",    0, NULL, 'h' },
        { "verbose", 0, NULL, 'v' },
        { "version", 0, NULL, 'V' },
        { "check",   0, NULL, 'c' },
        { "nthread", 1, NULL, 'n' },
        { "debug",   0, NULL, 'd' },
        { "wipe",    0, NULL, 'w' },
        { "output",  1, NULL, 'o' },
        { "prefix",  1, NULL, 'p' },
        { "input",   1, NULL, 'i' },
        { "psidatadir", 1, NULL, 'l' },
        { "scratch", 1, NULL, 's' },
        { "restart", 1, NULL, 'r' },
        { "messy",   0, NULL, 'm' },
        { "skip-preprocessor",  0, NULL, 'k' },
        { "interactive", 0, NULL, 't' },
        { "new-plugin", 1, NULL, 1 },  // requires 1 argument
        { NULL,      0, NULL,  0  }
    };

    // Check the command line arguments
    if (argc) {
        do {
            next_option = getopt_long(argc, argv, short_options, long_options, NULL);

            switch (next_option) {
            case 1: // --new-plugin
                // Let's check the next arugument to see if the user specified a +templatename
                pluginname = optarg;
                if (optind < argc) {
                    // check to see if the next arguments starts with + and that the user gave more
                    if (argv[optind][0] == '+' && strlen(argv[optind]) > 1) {
                        // yes, take it, skipping the +
                        templatename = &(argv[optind][1]);
                        optind++;
                    }
                }
                create_new_plugin(pluginname, templatename);
                exit(EXIT_SUCCESS);
                break;

            case 'a': // -a or --append
                append = true;
                break;

            case 'n': // -n or --nthread
                Process::environment.set_n_threads(atoi(optarg));
                break;

            case 'd': // -d or --debug
                debug = true;
                break;

            case 'h': // -h or --help
                print_usage();
                exit(EXIT_SUCCESS);
                break;

            case 'i': // -i or --input
                ifname = optarg;
                break;

            case 'l': // -l or --psidatadir
                // This is ok so long as all reads of this var are after here
                Process::environment.set("PSIDATADIR",optarg);
                break;

            case 's': // -s or --scratch
                Process::environment.set("PSI_SCRATCH",optarg);
                break;

            case 'o': // -o or --output
                ofname = optarg;
                break;

            case 'p': // -p or --prefix
                fprefix = optarg;
                break;

            case 'r': // -r or --restart
                restart_id = optarg;
                break;

            case 'v': // -v or --verbose
                verbose = true;
                break;

            case 'm': // -m or --messy
                messy = true;
                break;

            case 'V': // -V or --version
                print_version("stdout");
                exit(EXIT_SUCCESS);
                break;

                /*
        case 'c': // -c or --check
            check_only = true;
            outfile = stdout;
            append = true;
            break;
*/

            case 'w': // -w or --wipe
                clean_only = true;
                append = true; // to avoid deletion of an existing output file
                break;

            case 'k': // -y or --skip-parser
                skip_input_preprocess = true;
                break;

            case 't': // -t or --interactive
                interactive_python = true;
                break;

            case -1: // done with options
                break;

            default: // Something else, unexpected
                printf("Unexpected argument given: -%c\n", next_option);
                exit(EXIT_FAILURE);
                break;
            }
        } while (next_option != -1);

        if (optind <= argc-2) {
            ifname = argv[optind];
            ofname = argv[optind+1];
        }
        else if (optind <= argc-1) {
            ifname = argv[optind];
        }
    }

    if (interactive_python)
        ofname = "stdout";

    /* if some args were not specified on command-line - check the environment */
    if (ifname.empty() && Process::environment("PSI_INPUT").size())
        ifname = Process::environment("PSI_INPUT");
    if (ofname.empty() && Process::environment("PSI_OUTPUT").size())
        ofname = Process::environment("PSI_OUTPUT");
    if (fprefix.empty() && Process::environment("PSI_PREFIX").size())
        fprefix = Process::environment("PSI_PREFIX");
    if (restart_id.empty() && Process::environment("PSI_RESTART").size())
        restart_id = Process::environment("PSI_RESTART");

    /* if some arguments still not defined - assign default values */
    if (ifname.empty()) ifname = "input.dat";

    // if infile name is input.dat, make outfile name output.dat.  Otherwise,
    // strip off any .in or .dat at the end of the input file, and use the
    // base_name for the output file, appending .out. ---CDS 4/13
    if (ofname.empty()) {
        if (ifname == "input.dat") ofname = "output.dat";
        else {
          std::string base_name;
          if (ifname.size() > 4 && (ifname.rfind(".dat")==ifname.size()-4)){
              base_name = ifname.substr(0, ifname.rfind(".dat"));
          }
          else if (ifname.size()>3 && (ifname.rfind(".in")==ifname.size()-3)){
              base_name = ifname.substr(0, ifname.rfind(".in"));
          }
          else {
              base_name = ifname;
          }
         ofname = base_name + ".out";
        }
    }

    /* default prefix is not assigned here yet because need to check input file first */

    /* open input and output files */
#if !defined(MAKE_PYTHON_MODULE)
    infile = NULL;
    if (!interactive_python) {
        if(ifname == "stdin") {
            infile=stdin;
            infile_directory = ".";
        }
        else {
            infile = fopen(ifname.c_str(), "r");

            boost::filesystem::path ipath = boost::filesystem::system_complete(ifname);
            infile_directory = ipath.parent_path().string();
        }

        if(infile == NULL) {
            outfile->Printf( "Error: could not open input file %s\n",ifname.c_str());
            return(PSI_RETURN_FAILURE);
        }
    }
#endif

#if defined(MAKE_PYTHON_MODULE)
    outfile = stdout;
#else
        if(ofname == "stdout"){
            outfile=boost::shared_ptr<PsiOutStream>(new PsiOutStream());
        }
        else{
           outfile=boost::shared_ptr<PsiOutStream>
              (new OutFile(ofname,(append?APPEND:TRUNCATE)));
        }
#endif

    //if(debug)
    //    setbuf(outfile,NULL);

    // Initialize Yeti's env
    //yetiEnv.init(WorldComm->me(), ofname.c_str());
    // This seems a bit daft, but it's necessary to make sure the compiler doesn't
    // nuke it in the mistaken belief that it's unused; plugins may use it.
    //yeti::Env::out0() << "";

    /* if prefix still NULL - check input file */
    if (fprefix.empty()) {
        fprefix = PSI_DEFAULT_FILE_PREFIX;
    }

    /* copy over file prefix, etc. into their appropriate variables */
    psi_file_prefix = strdup(fprefix.c_str());

    // If check_only, force output to stdout because we don't need anything more
    //if (check_only)
    //    outfile = stdout;
    outfile_name = ofname;

    return(PSI_RETURN_SUCCESS);
}

/*! Print command-line usage information. */
void print_usage(void)
{
    printf("Usage:  psi4 [options] inputfile\n");
    printf(" -a  --append             Append results to output file. Default: Truncate first\n");
    printf(" -c  --check              Run input checks (not implemented).\n");
    printf(" -d  --debug              Flush the outfile at every psi::fprintf.\n"
           "                          Default: true iff --with-debug.\n");
    printf(" -h  --help               Display this usage information.\n");
    printf(" -i  --input filename     Input file name. Default: filename\n");
    printf(" -k  --skip-preprocessor  Skips input preprocessing. Expert mode.\n");
    printf(" -o  --output filename    Redirect output elsewhere. Default:\n");
    printf("                          filename.out if input is filename\n");
    printf("                          filename.out if input is filename.in\n");
    printf("                          output.dat if input is input.dat\n");
    printf(" -l  --psidatadir         Specify where to look for the Psi data directory.\n");
    printf("                          Overrides PSIDATADIR.\n");
    printf(" -m  --messy              Leave temporary files after the run is completed.\n");
    printf(" -r  --restart            Number to be used instead of process id.\n");
    printf(" -n  --nthread            Number of threads to use (overrides OMP_NUM_THREADS)\n");
    printf("     --new-plugin name    Creates a new directory with files for writing a\n"
           "                          new plugin. You can specify an additional argument\n"
           "                          that specifies a template to use, for example:\n"
           "                              --new-plugin name +mointegrals\n");
    printf(" -p  --prefix prefix      Prefix name for psi files. Default: psi\n");
    printf(" -t  --interactive        Start interactive Python session. Expert mode.\n");
    printf(" -v  --verbose            Print a lot of information.\n");
    printf(" -V  --version            Print version information.\n");
    printf(" -w  --wipe               Clean out your scratch area.\n");
    printf("\n");
    printf("Environment Variables\n");
    printf("     PSI_SCRATCH          Directory where scratch files are written. Default: /tmp/\n");
    printf("                          This should be a local, not network, disk\n");

    exit(EXIT_FAILURE);
}

}

