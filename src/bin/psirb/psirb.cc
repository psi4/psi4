/*! \defgroup psirb psirb: PSI3 Ruby Input Driver */

/*! 
** \file
** \ingroup psirb
** \brief PSI3 Ruby driver module
*/

#include <getopt.h>
#include <cstdio>
#include <ruby.h>
#define MAIN
#include "psirb.h"

namespace psi { namespace psirb {
// Functions defined in ruby.c that are only needed here
extern bool initialize_ruby();
extern void load_input_file_into_ruby();
extern void process_input_file();
extern void finalize_ruby();
extern void run_interactive_ruby();

// Functions defined here.
void parse_command_line(int, char**);
void redirect_output(char *szFilename, bool append=false);
void print_version();
void print_usage();
}}

// The following are completely ignored by psirb and psi
extern "C" {
	const char *gprgid() {
		const char *prog = "psirb";
		return prog;
	}
	char *psi_file_prefix;
};

int main(int argc, char *argv[])
{
	using namespace psi::psirb;
	
	// Parse the command line arguments
	parse_command_line(argc, argv);
	
	// Print the version of Psi first
	print_version();
	
	// Initialize the interpreter
	initialize_ruby();
	
	// At this point the following has happened:
	//  1. Output is being directed to where the user wants
	//  2. Input filename has been determined
	//  3. Ruby has been initialized
	//  4. So far, some basic commands have been registered with Ruby
	
	// Now we need to branch off. If the user did not specify -I or --irb then process the
	// input file.
	if (Globals::g_bIRB == false) {
		// Now have Ruby load in the input file and begin processing
		load_input_file_into_ruby();
		
		// Actually exist the input file.
		process_input_file();
	}
	// User requested irb (interactive Ruby) this requires other procedures.
	else {
		run_interactive_ruby();
	}
	
	// Close the interpreter
	finalize_ruby();
	
	return EXIT_SUCCESS;
}

namespace psi { namespace psirb {
/*! Handles the command line arguments and stores them in global variables. In some cases
	default globals variables are set.
	\param argc Number of command-line arguments received.
	\param argv Value of the command-line arguments.
*/
void parse_command_line(int argc, char *argv[])
{
	int next_option;
	char *output_filename = NULL;
	
	// Set the default verbosity value
	Globals::g_bVerbose = false;
	
	// A string listing of valid short option letters
	const char* const short_options = "IhvVo:";
	const struct option long_options[] = {
		{ "irb",     0, NULL, 'I' },
		{ "help",    0, NULL, 'h' },
		{ "verbose", 0, NULL, 'v' },
		{ "version", 0, NULL, 'V' },
		{ "output",  1, NULL, 'o' },
		{ NULL,      0, NULL,  0  }
	};
	
	// Check the command line arguments
	do {
		next_option = getopt_long(argc, argv, short_options, long_options, NULL);
		
		switch (next_option) {			
			case 'h': // -h or --help
			print_usage();
			break;
			
			case 'I': // -I or --irb
			Globals::g_bIRB = true;
			break;
			
			case 'o': // -o or --output
			output_filename = optarg;
			break;
			
			case 'v': // -v or --verbose
			Globals::g_bVerbose = true;
			break;
			
			case 'V': // -V or --version
			print_version();
			break;
			
			case -1: // done with options
			break;
			
			default: // Something else, unexpected
			fprintf(stderr, "Unknown command line option encountered: %c\n", next_option);
			exit(EXIT_FAILURE);
			break;
		}
	} while (next_option != -1);
	
	// If the output needs to be redirected do so now.
	redirect_output(output_filename);
	if (output_filename)
		Globals::g_szOutputFilename = output_filename;
	
	// Check the input file name now
	// Default to input.rb if not given
	if (optind == argc) {
		Globals::g_szInputFile = "input.rb";
	}
	else {
		Globals::g_szInputFile = argv[optind];
	}
}

/*! Redirects the screen output from the input file and this module to the file given 
	in the argument otherwise it is sent to the screen. This doesn't quite work right 
	in practice.
	\param szFilename File to redirect to.
	\param append Append to it?
*/
void redirect_output(char *szFilename, bool append)
{
	const char *szAppend;
	
	if (append == false)
		szAppend = "w+";
	else
		szAppend = "w";
		
	if (szFilename) {
		// Create the file, truncate if necessary
		Globals::g_fOutput = fopen(szFilename, szAppend);
		
		// Error?
		if (Globals::g_fOutput == NULL) {
			fprintf(stderr, "Unable to open: %s\n", szFilename);
			exit(EXIT_FAILURE);
		}
	}
	else {
		// Redirect to the screen
		Globals::g_fOutput = stdout;
	}
}

/*! Print PSI version information that is was set in configure.ac */
void print_version()
{
	fprintf(Globals::g_fOutput, "PSI %d.%d Ruby Driver\n", PSI_VERSION_MAJOR, PSI_VERSION_MINOR);
}

/*! Print command-line usage information. */
void print_usage()
{
	printf("Usage:  psirb [options] inputfile\n");
	printf(" -h  --help               Display this usage information.\n");
	printf(" -v  --verbose            Print alot of information.\n");
	printf(" -V  --version            Print version information.\n");
	printf(" -o  --output filename    Redirect output elsewhere.\n");
	printf(" -I  --irb                Run psirb in interactive mode.\n");
	
	exit(EXIT_FAILURE);
}

}} // namespace psi::psirb
