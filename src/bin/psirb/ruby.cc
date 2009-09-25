/*! \file
	\ingroup psirb
	\brief Contains Ruby initialization and some global functions.
*/
#include <ruby.h>
#include "psirb.h"

namespace psi { namespace psirb {
	
// Function declarations
bool initialize_ruby();
void load_input_file_into_ruby();
void process_input_file();
void finalize_ruby();
VALUE ruby_psi_print_version();
VALUE ruby_psi_get_version_major();
VALUE ruby_psi_get_version_minor();
bool create_ruby_psi_module();

// Defined elsewhere
extern void print_version();

/*! Initializes the Ruby interpreter, sets the Ruby $: and $0 variables, 
	adds the appropriate PSIDATADIR path to Ruby's search path,
    and also creates the Ruby-ized Psi objects. 
	\return true on success, false on failure
*/
bool initialize_ruby()
{
	char *psishare_dirname;
	
	// Setup and initialize interpreter
	ruby_init();
	
	// Initialize the $: (load path) variable; necessary if the input file loads any other modules
	ruby_init_loadpath();
	
	// Need to modify the load path of Ruby to include Psi's specific locations
	//  Get the global variable $:
	VALUE loadpath = rb_gv_get("$:");
	//  pop off the "." directory
	VALUE popped_value = rb_ary_pop(loadpath);
	
	// check to see what we need to add
	psishare_dirname = getenv("PSIDATADIR");
	if (psishare_dirname != NULL) {
		char *tmpstr = (char*)malloc(sizeof(char)*(strlen(psishare_dirname)+6));
		sprintf(tmpstr, "%s/ruby", psishare_dirname);
		psishare_dirname = tmpstr;
	}
	else {
		psishare_dirname = strdup(INSTALLEDPSIDATADIR "/ruby");
	}
	// Convert char* to Ruby string
	VALUE newpath = rb_str_new2(psishare_dirname);
	free(psishare_dirname);
	
	//  push on new location
	rb_ary_push(loadpath, newpath);
	//  push the "." location back on
	rb_ary_push(loadpath, popped_value);
	
	// Set the name of the Ruby script (and $0) to Psi
	ruby_script("Psi");
	
	// Add some basic functionality to Ruby
	create_ruby_psi_module();
	Task::create_ruby_class();
	ZEntry::create_ruby_class();
	Matrix::create_ruby_class();
	
	// Done
	return true;
}

/*! Loads the Ruby input file into the interpreter. Does not perform syntax checking, nor does
    it begin executing the code. Just checks to see if the file exists and loads it. */
void load_input_file_into_ruby()
{
	// Have Ruby do it
	rb_load_file(Globals::g_szInputFile.c_str());
}

/*! Run the input file. */
void process_input_file()
{
	// Process the input file
	// Since I do not want Ruby to take complete control of the system this is a 
	// hack version from Ruby's eval.c file function ruby_eval and ruby_stop
	// Typically you would call ruby_run but this function does not return at 
	// all from ruby_stop, it invokes an exit function call internally.
	int state;
    static int ex;

    if (ruby_nerrs > 0) exit(EXIT_FAILURE);
    state = ruby_exec();
    if (state && !ex) ex = state;
	ruby_cleanup(ex);
}

/*! Shutdown the interpreter */
void finalize_ruby()
{
	ruby_finalize();
}

//! Ruby function: Psi::Version::print_version 
/*! Prints out version information from the C-function print_version. 
	\return Qnil, the Ruby equivalent to NULL.
*/
VALUE ruby_psi_print_version()
{
	print_version();
	return Qnil;
}

//! Ruby function: indent_puts and Psi::indent_puts
/*! Increments the number of spaces that will appear at the beginning of a line printed
    with ''puts'' by 2. 
	\return Qnil, the Ruby equivalent to NULL.
*/
VALUE ruby_psi_indent_puts()
{
	Globals::g_iPutsIndent += 2;
	return Qnil;
}

//! Ruby funcation: unindent_puts and Psi::unindent_puts
/*! Decrements the number of spaces that will appear at the beginning of a line printed
    with ''puts' by 2. 
	\return Qnil, the Ruby equivalent to NULL.
*/
VALUE ruby_psi_unindent_puts()
{
	Globals::g_iPutsIndent -= 2;
	if (Globals::g_iPutsIndent < 0)
		Globals::g_iPutsIndent = 0;
	return Qnil;
}

//! Ruby function: Psi::Version::get_version_major
/*! \return The major version value for Psi */
VALUE ruby_psi_get_version_major()
{
	return INT2FIX(PSI_VERSION_MAJOR);
}

//! Ruby function: Psi::Version::get_version_minor
/*! \return The minor version value for Psi */
VALUE ruby_psi_get_version_minor()
{
	return INT2FIX(PSI_VERSION_MINOR);
}

//! Ruby function: puts
/*! This is an override of Ruby's standard puts function. This allows for indentation
    and redirection of the puts.
	\param argc Number of Ruby objects to print.
	\param argv C-array of Ruby objects to print.
	\return Qnil, the Ruby equvalent to NULL.
*/
VALUE ruby_psi_puts(int argc, VALUE *argv)
{
	int i;
	VALUE str;
	
	// If we are supposed to be quiet don't print anything
	if (Globals::g_bQuietRuby == false) {
		for (i=0; i<Globals::g_iPutsIndent; i++)
			fprintf(Globals::g_fOutput, " ");

		for (i=0; i<argc; ++i) {
		// StringValue calls to_str on the object, if needed, to convert the object to a string.
			VALUE str = StringValue(argv[i]);
		// Print the converted objected.
			fprintf(Globals::g_fOutput, RSTRING(str)->ptr);
		}
		fprintf(Globals::g_fOutput, "\n");
		fflush(Globals::g_fOutput);
	}
	
	return Qnil;
}

//! Ruby function: Psi::quiet=
/*! If quiet=true then puts will not print anything when called */
VALUE ruby_psi_quiet_set(VALUE self, VALUE value)
{
	if (value == Qtrue)
		Globals::g_bQuietRuby = true;
	else
		Globals::g_bQuietRuby = false;
		
	return value;
}

//! Ruby function: Psi::quiet
/*! Returns the quiet value */
VALUE ruby_psi_quiet_get(VALUE self, VALUE value)
{
	if (Globals::g_bQuietRuby == true)
		return Qtrue;
	else
		return Qfalse;
}

//! Creates a module in the Ruby address space named Psi.
/*! All Psi objects/functions reside in this address space. Each function that is to be
	available to a user in the input file must be registered with Ruby.
	\returns true on success, false on failure.
*/
bool create_ruby_psi_module()
{
	// Handle for Version module
	VALUE rubyVersion;
	
	// Define a Ruby module named Psi
	Globals::g_rbPsi = rb_define_module("Psi");
	
	// Override some commands to ensure output is redirected to where we want
	// This redefines the global puts function to use our's
	rb_define_global_function("puts",          RUBYCAST(ruby_psi_puts), -1);
	rb_define_global_function("indent_puts",   RUBYCAST(ruby_psi_indent_puts), 0);
	rb_define_global_function("unindent_puts", RUBYCAST(ruby_psi_unindent_puts), 0);
	
	// Define methods that belong to Psi
	rb_define_module_function(Globals::g_rbPsi, "indent_puts",   RUBYCAST(ruby_psi_indent_puts), 0);
	rb_define_module_function(Globals::g_rbPsi, "unindent_puts", RUBYCAST(ruby_psi_unindent_puts), 0);
	rb_define_module_function(Globals::g_rbPsi, "quiet=",        RUBYCAST(ruby_psi_quiet_set), 1);
	rb_define_module_function(Globals::g_rbPsi, "quiet",         RUBYCAST(ruby_psi_quiet_get), 0);
	
	// Define a sub-module of Psi named Version
	rubyVersion = rb_define_module_under(Globals::g_rbPsi, "Version");
	
	// Add methods to the new module: Version
	rb_define_module_function(rubyVersion, "print_version", RUBYCAST(ruby_psi_print_version), 0);
	rb_define_module_function(rubyVersion, "get_major",     RUBYCAST(ruby_psi_get_version_major), 0);
	rb_define_module_function(rubyVersion, "get_minor",     RUBYCAST(ruby_psi_get_version_minor), 0);
		
	// Done
	return true;
}

}} // namespace psi::psirb
