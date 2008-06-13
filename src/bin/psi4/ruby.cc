/*! \file ruby.cc
  \ingroup (psi4)
  \brief Contains Ruby initialization and some global functions.
*/

#include <ruby.h>
#include "psi4.h"

namespace psi { namespace psi4 {
  
  // Function declarations
  bool initialize_ruby();
  void load_input_file_into_ruby();
  void process_input_file();
  void finalize_ruby();
  int run_interactive_ruby();
  
  // Defined elsewhere
  extern void print_version();
  
  /*! Initializes the Ruby interpreter, sets the Ruby $: and $0 variables,
      adds the appropriate PSIDATADIR path to Ruby's search path,
      and also creates the Ruby-ized Psi objects.
      \return true on success, false on failure.
  */
  bool initialize_ruby()
  {
    char *psishare_dirname;
    
    // Setup the initialize the interpreter
    ruby_init();
    
    // Initialize the $: (load path) variable; necessary if the input file
    // loads any other modules
    ruby_init_loadpath();
    
    // Need to modify the load path to include Psi's specific locations
    //   Get the global variable $:
    VALUE loadpath = rb_gv_get("$:");
    //   Pop off the "." directory
    VALUE popped_value = rb_ary_pop(loadpath);
    
    // Check to see what we need to add
    psishare_dirname = getenv("PSIDATADIR");
    if (psishare_dirname != NULL) {
      char *tmpstr = (char*)malloc(sizeof(char)*(strlen(psishare_dirname)+6));
      sprintf(tmpstr, "%s/ruby", psishare_dirname);
      psishare_dirname = tmpstr;
    }
    else {
      psishare_dirname = strdup(INSTALLEDPSIDATADIR "/ruby");
    }
    // Convert the char* to a Ruby string
    VALUE newpath = rb_str_new2(psishare_dirname);
    free(psishare_dirname);
    
    //  push on new location
    rb_ary_push(loadpath, newpath);
    //  push the "." location back on
    rb_ary_push(loadpath, popped_value);
    
    // Set the name of the Ruby script (and $0) to Psi
    ruby_script("PSI");
    
    // Add some basic functionality to Ruby
    //      TODO: Need to work out the best way to do this.
    
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
  
  /*! This is used to make the call to rb_protect easier. */
  VALUE load_modules_for_irb(VALUE)
  {
    rb_require("psi4");
    rb_require("irb");
    return Qnil;
  }
  
  /*! Handles running Ruby interactively. This version makes use of the
      irb module found in Ruby to handle Ruby input parsing. */
  int run_interactive_ruby()
  {
    int error;
    
    printf("Starting interactive PSI4 driver:\n");
    // rb_protect is used to protect PSI from Ruby code exceptions
    // without it PSI would crash if a script did something nasty.
    rb_protect(load_modules_for_irb, 0, &error);
    if (error) {
      // If error is nonzero then an exception was thrown in the Ruby code.
      // Get the message from the exception
      VALUE err = rb_inspect(ruby_errinfo);
      // Print a backtrace to the screen, hopefully the user can make sense of it.
      rb_backtrace();
      // Print what the error message was.
      printf("ERROR %s\n", StringValuePtr(err));
      return -1;
    }
    
    // Run the Irb. This should probably be wrapped in a rb_protect as well.
    rb_eval_string("IRB.start");
    
    return true;
  }
}}

