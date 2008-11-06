/*! \file ruby.cc
  \ingroup (psi4)
  \brief Contains Ruby initialization and some global functions.
*/

#include <ruby.h>
#include "psi4.h"
#include "task.h"

namespace psi { //namespace psi4 {
  
  // Function declarations
  bool initialize_ruby();
  void load_input_file_into_ruby();
  void process_input_file();
  void finalize_ruby();
  int run_interactive_ruby();
  void handle_ruby_exception();
  bool create_global_task();
  void create_psi_module();
  
  // Defined elsewhere
  extern void print_version();
  
  /*! Initializes the Ruby interpreter, sets the Ruby $: and $0 variables,
      adds the appropriate PSIDATADIR path to Ruby's search path,
      and also creates the Ruby-ized Psi objects.
      \return true on success, false on failure.
  */
  bool initialize_ruby()
  {
    WHEREAMI();
    char *psishare_dirname;
    
    #if RUBY_MAJOR == 1 && RUBY_MINOR >= 9
    ruby_sysinit(0, NULL);
    RUBY_INIT_STACK;
    #endif
    
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
      char *tmpstr = NULL;
      asprintf(&tmpstr, "%s/ruby", psishare_dirname);
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
    
    // Set the name of the Ruby script (and $0) to PSI
    ruby_script("PSI");
    
    // Add some basic functionality to Ruby
    //      TODO: Need to work out the best way to do this.
    create_psi_module();
    Task::Ruby::create_ruby_class();
    
    // Done
    return true;
  }
  
  VALUE process_input_file_protected(VALUE)
  {
    WHEREAMI();
    
    // Run via the global task object.
    rb_funcall(g_rbTask, rb_intern("run"), 0);
    
    return Qnil;
  }
  
  /*! Loads the Ruby input file into the interpreter. Does not perform syntax checking, nor does
      it begin executing the code. Just checks to see if the file exists and loads it. */
  void load_input_file_into_ruby()
  {
    WHEREAMI();
  	// Have Ruby do it
  	#if RUBY_MAJOR == 1 && RUBY_MINOR < 9
  	// In pre-Ruby 1.9 the input file was loaded into a single global space.
  	rb_load_file(g_szInputFile.c_str());
  	#else
  	// However, in Ruby 1.9 it seems there is the possibility of loading multiple input files.
  	// So rb_load_file loads the input file and returns an execution node that we tell Ruby
  	// to execute when ready.
    g_rbExecNode = rb_load_file(g_szInputFile.c_str());
    #endif
  }
  
  /*! Run the input file. */
  void process_input_file()
  {
    WHEREAMI();
  	// Process the input file
    int error;
    
    // Call the protected version of this function. This catches exceptions (Ruby/C++)
    // and doesn't allow them to crash the program.
    rb_protect(process_input_file_protected, 0, &error);
    if (error) {
      // If error is nonzero then an exception was thrown in the code.
      handle_ruby_exception();
    }
  }
  
  VALUE create_global_task_protected(VALUE)
  {
    WHEREAMI();
    // If this fails it will throw a Ruby exception
    // Ruby is used to create a new Task instance, fully wrapped in Ruby.
    g_rbTask = rb_class_new_instance(0, 0, Task::Ruby::rbTask_);
    // Get the C++ object from Ruby to use.
    Data_Get_Struct(g_rbTask, Task, g_cTask);
    // Prevent Ruby from garbage collecting this new object when this
    // function returns by making it a global variable. I'm not sure
    // how you would access this global variable in Ruby.
    // Make sure the global task object know what input file to use.
    g_cTask->input_file(g_szInputFile);
    
    // This function doesn't need to return anything of value.
    // But ALL functions that Ruby calls MUST return something.
    return Qnil;
  }
  
  /*! Create the global Task object */
  bool create_global_task()
  {
    WHEREAMI();
    int error;
    
    // Call the protected version of this function.
    rb_protect(create_global_task_protected, 0, &error);
    if (error) {
      handle_ruby_exception();
      return false;
    }
    return true;
  }
  
  /*! Close out the global task object */
  void destroy_global_task()
  {
    WHEREAMI();
    // Ruby is responsible for freeing the global task when it is ready.
    g_rbTask = Qnil;
    g_cTask = NULL;
  }
  
  /*! Shutdown the interpreter */
  void finalize_ruby()
  {
    WHEREAMI();
  	ruby_finalize();
  }
  
  /*! This is used to make the call to rb_protect easier. */
  VALUE load_modules_for_irb(VALUE)
  {
    WHEREAMI();
    rb_require("psi4");
    rb_require("irb");
    return Qnil;
  }
  
  VALUE load_psi_module(VALUE)
  {
    WHEREAMI();
    rb_require("psi4");
    return Qnil;
  }
  
  void handle_ruby_exception()
  {
    WHEREAMI();
    VALUE err;
    
    // get the message from the exception.
    #if RUBY_MAJOR == 1 && RUBY_MINOR < 9
    err = rb_inspect(ruby_errinfo);
    #else
    err = rb_inspect(rb_errinfo());
    #endif
    // Print a backtrace to the screen, hopefully the user can make sense of it.
    rb_backtrace();
    // Print what the error message was.
    fprintf(stderr, "Ruby Exception caught: %s\n", StringValuePtr(err));
    exit(EXIT_FAILURE);
  }
  
  /*! Handles running Ruby interactively. This version makes use of the
      irb module found in Ruby to handle Ruby input parsing. */
  int run_interactive_ruby()
  {
    WHEREAMI();
    int error;
    
    printf("Starting interactive PSI4 driver:\n");
    // rb_protect is used to protect PSI from Ruby code exceptions
    // without it PSI would crash if a script did something nasty.
    rb_protect(load_modules_for_irb, 0, &error);
    if (error) {
      handle_ruby_exception();
      return -1;
    }
    
    // Run the Irb. This should probably be wrapped in a rb_protect as well.
    rb_eval_string("IRB.start");
    
    return true;
  }
  
  void create_psi_module()
  {
    WHEREAMI();
    g_rbPsi = rb_define_module("Psi");
    
    // Load in the psi4.rb file found in either lib/ruby or share/ruby.
    int error;
    
    // rb_protect is used to protect Psi from Ruby code exceptions
    // without it Psi would crash if a script did something nasty.
    rb_protect(load_psi_module, 0, &error);
    if (error) {
      handle_ruby_exception();
    }
  }
/*}*/}
