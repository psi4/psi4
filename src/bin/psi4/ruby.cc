/*! \file ruby.cc
  \ingroup (psi4)
  \brief Contains Ruby initialization and some global functions.
*/

#include <ruby.h>
#include "psi4.h"
#include "task.h"

namespace psi { namespace psi4 {
  
  // Function declarations
  bool initialize_ruby();
  void load_input_file_into_ruby();
  void process_input_file();
  void finalize_ruby();
  int run_interactive_ruby();
  void handle_ruby_exception();
  
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
    
    // Done
    return true;
  }
  
  VALUE process_input_file_protected(VALUE)
  {
    // Run via the global task object.
    rb_funcall(Globals::g_rTask, rb_intern("run"), 0);
  }
  
  /*! Loads the Ruby input file into the interpreter. Does not perform syntax checking, nor does
      it begin executing the code. Just checks to see if the file exists and loads it. */
  #if 0
  void load_input_file_into_ruby()
  {
  	// Have Ruby do it
  	#if RUBY_MAJOR == 1 && RUBY_MINOR < 9
  	// In pre-Ruby 1.9 the input file was loaded into a single global space.
  	rb_load_file(Globals::g_szInputFile.c_str());
  	#else
  	// However, in Ruby 1.9 it seems there is the possibility of loading multiple input files.
  	// So rb_load_file loads the input file and returns an execution node that we tell Ruby
  	// to execute when ready.
    Globals::g_rbExecNode = rb_load_file(Globals::g_szInputFile.c_str());
    #endif
  }
  #endif
  
  /*! Run the input file. */
  void process_input_file()
  {
    #if 0
  	// Process the input file
  	// Since I do not want Ruby to take complete control of the system this is a 
  	// hack version from Ruby's eval.c file function ruby_eval and ruby_stop
  	// Typically you would call ruby_run but this function does not return at 
  	// all from ruby_stop, it invokes an exit function call internally.
  	#if RUBY_MAJOR == 1 && RUBY_MINOR < 9
  	int state;
    static int ex;

    if (ruby_nerrs > 0) exit(EXIT_FAILURE);
    state = ruby_exec();
    if (state && !ex) ex = state;
  	ruby_cleanup(ex);
  	#else
    if (Globals::g_rbExecNode == NULL) exit(EXIT_FAILURE);
    ruby_run_node(Globals::g_rbExecNode);
  	#endif
  	#endif // 0
  	
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
    // If this fails it will throw a Ruby exception
    // Ruby is used to create a new Task instance, fully wrapped in Ruby.
    Globals::g_rbTask = rb_class_new_instance(0, 0, Task::m_rbTask);
    // Get the C++ object from Ruby to use.
    Data_Get_Struct(Globals::g_rbTask, Task, Globals::g_cTask);
    // Prevent Ruby from garbage collecting this new object when this
    // function returns by making it a global variable. I'm not sure
    // how you would access this global variable in Ruby.
    // Make sure the global task object know what input file to use.
    Globals::g_cTask->input_task(Globals::g_szInputFile);
    
    // This function doesn't need to return anything of value.
    // But ALL functions that Ruby calls MUST return something.
    return Qnil;
  }
  
  /*! Create the global Task object */
  bool create_global_task()
  {
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
    // Ruby is responsible for freeing the global task when it is ready.
    Globals::g_rTask = Qnil;
    Globals::g_cTask = NULL;
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
  
  VALUE load_psi_module(VALUE)
  {
    rb_require("psi4");
    return Qnil;
  }
  
  void handle_ruby_exception()
  {
    VALUE err;
    
    // get the message from the exception.
    #if RUBY_MAJOR == 1 && RUBY_MINOR < 9
    err = rb_inspect(ruby_errinfo);
    #else
    err = rb_inspect(rb_errinfo());
    #endif
    // Prin a backtrace to the screen, hopefully the user can make sense of it.
    rb_backtrace();
    // Print what the error message was.
    printf("ERROR: %s\n", StringValuePtr(err));
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
      handle_ruby_exception();
      return -1;
    }
    
    // Run the Irb. This should probably be wrapped in a rb_protect as well.
    rb_eval_string("IRB.start");
    
    return true;
  }
  
  void create_psi_module()
  {
    Globals::g_rbPsi = rb_define_module("Psi");
    
    // Load in the psi4.rb file found in either lib/ruby or share/ruby.
    int error;
    
    // rb_protect is used to protect Psi from Ruby code exceptions
    // without it Psi would crash if a script did something nasty.
    rb_protect(load_psi_module, 0, &error);
    if (error) {
      handle_ruby_exception();
    }
  }
}}
