#ifndef __psi_psi4_task_h_
#define __psi_psi4_task_h_

#include <stdio.h>
#include <string>
#include <vector>
#include <ruby.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

#ifndef NDEBUG
  #ifdef __GNUC__
    #define WHEREAMI() if (g_bVerbose) { fprintf(stderr, "@:%s:%-6d\t-> %s\n", __FILE__, __LINE__, __PRETTY_FUNCTION__ ); }
  #else
    #define WHEREAMI() if (g_bVerbose) { fprintf(stderr, "@:%s:%-6d\n"); }
  #endif
#else
  #define WHEREAMI(func)
#endif

namespace psi { namespace psi4 {
  
  /* A calculation is group into a Task object. Tasks provide a simple
     way of having several prefixes in a single input file. There is a 
     global Task object that is used by default, creaed by the driver.
     All commands given in a user's input file correlate to the global
     task unless a new Task object is created and used explicitly. */
     
  // Might want to derive this from another class which can be made
  // available to other modules so that they need not include Ruby
  // header files.
  class Task
  {
      /*! The libpsio object for this task. */
      psi::PSIO psiPSIO_;
    
      /*! Equivalent to the psi_file_prefix */
      std::string szFilePrefix_;
    
      /*! Default scratch space that this Task is to use */
      std::string szScratchPath_;
    
      /*! Name of the input file to use */
      std::string szInputFile_;
    
      /*! Array of module entry points */
      // TODO: Array of module entry points.

    public:
      /*! Default constructor; sets input to "input.dat"; 
          prefix to "psi"; and scratch to "/tmp/" */
      Task(std::string input = "input.dat", std::string prefix = "psi",
           std::string scratch = "/tmp/");
      
      /*! Destructor */
      virtual ~Task();
      
      /*! Register a new module with the driver */
      
      /*! Enable all modules; this function makes the necessary connection
          into Ruby allowing scripts to issue calls to psi modules */
      
      /*! Acessor function for the prefix
      \return current prefix */
      const std::string& prefix();
      
      /*! Accessor function for the prefix
      \param new_prefix new prefix */
      void prefix(std::string new_prefix);

      /*! Acessor function for the scratch
      \return current scratch location */
      const std::string& scratch();
      
      /*! Accessor function for the scratch
      \param new_scratch new scratch location */
      void scratch(std::string new_scratch);
      
      /*! Acessor function for the input file
      \return current input file */
      const std::string& input_file();
      
      /*! Accessor function for the input file
      \param new_input new input file */
      void input_file(std::string new_input);
      
      /*! Accessor function to retrieve the libpsio object associated
      with this Task. Do not delete this pointer. Should probably wrap
      this in a Ref<> object.
      \return libpsio object pointer */
      psi::PSIO *libpsio();
      
      /*! Accessor function to the file handle to send output to. */
      FILE *output();
      
      class Ruby {
          // The following commands work on the Task object that is sent to it.
          static VALUE prefix_set(Task *task, VALUE newPrefix);
          static VALUE prefix_get(Task *task);
          static VALUE scratch_set(Task *task, VALUE newScratch);
          static VALUE scratch_get(Task *task);
          static VALUE input_file_set(Task *task, VALUE newInput);
          static VALUE input_file_get(Task *task);
        
          static VALUE puts(Task *task, int argc, VALUE *argv);
          static VALUE run_module(Task *task, int argc, VALUE *argv);
        public:
          /*! Ruby reference to the Task class descriptor */
          static VALUE rbTask_;
          
          // Ruby framework for Task
          /*! Create the Ruby class framework */
          static void create_ruby_class();
          
          /*! Called by Ruby when it needs to delete a class */
          static void rb_free(void *p);
          
          /*! Called by Ruby during object creation */
          static VALUE rb_alloc(VALUE klass);
          
          /*! Called by Ruby during object creation */
          static VALUE rb_init(int argc, VALUE* argv, VALUE self);

          /*! Caled by Ruby during object copy creation */
          static VALUE rb_init_copy(VALUE copy, VALUE orig);

          /*! Called by Ruby if the user tries to print a Task object */
          static VALUE rb_to_s(VALUE self);

          /*! Ruby function: Task.prefix= or Task.set_prefix */
          static VALUE rb_prefix_set(VALUE self, VALUE newPrefix);
          static VALUE rb_global_prefix_set(VALUE, VALUE newPrefix);

          /*! Ruby function: Task.prefix or Task.get_prefix */
          static VALUE rb_prefix_get(VALUE self);
          static VALUE rb_global_prefix_get(VALUE);

          /*! Ruby function: Task.scratch= or Task.set_scratch */
          static VALUE rb_scratch_set(VALUE self, VALUE newScratch);
          static VALUE rb_global_scratch_set(VALUE, VALUE newScratch);

          /*! Ruby function: Task.scratch or Task.get_scratch */
          static VALUE rb_scratch_get(VALUE self);
          static VALUE rb_global_scratch_get(VALUE);

          /*! Ruby function: Task.input_file= or Task.set_input_file */
          static VALUE rb_input_file_set(VALUE self, VALUE newInput);
          static VALUE rb_global_input_file_set(VALUE, VALUE newInput);

          /*! Ruby function: Task.input_file or Task.get_input_file */
          static VALUE rb_input_file_get(VALUE self);
          static VALUE rb_global_input_file_get(VALUE);

          /*! Ruby function: Task.run 
          This function is the one responsible for loading and executing a user's input file.*/
          static VALUE rb_run_input_file(VALUE self);

          /*! Common entry point to run modules. It handles figuring out which module the
          user wanted to run. */
          static VALUE rb_run_module(int argc, VALUE* argv, VALUE self);
          static VALUE rb_global_run_module(int argc, VALUE* argv, VALUE);

          /*! Ruby print functions */
          static VALUE rb_puts(int argc, VALUE *argv, VALUE self);
          static VALUE rb_global_puts(int argc, VALUE *argv, VALUE);

          static VALUE rb_test(VALUE self);
      };

    friend class Task::Ruby;
  };
}}
#endif // __psi_psi4_task_h_
