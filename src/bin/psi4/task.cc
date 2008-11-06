#include <ruby.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include "psi4.h"
#include "task.h"

namespace psi { //namespace psi4 {
  VALUE Task::Ruby::rbTask_ = Qnil;
  // TODO: Module arrays.
  
  Task::Task(std::string input, std::string p, std::string s)
  {
    // Set the psiPSIO_ object to the values given.
    // Default: file32 is kept in the current working directory.
    // Default: only one volume is used.
    // The ability to change individual files will be added.
    psiPSIO_.filecfg_kwd("DEFAULT", "NAME",    -1, p.c_str());
    psiPSIO_.filecfg_kwd("DEFAULT", "NVOLUME", -1, "1");
    psiPSIO_.filecfg_kwd("DEFAULT", "VOLUME1", -1, s.c_str());
    psiPSIO_.filecfg_kwd("DEFAULT", "VOLUME1", 32, "./");
    
    // Save the variables for quick access.
    szFilePrefix_ = p;
    szScratchPath_ = s;
    szInputFile_ = input;
  }
  
  Task::~Task()
  {
  }
  
  // TODO: Add module registration functions.
  
  
  const std::string& Task::prefix()
  {
    return szFilePrefix_;
  }
  
  void Task::prefix(std::string new_prefix)
  {
    szFilePrefix_ = new_prefix;
    
    // change the setting in psiPSIO_
    psiPSIO_.filecfg_kwd("DEFAULT", "NAME", -1, new_prefix.c_str());
  }
  
  const std::string& Task::scratch()
  {
    return szScratchPath_;
  }
  
  void Task::scratch(std::string new_scratch)
  {
    szScratchPath_ = new_scratch;
    
    // change the setting in psiPSIO_
    psiPSIO_.filecfg_kwd("DEFAULT", "VOLUME1", -1, new_scratch.c_str());
  }
  
  const std::string& Task::input_file()
  {
    return szInputFile_;
  }
  
  void Task::input_file(std::string new_input)
  {
    szInputFile_ = new_input;
  }
  
  psi::PSIO* Task::libpsio()
  {
    return &psiPSIO_;
  }
  
  FILE *Task::output()
  {
    return outfile;
  }
  
  // Ruby framework for Task
  void Task::Ruby::create_ruby_class()
  {
    WHEREAMI();
    
    // Create the task class under the Psi module
    Task::Ruby::rbTask_ = rb_define_class_under(g_rbPsi, "Task", rb_cObject);
    
    // Register the allocation function with Ruby
    rb_define_alloc_func(Task::Ruby::rbTask_, Task::Ruby::rb_alloc);
    
    // Register the initialization function with Ruby
    //
    // The structure of a call to rb_define_method is the following. If # of arguments is -1 then a
    // C-array of Ruby objects is sent. Note every single function has a "VALUE self" parameter that
    // is not counted in the number listed here.
    //
    //
    rb_define_method(Task::Ruby::rbTask_, "initialize", RUBYCAST(Task::Ruby::rb_init), -1);
    
    // Register the initialize copy function with Ruby
    rb_define_method(Task::Ruby::rbTask_, "initialize_copy", RUBYCAST(Task::Ruby::rb_init_copy), 1);
    rb_define_method(Task::Ruby::rbTask_, "to_s", RUBYCAST(Task::Ruby::rb_to_s), 0);
    
    // Psi specific functions
    rb_define_method(Task::Ruby::rbTask_, "prefix=",  RUBYCAST(Task::Ruby::rb_prefix_set), 1);
    rb_define_method(Task::Ruby::rbTask_, "prefix",   RUBYCAST(Task::Ruby::rb_prefix_get), 0);
    rb_define_method(Task::Ruby::rbTask_, "scratch=", RUBYCAST(Task::Ruby::rb_scratch_set), 1);
    rb_define_method(Task::Ruby::rbTask_, "scratch",  RUBYCAST(Task::Ruby::rb_scratch_get), 0);
    rb_define_method(Task::Ruby::rbTask_, "input=",   RUBYCAST(Task::Ruby::rb_input_file_set), 1);
    rb_define_method(Task::Ruby::rbTask_, "input",    RUBYCAST(Task::Ruby::rb_input_file_get), 0);
    
    // Run input file
    // rb_define_method(Task::Ruby::rbTask_, "run", RUBYCAST(Task::Ruby::rb_run_input_file), 0);
    rb_require("task");
    
    // Global Ruby methods are only needed when running in interactive mode. When running normally
    // the user's input file is "load"ed into the "run" function task, so that whichever instance
    // of task the run is called on is in scope of that instance. When running interactively this
    // does not happen; the commands do not execute in the instance of the global task. These 
    // functions provide most of the functions of task to the user. Some functions do not make sense
    // to include in the global namespace; for instance, to_s and run.
    if (g_bIRB) {
      printf("Registering access functions to global task object.\n");
      rb_define_global_function("prefix=",     RUBYCAST(Task::Ruby::rb_global_prefix_set),                   1);
      rb_define_global_function("prefix",      RUBYCAST(Task::Ruby::rb_global_prefix_get),                   0);
      rb_define_global_function("scratch=",    RUBYCAST(Task::Ruby::rb_global_scratch_set),                  1);
      rb_define_global_function("scratch",     RUBYCAST(Task::Ruby::rb_global_scratch_get),                  0);
      rb_define_global_function("input_file=", RUBYCAST(Task::Ruby::rb_global_input_file_set),               1);
      rb_define_global_function("input_file",  RUBYCAST(Task::Ruby::rb_global_input_file_get),               0);
    }
  }

  void Task::Ruby::rb_free(void *p)
  {
    WHEREAMI();
    
    Task *pTask = (Task*)p;
    if (pTask) delete pTask;
  }
  
  VALUE Task::Ruby::rb_alloc(VALUE klass)
  {
    WHEREAMI();
    
    Task *newTask = new Task;
    VALUE newObj;
    
    // Wrap the newly created Task inside a Ruby object, and tell Ruby how to destory it
    // when done
    newObj = Data_Wrap_Struct(klass, 0, Task::Ruby::rb_free, newTask);
    // Ruby is now responsible for the object...not us.
    
    return newObj;
  }
  
  VALUE Task::Ruby::rb_init(int argc, VALUE* argv, VALUE self)
  {
    WHEREAMI();
    
    Task *task;
    VALUE rbObject;
    
    // Get the Task object from Ruby
    Data_Get_Struct(self, Task, task);
    
    // What Ruby does is if the user uses a Hash for their function argument
    // is that argv[0] is the Hash. If the user did a simple comma separated
    // array then argv[...] are the individual elements.
    // This function assumes a Hash. I personally think hashes for function
    // parameters are better as it allows for named parameters and order
    // doesn't matter.
    if (argc > 1) {
      // Throw an exception in Ruby
      rb_raise(rb_eArgError, "invalid function parameters");
      return Qnil;
    }
    
    if (argc == 1) {
      // Work through the Hash
      
      // Did the user give a prefix?
      rbObject = rb_hash_aref(argv[0], ID2SYM(rb_intern("prefix")));
      if (rbObject != Qnil) {
        // Get the string and tell task about it
        // StringValue calls to_s on the object, if needed
        VALUE str = StringValue(rbObject);
        // Set the prefix variable to it
        task->prefix(RSTRING_PTR(str));
      }
      
      // Did the user give a scratch?
      rbObject = rb_hash_aref(argv[0], ID2SYM(rb_intern("scratch")));
      if (rbObject != Qnil) {
        VALUE str = StringValue(rbObject);
        task->scratch(RSTRING_PTR(str));
      }
    }
    
    // If the user did not provide any arguments, that's ok.
    return self;
  }
  
  VALUE Task::Ruby::rb_init_copy(VALUE copy, VALUE orig)
  {
    WHEREAMI();
    
    Task *o, *c;
    
    // Do not self copy
    if (copy == orig) return copy;
    
    // We can only copy a Task to a Task. The only way to really tell
    // is what function will be used to delete the object.
    if (TYPE(orig) != T_DATA || RDATA(orig)->dfree != (RUBY_DATA_FUNC)Task::Ruby::rb_free) {
      rb_raise(rb_eTypeError, "wrong argument type; expecting Task object");
      return copy;
    }
    
    // Get the objects and copy
    Data_Get_Struct(copy, Task, c);
    Data_Get_Struct(orig, Task, o);
    c->prefix(o->prefix());
    c->scratch(o->scratch());
    c->input_file(o->input_file());
    
    return copy;
  }
  
  VALUE Task::Ruby::rb_to_s(VALUE self)
  {
    WHEREAMI();
    
    Task *s;
    std::string text("Task: prefix = ");
    
    Data_Get_Struct(self, Task, s);
    text += s->prefix() + "; scratch = " + s->scratch() + "; input file = " + s->input_file();
    
    return rb_str_new2(text.c_str());
  }
  
  VALUE Task::Ruby::prefix_set(Task *task, VALUE newPrefix)
  {
    WHEREAMI();
    
    VALUE str = StringValue(newPrefix);
    task->prefix(RSTRING_PTR(str));
    
    return newPrefix;
  }
  
  VALUE Task::Ruby::rb_prefix_set(VALUE self, VALUE newPrefix)
  {
    WHEREAMI();
    
    Task *task;
    Data_Get_Struct(self, Task, task);
    return prefix_set(task, newPrefix);
  }
  
  VALUE Task::Ruby::rb_global_prefix_set(VALUE, VALUE newPrefix)
  {
    WHEREAMI();
    return prefix_set(g_cTask, newPrefix);
  }
  
  VALUE Task::Ruby::prefix_get(Task *task)
  {
    WHEREAMI();
    return rb_str_new2(task->prefix().c_str());
  }
  
  VALUE Task::Ruby::rb_prefix_get(VALUE self)
  {
    WHEREAMI();
    Task *task;
    Data_Get_Struct(self, Task, task);
    return prefix_get(task);
  }
  
  VALUE Task::Ruby::rb_global_prefix_get(VALUE)
  {
    WHEREAMI();
    return prefix_get(g_cTask);
  }
  
  VALUE Task::Ruby::scratch_set(Task *task, VALUE newScratch)
  {
    WHEREAMI();
    VALUE str = StringValue(newScratch);
    task->scratch(RSTRING_PTR(str));
    return newScratch;
  }
  
  VALUE Task::Ruby::rb_scratch_set(VALUE self, VALUE newScratch)
  {
    WHEREAMI();
    Task *task;
    Data_Get_Struct(self, Task, task);
    return scratch_set(task, newScratch);
  }
  
  VALUE Task::Ruby::rb_global_scratch_set(VALUE, VALUE newScratch)
  {
    WHEREAMI();
    return scratch_set(g_cTask, newScratch);
  }
  
  VALUE Task::Ruby::scratch_get(Task *task)
  {
    WHEREAMI();
    return rb_str_new2(task->scratch().c_str());
  }
  
  VALUE Task::Ruby::rb_scratch_get(VALUE self)
  {
    WHEREAMI();
    Task *task;
    Data_Get_Struct(self, Task, task);
    return scratch_get(task);
  }
  
  VALUE Task::Ruby::rb_global_scratch_get(VALUE)
  {
    WHEREAMI();
    return scratch_get(g_cTask);
  }
  
  VALUE Task::Ruby::input_file_set(Task *task, VALUE newInput)
  {
    WHEREAMI();
    VALUE str = StringValue(newInput);
    task->input_file(RSTRING_PTR(str));
    return newInput;
  }
  
  VALUE Task::Ruby::rb_input_file_set(VALUE self, VALUE newInput)
  {
    WHEREAMI();
    Task *task;
    Data_Get_Struct(self, Task, task);
    return input_file_set(task, newInput);
  }
  
  VALUE Task::Ruby::rb_global_input_file_set(VALUE, VALUE newInput)
  {
    WHEREAMI();
    return input_file_set(g_cTask, newInput);
  }
  
  VALUE Task::Ruby::input_file_get(Task *task)
  {
    WHEREAMI();
    return rb_str_new2(task->input_file().c_str());
  }
  
  VALUE Task::Ruby::rb_input_file_get(VALUE self)
  {
    WHEREAMI();
    Task *task;
    Data_Get_Struct(self, Task, task);
    return input_file_get(task);
  }
  
  VALUE Task::Ruby::rb_global_input_file_get(VALUE)
  {
    WHEREAMI();
    return input_file_get(g_cTask);
  }
  
  VALUE Task::Ruby::rb_run_input_file(VALUE self)
  {
    WHEREAMI();
    Task *task;
    Data_Get_Struct(self, Task, task);
    
    // Get access to the IO class
    VALUE io = rb_path2class("IO");
    
    // Have IO load in the input file
    VALUE input = rb_str_new2(task->input_file().c_str());
    
    // The first argument to rb_funcall is the object scope to call into.
    VALUE file = rb_funcall(io, rb_intern("read"), 1, input);
    
    rb_funcall(self, rb_intern("puts"), 1, file);
    
    // file now contains the input file
    // Evaluate it using self (the task) as the scope
    VALUE retVal = rb_funcall(self, rb_intern("class_eval"), 1, file);
    
    // Return the value;
    return retVal;
  }
  
  VALUE Task::Ruby::run_module(Task *task, int argc, VALUE *argv)
  {
    WHEREAMI();
    
    // How did the user call this function in Ruby?
    #if RUBY_MAJOR == 1 && RUBY_MINOR < 9
    ID last_func = rb_frame_last_func();
    #else
    ID last_func = rb_frame_callee();
    #endif
    std::string function_name(rb_id2name(last_func));
    
    // Go through the list of known modules until we find  match
    // then issue a call to it.
    return Qnil;
  }
  
  VALUE Task::Ruby::rb_run_module(int argc, VALUE* argv, VALUE self)
  {
    WHEREAMI();
    Task *task;
    Data_Get_Struct(self, Task, task);
    return run_module(task, argc, argv);
  }
  
  VALUE Task::Ruby::rb_global_run_module(int argc, VALUE* argv, VALUE)
  {
    WHEREAMI();
    return run_module(g_cTask, argc, argv);
  }
  
  VALUE Task::Ruby::puts(Task *task, int argc, VALUE *argv)
  {
    WHEREAMI();
    FILE *file = task->output();
    for (int i=0; i < argc; ++i) {
      VALUE str = StringValue(argv[i]);
      fprintf(file, "%s\n", RSTRING_PTR(str));
    }
    fprintf(file, "\n");
    fflush(file);
    
    return Qnil;
  }
  
  VALUE Task::Ruby::rb_puts(int argc, VALUE *argv, VALUE self)
  {
    WHEREAMI();
    Task *task;
    Data_Get_Struct(self, Task, task);
    return puts(task, argc, argv);
  }
  
  VALUE Task::Ruby::rb_global_puts(int argc, VALUE *argv, VALUE)
  {
    WHEREAMI();
    return puts(g_cTask, argc, argv);
  }
/*}*/}