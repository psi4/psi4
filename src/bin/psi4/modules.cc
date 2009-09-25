#include <stdio.h>
#include <ruby.h>
#include <psiconfig.h>
#include <liboptions/liboptions.hpppp>

#include <psi4-dec.h>
#include "task.h"

  
namespace psi { 
  namespace input    { int input(); void register_input_options(); }
  // namespace CINTS    { int cints(Options &, int argc, char *argv[]); }
  // namespace cscf     { int cscf(int argc, char *argv[]); }
  // namespace psiclean { int psiclean(int argc, char *argv[]); }
  
  ModuleInformation modules_[] = {
    { "input", &input::input, &input::register_input_options },
    { "",       NULL,          NULL}                       // Last module
  };
  
  // Could be called multiple times.
  void register_modules()
  {
    WHEREAMI();
    unsigned int i = 0;
    
    while (modules_[i].entryPoint != NULL) {
      
      // Have the module register its options
      if (modules_[i].registerOptions)
        modules_[i].registerOptions();
      
      // Add the 
      i += 1;
    }    
  }
  
  // Will only be called once.
  void enable_modules()
  {
    WHEREAMI();
    unsigned int i = 0;
    
    while (modules_[i].entryPoint != NULL && modules_[i].rubyFunctionName.size() > 0) {
      rb_define_method(Task::Ruby::rbTask_, modules_[i].rubyFunctionName.c_str(), RUBYCAST(Task::Ruby::rb_run_module), -1);
      if (g_bIRB)
        rb_define_global_function(modules_[i].rubyFunctionName.c_str(), RUBYCAST(Task::Ruby::rb_run_module), -1);
    }
  }
}
