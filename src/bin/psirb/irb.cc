#include <ruby.h>
#include "psirb.h"

namespace psi { namespace psirb { 
/*! Handles running Ruby interactively. This version makes use of the
    irb module found in Ruby to handle Ruby input parsing. */
void run_interactive_ruby()
{
	printf("Including Psi3's namespace ... ");
	// Run the equivalent of "require 'psi3'"
	rb_require("psi3");
	printf("done\nIncluding Ruby's IRB namespace ... ");
	// Run the equivalent of "require 'irb'"
	rb_require("irb");
	printf("done\n");
	
	// Run the Irb
	printf("Starting interactive Psi3 Ruby driver.\n");
	rb_eval_string("IRB.start");
}
}} // namespace psi::psirb
