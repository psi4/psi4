#include <ruby.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include "psirb.h"

namespace psi { namespace psirb {
	
using namespace psi;

extern "C" {
extern void free_block(double **);
};

//
// class Task
VALUE Task::m_rbTask = Qnil;

Task::Task(std::string prefix, std::string scratch)
{
	// Set the m_psiPSIO object to the values given:
	// By default file32 is kept in the current working directory.
	// By default only one volume is used.
	// The ability to change individual files will be added.
	m_psiPSIO.filecfg_kwd("DEFAULT", "NAME",    -1, prefix.c_str());
	m_psiPSIO.filecfg_kwd("DEFAULT", "NVOLUME", -1, "1");
	m_psiPSIO.filecfg_kwd("DEFAULT", "VOLUME1", -1, scratch.c_str());
	m_psiPSIO.filecfg_kwd("DEFAULT", "VOLUME1", 32, "./");
	
	// Save the variables for quick access
	m_szFilePrefix = prefix;
	m_szScratchPath = scratch;
}

Task::~Task()
{
	// Currently nothing to do.
}

const std::string& Task::prefix()
{
	return m_szFilePrefix;
}

void Task::prefix(std::string new_prefix)
{
	m_szFilePrefix = new_prefix;
	
	// Change the setting in m_psiPSIO
	m_psiPSIO.filecfg_kwd("DEFAULT", "NAME", -1, new_prefix.c_str());
}

const std::string& Task::scratch()
{
	return m_szScratchPath;
}

void Task::scratch(std::string new_scratch)
{
	m_szScratchPath = new_scratch;
	
	// Change the setting in m_psiPSIO;
	m_psiPSIO.filecfg_kwd("DEFAULT", "VOLUME1", -1, new_scratch.c_str());
}

psi::PSIO *Task::libpsio()
{
	return &m_psiPSIO;
}

//
// Ruby framework for Task
//

void Task::create_ruby_class()
{
	// Create the Task class under the Psi module
	Task::m_rbTask = rb_define_class_under(Globals::g_rbPsi, "Task", rb_cObject);
	
	// Register the allocation functioin with Ruby
	rb_define_alloc_func(Task::m_rbTask, Task::rb_alloc);
	
	
	// Register the initialization function with Ruby
	//
	// The structure of a call to rb_define_method is the following. If # of arguments is -1 then an C-array of
	// Ruby objects is sent. Note every single function has a "VALUE self" parameter that is not counted in the
	// number listed here.
	//
	//               Ruby object      Ruby function name    C/++ function to be called                  # of arguments
	rb_define_method(Task::m_rbTask, "initialize",          RUBYCAST(Task::rb_init),                        -1);
	
	// Register the initialization copy function with Ruby
	rb_define_method(Task::m_rbTask, "initialize_copy",     RUBYCAST(Task::rb_init_copy),                    1);
	
	// Register the the to_s function.
	rb_define_method(Task::m_rbTask, "to_s",                RUBYCAST(Task::rb_to_s),                         0);
	
	// Add other functions starting here.
	rb_define_method(Task::m_rbTask, "prefix=",             RUBYCAST(Task::rb_prefix_set),                   1);
	rb_define_method(Task::m_rbTask, "prefix",              RUBYCAST(Task::rb_prefix_get),                   0);
	rb_define_method(Task::m_rbTask, "scratch=",            RUBYCAST(Task::rb_scratch_set),                  1);
	rb_define_method(Task::m_rbTask, "scratch",             RUBYCAST(Task::rb_scratch_get),                  0);
	rb_define_method(Task::m_rbTask, "load_zmat",           RUBYCAST(Task::rb_load_zmat),                    0);
	rb_define_method(Task::m_rbTask, "load_cartesian",      RUBYCAST(Task::rb_load_cartesian),               0);
	
	// Interface to libpsio++
	rb_define_method(Task::m_rbTask, "print_toc",           RUBYCAST(Task::rb_print_toc),                    1);
	
	// Checkpoint interface
	rb_define_method(Task::m_rbTask, "exist?",              RUBYCAST(Task::rb_chkpt_exist),                  1);
	rb_define_method(Task::m_rbTask, "exists?",             RUBYCAST(Task::rb_chkpt_exist),                  1);
	rb_define_method(Task::m_rbTask, "label",               RUBYCAST(Task::rb_chkpt_label_get),              0);
	rb_define_method(Task::m_rbTask, "label=",              RUBYCAST(Task::rb_chkpt_label_set),              1); 
	rb_define_method(Task::m_rbTask, "escf",                RUBYCAST(Task::rb_chkpt_escf_get),               0);
	rb_define_method(Task::m_rbTask, "escf=",               RUBYCAST(Task::rb_chkpt_escf_set),               1);
	rb_define_method(Task::m_rbTask, "eref",                RUBYCAST(Task::rb_chkpt_eref_get),               0);
	rb_define_method(Task::m_rbTask, "eref=",               RUBYCAST(Task::rb_chkpt_eref_set),               1);
	rb_define_method(Task::m_rbTask, "ecorr",               RUBYCAST(Task::rb_chkpt_ecorr_get),              0);
	rb_define_method(Task::m_rbTask, "ecorr=",              RUBYCAST(Task::rb_chkpt_ecorr_set),              1);
	rb_define_method(Task::m_rbTask, "enuc",                RUBYCAST(Task::rb_chkpt_enuc_get),               0);
	rb_define_method(Task::m_rbTask, "enuc=",               RUBYCAST(Task::rb_chkpt_enuc_set),               1);
	rb_define_method(Task::m_rbTask, "efzc",                RUBYCAST(Task::rb_chkpt_efzc_get),               0);
	rb_define_method(Task::m_rbTask, "efzc=",               RUBYCAST(Task::rb_chkpt_efzc_set),               1);
	rb_define_method(Task::m_rbTask, "etot",                RUBYCAST(Task::rb_chkpt_etot_get),               0);
	rb_define_method(Task::m_rbTask, "etot=",               RUBYCAST(Task::rb_chkpt_etot_set),               1);
	rb_define_method(Task::m_rbTask, "disp",                RUBYCAST(Task::rb_chkpt_disp_get),               0);
	rb_define_method(Task::m_rbTask, "disp=", 	            RUBYCAST(Task::rb_chkpt_disp_set),               1);
	rb_define_method(Task::m_rbTask, "eccsd",               RUBYCAST(Task::rb_chkpt_eccsd_get),              0);
	rb_define_method(Task::m_rbTask, "e_t",                 RUBYCAST(Task::rb_chkpt_e_t_get),                0);
	rb_define_method(Task::m_rbTask, "emp2",                RUBYCAST(Task::rb_chkpt_emp2_get),               0);
	rb_define_method(Task::m_rbTask, "eom_states_energy",   RUBYCAST(Task::rb_chkpt_eom_state_energies_get), 0);
	rb_define_method(Task::m_rbTask, "num_irreps",          RUBYCAST(Task::rb_chkpt_num_irreps_get),         0);
	rb_define_method(Task::m_rbTask, "clsdpi=",             RUBYCAST(Task::rb_chkpt_clsdpi_set),             1);
	rb_define_method(Task::m_rbTask, "clsdpi",              RUBYCAST(Task::rb_chkpt_clsdpi_get),             0);
	rb_define_method(Task::m_rbTask, "frzcpi=",             RUBYCAST(Task::rb_chkpt_frzcpi_set),             1);
	rb_define_method(Task::m_rbTask, "frzcpi",              RUBYCAST(Task::rb_chkpt_frzcpi_get),             0);
	rb_define_method(Task::m_rbTask, "frzvpi=",             RUBYCAST(Task::rb_chkpt_frzvpi_set),             1);
	rb_define_method(Task::m_rbTask, "frzvpi",              RUBYCAST(Task::rb_chkpt_frzvpi_get),             0);
	rb_define_method(Task::m_rbTask, "evals",               RUBYCAST(Task::rb_chkpt_evals_get),              0);
	rb_define_method(Task::m_rbTask, "evals=",              RUBYCAST(Task::rb_chkpt_evals_set),              1);
	rb_define_method(Task::m_rbTask, "alpha_evals",         RUBYCAST(Task::rb_chkpt_alpha_evals_get),        0);
	rb_define_method(Task::m_rbTask, "alpha_evals=",        RUBYCAST(Task::rb_chkpt_alpha_evals_set),        1);
	rb_define_method(Task::m_rbTask, "beta_evals",          RUBYCAST(Task::rb_chkpt_beta_evals_get),         0);
	rb_define_method(Task::m_rbTask, "beta_evals=",         RUBYCAST(Task::rb_chkpt_beta_evals_set),         1);
	rb_define_method(Task::m_rbTask, "exps",                RUBYCAST(Task::rb_chkpt_exps_get),               0);
	rb_define_method(Task::m_rbTask, "exps=",               RUBYCAST(Task::rb_chkpt_exps_set),               1);
	
	// A checkpoint function that uses Psi::Matrix
//	rb_define_method(Task::m_rbTask, "fgeom=",              RUBYCAST(Task::rb_chkpt_fgeom_get),              0);
}

void Task::rb_free(void *p)
{
	Task *pTask = (Task*)p;
	if (pTask)
		delete pTask;
}

VALUE Task::rb_alloc(VALUE klass)
{
	Task *newTask = new Task;
	VALUE newObj;
	
	// Wrap the newly created Task inside a Ruby object.
	newObj = Data_Wrap_Struct(klass, 0, Task::rb_free, newTask);
	// Ruby is now responsible for the object, not us.
	
	// Return the new object.
	return newObj;
}

VALUE Task::rb_init(int argc, VALUE* argv, VALUE self)
{
	Task *task;
	VALUE rbObject;
	
	// Get the Task object from Ruby
	Data_Get_Struct(self, Task, task);
	
	// What Ruby does is that if the user uses a Hash for their function
	// arguments argv[0] is the Hash. If the user did a simple comma
	// separated array then argv[...] are the individual elements.
	// This function assumes a Hash. I think hashes for function parameters are
	// better as it allows for any parameter ordering.
	
	if (argc > 1) {
		// Throw an exception in Ruby
		rb_raise(rb_eArgError, "must use a Hash as your function parameters");
		return Qnil;
	}
	
	if (argc == 1) {
		// Work through the Hash
		
		// Did the user give a prefix?
		// Since I expect the user to use :prefix need to convert ID prefix
		// to Symbol :prefix.
		rbObject = rb_hash_aref(argv[0], ID2SYM(rb_intern("prefix")));
		if (rbObject != Qnil) {
			// Get the string and tell Task about it
			// StringValue calls to_str on the object, if needed
			VALUE str = StringValue(rbObject);
			// Set the prefix variable to this
			task->prefix(RSTRING(str)->ptr);
		}
		
		// Did the user give a scratch?
		rbObject = rb_hash_aref(argv[0], ID2SYM(rb_intern("scratch")));
		if (rbObject != Qnil) {
			// Get the string and tell Task about it
			// StringValue calls to_str on the object, if needed
			VALUE str = StringValue(rbObject);
			// Set the prefix variable to this
			task->scratch(RSTRING(str)->ptr);
		}
	}
	
	// If the user did not provide any arguments, that's ok
	// Task sets some values by default.
	return self;
}

VALUE Task::rb_init_copy(VALUE copy, VALUE orig)
{
	Task *o, *c;
	
	// Do not self copy
	if (copy == orig)
		return copy;
		
	// We can only copy a Task to a Task. This function is called on the copy
	// object by Ruby, so we must make sure orig is a Task. Simple thing to
	// do is check if orig is freed by Task::free.
	if (TYPE(orig) != T_DATA || RDATA(orig)->dfree != (RUBY_DATA_FUNC)Task::rb_free) {
		rb_raise(rb_eTypeError, "wrong argument type");
		return copy;
	}

	// Get the objects and copy
	Data_Get_Struct(copy, Task, c);
	Data_Get_Struct(orig, Task, o);
	
	c->prefix(o->prefix());
	c->scratch(o->scratch());
	
	return copy;
}

VALUE Task::rb_to_s(VALUE self)
{
	Task *s;
	std::string text("Task: prefix = ");
	
	Data_Get_Struct(self, Task, s);
	text += s->prefix() + "; scratch = " + s->scratch();
	
	return rb_str_new2(text.c_str());
}

VALUE Task::rb_prefix_set(VALUE self, VALUE newPrefix)
{
	Task *task;
	
	Data_Get_Struct(self, Task, task);
	
	// Get the string and tell Task about it
	// StringValue calls to_str on the object, if needed
	VALUE str = StringValue(newPrefix);
	// Set the prefix variable to this
	task->prefix(RSTRING(str)->ptr);

	return self;
}

VALUE Task::rb_prefix_get(VALUE self)
{
	Task *task;
	
	Data_Get_Struct(self, Task, task);
	
	// Get the string from Task and return it.
	return rb_str_new2(task->prefix().c_str());
}

VALUE Task::rb_scratch_set(VALUE self, VALUE newScratch)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	
	// Get the string and tell Task about it
	// StringValue calls to_str on the object, if needed
	VALUE str = StringValue(newScratch);
	// Set the prefix variable to this
	task->scratch(RSTRING(str)->ptr);

	return self;
}

VALUE Task::rb_scratch_get(VALUE self)
{
	Task *task;
	
	Data_Get_Struct(self, Task, task);
	
	// Get the string from Task and return it.
	return rb_str_new2(task->scratch().c_str());
}

//! Ruby functions: Psi::Task.exist? and Psi::Task.exists? 
/*! Ruby interface to chkpt_exist. Checks to see if the requested keyword exists in the
	checkpoint file.
	\param self Ruby object calling this function.
	\param keyword Keyword to check for.
	\return Qtrue (C) / true (Ruby) if it exists, or Qfalse (C) / false (Ruby) if it does not.
*/
VALUE Task::rb_chkpt_exist(VALUE self, VALUE keyword)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	
	// Convert the given keyword to a C-string
	VALUE str = StringValue(keyword);
	char *p = RSTRING(str)->ptr;
	char *keyw = NULL;
	
	if (p == NULL) {
		rb_raise(rb_eArgError, "wrong argument, expected a string");
	}
	
	// Call the Psi chkpt function
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	keyw = chkpt.build_keyword(p);
	int result = chkpt.exist(keyw);
	free(keyw);
	
	if (result) return Qtrue;
	else        return Qfalse;
}

//! Ruby function: Psi::Task.label
/*! Ruby interface to chkpt_rd_label. Reads the label from checkpoint.
	\param self Ruby object calling this function.
	\return Label as a Ruby string.
*/
VALUE Task::rb_chkpt_label_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	// Read in the label from Chkpt
	char *label = NULL;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	label = chkpt.rd_label();
	
	// Create ruby string
	VALUE result = rb_str_new2(label);
	
	return result;
}

//! Ruby function: Psi::Task.label
/*! Ruby interface to chkpt_wt_label. Writes the label to checkpoint.
	\param self Ruby object calling this function.
	\param label New label.
*/
VALUE Task::rb_chkpt_label_set(VALUE self, VALUE label)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	
	// Convert the given keyword to a C-string
	VALUE str = StringValue(label);
	char *p = RSTRING(str)->ptr;	
	if (p == NULL) {
		rb_raise(rb_eArgError, "wrong argument, expected a string");
	}
	
	// Call the Psi chkpt function
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	chkpt.wt_label(p);
	
	return self;
}

//! Ruby function: Psi::Task.escf
/*! Ruby interface to chkpt_rd_escf. Reads SCF energy from checkpoint.
	\param self Ruby object calling this function.
	\return SCF energy in a Ruby object.
*/
VALUE Task::rb_chkpt_escf_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	// Read in the scf
	double escf;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	escf = chkpt.rd_escf();
	
	VALUE result = rb_float_new(escf);
	return result;
}

//! Ruby function: Psi::Task.escf=
/*! Ruby interface to chkpt_wt_escf. Writes the SCF energy to checkpoint.
	\param self Ruby object calling this function.
	\param vescf New SCF energy.
*/
VALUE Task::rb_chkpt_escf_set(VALUE self, VALUE vescf)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	// Read in the scf
	double escf = NUM2DBL(vescf);
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	chkpt.wt_escf(escf);
	
	return self;
}

//! Ruby function: Psi::Task.eref
/*! Ruby interface to chkpt_rd_eref. Read the reference energy from checkpoint.
	\param self Ruby object calling this function.
	\return Reference energy in a Ruby object.
*/
VALUE Task::rb_chkpt_eref_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	// Read in the scf
	double escf;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	escf = chkpt.rd_escf();
	
	VALUE result = rb_float_new(escf);
	return result;
}

//! Ruby function: Psi::Task.eref=
/*! Ruby interface to chkpt_wt_ref. Writes the reference energy to checkpoint.
	\param self Ruby object calling this function.
	\param veref New reference energy as a Ruby object.
*/
VALUE Task::rb_chkpt_eref_set(VALUE self, VALUE veref)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	// Read in the scf
	double eref = NUM2DBL(veref);
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	chkpt.wt_eref(eref);
	
	return self;
}

//! Ruby function: Psi::Task.ecorr
/*! Ruby interface to chkpt_rd_ecorr. Read the correlation energy from checkpoint.
	\param self Ruby object calling this function.
	\return Correlation energy as a Ruby object.
*/
VALUE Task::rb_chkpt_ecorr_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	double ecorr;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	ecorr = chkpt.rd_ecorr();
	
	VALUE result = rb_float_new(ecorr);
	return result;
}

//! Ruby function: Psi::Task.ecorr=
/*! Ruby interface to chkpt_wt_ecorr. Writes the new correlation energy to checkpoint.
	\param self Ruby object calling this function.
	\param vecorr New correlation energy as a Ruby object.
*/
VALUE Task::rb_chkpt_ecorr_set(VALUE self, VALUE vecorr)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	double ecorr = NUM2DBL(vecorr);
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	chkpt.wt_ecorr(ecorr);
	
	return self;
}

//! Ruby function: Psi::Task.enuc
/*! Ruby interface to chkpt_rd_enuc. Read the nuclear repulsion energy to checkpoint.
	\param self Ruby object calling this function.
	\return Nuclear repulsion energy as a Ruby object.
*/
VALUE Task::rb_chkpt_enuc_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	double enuc;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	enuc = chkpt.rd_enuc();
	
	VALUE result = rb_float_new(enuc);
	return result;
}

//! Ruby function: Psi::Chkpt:enuc=
/*! Ruby interface to chkpt_wt_enuc. Writes the nuclear repulsion energy to checkpoint.
	\param self Ruby object calling this function.
	\param venuc New nuclear repulsion energy as a Ruby object.
*/
VALUE Task::rb_chkpt_enuc_set(VALUE self, VALUE venuc)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);

	double enuc = NUM2DBL(venuc);
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	chkpt.wt_enuc(enuc);
	
	return self;
}

//! Ruby function: Psi::Task.efzc
/*! Ruby interface to chkpt_rd_efzc. Reads the frozen core energy from checkpoint.
	\param self Ruby object that is calling this function.
	\return Frozen core energy in a Ruby object.
*/
VALUE Task::rb_chkpt_efzc_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	double efzc;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	efzc = chkpt.rd_efzc();
	
	VALUE result = rb_float_new(efzc);
	return result;
}

//! Ruby function: Psi::Task.efzc=
/*! Ruby interface to chkpt_wt_efzc. Writes the frozen core energy to checkpoint.
	\param self Ruby object that is calling this function.
	\param vefzc New frozen core energy.
*/
VALUE Task::rb_chkpt_efzc_set(VALUE self, VALUE vefzc)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	double efzc = NUM2DBL(vefzc);
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	chkpt.wt_efzc(efzc);
	
	return self;
}

//! Ruby function: Psi::Task.etot
/*! Ruby interface to chkpt_rd_etot. Reads the total energy from checkpoint file.
	\param self Ruby object that is calling this function.
	\return Total energy in a Ruby object.
*/
VALUE Task::rb_chkpt_etot_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	double etot;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	etot = chkpt.rd_etot();
	
	VALUE result = rb_float_new(etot);
	return result;
}

//! Ruby function: Psi::Task.etot=
/*! Ruby interface to chkpt_wt_etot. Write the new total energy to checkpoint file.
	\param self Ruby object that is calling this function.
	\param vetot New total energy.
*/
VALUE Task::rb_chkpt_etot_set(VALUE self, VALUE vetot)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	double etot = NUM2DBL(vetot);
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	chkpt.wt_etot(etot);
	
	return self;
}

//! Ruby function: Psi::Task.num_irreps
/*! Ruby interface to chkpt_rd_nirreps. */
VALUE Task::rb_chkpt_num_irreps_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	int num_irreps;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	num_irreps = chkpt.rd_nirreps();
	
	return INT2FIX(num_irreps);
}

//! Ruby function: Psi::Task.disp
/*! Ruby interface to chkpt_rd_disp. Reads the current geometry displacement number.
	\param self Ruby object that is calling this function.
	\return Current displacement number.
*/
VALUE Task::rb_chkpt_disp_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	int disp;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	disp = chkpt.rd_disp();
	
	VALUE result = INT2FIX(disp);
	return result;
}

//! Ruby function: Psi::Task.disp=
/*! Ruby interface to chkpt_wt_disp. Writes out the current geometry displacement number.
	\param self Ruby object that is calling this function.
	\param ndisp New displacement number.
*/
VALUE Task::rb_chkpt_disp_set(VALUE self, VALUE ndisp)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	int disp = NUM2INT(ndisp);
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	chkpt.wt_disp(disp);
	
	return self;
}

//! Ruby function: Psi::Task.eccsd
/*! Returns the CCSD contribution to the total energy.
	\param self The Ruby object that is calling this function.
	\return CCSD energy contribution as a Ruby object.
*/
VALUE Task::rb_chkpt_eccsd_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	double energy;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	energy = chkpt.rd_eccsd();
	
	// Return the value to the user
	return rb_float_new(energy);
}

//! Ruby function: Psi::Task.e_t
/*! Returns the (T) contribution to the total energy.
	\param self The Ruby object that is calling this function.
	\return (T) energy contribution as a Ruby object.
*/
VALUE Task::rb_chkpt_e_t_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	double energy;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	energy = chkpt.rd_e_t();
	
	// Return the value to the user
	return rb_float_new(energy);
}

//! Ruby function: Psi::Task.emp2
/*! Returns the MP2 contribution to the total energy.
	\param self The Ruby object that is calling this function.
	\return MP2 energy contribution as a Ruby object.
*/
VALUE Task::rb_chkpt_emp2_get(VALUE self)
{
	Task *task;	
	Data_Get_Struct(self, Task, task);
	double energy;

	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	energy = chkpt.rd_emp2();

	// Return the value to the user
	return rb_float_new(energy);
}

//! Ruby function: Psi::Task.eom_states_energy
/*! Returns an array containing the EOM energies for the requested states.
    \param self The Ruby object that is calling this function.
    \return An array containing the energies.
*/
VALUE Task::rb_chkpt_eom_state_energies_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	int count, i;
	double *energies;
	char *keyword;
	
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	keyword = chkpt.build_keyword("EOM Number of States");
	
	// Read in the number
	task->m_psiPSIO.read_entry(32, keyword, (char*)&count, sizeof(int));
	free(keyword);
	
	// Allocate memory
	energies = (double*)malloc(sizeof(double) * count);
	
	// Read in the energies
	keyword = chkpt.build_keyword("EOM State Energies");
	task->m_psiPSIO.read_entry(32, keyword, (char*)energies, sizeof(double) * count);
	free(keyword);

	// Convert the energy array to a Ruby array
	VALUE array = rb_ary_new();
	for (i=0; i<count; i++)
		rb_ary_push(array, rb_float_new(energies[i]));
		
	psi::Chkpt::free(energies);
	
	return array;
}

//! Ruby function: Psi::Task.print_toc(unit)
/*! Prints the TOC entries for unit to Globals::g_fOutput */
VALUE Task::rb_print_toc(VALUE self, VALUE rUnit)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	unsigned int unit = NUM2UINT(rUnit);
	int bAlreadyOpen;
	
	// Is the file already open?
	bAlreadyOpen = task->m_psiPSIO.open_check(unit);
	if (bAlreadyOpen == false)
		task->m_psiPSIO.open(unit, PSIO_OPEN_OLD);
		
	// Print the toc
	task->m_psiPSIO.tocprint(unit, Globals::g_fOutput);
	
	// If it was already open, do not close it.
	if (bAlreadyOpen == false)
		task->m_psiPSIO.close(unit, 1);
		
	return self;
}

//! Ruby function: Psi::Task.clsdpi=(array)
/*! Saves the clsdpi array to checkpoint */
VALUE Task::rb_chkpt_clsdpi_set(VALUE self, VALUE arr)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	int *clsdpi;
	int count, expectedCount;
	
	// Convert the arr to a C-array. The function call handles type checking
	count = create_array(arr, &clsdpi);
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	expectedCount = chkpt.rd_nirreps();
	
	if (count != expectedCount) {
		rb_raise(rb_eArgError, "array is of length %d, expecting %d", count, expectedCount);
		return self;
	}
	
	// Save a new clsdpi
	chkpt.wt_clsdpi(clsdpi);
	
	// Free memory
	psi::Chkpt::free(clsdpi);
	
	return self;
}

//! Ruby function: Psi::Task.clsdpi
/*! Reads in the clsdpi array from checkpoint */
VALUE Task::rb_chkpt_clsdpi_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	int *clsdpi;
	int count;
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	count = chkpt.rd_nirreps();
	clsdpi = chkpt.rd_clsdpi();
	
	VALUE arr = create_array(count, clsdpi);
	psi::Chkpt::free(clsdpi);
	
	return arr;
}

//! Ruby function: Psi::Task.frzcpi=(array)
/*! Saves the frzcpi array to checkpoint */
VALUE Task::rb_chkpt_frzcpi_set(VALUE self, VALUE arr)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	int *frzcpi;
	int count, expectedCount;
	
	// Convert the arr to a C-array. The function call handles type checking
	count = create_array(arr, &frzcpi);
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	expectedCount = chkpt.rd_nirreps();
	
	if (count != expectedCount) {
		rb_raise(rb_eArgError, "array is of length %d, expecting %d", count, expectedCount);
		return self;
	}
	
	// Save a new clsdpi
	chkpt.wt_frzcpi(frzcpi);
	
	// Free memory
	psi::Chkpt::free(frzcpi);
	
	return self;
}

//! Ruby function: Psi::Task.frzcpi
/*! Reads in the frzcpi array from checkpoint */
VALUE Task::rb_chkpt_frzcpi_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	int *frzcpi;
	int count;
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	count = chkpt.rd_nirreps();
	frzcpi = chkpt.rd_frzcpi();
	
	VALUE arr = create_array(count, frzcpi);
	psi::Chkpt::free(frzcpi);
	
	return arr;
}

//! Ruby function: Psi::Task.frzvpi=(array)
/*! Saves the frzvpi array to checkpoint */
VALUE Task::rb_chkpt_frzvpi_set(VALUE self, VALUE arr)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	int *frzvpi;
	int count, expectedCount;
	
	// Convert the arr to a C-array. The function call handles type checking
	count = create_array(arr, &frzvpi);
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	expectedCount = chkpt.rd_nirreps();
	
	if (count != expectedCount) {
		rb_raise(rb_eArgError, "array is of length %d, expecting %d", count, expectedCount);
		return self;
	}
	
	// Save a new clsdpi
	chkpt.wt_frzvpi(frzvpi);
	
	// Free memory
	psi::Chkpt::free(frzvpi);
	
	return self;
}

//! Ruby function: Psi::Task.frzcpi
/*! Reads in the frzcpi array from checkpoint */
VALUE Task::rb_chkpt_frzvpi_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	int *frzvpi;
	int count;
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	count = chkpt.rd_nirreps();
	frzvpi = chkpt.rd_frzvpi();
	
	VALUE arr = create_array(count, frzvpi);
	psi::Chkpt::free(frzvpi);
	
	return arr;
}

VALUE Task::rb_chkpt_evals_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	double *evals;
	int count;
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	count = chkpt.rd_nmo();
	evals = chkpt.rd_evals();
	
	VALUE arr = create_array(count, evals);
	psi::Chkpt::free(evals);
	
	return arr;
}

VALUE Task::rb_chkpt_alpha_evals_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	double *evals;
	int count;
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	count = chkpt.rd_nmo();
	evals = chkpt.rd_alpha_evals();
	
	VALUE arr = create_array(count, evals);
	psi::Chkpt::free(evals);
	
	return arr;
}

VALUE Task::rb_chkpt_beta_evals_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	double *evals;
	int count;
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	count = chkpt.rd_nmo();
	evals = chkpt.rd_beta_evals();
	
	VALUE arr = create_array(count, evals);
	psi::Chkpt::free(evals);
	
	return arr;
}

VALUE Task::rb_chkpt_evals_set(VALUE self, VALUE arr)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	double *evals;
	int count, expectedCount;
	
	// Convert the arr to a C-array. The function call handles type checking
	count = create_array(arr, &evals);
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	expectedCount = chkpt.rd_nmo();
	
	if (count != expectedCount) {
		rb_raise(rb_eArgError, "array is of length %d, expecting %d", count, expectedCount);
		return self;
	}
	
	// Save a new evals
	chkpt.wt_evals(evals);
	
	// Free memory
	psi::Chkpt::free(evals);
	
	return self;
}

VALUE Task::rb_chkpt_alpha_evals_set(VALUE self, VALUE arr)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	double *evals;
	int count, expectedCount;
	
	// Convert the arr to a C-array. The function call handles type checking
	count = create_array(arr, &evals);
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	expectedCount = chkpt.rd_nmo();
	
	if (count != expectedCount) {
		rb_raise(rb_eArgError, "array is of length %d, expecting %d", count, expectedCount);
		return self;
	}
	
	// Save a new evals
	chkpt.wt_alpha_evals(evals);
	
	// Free memory
	psi::Chkpt::free(evals);
	
	return self;
}

VALUE Task::rb_chkpt_beta_evals_set(VALUE self, VALUE arr)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	double *evals;
	int count, expectedCount;

	// Convert the arr to a C-array. The function call handles type checking
	count = create_array(arr, &evals);

	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	expectedCount = chkpt.rd_nmo();

	if (count != expectedCount) {
		rb_raise(rb_eArgError, "array is of length %d, expecting %d", count, expectedCount);
		return self;
	}

	// Save a new evals
	chkpt.wt_beta_evals(evals);

	// Free memory
	psi::Chkpt::free(evals);

	return self;
}

VALUE Task::rb_chkpt_exps_get(VALUE self)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	double *exps;
	int count;
	
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	count = chkpt.rd_nprim();
	exps = chkpt.rd_exps();
	
	VALUE arr = create_array(count, exps);
	psi::Chkpt::free(exps);
	
	return arr;
}

VALUE Task::rb_chkpt_exps_set(VALUE self, VALUE arr)
{
	Task *task;
	Data_Get_Struct(self, Task, task);
	double *exps;
	int count, expectedCount;

	// Convert the arr to a C-array. The function call handles type checking
	count = create_array(arr, &exps);

	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	expectedCount = chkpt.rd_nprim();

	if (count != expectedCount) {
		rb_raise(rb_eArgError, "array is of length %d, expecting %d", count, expectedCount);
		return self;
	}

	// Save a new evals
	chkpt.wt_exps(exps);

	// Free memory
	psi::Chkpt::free(exps);

	return self;
}

VALUE Task::rb_load_zmat(VALUE self)
{
	Task *task;
	ZEntry z;
	VALUE zmat;
	
	Data_Get_Struct(self, Task, task);
	z.attach_to(task);
	zmat = z.to_a();
	
	return zmat;
}

VALUE Task::rb_load_cartesian(VALUE self)
{
	Task *task;
	VALUE xyz;
	char **felements;
	double **fgeom;
	int nallatoms;
	int i;
	
	Data_Get_Struct(self, Task, task);
	// Gain access to checkpoint
	Chkpt chkpt(&task->m_psiPSIO, PSIO_OPEN_OLD);
	felements = chkpt.rd_felement();
	fgeom = chkpt.rd_fgeom();
	nallatoms = chkpt.rd_nallatom();
	
	// Go through and construct the cartesian geometry matrix for psirb
	xyz = rb_ary_new();
	for (i=0; i<nallatoms; ++i) {
		VALUE entry = rb_ary_new();
		
		// Add the element name to the entry
		rb_ary_push(entry, rb_cvar_get(rb_cObject, rb_intern(felements[i])));
		
		// Push the x, y, z values
		rb_ary_push(entry, rb_float_new(fgeom[i][0]));
		rb_ary_push(entry, rb_float_new(fgeom[i][1]));
		rb_ary_push(entry, rb_float_new(fgeom[i][2]));
		
		// Push the entry to the main array
		rb_ary_push(xyz, entry);
	} 
	
	// Clean up
	psi::Chkpt::free(fgeom);
	psi::Chkpt::free(felements);
	
	return xyz;
}

//VALUE Task::

}} // namespace psi::psirb
