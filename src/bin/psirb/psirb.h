#ifndef __PSIRB_H__
#define __PSIRB_H__

#ifdef MAIN
#define EXT
#else
#define EXT extern
#endif

#include <cstdio>
#include <string>
#include <ruby.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

namespace psi { namespace psirb {
	
#define PSI_VERSION_MAJOR 3
#define PSI_VERSION_MINOR 3

//
// Useful functions
extern VALUE create_array(unsigned int count, int *array);
extern VALUE create_array(unsigned int count, double *array);
extern int create_array(VALUE arr, int **array);
extern int create_array(VALUE arr, double **array);

/*! Namespace for containing global variables */
namespace Globals {
/*! All output should be sent to this variable */
	EXT FILE* g_fOutput;

/*! The name of the input file, could be useful for something */
	EXT std::string g_szInputFile;
	
/*! Name of the output filename */
	EXT std::string g_szOutputFilename;
	
/*! Verbosity */
	EXT bool g_bVerbose;
		
#ifdef MAIN
/*! All classes that provide a Ruby interface need this to say which module they belong to */
	EXT VALUE g_rbPsi = Qnil;
/*! How much to indent the Ruby puts output by. */
	EXT int g_iPutsIndent = 0;
	EXT bool g_bQuietRuby = false;
	EXT bool g_bIRB = false;
#else
	EXT VALUE g_rbPsi;
	EXT int g_iPutsIndent;
	EXT bool g_bQuietRuby;
	EXT bool g_bIRB;
#endif
};

// Helper macros
/*! Cast a C/++ function to a Ruby acceptable form */
#define RUBYCAST(x) (VALUE (*)(...))(x)

/*! Help the user get the Psi object from Ruby. Psi object MUST have an assign function */
#define RUBYPSIDATA(RObj, OType, OPtr) \
	Data_Get_Object(RObj, OType, OPtr);

/*! A calculation will now be grouped into a Task object. Tasks provide a 
	simple way of having several prefixes in a single input file.
	There is a global Task object that is used by default. */
class Task {
	/*! The libpsio object for this task */
	psi::PSIO m_psiPSIO;
	
	/*! Equivalent to the psi_file_prefix */
	std::string m_szFilePrefix;

	/*! Default scratch space that Psi is supposed to use. */
	std::string m_szScratchPath;
	
	/*! Ruby reference to the Task class descriptor */
	static VALUE m_rbTask;
	
public:
	/*! Default constructor; sets prefix to "psi" and scratch to "/tmp/" */
	Task(std::string prefix = "psi", std::string scratch = "/tmp/");
	
	/*! Destructor */
	~Task();
	
	/*! Accessor function get for the prefix 
		\return current prefix */
	const std::string& prefix();
	/*! Accessor function put for the prefix
		\param new_prefix what to change the prefix to */
	void prefix(std::string new_prefix);
	
	/*! Accessor function get for the scratch location
		\return current location */
	const std::string& scratch();
	/*! Accessor function put for the scratch location
		\param new_scratch what to change the scratch to */
	void scratch(std::string new_scratch);
	
	/*! Accessor get function for the libpsio object associated with this Task
		\return libpsio object pointer */
	psi::PSIO *libpsio();

	//
	// Ruby framework for Task
	//
	
	/*! Creates the Ruby class framework */
	static void create_ruby_class();

	/*! Called by Ruby when it needs to delete a class. */
	static void rb_free(void *p);
	
	/*! Called by Ruby during object creation */
	static VALUE rb_alloc(VALUE klass);
	
	/*! Called by Ruby during object creation */
	static VALUE rb_init(int argc, VALUE* argv, VALUE self);
	
	/*! Called by Ruby during object copy creation */
	static VALUE rb_init_copy(VALUE copy, VALUE orig);
	
	/*! Called by Ruby if the user tries to print a Task object */
	static VALUE rb_to_s(VALUE self);
	
	/*! Ruby function: Task.prefix= */
	static VALUE rb_prefix_set(VALUE self, VALUE newPrefix);
	/*! Ruby function: Task.prefix */
	static VALUE rb_prefix_get(VALUE self);

	/*! Ruby function: Task.scratch= */
	static VALUE rb_scratch_set(VALUE self, VALUE newsScratch);
	/*! Ruby function: Task.scratch */
	static VALUE rb_scratch_get(VALUE self);
	
	static VALUE rb_load_zmat(VALUE self);
	static VALUE rb_load_cartesian(VALUE self);
	
	//
	// Libpsio++ interface
	static VALUE rb_print_toc(VALUE, VALUE);
	
	//
	// Checkpoint interface
	static VALUE rb_chkpt_exist(VALUE, VALUE);
	static VALUE rb_chkpt_label_get(VALUE);
	static VALUE rb_chkpt_label_set(VALUE self, VALUE label);
	static VALUE rb_chkpt_escf_get(VALUE);
	static VALUE rb_chkpt_escf_set(VALUE, VALUE);
	static VALUE rb_chkpt_eref_get(VALUE);
	static VALUE rb_chkpt_eref_set(VALUE, VALUE);
	static VALUE rb_chkpt_ecorr_get(VALUE);
	static VALUE rb_chkpt_ecorr_set(VALUE, VALUE);
	static VALUE rb_chkpt_enuc_get(VALUE);
	static VALUE rb_chkpt_enuc_set(VALUE, VALUE);
	static VALUE rb_chkpt_efzc_get(VALUE);
	static VALUE rb_chkpt_efzc_set(VALUE, VALUE);
	static VALUE rb_chkpt_etot_get(VALUE);
	static VALUE rb_chkpt_etot_set(VALUE, VALUE);
	static VALUE rb_chkpt_disp_get(VALUE);
	static VALUE rb_chkpt_disp_set(VALUE, VALUE);
	static VALUE rb_chkpt_eccsd_get(VALUE);
	static VALUE rb_chkpt_e_t_get(VALUE);
	static VALUE rb_chkpt_emp2_get(VALUE);	
	static VALUE rb_chkpt_eom_state_energies_get(VALUE);
	static VALUE rb_chkpt_num_irreps_get(VALUE);
	static VALUE rb_chkpt_clsdpi_set(VALUE, VALUE);
	static VALUE rb_chkpt_clsdpi_get(VALUE);
	static VALUE rb_chkpt_frzcpi_set(VALUE, VALUE);
	static VALUE rb_chkpt_frzcpi_get(VALUE);
	static VALUE rb_chkpt_frzvpi_set(VALUE, VALUE);
	static VALUE rb_chkpt_frzvpi_get(VALUE);
	static VALUE rb_chkpt_evals_get(VALUE);
	static VALUE rb_chkpt_alpha_evals_get(VALUE);
	static VALUE rb_chkpt_beta_evals_get(VALUE);
	static VALUE rb_chkpt_evals_set(VALUE, VALUE);
	static VALUE rb_chkpt_alpha_evals_set(VALUE, VALUE);
	static VALUE rb_chkpt_beta_evals_set(VALUE, VALUE);
	static VALUE rb_chkpt_exps_set(VALUE, VALUE);
	static VALUE rb_chkpt_exps_get(VALUE);
};

//
// Used to read in the z-matrix from checkpoint file.
class ZEntry {
private:
	//! Z-Matrix entries read from checkpoint
	z_entry *m_zEntry;

	//! How many z_entry are there in m_zEntry
	unsigned int m_numZEntries;
	
	//! Full geometry elements (includes dummies)
	char **m_szFElement;
	
	//! Ruby reference to the Z-Matrix class descriptor
	static VALUE m_rbZEntry;
	
	//! Task to use for libchkpt interface.
	Task *m_pTask;
	
public:
	//! Default constructor
	ZEntry();
	~ZEntry();
	
	//! Attach to a specific Task to gain access to a checkpoint file.
	void attach_to(Task* task);
	
	//! Read in the z_entry from checkpoint. Must be attached to a Task first.
	void read();
	
	VALUE to_a();
	
	//
	// Ruby framework for ZEntry
	//
	
	//! Creates the Ruby class framework
	static void create_ruby_class();
		
	/*! Called by Ruby when it needs to delete a class. */
	static void rb_free(void *p);
	
	/*! Called by Ruby during object creation */
	static VALUE rb_alloc(VALUE klass);
	
	/*! Called by Ruby during object creation */
	static VALUE rb_init(VALUE self, VALUE arg);
	
	/*! Called by Ruby during object copy creation */
	static VALUE rb_init_copy(VALUE copy, VALUE orig);
	
	/*! Called by Ruby if the user try to print a Task object */
	static VALUE rb_to_s(VALUE self);

	/*! Called by Ruby if the user try to convert to an array */
	static VALUE rb_to_a(VALUE self);
};

//
// Ruby framework for handling BLAS compatible matrices
class Matrix {
	double **m_pMatrix;
	size_t m_nRows, m_nCols;

	//! Ruby reference to the Z-Matrix class descriptor
	static VALUE m_rbMatrix;
	
public:
	Matrix();
	~Matrix();
	
	bool allocate(size_t rows, size_t cols);
	void release();
	void copy(Matrix *);
	
	void set(size_t x, size_t y, double value) {
		m_pMatrix[x][y] = value;
	}
	double get(size_t x, size_t y) {
		return m_pMatrix[x][y];
	}
	
	//
	// Ruby framework for Matrix
	//
	//! Creates the Ruby class framework
	static void create_ruby_class();
	
	//! Called by Ruby when it needs to delete a class
	static void rb_free(void *p);
	
	//! Called by Ruby during object creation
	static VALUE rb_alloc(VALUE klass);
	
	//! Called by Ruby during object creation
	static VALUE rb_init(int argc, VALUE *argv, VALUE self);
	
	//! Called by Ruby during object creation
	static VALUE rb_init_copy(VALUE copy, VALUE orig);
	
	//! Accessor methods
	static VALUE rb_element_get(VALUE self, VALUE i, VALUE j);
	static VALUE rb_element_set(VALUE self, VALUE i, VALUE j, VALUE val);
	
	//! Conversion routines
	static VALUE rb_to_s(VALUE self);
};

}} // namespace psi::psirb

#endif // __PSIRB_H__
