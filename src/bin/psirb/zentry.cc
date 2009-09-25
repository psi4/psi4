#include <ruby.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "psirb.h"

namespace psi { namespace psirb {
	
using namespace psi;

//
// class ZEntry
VALUE ZEntry::m_rbZEntry = Qnil;

ZEntry::ZEntry() : m_zEntry(NULL), m_pTask(NULL), m_numZEntries(0)
{
}

ZEntry::~ZEntry()
{
	if (m_zEntry)
		free(m_zEntry);
	if (m_szFElement) {
		for (int i = 0; i < m_numZEntries; ++i)
			free(m_szFElement[i]);
		free(m_szFElement);
	}
	m_szFElement = NULL;
	m_zEntry = NULL;
	m_pTask = NULL;
}

void ZEntry::create_ruby_class()
{
	// Create the ZEntry class under the Psi module
	ZEntry::m_rbZEntry = rb_define_class_under(Globals::g_rbPsi, "ZEntry", rb_cObject);
	
	// Register the allocation function with Ruby
	rb_define_alloc_func(ZEntry::m_rbZEntry, ZEntry::rb_alloc);
	
	// Register the initialization function with Ruby
	rb_define_method(ZEntry::m_rbZEntry, "initialize", RUBYCAST(ZEntry::rb_init), 1);
	rb_define_method(ZEntry::m_rbZEntry, "initialize_copy", RUBYCAST(ZEntry::rb_init_copy), 1);
	
	// Register the to_s function
	rb_define_method(ZEntry::m_rbZEntry, "to_s", RUBYCAST(ZEntry::rb_to_s), 0);
	// Register the to_a function
	rb_define_method(ZEntry::m_rbZEntry, "to_a", RUBYCAST(ZEntry::rb_to_a), 0);
}

void ZEntry::rb_free(void *p)
{
	ZEntry *pZEntry = (ZEntry*)p;
	if (pZEntry)
		delete pZEntry;
}

VALUE ZEntry::rb_alloc(VALUE klass)
{
	ZEntry *newZEntry = new ZEntry;
	VALUE newObj;
	
	// Wrap the newly created ZEntry inside a Ruby object
	newObj = Data_Wrap_Struct(klass, 0, ZEntry::rb_free, newZEntry);
	// Ruby is now responsible for the object, not us.
	
	// Return the new object.
	return newObj;
}

VALUE ZEntry::rb_init(VALUE self, VALUE arg)
{
	// We are expecting a Task object in arg
	Task *task;
	ZEntry *zentry;
	
	// Get the objects
	Data_Get_Struct(self, ZEntry, zentry);
	
	// Make sure the user passed a Task object
	if (TYPE(arg) != T_DATA || RDATA(arg)->dfree != (RUBY_DATA_FUNC)Task::rb_free) {
		rb_raise(rb_eTypeError, "wrong argument type, expecting a Task object");
		return self;
	}
	Data_Get_Struct(arg, Task, task);
	
	// Get the libpsio object from task
	zentry->attach_to(task);
	
	// At this point the information has been read in, the attach_to method triggers a read call.
	return self;
}

VALUE ZEntry::rb_init_copy(VALUE copy, VALUE orig)
{
	ZEntry *o, *c;
	
	// Do not self copy
	if (copy == orig)
		return copy;
		
	// We can only copy a ZEntry to a ZEntry.
	if (TYPE(orig) != T_DATA || RDATA(orig)->dfree != (RUBY_DATA_FUNC)ZEntry::rb_free) {
		rb_raise(rb_eTypeError, "wrong argument type");
		return copy;
	}
	
	// Get the objects and copy
	Data_Get_Struct(copy, ZEntry, c);
	Data_Get_Struct(orig, ZEntry, o);
	
	// Attach the copy to the original's Task
	c->attach_to(o->m_pTask);
	
	return copy;
}

void ZEntry::attach_to(Task* task)
{
	m_pTask = task;
	read();
}

void ZEntry::read()
{
	if (m_pTask == NULL)
		return;
	
	// If m_zEntry already contains data free it
	if (m_zEntry)
		free(m_zEntry);
	m_numZEntries = 0;
	
	// Get the libpsio object associated with the assigned task
	psi::PSIO *libpsio = m_pTask->libpsio();
	
	// Open checkpoint
	psi::Chkpt chkpt(libpsio, PSIO_OPEN_OLD);
	
	// How many are we going to get?
	m_numZEntries = chkpt.rd_nallatom();
	
	// Read the full element names
	m_szFElement = chkpt.rd_felement();
  
	// Read in the entries
	m_zEntry = chkpt.rd_zmat();
}

VALUE ZEntry::rb_to_s(VALUE self)
{
	ZEntry *s;
	z_entry *zmat;
	
	Data_Get_Struct(self, ZEntry, s);
	std::string text;
	std::stringstream stream;
	
	zmat = s->m_zEntry;
	
	text += "zmat = [\n";
	for (int i=0; i < s->m_numZEntries; ++i) {
		text += "  [ ";
		text += s->m_szFElement[i];
		stream.str(std::string());
		
		if (i > 0) {
			stream << " " << zmat[i].bond_atom;
			stream << " " << std::setprecision(5) << std::setw(10) << zmat[i].bond_val;
		}
		if (i > 1) {
			stream << " " << zmat[i].angle_atom;
			stream << " " << std::setprecision(5) << std::setw(10) << zmat[i].angle_val;
		}
		if (i > 2) {
			stream << " " << zmat[i].tors_atom;
			stream << " " << std::setprecision(5) << std::setw(10) << zmat[i].tors_val;
		}
		text += stream.str();
		text += " ]\n";
	}
	text += "]\n";
	
	return rb_str_new2(text.c_str());
}

VALUE ZEntry::to_a()
{
	VALUE final_array = rb_ary_new();
	
	for (int i=0; i < m_numZEntries; ++i) {
		VALUE entry = rb_ary_new();
		
		// Add the element name to the entry
		rb_ary_push(entry, rb_cvar_get(rb_cObject, rb_intern(m_szFElement[i])));
		
		if (i > 0) {
			rb_ary_push(entry, INT2FIX(m_zEntry[i].bond_atom));
			rb_ary_push(entry, rb_float_new(m_zEntry[i].bond_val));
		}
		if (i > 1) {
			rb_ary_push(entry, INT2FIX(m_zEntry[i].angle_atom));
			rb_ary_push(entry, rb_float_new(m_zEntry[i].angle_val));
		}
		if (i > 2) {
			rb_ary_push(entry, INT2FIX(m_zEntry[i].tors_atom));
			rb_ary_push(entry, rb_float_new(m_zEntry[i].tors_val));
		}
		rb_ary_push(final_array, entry);
	}
	
	return final_array;
}

VALUE ZEntry::rb_to_a(VALUE self)
{
	ZEntry *s;
	
	Data_Get_Struct(self, ZEntry, s);
	VALUE final_array = s->to_a();
	
	return final_array;
}

}} // namespace psi::psirb
