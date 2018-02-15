/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_liboptions_liboptions_hpp
#define _psi_src_lib_liboptions_liboptions_hpp

#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif

#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/libpsi4util.h" // Needed for Ref counting, string splitting, and conversions

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace psi {

class DataTypeException : public PsiException {
   public:
    DataTypeException(const std::string& message) : PSIEXCEPTION(message) {}
};

class IndexException : public PsiException {
   public:
    IndexException(const std::string& message)
        : PSIEXCEPTION(message + " is not a valid option.") {}
    IndexException(const std::string& message, const std::string& module)
        : PSIEXCEPTION(message + " is not a valid option for module " + module) {}
};

class DuplicateKeyException : public PsiException {
   public:
    DuplicateKeyException(const std::string& key, const std::string& type1,
                          const std::string& type2, const char* lfile, int lline)
        : PsiException("Option " + key + " has been declared as a " + type1 + " and a " + type2,
                       lfile, lline) {}
};

class OptionsException : public PsiException {
   public:
    OptionsException(const std::string& message) : PSIEXCEPTION("Options Exception: " + message) {}
};

class Data;
class DataType {
    bool changed_;

   public:
    DataType();
    virtual ~DataType();

    bool has_changed() const;
    void changed();
    void dechanged();

    void to_upper(std::string& str);

    virtual void add_choices(std::string str);
    virtual std::string type() const;

    virtual bool is_array() const;
    virtual size_t size() const;
    virtual void add(DataType*);
    virtual void add(std::string, DataType*);
    virtual void add(bool);
    virtual void add(int);
    virtual void add(double);
    virtual void add(std::string, bool);
    virtual void add(std::string, std::string);
    virtual void add(std::string, int);
    virtual void add(std::string, double);
    virtual void add(std::string, std::string, std::string);

    virtual bool exists(std::string);

    virtual std::string to_string() const;
    virtual int to_integer() const;
    virtual double to_double() const;

    virtual void assign(DataType*);
    virtual void assign(bool);
    virtual void assign(int);
    virtual void assign(double);
    virtual void assign(std::string);

    virtual void reset();

    virtual Data& operator[](std::string);
    virtual Data& operator[](size_t);
};

#ifdef __INTEL_COMPILER
#pragma warning disable 654
#endif
class BooleanDataType : public DataType {
    bool boolean_;

   public:
    BooleanDataType();
    BooleanDataType(bool b);
    virtual ~BooleanDataType();

    virtual std::string type() const;

    virtual std::string to_string() const;
    virtual int to_integer() const;
    virtual double to_double() const;

    virtual void assign(bool b);
    virtual void assign(int i);
    virtual void assign(double d);
    virtual void assign(std::string s);
};

#ifdef __INTEL_COMPILER
#pragma warning disable 654
#endif
class IntDataType : public DataType {
    int integer_;

   public:
    IntDataType();
    IntDataType(int i);
    virtual ~IntDataType();

    virtual std::string type() const;

    virtual std::string to_string() const;
    virtual int to_integer() const;
    virtual double to_double() const;

    virtual void assign(bool b);
    virtual void assign(int i);
    virtual void assign(double d);
    virtual void assign(std::string s);
};

#ifdef __INTEL_COMPILER
#pragma warning disable 654
#endif
class DoubleDataType : public DataType {
    double double_;

   public:
    DoubleDataType();
    DoubleDataType(double d);
    virtual ~DoubleDataType();

    virtual std::string type() const;

    virtual std::string to_string() const;
    virtual int to_integer() const;
    virtual double to_double() const;

    virtual void assign(bool b);
    virtual void assign(int i);
    virtual void assign(double d);
    virtual void assign(std::string s);
};

#ifdef __INTEL_COMPILER
#pragma warning disable 654
#endif
class StringDataType : public DataType {
    std::string str_;
    std::vector<std::string> choices_;

   public:
    StringDataType();
    StringDataType(std::string s);
    StringDataType(std::string s, std::string c);

    virtual ~StringDataType();

    virtual void add_choices(std::string str);

    virtual std::string type() const;

    virtual std::string to_string() const;
    virtual int to_integer() const;
    virtual double to_double() const;

    virtual void assign(bool b);
    virtual void assign(int i);
    virtual void assign(double d);
    virtual void assign(std::string s);
};

#ifdef __INTEL_COMPILER
#pragma warning disable 654
#endif
class IStringDataType : public DataType {
    std::string str_;
    std::vector<std::string> choices_;

   public:
    IStringDataType();
    IStringDataType(std::string s);
    IStringDataType(std::string s, std::string c);
    virtual ~IStringDataType();

    virtual void add_choices(std::string str);

    virtual std::string type() const;

    virtual std::string to_string() const;
    virtual int to_integer() const;
    virtual double to_double() const;

    virtual void assign(bool b);
    virtual void assign(int i);
    virtual void assign(double d);
    virtual void assign(std::string s);
};

class Data {
    std::shared_ptr<DataType> ptr_;

   public:
    Data();
    Data(DataType* t);
    Data(const Data& copy);

    std::string to_string() const;
    int to_integer() const;
    double to_double() const;

    bool is_array() const;
    size_t size() const;

    bool has_changed() const;

    void changed();
    void dechanged();

    void add_choices(std::string str);

    std::string type() const;

    void add(DataType* data);
    void add(std::string s, DataType* data);
    void add(bool b);
    void add(int i);
    void add(double d);
    void add(std::string s, std::string c);
    void add(std::string key, bool b);
    void add(std::string key, int i);
    void add(std::string key, double d);
    void add(std::string key, std::string s, std::string c);

    void assign(DataType* data);
    void assign(bool b);
    void assign(int i);
    void assign(double d);
    void assign(std::string s);

    void reset();

    DataType* get() const;

    Data& operator[](int i);
    Data& operator[](std::string s);
};

#ifdef __INTEL_COMPILER
#pragma warning disable 654
#endif
class ArrayType : public DataType {
    std::vector<Data> array_;

   public:
    ArrayType();

    virtual std::string type() const;

    virtual void add(DataType* data);
    virtual void add(bool b);
    virtual void add(int i);
    virtual void add(double d);
    virtual void add(std::string s, std::string c = "");
    virtual void assign(DataType* data);

    virtual Data& operator[](size_t i);
    virtual Data& operator[](std::string s);
    virtual bool is_array() const;

    virtual size_t size() const;

    virtual std::string to_string() const;

    virtual void reset();
    std::vector<Data> data() { return array_; }
};

#ifdef __INTEL_COMPILER
#pragma warning disable 654
#endif
class MapType : public DataType {
    std::map<std::string, Data> keyvals_;
    typedef std::map<std::string, Data>::iterator iterator;
    typedef std::map<std::string, Data>::const_iterator const_iterator;

   public:
    MapType();

    virtual std::string type() const;

    virtual void add(std::string key, DataType* data);
    virtual void add(std::string key, bool b);
    virtual void add(std::string key, int i);
    virtual void add(std::string key, double d);
    virtual void add(std::string key, std::string s, std::string c = "");

    virtual bool exists(std::string key);

    virtual Data& operator[](std::string s);
    virtual bool is_array() const;

    virtual size_t size() const;

    virtual std::string to_string() const;
};

class PSI_API Options {
    bool edit_globals_;

    /// A temporary map used for validation of local options
    std::map<std::string, Data> all_local_options_;
    /// The module that's active right now
    std::string current_module_;

    /// "Active" set of options
    std::map<std::string, std::map<std::string, Data> > locals_;

    /// "Global" set of options
    std::map<std::string, Data> globals_;

    typedef std::map<std::string, Data>::iterator iterator;
    typedef std::map<std::string, Data>::const_iterator const_iterator;
    typedef std::map<std::string, std::map<std::string, Data> >::const_iterator const_mod_iterator;

   public:
    Options();

    Options& operator=(const Options& rhs);
    bool read_globals() const;
    void set_read_globals(bool _b);
    void set_current_module(const std::string s);
    std::string get_current_module() const { return current_module_; }

    void to_upper(std::string& str);

    void validate_options();

    void add(std::string key, DataType* data);
    void add(std::string key, bool b);
    void add(std::string key, int i);
    void add(std::string key, double d);
    void add(std::string key, std::string s, std::string c = "");
    void add_i(std::string key, std::string s, std::string c = "");

    void add_bool(std::string key, bool b);
    void add_int(std::string key, int i);
    void add_double(std::string key, double d);
    void add_str(std::string key, std::string s, std::string c = "");
    void add_str_i(std::string key, std::string s, std::string c = "");
    void add_array(std::string key);
    void set_bool(const std::string& module, const std::string& key, bool b);
    void set_int(const std::string& module, const std::string& key, int i);
    void set_double(const std::string& module, const std::string& key, double d);
    void set_str(const std::string& module, const std::string& key, std::string s);
    void set_str_i(const std::string& module, const std::string& key, std::string s);
    void set_array(const std::string& module, const std::string& key);

    void set_global_bool(const std::string& key, bool b);
    void set_global_int(const std::string& key, int i);
    void set_global_double(const std::string& key, double d);
    void set_global_str(const std::string& key, const std::string& s);
    void set_global_array(const std::string& key);

    DataType* set_global_array_entry(const std::string& key, DataType* entry, DataType* loc);
    void set_global_array_double(std::string key, double val, DataType* entry);
    void set_global_array_string(std::string key, std::string val, DataType* entry);
    void set_global_array_int(std::string key, int val, DataType* entry);
    DataType* set_global_array_array(std::string key, DataType* entry);
    DataType* set_local_array_entry(const std::string& module, const std::string& key,
                                    DataType* entry, DataType* loc);
    void set_local_array_double(const std::string& module, const std::string& key, double val,
                                DataType* entry);
    void set_local_array_string(const std::string& module, const std::string& key, std::string val,
                                DataType* entry);
    void set_local_array_int(const std::string& module, const std::string& key, int val,
                             DataType* entry);
    DataType* set_local_array_array(const std::string& module, const std::string& key,
                                    DataType* entry);

    void clear(void);

    bool exists_in_active(std::string key);

    bool exists_in_global(std::string key);
    bool exists(std::string key);

    Data& get(std::string key);

    Data& get(std::map<std::string, Data>& m, std::string& key);

    Data& get_global(std::string key);

    Data& get_local(std::string& key);

    Data& use(std::string& key);

    Data& use_local(std::string& key);

    bool get_bool(std::string key);
    int get_int(std::string key);
    double get_double(std::string key);
    std::string get_str(std::string key);
    int* get_int_array(std::string key);
    void fill_int_array(std::string key, int* empty_array);
    std::vector<int> get_int_vector(std::string key);
    double* get_double_array(std::string key);

    std::vector<double> get_double_vector(std::string key);

    const char* get_cstr(std::string key);

    Data& operator[](std::string key);

    std::string to_string() const;
    std::string globals_to_string() const;

    void print();
    void print_globals();
    std::vector<std::string> list_globals();
};
}
#endif
