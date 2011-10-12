#include <iostream>
#include <vector>
#include <map>
#include <cstddef>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <assert.h>

#include <exception.h>
#include <libutil/libutil.h> // Needed for Ref counting, string splitting, and conversions
#include <libutil/ref.h> // Needed for Ref counting, string splitting, and conversions
#include <boost/shared_ptr.hpp>

#include "liboptions.h"

namespace psi {

// DataType base
DataType::DataType()
    : changed_(false)
{ }

DataType::~DataType()
{ }

bool DataType::has_changed() const
{
    return changed_;
}

void DataType::changed()
{
    changed_ = true;
}

void DataType::to_upper(std::string& str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

void DataType::add_choices(std::string str)
{
    throw NotImplementedException("add(bool)");
}

std::string DataType::type() const
{
    return std::string("unknown");
}

bool DataType::is_array() const
{
    return false;
}

unsigned int DataType::size() const
{
    throw NotImplementedException("size()");
}

void DataType::add(DataType *)
{
    throw NotImplementedException("add(DataType*)");
}

void DataType::add(std::string, DataType*)
{
    throw NotImplementedException("add(std::string, DataType*)");
}

void DataType::add(bool)
{
    throw NotImplementedException("add(bool)");
}

void DataType::add(int)
{
    throw NotImplementedException("add(int)");
}

void DataType::add(double)
{
    throw NotImplementedException("add(double)");
}

void DataType::add(std::string, bool)
{
    throw NotImplementedException("add(std::string, bool)");
}

void DataType::add(std::string, std::string)
{
    throw NotImplementedException("add(std::string, std::string)");
}

void DataType::add(std::string, int)
{
    throw NotImplementedException("add(std::string, int)");
}

void DataType::add(std::string, double)
{
    throw NotImplementedException("add(std::string, double)");
}

void DataType::add(std::string, std::string, std::string)
{
    throw NotImplementedException("add(std::string, std::string, std::string)");
}

bool DataType::exists(std::string)
{
    throw NotImplementedException("exists(std::string)");
}

std::string DataType::to_string() const
{
    throw DataTypeException("don't know how to convert to a string");
}

int DataType::to_integer() const
{
    throw DataTypeException("don't know how to convert to an integer");
}

double DataType::to_double() const
{
    throw DataTypeException("don't know how to convert to a double");
}

void DataType::assign(DataType*)
{
    throw DataTypeException("assign(DataType*) failure");
}

void DataType::assign(bool)
{
    throw DataTypeException("assign(bool) failure");
}

void DataType::assign(int)
{
    throw DataTypeException("assign(int) failure");
}

void DataType::assign(double)
{
    throw DataTypeException("assign(double) failure");
}

void DataType::assign(std::string)
{
    throw DataTypeException("assign(std:string) failure");
}

void DataType::reset()
{
    throw DataTypeException("reset() failure");
}

Data& DataType::operator[](std::string)
{
    throw NotImplementedException("Data& [string]");
}

Data& DataType::operator[](unsigned int)
{
    throw NotImplementedException("Data& [unsigned int]");
}

// BooleanDataType
BooleanDataType::BooleanDataType()
    : DataType(), boolean_(false)
{ }

BooleanDataType::BooleanDataType(bool b)
    : DataType(), boolean_(b)
{ }

BooleanDataType::~BooleanDataType()
{ }

std::string BooleanDataType::type() const
{
    return std::string("boolean");
}

std::string BooleanDataType::to_string() const
{
    std::string ret;
    if (boolean_)
        ret = "TRUE";
    else
        ret = "FALSE";
    return ret;
}

int BooleanDataType::to_integer() const
{
    return static_cast<int>(boolean_);
}

double BooleanDataType::to_double() const
{
    return static_cast<double>(boolean_);
}

void BooleanDataType::assign(bool b)
{
    changed();
    boolean_ = b;
}

void BooleanDataType::assign(int i)
{
    assign(static_cast<bool>(i));
}

void BooleanDataType::assign(double d)
{
    assign(static_cast<bool>(d));
}

void BooleanDataType::assign(std::string s)
{
    assign(static_cast<bool>(std::strtod(s.c_str(), NULL)));
}

// IntDataType
IntDataType::IntDataType()
    : DataType(), integer_(0)
{ }

IntDataType::IntDataType(int i)
    : DataType(), integer_(i)
{ }

IntDataType::~IntDataType()
{ }

std::string IntDataType::type() const
{
    return std::string("int");
}

std::string IntDataType::to_string() const
{
    std::stringstream strm;
    strm << integer_;
    return strm.str();
}

int IntDataType::to_integer() const
{
    return integer_;
}

double IntDataType::to_double() const
{
    return static_cast<double>(integer_);
}

void IntDataType::assign(bool b)
{
    assign(static_cast<int>(b));
}

void IntDataType::assign(int i)
{
    changed();
    integer_ = i;
}

void IntDataType::assign(double d)
{
    assign(static_cast<int>(d));
}

void IntDataType::assign(std::string s)
{
    assign(static_cast<int>(std::strtod(s.c_str(), NULL)));
}

// DoubleDataType
DoubleDataType::DoubleDataType()
    : DataType(), double_(0.0)
{ }

DoubleDataType::DoubleDataType(double d)
    : DataType(), double_(d)
{ }

DoubleDataType::~DoubleDataType()
{ }

std::string DoubleDataType::type() const
{
    return std::string("double");
}

std::string DoubleDataType::to_string() const
{
    std::stringstream strm;
    strm << double_;
    return strm.str();
}

int DoubleDataType::to_integer() const
{
    return static_cast<int>(double_);
}

double DoubleDataType::to_double() const
{
    return double_;
}

void DoubleDataType::assign(bool b)
{
    assign(static_cast<double>(b));
}

void DoubleDataType::assign(int i)
{
    assign(static_cast<double>(i));
}

void DoubleDataType::assign(double d)
{
    changed();
    double_ = d;
}

void DoubleDataType::assign(std::string s)
{
    assign(std::strtod(s.c_str(), NULL));
}

// StringDataType
StringDataType::StringDataType()
    : DataType(), str_(), choices_()
{ }

StringDataType::StringDataType(std::string s)
    : DataType(), str_(s), choices_()
{
    to_upper(str_);
}

StringDataType::StringDataType(std::string s, std::string c)
    : DataType(), str_(s), choices_()
{
    to_upper(str_);
    to_upper(c);
    choices_ = split(c);
}

StringDataType::~StringDataType()
{ }

void StringDataType::add_choices(std::string str)
{
    to_upper(str);
    std::vector<std::string> temp = split(str);
    for(int i = 0; i < temp.size(); ++i)
        choices_.push_back(temp[i]);
}

std::string StringDataType::type() const
{
    return std::string("string");
}

std::string StringDataType::to_string() const
{
    return str_;
}

int StringDataType::to_integer() const
{
    return static_cast<int>(std::strtod(str_.c_str(), NULL));
}

double StringDataType::to_double() const
{
    return std::strtod(str_.c_str(), NULL);
}

void StringDataType::assign(bool b)
{
    if (b)
        assign("TRUE");
    else
        assign("FALSE");
}

void StringDataType::assign(int i)
{
    std::stringstream strm;
    strm << i;
    assign(strm.str());
}

void StringDataType::assign(double d)
{
    std::stringstream strm;
    strm << d;
    assign(strm.str());
}

void StringDataType::assign(std::string s)
{
    to_upper(s);
    if (choices_.size() > 0) {
        bool wrong_input = true;
        for (unsigned int i=0; i<choices_.size(); ++i)
            if (s == choices_[i])
                wrong_input = false;
        if (wrong_input)
            throw DataTypeException(s + " is not a valid choice");
        changed();
        str_ = s;
    }
    else {
        changed();
        str_ = s;
    }
}

// IStringDataType
IStringDataType::IStringDataType()
    : DataType(), str_(), choices_()
{ }

IStringDataType::IStringDataType(std::string s)
    : DataType(), str_(s), choices_()
{ }

IStringDataType::IStringDataType(std::string s, std::string c)
    : DataType(), str_(s), choices_()
{
    choices_ = split(c);
}

IStringDataType::~IStringDataType()
{ }

void IStringDataType::add_choices(std::string str)
{
    std::vector<std::string> temp = split(str);
    for(int i = 0; i < temp.size(); ++i)
        choices_.push_back(temp[i]);
}

std::string IStringDataType::type() const
{
    return std::string("string");
}

std::string IStringDataType::to_string() const
{
    return str_;
}

int IStringDataType::to_integer() const
{
    return static_cast<int>(std::strtod(str_.c_str(), NULL));
}

double IStringDataType::to_double() const
{
    return std::strtod(str_.c_str(), NULL);
}

void IStringDataType::assign(bool b)
{
    if (b)
        assign("TRUE");
    else
        assign("FALSE");
}

void IStringDataType::assign(int i)
{
    std::stringstream strm;
    strm << i;
    assign(strm.str());
}

void IStringDataType::assign(double d)
{
    std::stringstream strm;
    strm << d;
    assign(strm.str());
}

void IStringDataType::assign(std::string s)
{
    if (choices_.size() > 0) {
        bool wrong_input = true;
        for (unsigned int i=0; i<choices_.size(); ++i)
            if (s == choices_[i])
                wrong_input = false;
        if (wrong_input)
            throw DataTypeException(s + " is not a valid choice");
        changed();
        str_ = s;
    }
    else {
        changed();
        str_ = s;
    }
}

// Data
Data::Data()
{ }

Data::Data(DataType *t)
    : ptr_(t)
{ }

Data::Data(const Data& copy)
{
    ptr_ = copy.ptr_;
}

std::string Data::to_string() const
{
    return ptr_->to_string();
}

int Data::to_integer() const
{
    return ptr_->to_integer();
}

double Data::to_double() const
{
    return ptr_->to_double();
}

bool Data::is_array() const
{
    return ptr_->is_array();
}

unsigned int Data::size() const
{
    return ptr_->size();
}

bool Data::has_changed() const
{
    return ptr_->has_changed();
}

void Data::changed()
{
    ptr_->changed();
}

void Data::add_choices(std::string str)
{
    ptr_->add_choices(str);
}

std::string Data::type() const
{
    return ptr_->type();
}

void Data::add(DataType *data)
{
    ptr_->add(data);
}

void Data::add(std::string s, DataType *data)
{
    ptr_->add(s, data);
}

void Data::add(bool b)
{
    ptr_->add(b);
}

void Data::add(int i)
{
    ptr_->add(i);
}

void Data::add(double d)
{
    ptr_->add(d);
}

void Data::add(std::string s, std::string c)
{
    ptr_->add(s, c);
}

void Data::add(std::string key, bool b)
{
    ptr_->add(key, b);
}

void Data::add(std::string key, int i)
{
    ptr_->add(key, i);
}

void Data::add(std::string key, double d)
{
    ptr_->add(key, d);
}

void Data::add(std::string key, std::string s, std::string c)
{
    ptr_->add(key, s, c);
}

void Data::assign(DataType *data)
{
    ptr_->assign(data);
}

void Data::assign(bool b)
{
    ptr_->assign(b);
}

void Data::assign(int i)
{
    ptr_->assign(i);
}

void Data::assign(double d)
{
    ptr_->assign(d);
}

void Data::assign(std::string s)
{
    ptr_->assign(s);
}

void Data::reset()
{
    ptr_->reset();
}

DataType* Data::get() const
{
    return ptr_.pointer();
}

Data& Data::operator[](int i)
{
    return (*(ptr_.pointer()))[i];
}

Data& Data::operator[](std::string s)
{
    return (*(ptr_.pointer()))[s];
}

// ArrayType
ArrayType::ArrayType()
{ }

std::string ArrayType::type() const
{
    return std::string("array");
}

void ArrayType::add(DataType *data)
{
    array_.push_back(Data(data));
}

void ArrayType::add(bool b)
{
    add(new BooleanDataType(b));
}

void ArrayType::add(int i)
{
    add(new IntDataType(i));
}

void ArrayType::add(double d)
{
    add(new DoubleDataType(d));
}

void ArrayType::add(std::string s, std::string c)
{
    add(new StringDataType(s, c));
}

void ArrayType::assign(DataType* data)
{
    changed();
    array_.push_back(Data(data));
}

Data& ArrayType::operator[](unsigned int i)
{
    if (i >= array_.size())
        throw IndexException("out of range");
    changed();
    return array_[i];
}

Data& ArrayType::operator[](std::string s)
{
    unsigned int i = static_cast<unsigned int>(std::strtod(s.c_str(), NULL));
    if (i >= array_.size())
        throw IndexException("out of range");
    changed();
    return array_[i];
}

bool ArrayType::is_array() const
{
    return true;
}

unsigned int ArrayType::size() const
{
    return array_.size();
}

std::string ArrayType::to_string() const
{
    std::string str = "[ ";
    for (unsigned int i=0; i<array_.size(); ++i) {
        str += array_[i].to_string();
        if (i != array_.size() - 1)
            str += ", ";
    }
    str += " ]";
    return str;
}

void ArrayType::reset()
{
    array_.clear();
}

MapType::MapType()
{ }

std::string MapType::type() const
{
    return std::string("map");
}

void MapType::add(std::string key, DataType *data)
{
    to_upper(key);

    iterator pos = keyvals_.find(key);
    if (pos != keyvals_.end())
        throw DuplicateKeyException(key, data->type(), pos->second.type(), __FILE__, __LINE__);
    keyvals_[key] = Data(data);
}

void MapType::add(std::string key, bool b)
{
    add(key, new BooleanDataType(b));
}

void MapType::add(std::string key, int i)
{
    add(key, new IntDataType(i));
}

void MapType::add(std::string key, double d)
{
    add(key, new DoubleDataType(d));
}

void MapType::add(std::string key, std::string s, std::string c)
{
    add(key, new StringDataType(s, c));
}

bool MapType::exists(std::string key)
{
    to_upper(key);
    iterator pos = keyvals_.find(key);
    if (pos != keyvals_.end())
        return true;
    return false;
}

Data& MapType::operator[](std::string s)
{
    to_upper(s);
    if (!exists(s))
        throw IndexException(s);
    return keyvals_[s];
}

bool MapType::is_array() const
{
    return true;
}

unsigned int MapType::size() const
{
    return keyvals_.size();
}

std::string MapType::to_string() const
{
    std::string str = "{ ";
    for (const_iterator pos = keyvals_.begin(); pos != keyvals_.end(); ++pos) {
        str += pos->first + " => " + pos->second.to_string() + ", ";
    }
    str += "}";
    return str;
}

Options::Options()
    : edit_globals_(false)
{ }

Options& Options::operator=(const Options& rhs)
{
    // Don't self copy
    if (this == &rhs)
        return *this;

    locals_ = rhs.locals_;
    globals_ = rhs.globals_;

    return *this;
}

bool Options::read_globals() const
{
    return edit_globals_;
}

void Options::set_read_globals(bool _b)
{
    edit_globals_ = _b;
}

void Options::set_current_module(const std::string s)
{
    current_module_ = s;
    all_local_options_.clear();
}

void Options::to_upper(std::string& str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

void Options::validate_options()
{
    std::map<std::string, Data>::const_iterator iter = locals_[current_module_].begin();
    std::map<std::string, Data>::const_iterator stop = locals_[current_module_].end();
    std::map<std::string, Data>::const_iterator not_found = all_local_options_.end();
    for(; iter != stop; ++iter){
        if(iter->second.has_changed()){
            if(all_local_options_.find(iter->first) == not_found)
                throw PSIEXCEPTION("Option " + iter->first +
                                   " is not recognized by the " + current_module_ + " module.");
        }
    }
    all_local_options_.clear();
}

void Options::add(std::string key, DataType *data)
{
    to_upper(key);

    std::map<std::string, Data> & local = edit_globals_ ? globals_ : locals_[current_module_];

    Data val(data);
    all_local_options_[key] = val;

    // Make sure the key isn't already there
    iterator pos = local.find(key);
    if (pos != local.end()) { // If it is there, make sure they are the same type
        if (pos->second.type() != data->type())
            throw DuplicateKeyException(key, data->type(), pos->second.type(), __FILE__, __LINE__);
        return;
    }
    local[key] = val;
}

void Options::add(std::string key, bool b)
{
    add(key, new BooleanDataType(b));
}

void Options::add(std::string key, int i)
{
    add(key, new IntDataType(i));
}

void Options::add(std::string key, double d)
{
    add(key, new DoubleDataType(d));
}

void Options::add(std::string key, std::string s, std::string c)
{
    if(edit_globals_  && globals_.count(key)){
        globals_[key].add_choices(c);
    }
    else{
        add(key, new StringDataType(s, c));
    }
}

void Options::add_i(std::string key, std::string s, std::string c)
{
    if(edit_globals_  && globals_.count(key)){
        globals_[key].add_choices(c);
    }
    else{
        add(key, new IStringDataType(s, c));
    }
}

void Options::add_bool(std::string key, bool b)
{
    add(key,b);
}

void Options::add_int(std::string key, int i)
{
    add(key,i);
}

void Options::add_double(std::string key, double d)
{
    add(key,d);
}

void Options::add_str(std::string key, std::string s, std::string c)
{
    add(key,s,c);
}

void Options::add_str_i(std::string key, std::string s, std::string c)
{
    add_i(key,s,c);
}

void Options::add_array(std::string key)
{
    add(key, new ArrayType());
}

void Options::set_bool(const std::string &module, const std::string &key, bool b)
{
    locals_[module][key] = new BooleanDataType(b);
    locals_[module][key].changed();
}

void Options::set_int(const std::string &module, const std::string &key, int i)
{
    locals_[module][key] = new IntDataType(i);
    locals_[module][key].changed();
}

void Options::set_double(const std::string & module, const std::string &key, double d)
{
    locals_[module][key] = new DoubleDataType(d);
    locals_[module][key].changed();
}

void Options::set_str(const std::string & module, const std::string &key, std::string s)
{
    locals_[module][key] = new StringDataType(s);
    locals_[module][key].changed();
}

void Options::set_array(const std::string &module, const std::string& key)
{
    locals_[module][key] = Data(new ArrayType);
    locals_[module][key].changed();
}

void Options::set_global_bool(const std::string &key, bool b)
{
    get_global(key).assign(b);
}

void Options::set_global_int(const std::string &key, int i)
{
    get_global(key).assign(i);
}

void Options::set_global_double(const std::string &key, double d)
{
    get_global(key).assign(d);
}

void Options::set_global_str(const std::string &key, const std::string &s)
{
    get_global(key).assign(s);
}

void Options::set_global_array(const std::string& key)
{
    globals_[key] = Data(new ArrayType());
    globals_[key].changed();
}

DataType* Options::set_global_array_entry(const std::string& key, DataType* entry, DataType* loc)
{
    if(loc == NULL){
        // This is the first entry to be added
        Data &data = get_global(key);
        data.assign(entry);
    }
    else{
        // We're adding to an existing entry
        ArrayType* arrptr(dynamic_cast<ArrayType*>(loc));
        arrptr->assign(entry);
    }
    return entry;
}

void Options::set_global_array_double(std::string key, double val, DataType *entry)
{
    set_global_array_entry(key, new DoubleDataType(val), entry);
}

void Options::set_global_array_string(std::string key, std::string val, DataType *entry)
{
    set_global_array_entry(key, new StringDataType(val), entry);
}

void Options::set_global_array_int(std::string key, int val, DataType *entry)
{
    set_global_array_entry(key, new IntDataType(val), entry);
}

DataType* Options::set_global_array_array(std::string key, DataType *entry)
{
    return set_global_array_entry(key, new ArrayType(), entry);
}

DataType* Options::set_local_array_entry(const std::string &module, const std::string& key, DataType* entry, DataType* loc)
{
    if(loc == NULL){
        // This is the first entry to be added
        locals_[module][key].assign(entry);
    }
    else{
        // We're adding to an existing entry
        ArrayType* arrptr(dynamic_cast<ArrayType*>(loc));
        arrptr->assign(entry);
    }
    return entry;
}

void Options::set_local_array_double(const std::string &module, const std::string &key, double val, DataType *entry)
{
    set_local_array_entry(module, key, new DoubleDataType(val), entry);
}

void Options::set_local_array_string(const std::string &module, const std::string &key, std::string val, DataType *entry)
{
    set_local_array_entry(module, key, new StringDataType(val), entry);
}

void Options::set_local_array_int(const std::string &module, const std::string &key, int val, DataType *entry)
{
    set_local_array_entry(module, key, new IntDataType(val), entry);
}

DataType* Options::set_local_array_array(const std::string &module, const std::string &key, DataType *entry)
{
    return set_local_array_entry(module, key, new ArrayType(), entry);
}

void Options::clear(void)
{
    locals_.clear();
}

bool Options::exists_in_active(std::string key)
{
    to_upper(key);

    if(!locals_.count(current_module_)) return false;
    return (locals_[current_module_].count(key));
}

bool Options::exists_in_global(std::string key)
{
    to_upper(key);

    iterator pos = globals_.find(key);
    if (pos != globals_.end())
        return true;
    return false;
}

bool Options::exists(std::string key)
{
    return exists_in_active(key) || exists_in_global(key);
}


Data& Options::get(std::string key)
{
    to_upper(key);
    if (!exists_in_active(key)) {
        // Key not found. Throw an error
        throw IndexException(key);
    }
    return locals_[current_module_][key];
}

Data& Options::get(std::map<std::string, Data>& m, std::string& key)
{
    to_upper(key);
    return m[key];
}

Data& Options::get_global(std::string key)
{
    to_upper(key);
    if (!exists_in_global(key)) {
        // Key not found. Throw an error
        throw IndexException(key);
    }
    return globals_[key];
}

Data& Options::use(std::string& key)
{
    to_upper(key);

    // edit globals being true overrides everything
    if (edit_globals_){
        return get(globals_, key);
    }

    if (!exists_in_active(key) && !exists_in_global(key))
        throw IndexException(key);
    else if (!exists_in_active(key) && exists_in_global(key))
        return get(globals_, key);
    else if (exists_in_active(key) && exists_in_global(key)) {
        Data& active = get(locals_[current_module_], key);
        Data& global = get(globals_, key);

        if (active.has_changed()) {
            // Pull from keyvals
            return active;
        }
        else if (global.has_changed()){
            // Pull from globals
            return global;
        }
        else{
            // No user input - the default should come from local vals
            return active;
        }
    }
    else
        return get(locals_[current_module_], key);
}

bool Options::get_bool(std::string key)
{
    return(static_cast<bool>(use(key).to_integer()));
}

int Options::get_int(std::string key)
{
    return(use(key).to_integer());
}

double Options::get_double(std::string key)
{
    return(use(key).to_double());
}

std::string Options::get_str(std::string key) {
    return(use(key).to_string());
}

int* Options::get_int_array(std::string key)
{
    int *array = new int[use(key).size()];
    for (unsigned int i=0; i<use(key).size(); ++i) {
        array[i] = use(key)[i].to_integer();
    }
    return array;
}

void Options::fill_int_array(std::string key, int* empty_array)
{
    for (unsigned int i=0; i<use(key).size(); ++i) {
        empty_array[i] = use(key)[i].to_integer();
    }
}

std::vector<int> Options::get_int_vector(std::string key)
{
    std::vector<int> array;
    for (unsigned int i=0; i<use(key).size(); ++i) {
        array.push_back(use(key)[i].to_integer());
    }
    return array;
}

double* Options::get_double_array(std::string key)
{
    double *array = new double[use(key).size()];
    for (unsigned int i=0; i<use(key).size(); ++i) {
        array[i] = use(key)[i].to_double();
    }
    return array;
}

std::vector<double> Options::get_double_vector(std::string key)
{
    std::vector<double> array;
    for (unsigned int i=0; i<use(key).size(); ++i) {
        array.push_back(use(key)[i].to_double());
    }
    return array;
}

const char* Options::get_cstr(std::string key)
{
    return(use(key).to_string().c_str());
}

Data& Options::operator[](std::string key)
{
    return use(key);
}

}
