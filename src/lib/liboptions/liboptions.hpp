#ifndef _psi_src_lib_liboptions_liboptions_h
#define _psi_src_lib_liboptions_liboptions_h

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

#include <libutil/libutil.h> // Needed for Ref counting, string splitting, and conversions

namespace psi {

  class DataTypeException : public std::runtime_error
  {
  public:
    DataTypeException(const std::string& message) : std::runtime_error(message) { }
  };

  class IndexException : public std::runtime_error
  {
  public:
    IndexException(const std::string& message) : std::runtime_error("unable to find index " + message) { }
  };

  class DuplicateKeyException : public std::runtime_error
  {
  public:
    DuplicateKeyException() : std::runtime_error("duplicate key found") { }
  };

  class NotImplementedException : public std::runtime_error
  {
  public:
    NotImplementedException(const std::string& message) : std::runtime_error(message + " function not implemented") { }
  };

  class Data;
  class DataType 
  {
    bool changed_;
  public:
    DataType() : changed_(false) { }
    virtual ~DataType() { }

    bool has_changed() const {
      return changed_;
    }
    void changed() {
      changed_ = true;
    }

    void to_upper(std::string& str) {
      std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    }

    virtual bool is_array() const {
      return false;
    }
    virtual unsigned int size() const {
      throw NotImplementedException("size()");
    }
    virtual void add(DataType *) {
      throw NotImplementedException("add(DataType*)");
    }
    virtual void add(std::string, DataType*) {
      throw NotImplementedException("add(std::string, DataType*)");
    }
    virtual void add(int i) {
      throw NotImplementedException("add(int)");
    }
    virtual void add(double d) {
      throw NotImplementedException("add(double)");
    }
    virtual void add(std::string s, std::string c) {
      throw NotImplementedException("add(std::string, std::string)");
    }    
    virtual void add(std::string key, int i) {
      throw NotImplementedException("add(std::string, int)");
    }
    virtual void add(std::string key, double d) {
      throw NotImplementedException("add(std::string, double)");
    }
    virtual void add(std::string key, std::string s, std::string) {
      throw NotImplementedException("add(std::string, std::string, std::string)");
    }

    virtual bool exists(std::string key) {
      throw NotImplementedException("exists(std::string)");
    }

    virtual std::string to_string() {
      throw DataTypeException("don't know how to convert to a string");
    };
    virtual int to_integer() const {
      throw DataTypeException("don't know how to convert to an integer");
    }
    virtual double to_double() const {
      throw DataTypeException("don't know how to convert to a double");
    }

    virtual void assign(int) {
      throw DataTypeException("assign(int) failure");
    }
    virtual void assign(double) {
      throw DataTypeException("assign(double) failure");
    }
    virtual void assign(std::string) {
      throw DataTypeException("assign(std:string) failure");
    }

    virtual Data& operator[](std::string key) {
      throw NotImplementedException("Data& [string]");
    }    
    virtual Data& operator[](unsigned int i) {
      throw NotImplementedException("Data& [unsigned int]");
    }
  };

  class IntDataType : public DataType
  {
    int integer_;
  public:
    IntDataType() : DataType() { }
    IntDataType(int i) : DataType(), integer_(i) { }
    virtual ~IntDataType() { }

    virtual std::string to_string() {
      std::stringstream strm;
      strm << integer_;
      return strm.str();
    }
    virtual int to_integer() const {
      return integer_;
    }
    virtual double to_double() const {
      return static_cast<double>(integer_);
    }

    virtual void assign(int i) {
      changed();
      integer_ = i;
    }
    virtual void assign(double d) {
      assign(static_cast<int>(d));
    }
    virtual void assign(std::string s) {
      assign(static_cast<int>(std::strtod(s.c_str(), NULL)));
    }
  };

  class DoubleDataType : public DataType
  {
    double double_;
  public:
    DoubleDataType() : DataType() { }
    DoubleDataType(double d) : DataType(), double_(d) { }
    virtual ~DoubleDataType() { }

    virtual std::string to_string() {
      std::stringstream strm;
      strm << double_;
      return strm.str();
    }
    virtual int to_integer() const {
      return static_cast<int>(double_);
    }
    virtual double to_double() const {
      return double_;
    }

    virtual void assign(int i) {
      assign(tatic_cast<double>(i));
    }
    virtual void assign(double d) {
      changed();
      double_ = d;
    }
    virtual void assign(std::string s) {
      assign(std::strtod(s.c_str(), NULL));
    }
  };

  class StringDataType : public DataType
  {
    std::string str_;
    std::vector<std::string> choices_;
  public:
    StringDataType() : DataType() { }
    StringDataType(std::string s) : DataType(), str_(s) { to_upper(str_); }
    StringDataType(std::string s, std::string c) : DataType(), str_(s) { to_upper(str_); to_upper(c); choices_ = split(c); }
    virtual ~StringDataType() { }

    std::vector<std::string> split(const std::string& str){
      // Split a string
      typedef std::string::const_iterator iter;
      std::vector<std::string> splitted_string;
      iter i = str.begin();
      while(i != str.end()){
        // Ignore leading blanks
        i = find_if(i,str.end(), not_space);
        // Find the end of next word
        iter j = find_if(i,str.end(),space);
        // Copy the characters in [i,j)
        if(i!=str.end())
          splitted_string.push_back(std::string(i,j));
        i = j;
      }
      return(splitted_string);
    }

    virtual std::string to_string() {
      return str_;
    }
    virtual int to_integer() const {
      return static_cast<int>(std::strtod(str_.c_str(), NULL));
    }
    virtual double to_double() const {
      return std::strtod(str_.c_str(), NULL);
    }

    virtual void assign(int i) {
      std::stringstream strm;
      strm << i;
      assign(strm.str());
    }
    virtual void assign(double d) {
      std::stringstream strm;
      strm << d;
      assign(strm.str());
    }
    virtual void assign(std::string s) {
      to_upper(s);
      if (choices_.size() > 0) {
        bool wrong_input = true;
        for (int i=0; i<choices_.size(); ++i)
          if (s == choices_[i])
          wrong_input = false;
        if (wrong_input)
          throw DataTypeException(s + " is not a valid choice");
        str_ = s;
      }
      else {
        changed();
        str_ = s;
      }
    }
  };

  class Data
  {
    Ref<DataType> ptr_;
  public:
    Data() { }
    Data(DataType *type) : ptr_(type) { }
    Data(const Data& copy) { ptr_ = copy.ptr_; }

    std::string to_string() {
      return ptr_->to_string();
    }
    int to_integer() const {
      return ptr_->to_integer();
    }
    double to_double() const {
      return ptr_->to_double();
    }

    bool is_array() const {
      return ptr_->is_array();
    }

    void add(DataType *data) {
      ptr_->add(data);
    }
    void add(std::string s, DataType *data) {
      ptr_->add(s, data);
    }
    void add(int i) {
      ptr_->add(i);
    }
    void add(double d) {
      ptr_->add(d);
    }
    void add(std::string s, std::string c) {
      ptr_->add(s, c);
    }    
    void add(std::string key, int i) {
      ptr_->add(key, i);
    }
    void add(std::string key, double d) {
      ptr_->add(key, d);
    }
    void add(std::string key, std::string s, std::string c) {
      ptr_->add(key, s, c);
    }

    void assign(int i) {
      ptr_->assign(i);
    }
    void assign(double d) {
      ptr_->assign(d);
    }
    void assign(std::string s) {
      ptr_->assign(s);
    }

    Data& operator[](int i) {
      return (*(ptr_.pointer()))[i];
    }
    Data& operator[](std::string s) {
      return (*(ptr_.pointer()))[s];
    }
  };

  class ArrayType : public DataType
  {
    std::vector<Data> array_;
  public:
    ArrayType() { }

    virtual void add(DataType *data) {
      array_.push_back(Data(data));
    }
    virtual void add(int i) {
      add(new IntDataType(i));
    }
    virtual void add(double d) {
      add(new DoubleDataType(d));
    }
    virtual void add(std::string s, std::string c = "") {
      add(new StringDataType(s, c));
    }

    virtual Data& operator[](unsigned int i) {
      if (i >= array_.size())
        throw IndexException("out of range");
      return array_[i];
    }
    virtual Data& operator[](std::string s) {
      int i = static_cast<int>(std::strtod(s.c_str(), NULL));
      if (i >= array_.size())
        throw IndexException("out of range");
      return array_[i];
    }
    virtual bool is_array() {
      return true;
    }

    virtual unsigned int size() const {
      return array_.size();
    }

    virtual std::string to_string() {
      std::string str = "[ ";
      for (int i=0; i<array_.size(); ++i) {
        str += array_[i].to_string();
        if (i != array_.size() - 1)
          str += ", ";
      }
      str += " ]";
      return str;
    }    
  };

  class MapType : public DataType
  {
    std::map<std::string, Data> keyvals_;
    typedef std::map<std::string, Data>::iterator iterator;
  public:
    MapType() { }

    virtual void add(std::string key, DataType *data) {
      to_upper(key);

      iterator pos = keyvals_.find(key);
      if (pos != keyvals_.end())
        throw DuplicateKeyException();
      keyvals_[key] = Data(data);
    }

    virtual void add(std::string key, int i) {
      add(key, new IntDataType(i));
    }
    virtual void add(std::string key, double d) {
      add(key, new DoubleDataType(d));
    }
    virtual void add(std::string key, std::string s, std::string c = "") {
      add(key, new StringDataType(s, c));
    }

    virtual bool exists(std::string key) {
      to_upper(key);
      iterator pos = keyvals_.find(key);
      if (pos != keyvals_.end())
        return true;
      return false;
    }

    virtual Data& operator[](std::string s) {
      to_upper(s);
      if (!exists(s))
        throw IndexException(s);
      return keyvals_[s];
    }
    virtual bool is_array() {
      return true;
    }

    virtual unsigned int size() const {
      return keyvals_.size();
    }

    virtual std::string to_string() {
      std::string str = "{ ";
      for (iterator pos = keyvals_.begin(); pos != keyvals_.end(); ++pos) {
        str += pos->first + " => " + pos->second.to_string() + ", ";
      }
      str += "}";
      return str;
    }    
  };

  class Options
  {
    std::map<std::string, Data> keyvals_;
    typedef std::map<std::string, Data>::iterator iterator;
  public:
    Options() { }

    void to_upper(std::string& str) {
      std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    }

    void add(std::string key, DataType *data) {
      to_upper(key);

      // Make sure the key isn't already there
      iterator pos = keyvals_.find(key);
      if (pos != keyvals_.end())
        throw DuplicateKeyException();
      keyvals_[key] = Data(data);
    }

    void add(std::string key, int i) {
      add(key, new IntDataType(i));
    }
    void add(std::string key, double d) {
      add(key, new DoubleDataType(d));
    }
    void add(std::string key, std::string s, std::string c = "") {
      add(key, new StringDataType(s, c));
    }

    bool exists(std::string key) {
      to_upper(key);
      iterator pos = keyvals_.find(key);
      if (pos != keyvals_.end())
        return true;
      return false;
    }

    Data& get(std::string key) {
      to_upper(key);      
      if (!exists(key)) {
        // Key not found. Throw an error
        throw IndexException(key);
      }
      return keyvals_[key];
    }
    Data& operator[](std::string key) {
      return get(key);
    }

    std::string to_string() {
      std::stringstream str;
      str << "Options: \n";
      for (iterator pos = keyvals_.begin(); pos != keyvals_.end(); ++pos) {
        str << "  " << std::setw(12) << pos->first << " => " << pos->second.to_string() << std::endl;
      }
      return str.str();
    }
  };

}
#endif
