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
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef THCE_H
#define THCE_H

#include "psi4/libmints/typedefs.h"
#include "psi4/libpsi4util/exception.h"

#include <map>

namespace psi {

class Tensor;
class CoreTensor;
class DiskTensor;

/**
 * Class THCE is a container and interface to
 * a collection of Tensor objects
 **/
class THCE {


protected:

    // => Object Data <= //

    /// List of currently declared dimensions (dimension classes)
    std::map<std::string, int> dimensions_;
    /// List of currently declated tensors
    std::map<std::string, std::shared_ptr<Tensor> > tensors_;

    // => Safety/Debugging <= //

    /// Throws if a dimension name is not in dimensions_
    void dimension_check(const std::string& name);
    /// Throws if a tensor name is not in tensors_
    void tensor_check(const std::string& name);

    /// Get a unique ID number for a tensor (prevents disk clash)
    static int unique_id();

public:

    // => Constructors <= //

    /// Instance Constructor
    THCE();
    /// Instance Destructor
    virtual ~THCE();

    // => Accessors <= //

    /// Accessor to the list of currently declared dimensions (dimension classes)
    std::map<std::string, int>& dimensions() { return dimensions_; }
    /// Accessor to the list of currently declared tensors
    std::map<std::string, std::shared_ptr<Tensor> >& tensors() { return tensors_; }
    /// Direct tensor element accessor
    std::shared_ptr<Tensor>& operator[] (const std::string& key) { return tensors_[key]; }

    /// Print a brief summary
    void print_header(std::string fh="outfile") const { print(fh,0); }
    /// Print a more detailed trace of the object
    void print(std::string fh="outfile", int print = 1) const;

    /// Return the current core memory utilization of this instance, in doubles
    size_t core_doubles() const;
    /// Return the current disk memory utilization of this instance, in doubles
    size_t disk_doubles() const;

    // => Front-end interface <= //

    /// Add or overwrite a dimension to this instance
    void new_dimension(const std::string& name, int size);
    /// Remove a dimension from this instance
    void delete_dimension(const std::string& name);
    /// Add or overwrite a new core tensor, given a name and comma-separated dimension list
    void new_core_tensor(const std::string& name, const std::string& dimensions, double* data = NULL, bool trust = false);
    /// Add or overwrite a new disk tensor, given a name and comma-separated dimension list
    void new_disk_tensor(const std::string& name, const std::string& dimensions, bool save = false, bool load = false);
    /// Add or alias an existing Tensor of any type
    void add_tensor(const std::string& name, std::shared_ptr<Tensor> tensor);
    /// Decrement a Tensor of any type (removes from tensors_)
    void delete_tensor(const std::string& name);

};

/**
 * Abstract class Tensor provides the high-level interface for
 * Tensor creation and manipulation in disk or core forms
 **/
class Tensor {

private:

    /// Unique ID to prevent disk clash
    static long int unique_id;

protected:

    // => Object Data <= //

    /// Name of this Tensor
    std::string name_;
    /// Unique filename to prevent disk clash
    std::string filename_;
    /// Total number of elements in this Tensor
    size_t numel_;
    /// Order (number of dimensions) in this Tensor
    int order_;
    /// Classes of dimensions
    std::vector<std::string> dimensions_;
    /// Sizes of dimensions (lda-style)
    std::vector<int> sizes_;
    /// Active sizes of dimensions (gimp-style)
    std::vector<int> active_sizes_;

    // => Helper Methods <= //

    /// Set the filename to scratch, PID, namespace, unique ID, name
    void set_filename();

public:

    // => Constructors <= //

    /// Master constructor
    Tensor(const std::string& name,
           std::vector<std::string>& dimensions,
           std::vector<int>& sizes);

    /// Master destructor
    virtual ~Tensor();

    // => Universal Accessors <= //

    /// Name of this Tensor
    std::string name() const { return name_; }
    /// Unique filename of this Tensor
    std::string filename() const { return filename_; }
    /// Total number of elements in this Tensor
    size_t numel() const { return numel_; }
    /// Order (number of dimensions) in this Tensor
    int order() const { return order_; }
    /// Classes of dimensions
    std::vector<std::string>& dimensions() { return dimensions_; }
    /// Allocated sizes of dimensions
    std::vector<int>& sizes() { return sizes_; }
    /// Active sizes of dimensions
    std::vector<int>& active_sizes() { return active_sizes_; }

    /// Is this a core tensor?
    virtual bool core() const = 0;
    /// Is this a disk tensor?
    virtual bool disk() const = 0;
    /// Is this a trust core tensor?
    virtual bool trust() const = 0;
    /// Is this a swapped-out core tensor?
    virtual bool swapped() const = 0;
    /// Throw an exception if the tensor is swapped out or disk based
    virtual void swap_check() const = 0;

    /// How many core doubles are currently allocated by this tensor?
    virtual size_t core_doubles() const = 0;
    /// How many disk doubles are currently allocated by this tensor?
    virtual size_t disk_doubles() const = 0;

    /// Print the full tensor data in scalar/vector/matrix/pages style if print >= 2
    virtual void print(std::string fh="outfile", int level = 2) const = 0;
    /// Print only name and sizing data
    void print_header(std::string fh="outfile") const { print(fh, 0); }


    // => Common Interface <= //

    /// Zero the tensor out (ignores active dims). Prestripes if DiskTensor
    virtual void zero() = 0;
    /// Slice assignment from A into C (respects active dims)
    /// Topology is a list following the joint rank of A and C, which must be coincidentally ordered (less singletons)
    /// in the form <inC?,startC,endC,inA?,startA,endA>, e.g.,
    /// C_{0,:,2:3,:,1:2} = A{1,:,:,2:3,4:5} =>
    /// <true,0,1,false,-1,-1>  % Singleton A (Singleton A trumps Singleton C)
    /// <false,-1,-1,true,1,2>  % Singleton C
    /// <true,-1,-1,true,-1,-1> % Full/Full
    /// <true,2,3,true,-1,-1>   % Slice/Full
    /// <true,-1,-1,true,2,3>   % Full/Slice
    /// <true,1,2,true,4,5>     % Slice/Slice
    /// Slice indices are 0-based, and of the form [start,end)
    virtual void slice(std::shared_ptr<Tensor> A, std::vector<std::tuple<bool,int,int,bool,int,int> >& topology);
    /// Access the data pointer
    virtual double* pointer() const { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }
    /// Access the file pointer
    virtual FILE* file_pointer() { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }

    // => CoreTensor Interface <= //

    // > Unary Operations < //

    /// Swap this tensor out to disk (does not write if existing disk mirror and not changed)
    virtual void swap_out(bool changed = true) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }
    /// Swap this tensor in to core (explicitly reads if true)
    virtual void swap_in(bool read = true) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }
    /// Scale the tensor by val (ignores active dims)
    virtual void scale(double val) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }
    /// Copy contents of data (ignores active dims)
    virtual void set_data(double* data) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }
    /// Update the trust pointer to data (throws if not trust)
    virtual void set_pointer(double* data) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }

    // > Binary Operations < //

    /// Direct elementwise generalized DAXPY: C = alpha * A + beta * C (ignores active dims)
    virtual void add(std::shared_ptr<Tensor> A, double alpha = 1.0, double beta = 0.0) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }
    /// Permute from A into C (ignores active dims)
    /// Topology (P) contains rank of index in C indexed by rank in A
    /// C_{P=0,P=1,...} = A_{0,1,...}, e.g.,
    /// C_{k,i,l,j} = A_{i,j,k,l} =>
    /// [1,3,0,2]
    virtual void permute(std::shared_ptr<Tensor> A, std::vector<int>& topology) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }

    // > Ternary Operations < //

    /// Main contraction operation method (respects active dims)
    /// C = alpha * Op(A, B) + beta * C
    /// Topology contains the specification of the contraction,
    /// in the form <name,rC,rA,rB>, e.g.,
    /// C_{j,i,l} = A_{j,i,k} B_{j,l,k} =>
    /// <j,0,0,0>  % Hadamard Dimension (H)
    /// <i,1,1,-1> % Outer Dimension (L)
    /// <l,2,-1,1> % Outer Dimension (R)
    /// <k,-1,2,2> % Inner Dimension (I)
    ///
    /// Active dimensions H, L, and R are inherited in C from A and B tensors.
    /// If an index is H or I both A and B active dimensions must agree
    virtual void contract(std::shared_ptr<Tensor> A, std::shared_ptr<Tensor> B, std::vector<std::tuple<std::string,int,int,int> >& topology, double alpha = 1.0, double beta = 0.0) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }

    // => DiskTensor Interface <= //

    /// Set save flag in DiskTensor
    virtual void set_save(bool save) { throw PSIEXCEPTION("Not implemented in this Tensor subclass."); }
};

/**
 * Class CoreTensor provides functionality for the
 * creation and manipulation of Tensor objects in
 * core memory
 **/
class CoreTensor : public Tensor {

protected:

    // => Object Data <= //

    /// Is this a trust pointer?
    bool trust_;
    /// Data underlying this Tensor
    double* data_;
    /// Is this tensor swapped?
    bool swapped_;
    /// File handle (valid if this tensor has ever been swapped out)
    FILE* fh_;

public:

    // => Constructors <= //

    /// Master constructor
    CoreTensor(const std::string& name,
        std::vector<std::string>& dimensions, std::vector<int>& sizes,
        double* data = NULL, bool trust = false);

    /// Master destructor
    virtual ~CoreTensor();

    /// Order-0 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        double* data = NULL, bool trust = false);
    /// Order-1 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        const std::string& dimension1, int size1,
        double* data = NULL, bool trust = false);
    /// Order-2 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        const std::string& dimension1, int size1,
        const std::string& dimension2, int size2,
        double* data = NULL, bool trust = false);
    /// Order-3 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        const std::string& dimension1, int size1,
        const std::string& dimension2, int size2,
        const std::string& dimension3, int size3,
        double* data = NULL, bool trust = false);
    /// Order-4 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        const std::string& dimension1, int size1,
        const std::string& dimension2, int size2,
        const std::string& dimension3, int size3,
        const std::string& dimension4, int size4,
        double* data = NULL, bool trust = false);

    // => Universal Accessors <= //

    /// Is this a core tensor?
    virtual bool core() const { return true; }
    /// Is this a disk tensor?
    virtual bool disk() const { return false; }
    /// Is this a trust core tensor?
    virtual bool trust() const { return trust_; }
    /// Is this a swapped-out core tensor?
    virtual bool swapped() const { return swapped_; }
    /// Throw an exception if the tensor is swapped out or disk based
    virtual void swap_check() const;

    /// How many core doubles are currently allocated by this tensor?
    virtual size_t core_doubles() const { return (swapped() | trust() ? 0L : numel_); }
    /// How many disk doubles are currently allocated by this tensor?
    virtual size_t disk_doubles() const { return 0L; }

    /// Print the full tensor data in scalar/vector/matrix/pages style if print >= 2
    virtual void print(std::string fh="outfile", int level = 2) const;

    // => Conditional Accessors <= //

    /// Data underlying this Tensor (careful, returns NULL if swapped)
    virtual double* pointer() const { return data_; }
    /// File pointer underlying this Tensor (careful, returns NULL if not swapped, mirror not guaranteed)
    virtual FILE* file_pointer() { return fh_; }

    // > Unary Operations < //

    /// Swap this tensor out to disk (if existing file, only write if changed)
    virtual void swap_out(bool changed = true);
    /// Swap this tensor in to core (if read is false, allocate with zeros)
    virtual void swap_in(bool read = true);
    /// Zero the tensor out
    virtual void zero();
    /// Scale the tensor by val
    virtual void scale(double val);
    /// Copy contents of data (ignores active dims)
    virtual void set_data(double* data);
    /// Update the trust pointer to data (throws if not trust)
    virtual void set_pointer(double* data);

    // > Binary Operations < //

    /// Direct elementwise generalized DAXPY: C = alpha * A + beta * C (ignores active dims)
    virtual void add(std::shared_ptr<Tensor> A, double alpha = 1.0, double beta = 0.0);
    /// Permute from A into C
    /// Topology (P) contains rank of index in C indexed by rank in A
    /// C_{P=0,P=1,...} = A_{0,1,...}, e.g.,
    /// C_{k,i,l,j} = A_{i,j,k,l} =>
    /// [1,3,0,2]
    virtual void permute(std::shared_ptr<Tensor> A, std::vector<int>& topology);

    // > Ternary Operations < //

    /// C = alpha * Op(A, B) + beta * C
    /// Topology contains the specification of the contraction,
    /// in the form <name,rC,rA,rB>, e.g.,
    /// C_{j,i,l} = A_{j,i,k} B_{j,l,k} =>
    /// <j,0,0,0>  % Hadamard Dimension (H)
    /// <i,1,1,-1> % Outer Dimension (L)
    /// <l,2,-1,1> % Outer Dimension (R)
    /// <k,-1,2,2> % Inner Dimension (I)
    ///
    /// Active dimensions H, L, and R are inherited in C from A and B tensors.
    /// If an index is H or I both A and B active dimensions must agree
    virtual void contract(std::shared_ptr<Tensor> A, std::shared_ptr<Tensor> B, std::vector<std::tuple<std::string,int,int,int> >& topology, double alpha = 1.0, double beta = 0.0);

};

/**
 * Class DiskTensor provides functionality for the
 * creation and manipulation of Tensor objects in
 * scratch disk space
 **/
class DiskTensor : public Tensor {

protected:

    // => Object Data <= //

    /// Save file upon destruct?
    bool save_;
    /// File pointer
    FILE* fh_;

public:

    // => Constructors <= //

    /// Master constructor
    /// Will create new file (w+) if load
    /// Will open existing file (r+) if !load
    /// Will save file upon destruct if save
    /// Will delete file upon destruct if !save
    /// Does not explicitly prestripe, call zero to accomplish this
    DiskTensor(const std::string& name,
        std::vector<std::string>& dimensions, std::vector<int>& sizes,
        bool save = false, bool load = false);

    /// Master destructor
    virtual ~DiskTensor();

    /// Order-0 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        bool save = false, bool load = false);
    /// Order-1 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        const std::string& dimension1, int size1,
        bool save = false, bool load = false);
    /// Order-2 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        const std::string& dimension1, int size1,
        const std::string& dimension2, int size2,
        bool save = false, bool load = false);
    /// Order-3 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        const std::string& dimension1, int size1,
        const std::string& dimension2, int size2,
        const std::string& dimension3, int size3,
        bool save = false, bool load = false);
    /// Order-4 Constructor
    static std::shared_ptr<Tensor> build(const std::string& name,
        const std::string& dimension1, int size1,
        const std::string& dimension2, int size2,
        const std::string& dimension3, int size3,
        const std::string& dimension4, int size4,
        bool save = false, bool load = false);

    // => Universal Accessors <= //

    /// Is this a core tensor?
    virtual bool core() const { return false; }
    /// Is this a disk tensor?
    virtual bool disk() const { return true; }
    /// Is this a trust core tensor?
    virtual bool trust() const { return false; }
    /// Is this a swapped-out core tensor?
    virtual bool swapped() const { return false; }
    /// Throw an exception if the tensor is swapped out or disk based
    virtual void swap_check() const;

    /// How many core doubles are currently allocated by this tensor?
    virtual size_t core_doubles() const { return 0L; }
    /// How many disk doubles are currently allocated by this tensor?
    virtual size_t disk_doubles() const { return numel_; }

    /// Print the available tensor data
    virtual void print(std::string fh="outfile", int level = 2) const;

    // => Conditional Accessors <= //

    /// File underlying this Tensor
    virtual FILE* file_pointer() { return fh_; }

    /// Set save flag for disk resuse
    virtual void set_save(bool save) { save_ = save; }

    /// Zero the tensor out and prestripe
    virtual void zero();

};

} // End namespace

#endif
