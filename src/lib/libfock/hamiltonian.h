#ifndef HAMILTONIAN_H 
#define HAMILTONIAN_H

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Matrix;
class Vector;
class JK;
class VBase;

// => BASE CLASSES <= //

class Hamiltonian {

protected:

    /// Print flag, defaults to 1
    int print_;	
    /// Debug flag, defaults to 0
    int debug_;	
    /// Bench flag, defaults to 0
    int bench_;
    /// jk object
    boost::shared_ptr<JK> jk_;  
    /// v object
    boost::shared_ptr<VBase> v_;  

    void common_init();
    
public:
    // => Constructors < = //

    Hamiltonian(boost::shared_ptr<JK> jk);
    Hamiltonian(boost::shared_ptr<JK> jk, boost::shared_ptr<VBase> v);
    /// Destructor
    virtual ~Hamiltonian();

    // => Accessors <= //

    /**
    * Pointer to the JK object
    * @return current JK object
    */
    boost::shared_ptr<JK> jk() const { return jk_; }
    /**
    * Pointer to the V object
    * @return current VBase object
    */
    boost::shared_ptr<VBase> v() const { return v_; }

    /**
    * Print header information regarding Hamiltonian
    * type on output file
    */
    virtual void print_header() const = 0;

    // => Knobs <= //

    /**
    * Knob to swap out a JK object
    * @param jk new JK object
    */
    void set_JK(boost::shared_ptr<JK> jk) { jk_ = jk; }
    /**
    * Knob to swap out a V object
    * @param v new V object
    */
    void set_V(boost::shared_ptr<VBase> v) { v_ = v; }
    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }
};

class RHamiltonian : public Hamiltonian {

public:
    // => Constructors < = //

    RHamiltonian(boost::shared_ptr<JK> jk);
    RHamiltonian(boost::shared_ptr<JK> jk, boost::shared_ptr<VBase> v);
    /// Destructor
    virtual ~RHamiltonian();

    // => Required Methods <= //

    /**
    * Return the approximate diagonal of the Hamiltonian
    * (typically orbital energy differences). This is used
    * for vector guess, eigenvector occupation, and preconditioning.
    * @return diagonal approximation, blocked by symmetry
    */
    virtual boost::shared_ptr<Vector> diagonal() = 0; 
    /**
    * Form the product Hx for each x found in the argument, placing in second argument 
    * @param x vector of state functions, blocked by symmetry 
    * @param b vector of product functions, blocked by symmetry (preallocated)
    */
    virtual void  product(const std::vector<boost::shared_ptr<Vector> >& x,
                                std::vector<boost::shared_ptr<Vector> >& b) = 0; 

    /**
    * Form the explicit hamiltonian for debugging purposes
    */
    SharedMatrix explicit_hamiltonian();

};

class UHamiltonian : public Hamiltonian {

public:
    // => Constructors < = //

    UHamiltonian(boost::shared_ptr<JK> jk);
    UHamiltonian(boost::shared_ptr<JK> jk, boost::shared_ptr<VBase> v);
    /// Destructor
    virtual ~UHamiltonian();

    // => Required Methods <= //

    /**
    * Return the approximate diagonal of the Hamiltonian
    * (typically orbital energy differences). This is used
    * for vector guess, eigenvector occupation, and preconditioning.
    * @return diagonal approximation, \alpha and \beta vectors,
    * blocked by symmetry
    */
    virtual std::pair<boost::shared_ptr<Vector>,
                      boost::shared_ptr<Vector> > diagonal() = 0; 
    /**
    * Form the product Hx for each x found in the argument, placing in second argument 
    * @param x vector of state functions, \alpha and \beta, blocked by symmetry 
    * @param b vector of product functions, \alpha and \beta, blocked by symmetry (preallocated)
    */
    virtual void product(const std::vector<std::pair<boost::shared_ptr<Vector>, boost::shared_ptr<Vector> > >& x,
                               std::vector<std::pair<boost::shared_ptr<Vector>, boost::shared_ptr<Vector> > >& b) = 0; 
    /**
    * Form the explicit hamiltonian for debugging purposes
    */
    std::pair<SharedMatrix, SharedMatrix > explicit_hamiltonian();
};

// => APPLIED CLASSES <= //

class MatrixRHamiltonian : public RHamiltonian {

protected:
    
    SharedMatrix M_;

public:

    MatrixRHamiltonian(SharedMatrix M);
    virtual ~MatrixRHamiltonian();

    virtual void print_header() const;
    virtual boost::shared_ptr<Vector> diagonal();
    virtual void product(const std::vector<boost::shared_ptr<Vector> >& x,
                               std::vector<boost::shared_ptr<Vector> >& b);

};

class MatrixUHamiltonian : public UHamiltonian {

protected:

    std::pair<SharedMatrix, SharedMatrix > M_;

public:

    MatrixUHamiltonian(std::pair<SharedMatrix, SharedMatrix > M);
    virtual ~MatrixUHamiltonian();

    virtual void print_header() const;
    virtual std::pair<boost::shared_ptr<Vector>,
                      boost::shared_ptr<Vector> > diagonal();
    virtual void product(const std::vector<std::pair<boost::shared_ptr<Vector>, boost::shared_ptr<Vector> > >& x,
                               std::vector<std::pair<boost::shared_ptr<Vector>, boost::shared_ptr<Vector> > >& b);

};

class CISRHamiltonian : public RHamiltonian {

protected:
    
    bool singlet_;
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    boost::shared_ptr<Vector> eps_aocc_;
    boost::shared_ptr<Vector> eps_avir_;

    bool rank_reduction_;
    int max_rank_;
    double rank_cutoff_;

public:
    CISRHamiltonian(boost::shared_ptr<JK> jk, 
                    SharedMatrix Caocc, 
                    SharedMatrix Cavir, 
                    boost::shared_ptr<Vector> eps_aocc,
                    boost::shared_ptr<Vector> eps_avir, 
                    boost::shared_ptr<VBase> v = boost::shared_ptr<VBase>());
    virtual ~CISRHamiltonian();

    virtual void print_header() const;
    virtual boost::shared_ptr<Vector> diagonal();
    virtual void product(const std::vector<boost::shared_ptr<Vector> >& x,
                               std::vector<boost::shared_ptr<Vector> >& b);

    virtual std::vector<SharedMatrix > unpack(const boost::shared_ptr<Vector>& x);

    void set_rank_cutoff(double cutoff) { rank_cutoff_ = cutoff; rank_reduction_ = (rank_cutoff_ != 0.0 || max_rank_ != 0);  }
    void set_max_rank(int rank) { max_rank_ = rank; rank_reduction_ = (rank_cutoff_ != 0.0 || max_rank_ != 0);  }

    void set_singlet(bool singlet) { singlet_ = singlet; }

};

class TDHFRHamiltonian : public RHamiltonian {

protected:

    bool singlet_;
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    boost::shared_ptr<Vector> eps_aocc_;
    boost::shared_ptr<Vector> eps_avir_;

public:
    TDHFRHamiltonian(boost::shared_ptr<JK> jk, 
                    SharedMatrix Caocc, 
                    SharedMatrix Cavir, 
                    boost::shared_ptr<Vector> eps_aocc,
                    boost::shared_ptr<Vector> eps_avir,
                    boost::shared_ptr<VBase> v = boost::shared_ptr<VBase>());
    virtual ~TDHFRHamiltonian();

    virtual void print_header() const;
    virtual boost::shared_ptr<Vector> diagonal();
    virtual void product(const std::vector<boost::shared_ptr<Vector> >& x,
                               std::vector<boost::shared_ptr<Vector> >& b);

    void set_singlet(bool singlet) { singlet_ = singlet; }
};

class CPHFRHamiltonian : public RHamiltonian {

protected:
    
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    boost::shared_ptr<Vector> eps_aocc_;
    boost::shared_ptr<Vector> eps_avir_;

public:
    CPHFRHamiltonian(boost::shared_ptr<JK> jk, 
                     SharedMatrix Caocc, 
                     SharedMatrix Cavir, 
                     boost::shared_ptr<Vector> eps_aocc,
                     boost::shared_ptr<Vector> eps_avir,
                     boost::shared_ptr<VBase> v = boost::shared_ptr<VBase>());
    virtual ~CPHFRHamiltonian();

    virtual void print_header() const;
    virtual boost::shared_ptr<Vector> diagonal();
    virtual void product(const std::vector<boost::shared_ptr<Vector> >& x,
                               std::vector<boost::shared_ptr<Vector> >& b);

    virtual std::map<std::string, SharedVector> pack(const std::map<std::string, boost::shared_ptr<Matrix> >& b);
    virtual std::vector<SharedMatrix > unpack(const std::vector<SharedVector >& x);
};

class TDARHamiltonian : public CISRHamiltonian {

protected:
    
    SharedMatrix Cocc_;

public:
    TDARHamiltonian(boost::shared_ptr<JK> jk, 
                    boost::shared_ptr<VBase> v,
                    SharedMatrix Cocc,
                    SharedMatrix Caocc, 
                    SharedMatrix Cavir, 
                    boost::shared_ptr<Vector> eps_aocc,
                    boost::shared_ptr<Vector> eps_avir);
    virtual ~TDARHamiltonian();

    virtual void print_header() const;
    virtual void product(const std::vector<boost::shared_ptr<Vector> >& x,
                               std::vector<boost::shared_ptr<Vector> >& b);
};

class TDDFTRHamiltonian : public TDHFRHamiltonian {

protected:
    
    SharedMatrix Cocc_;

public:
    TDDFTRHamiltonian(boost::shared_ptr<JK> jk, 
                    boost::shared_ptr<VBase> v,
                    SharedMatrix Cocc,
                    SharedMatrix Caocc, 
                    SharedMatrix Cavir, 
                    boost::shared_ptr<Vector> eps_aocc,
                    boost::shared_ptr<Vector> eps_avir);
    virtual ~TDDFTRHamiltonian();

    virtual void print_header() const;
    virtual void product(const std::vector<boost::shared_ptr<Vector> >& x,
                               std::vector<boost::shared_ptr<Vector> >& b);
};

class CPKSRHamiltonian : public CPHFRHamiltonian {

protected:
    
    SharedMatrix Cocc_;

public:
    CPKSRHamiltonian(boost::shared_ptr<JK> jk, 
                    boost::shared_ptr<VBase> v,
                    SharedMatrix Cocc,
                    SharedMatrix Caocc, 
                    SharedMatrix Cavir, 
                    boost::shared_ptr<Vector> eps_aocc,
                    boost::shared_ptr<Vector> eps_avir);
    virtual ~CPKSRHamiltonian();

    virtual void print_header() const;
    virtual void product(const std::vector<boost::shared_ptr<Vector> >& x,
                               std::vector<boost::shared_ptr<Vector> >& b);
};

}
#endif
