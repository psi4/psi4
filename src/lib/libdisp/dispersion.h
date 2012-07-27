#ifndef dispersion_h
#define dispersion_h

/**********************************************************
* dispersion.h: declarations -D(1-3) for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include <psi4-dec.h>
#include <string>

namespace psi {

class Molecule;
class Matrix;

class Dispersion {

public:
    enum C6_type { C6_arit, C6_geom };
    enum C8_type { C8_geom };
    enum Damping_type { Damping_D1, Damping_CHG, Damping_TT };
    enum Spherical_type { Spherical_Das, Spherical_zero };

protected:

    std::string name_;
    std::string description_;
    std::string citation_;

    C6_type C6_type_;
    C8_type C8_type_;
    Damping_type Damping_type_;
    Spherical_type Spherical_type_;

    double s6_;
    double d_;
    const double *RvdW_;
    const double *C6_;
    const double *C8_;
    const double *A_;
    const double *Beta_;

public:

    Dispersion();
    virtual ~Dispersion();

    static boost::shared_ptr<Dispersion> build(const std::string & type, double s6 = 0.0);

    std::string name() const { return name_; }
    std::string description() const { return description_; }
    std::string citation() const { return citation_; }
    void set_name(const std::string & name) { name_ = name; }
    void set_description(const std::string & description) { description_ = description; }
    void set_citation(const std::string & citation) { citation_ = citation; }

    boost::shared_ptr<Vector> set_atom_list(boost::shared_ptr<Molecule> mol);

    double get_d() const { return d_; }
    double get_s6() const { return s6_; }

    void set_d(double d) { d_ = d; }
    void set_s6(double s6) { s6_ = s6; }

    std::string print_energy(boost::shared_ptr<Molecule> m);
    std::string print_gradient(boost::shared_ptr<Molecule> m);
    std::string print_hessian(boost::shared_ptr<Molecule> m);

    virtual double compute_energy(boost::shared_ptr<Molecule> m);
    virtual SharedMatrix compute_gradient(boost::shared_ptr<Molecule> m);
    virtual SharedMatrix compute_hessian(boost::shared_ptr<Molecule> m);

    virtual void print(FILE* out = outfile, int level = 1) const;
    void py_print() const { print(outfile, 1); }
};

}

#endif
