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

    protected:

        std::string name_;
        std::string description_;
        std::string citation_;
        double s6_;
        double s8_;

        // -D1 and -D2 values
        double d_;
        const double *RvdW_;
        const double *C6_;

        // -D2 values
        double sr6_;
        double sr8_;
        double alpha6_;
        double alpha8_;
        double k1_;
        double k2_;
        double k3_;

    public:

        Dispersion();
        ~Dispersion();

        static boost::shared_ptr<Dispersion> build(const std::string & type, double s6 = 0.0, double s8 = 0.0, \
            double sr6 = 0.0, double sr8 = 0.0);

        std::string name() const { return name_; }
        std::string description() const { return description_; }
        std::string citation() const { return citation_; }
        void set_name(const std::string & name) { name_ = name; }
        void set_description(const std::string & description) { description_ = description; }
        void set_citation(const std::string & citation) { citation_ = citation; }

        double get_DampingD() const { return d_; }
        double get_DampingAlpha6() const { return alpha6_; }
        double get_DampingAlpha8() const { return alpha8_; }
        double get_S6() const { return s6_; }
        double get_S8() const { return s8_; }
        double get_SR6() const { return sr6_; }
        double get_SR8() const { return sr8_; }
        double get_K1() const { return k1_; }
        double get_K2() const { return k2_; }
        double get_K3() const { return k3_; }

        void set_DampingD(double d) { d_ = d; }
        void set_DampingAlpha6(double d) { alpha6_ = d; }
        void set_DampingAlpha8(double d) { alpha8_ = d; }
        void set_S6(double s6) { s6_ = s6; }
        void set_S8(double s8) { s8_ = s8; }
        void set_SR6(double sr6) { sr6_ = sr6; }
        void set_SR8(double sr8) { sr8_ = sr8; }
        void set_K1(double k1) { k1_ = k1; }
        void set_K2(double k2) { k2_ = k2; }
        void set_K3(double k3) { k3_ = k3; }

        std::string print_energy(boost::shared_ptr<Molecule> m);
        std::string print_gradient(boost::shared_ptr<Molecule> m);
        std::string print_hessian(boost::shared_ptr<Molecule> m);

        virtual double compute_energy(boost::shared_ptr<Molecule> m) = 0; 
        virtual SharedMatrix compute_gradient(boost::shared_ptr<Molecule> m) = 0; 
        virtual SharedMatrix compute_hessian(boost::shared_ptr<Molecule> m) = 0; 

        virtual void print(FILE* out = outfile, int level = 1) const;
};

class D1 : public Dispersion {
    protected:
    public:
        D1(double s6 = 1.0);
        ~D1();
        double compute_energy(boost::shared_ptr<Molecule> m);
        SharedMatrix compute_gradient(boost::shared_ptr<Molecule> m);
        SharedMatrix compute_hessian(boost::shared_ptr<Molecule> m);
};
class D2 : public Dispersion {
    protected:
    public:
        D2(double s6 = 1.0);
        ~D2();
        double compute_energy(boost::shared_ptr<Molecule> m);
        SharedMatrix compute_gradient(boost::shared_ptr<Molecule> m);
        SharedMatrix compute_hessian(boost::shared_ptr<Molecule> m);

};
class D3 : public Dispersion {
    protected:
    public:
        D3(double s6 = 1.0, double s8 = 1.0, double sr6 = 1.0, double sr8 = 1.0);
        ~D3();
        double compute_energy(boost::shared_ptr<Molecule> m);
        SharedMatrix compute_gradient(boost::shared_ptr<Molecule> m);
        SharedMatrix compute_hessian(boost::shared_ptr<Molecule> m);
};

}

#endif
