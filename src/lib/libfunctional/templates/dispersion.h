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

namespace functional {

class Dispersion {
    public:

        static boost::shared_ptr<Dispersion> createDispersion(const std::string & type, double s6 = 0.0, double s8 = 0.0, \
            double sr6 = 0.0, double sr8 = 0.0);

        Dispersion();
        ~Dispersion();

        std::string getName() const { return name_; }
        std::string getDescription() const { return description_; }
        std::string getCitation() const { return citation_; }
        void setName(const std::string & name) { name_ = name; }
        void setDescription(const std::string & description) { description_ = description; }
        void setCitation(const std::string & citation) { citation_ = citation; }

        double getDampingD() const { return d_; }
        double getDampingAlpha6() const { return alpha6_; }
        double getDampingAlpha8() const { return alpha8_; }
        double getS6() const { return s6_; }
        double getS8() const { return s8_; }
        double getSR6() const { return sr6_; }
        double getSR8() const { return sr8_; }
        double getK1() const { return k1_; }
        double getK2() const { return k2_; }
        double getK3() const { return k3_; }

        void setDampingD(double d) { d_ = d; }
        void setDampingAlpha6(double d) { alpha6_ = d; }
        void setDampingAlpha8(double d) { alpha8_ = d; }
        void setS6(double s6) { s6_ = s6; }
        void setS8(double s8) { s8_ = s8; }
        void setSR6(double sr6) { sr6_ = sr6; }
        void setSR8(double sr8) { sr8_ = sr8; }
        void setK1(double k1) { k1_ = k1; }
        void setK2(double k2) { k2_ = k2; }
        void setK3(double k3) { k3_ = k3; }

        std::string printEnergy(boost::shared_ptr<Molecule> m);
        std::string printGradient(boost::shared_ptr<Molecule> m);
        std::string printHessian(boost::shared_ptr<Molecule> m);

        virtual double computeEnergy(boost::shared_ptr<Molecule> m) { return 0.0; }
        virtual double* computeGradient(boost::shared_ptr<Molecule> m) { double *blank = 0; return blank; }
        virtual double** computeHessian(boost::shared_ptr<Molecule> m) { double **blank = 0; return blank; }

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
};

class D1 : public Dispersion {
    public:
        D1(double s6 = 1.0);
        ~D1();
        double computeEnergy(boost::shared_ptr<Molecule> m);
        double* computeGradient(boost::shared_ptr<Molecule> m);
        double** computeHessian(boost::shared_ptr<Molecule> m);

    protected:
};
class D2 : public Dispersion {
    public:
        D2(double s6 = 1.0);
        ~D2();
        double computeEnergy(boost::shared_ptr<Molecule> m);
        double* computeGradient(boost::shared_ptr<Molecule> m);
        double** computeHessian(boost::shared_ptr<Molecule> m);

    protected:


};
class D3 : public Dispersion {
    public:
        D3(double s6 = 1.0, double s8 = 1.0, double sr6 = 1.0, double sr8 = 1.0);
        ~D3();
        double computeEnergy(boost::shared_ptr<Molecule> m);
        double* computeGradient(boost::shared_ptr<Molecule> m);
        double** computeHessian(boost::shared_ptr<Molecule> m);

    protected:


};

}}

#endif
