#ifndef dispersion_h
#define dispersion_h

/**********************************************************
* dispersion.h: declarations -D(1-3) for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include <psi4-dec.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

class Molecule;
class Matrix;

namespace functional {

class Dispersion {
    public:

        static boost::shared_ptr<Dispersion> createDispersion(const std::string & type);

        Dispersion();
        virtual ~Dispersion();

        std::string getName() const { return name_; }        
        std::string getDescription() const { return description_; }        
        std::string getCitation() const { return citation_; }        
        void setName(const std::string & name) { name_ = name; }
        void setDescription(const std::string & description) { description_ = description; }
        void setCitation(const std::string & citation) { citation_ = citation; }

        virtual double computeEnergy(boost::shared_ptr<Molecule> m) = 0; 
        virtual boost::shared_ptr<Matrix> computeGradient(boost::shared_ptr<Molecule> m) = 0; 
        virtual boost::shared_ptr<Matrix> computeHessian(boost::shared_ptr<Molecule> m) = 0; 
     
    protected:

        std::string name_;
        std::string description_;
        std::string citation_;
};

class D1 : public Dispersion {
    public:
        D1(); 
        virtual ~D1(); 
        double computeEnergy(boost::shared_ptr<Molecule> m); 
        boost::shared_ptr<Matrix> computeGradient(boost::shared_ptr<Molecule> m); 
        boost::shared_ptr<Matrix> computeHessian(boost::shared_ptr<Molecule> m); 
         
    protected:
       
         
};
class D2 : public Dispersion {
    public:
        D2(); 
        virtual ~D2(); 
        double computeEnergy(boost::shared_ptr<Molecule> m); 
        boost::shared_ptr<Matrix> computeGradient(boost::shared_ptr<Molecule> m); 
        boost::shared_ptr<Matrix> computeHessian(boost::shared_ptr<Molecule> m); 
         
    protected:
       
         
};
class D3 : public Dispersion {
    public:
        D3(); 
        virtual ~D3(); 
        double computeEnergy(boost::shared_ptr<Molecule> m); 
        boost::shared_ptr<Matrix> computeGradient(boost::shared_ptr<Molecule> m); 
        boost::shared_ptr<Matrix> computeHessian(boost::shared_ptr<Molecule> m); 
         
    protected:
       
         
};

}}

#endif 
