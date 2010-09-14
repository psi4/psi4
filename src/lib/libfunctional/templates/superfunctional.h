#ifndef super_functional_h
#define super_functional_h

/**********************************************************
* superfunctional.h: declarations for super (aliased) functionals for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/
#include "functional.h"
#include "dispersion.h"
#include <psi4-dec.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

class Properties;

namespace functional {

class SuperFunctional {
    public:

        static boost::shared_ptr<SuperFunctional> createSuperFunctional(int npoints = 5000, int deriv = 1);
        static boost::shared_ptr<SuperFunctional> buildSuperFunctional(const std::string & build, int npoints = 5000, int deriv = 1);        

        SuperFunctional(int npoints = 5000, int deriv = 1);
        virtual ~SuperFunctional();

        std::string getName() const { return name_; }        
        std::string getDescription() const { return description_; }        
        std::string getCitation() const { return citation_; }        
        void setName(const std::string & name) { name_ = name; }
        void setDescription(const std::string & description) { description_ = description; }
        void setCitation(const std::string & citation) { citation_ = citation; }

        std::string getComposition();        
        
        bool isGGA();    
        bool isMeta();    
        
        bool isHybrid() { return  exact_exchange_ > 0.0; }  
        bool isDoubleHybrid() { return pt2_ > 0.0; }  
        bool isRangeCorrected() { return omega_ > 0.0; }  
        bool isDashD() { return dashD_weight_ > 0.0; }  
        void setParameter(const std::string & functional, const std::string & param, double val);

        double getExactExchange() { return exact_exchange_; } 
        void setExactExchange(double exch) {exact_exchange_ = exch; }
     
        double getPT2() { return pt2_; } 
        void setPT2(double pt2) {pt2_ = pt2; }

        double getOmega() { return omega_; } 
        void setOmega(double omega) {omega_ = omega; }

        double getDashDWeight() { return dashD_weight_; }
        boost::shared_ptr<Dispersion> getDashD() { return dashD_; }
        void setDashD(shared_ptr<Dispersion> disp, double weight);    

        void setNPoints(int npoints) { npoints_ = npoints; reallocate(); }
        void setDeriv(int deriv) { deriv_ = deriv; reallocate(); }
        int getNPoints() const { return npoints_; }
        int getDeriv() const { return deriv_; }

        void computeRKSFunctional(boost::shared_ptr<Properties> props);
        void computeUKSFunctional(boost::shared_ptr<Properties> props);

        void addFunctional(const boost::shared_ptr<Functional> & f, double weight);
        void addFunctional(const std::string & name, double weight);
        boost::shared_ptr<Functional>& operator[](int index) { return functionals_[index].first; }
        boost::shared_ptr<Functional> getFunctional(int index) { return functionals_[index].first; }
        double getWeight(int index) { return functionals_[index].second; }
        int size() { return functionals_.size(); }
    
     
    protected:
        void reallocate();
 
        std::vector<std::pair<boost::shared_ptr<Functional>, double> > functionals_;

        std::string name_;
        std::string description_;
        std::string citation_;

        double exact_exchange_;
        double pt2_;
        double omega_;
        double dashD_weight_;
        boost::shared_ptr<Dispersion> dashD_;

        int npoints_;
        int deriv_;
};

}}

#endif 
