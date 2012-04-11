#ifndef SRC_LIB_KS_H
#define SRC_LIB_KS_H    

#include <libscf_solver/ks.h>
#include <libmints/typedefs.h>

namespace psi{

class Options;

namespace scf{


class LHF : public RKS {

protected:

    // Current LHF EXX energy
    double EXX_;

    // => Potentials (On Grid) <= //

    // Current Local Exchange potential on grid
    SharedVector V_X_;
    // Current Slater Exchange potential on grid
    SharedVector V_S_;
    // Current Correction potential on grid
    SharedVector V_C_;

    // => Density Fitting <= //

    // Auxiliary basis set
    boost::shared_ptr<BasisSet> auxiliary_;
    // D_ij^A fitted three-index tensor
    SharedMatrix build_Dij(); 

    // => Methods to build potentials <= //
    
    virtual void setup_V();    
    virtual void build_V_S();
    virtual void build_V_C();
    virtual void build_V_X();

    // => Overloaded Methods <= //

    virtual void form_G();
    virtual void form_V();
    virtual double compute_E();
    
    void common_init();
public:
    LHF(Options &options);
    virtual ~LHF();
};


}} // Namespaces

#endif // Header guard
