/* class ECP contains the contracted expansion of primitive GaussianECPs that define a particular ECP
   class GaussianECP is simply a data structure for the primitive gaussian parameters.  
   class ECPBasis is just a glorified list of all ECPs being used
   These are just skeletons to be replaced at a later date. 

	Robert A. Shaw 2016

	REFERENCES:
	L.E. McMurchie, E.R. Davidson, J. Comput. Phys. 44 (1981), 289 - 301
	R. Flores-Moreno et al., J. Comp. Chem. 27 (2006), 1009 - 1019
 */

#ifndef ECP_HEAD
#define ECP_HEAD

#include <vector>
#include <string>
#include "psi4/libpsi4util/exception.h"
#include "basisset.h"
#include "basisset_parser.h"
#include "gshell.h"

namespace psi {

    class ECPShellInfo : public ShellInfo {
    private:
        std::vector<int> n_;
        std::vector<int> sub_l_;
    public:
        ECPShellInfo(int am,
                     const std::vector<double>& c,
                     const std::vector<double>& e,
                     const std::vector<int>& n,
                     const std::vector<int>& sub_l,
                     int nc,
                     const Vector3& center,
                     int start);
        
        ECPShellInfo copy();
        ECPShellInfo copy(int nc, const Vector3& c);
        
        int n(int prim) const { return n_[prim]; }
        int subl(int prim) const { return sub_l_[prim]; }
        
        bool operator==(const ECPShellInfo& RHS) const;
        bool operator!=(const ECPShellInfo& RHS) const {
            return !(*this==RHS);
        }
    };

class ECPBasisSetParser : public BasisSetParser {
public:
    ECPBasisSetParser() : BasisSetParser() { }
    virtual ~ECPBasisSetParser() { }
    
    std::vector<ECPShellInfo> ecp_parse(const std::string& symbol, const std::string& dataset) {
        return ecp_parse(symbol, string_to_vector(dataset));
    }
    
    virtual std::vector<ECPShellInfo> ecp_parse(const std::string& symbol, const std::vector<std::string>& dataset) = 0;

};
    
class Gaussian94ECPBasisSetParser : public ECPBasisSetParser {
public:
    Gaussian94ECPBasisSetParser() : ECPBasisSetParser() {}
    virtual std::vector<ECPShellInfo> ecp_parse(const std::string& symbol, const std::vector<std::string>& dataset);
};

    
class GaussianECPShell : public GaussianShell {
private:
    const int* n_;
    const int* sub_l_;
	
public:
	GaussianECPShell(int am,
                     int nprimitive,
                     const double *oc,
                     const double *e,
                     const int *n,
                     const int *subl,
                     int nc,
                     const double* center,
                     int start);
    GaussianECPShell() { }
	
    int n(int prim) const { return n_[prim]; }
    int subl(int prim) const { return sub_l_[prim]; }
	
	// Evaluate U_l(r)
	double evaluate(double r, int l) const;
	
};

class ECPBasisSet : public BasisSet {
    friend class ECPBasisSetParser;
    
    int *uns_;
    int *usubls_;
    
public:
    ECPBasisSet();
    ECPBasisSet(const std::string &basistype, SharedMolecule mol,
                std::map<std::string, std::map<std::string, std::vector<ECPShellInfo> > > &shell_map);

    /** Returns a new basis set object
     * Constructs an ECP basis set from the parsed information
     *
     * @param mol           Psi4 molecule.  WARNING: The nuclear charges are modified by this routine
     * @param py::dict      Python dictionary containing the basis information
     * @param forced_puream Force puream or not
    **/
    static std::shared_ptr<ECPBasisSet> construct_ecp_from_pydict(std::shared_ptr <Molecule> mol, py::dict pybs, const int forced_puream);

	
    
};

}

#endif
