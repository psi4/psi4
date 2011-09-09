#ifndef three_index_ps_H
#define three_index_ps_H

#include <psi4-dec.h>

namespace boost {
template <class T>
class shared_ptr;
}

namespace psi {

class PSIO;
class BasisSet;
class Matrix;
class Vector;
class IntVector;
class Vector3;
class PseudospectralInt;
class PseudospectralGrid;
class GridBlock;

class DealiasBasisSet {

protected:
    /// The options object
    Options& options_; 
    /// The primary basis set
    boost::shared_ptr<BasisSet> primary_;
    /// The effective alphas of the primary set [center][am][index]
    std::vector<std::vector<std::vector<double> > > primary_alpha_;
    /// The alphas of the dealias set [center][am][index]
    std::vector<std::vector<std::vector<double> > > dealias_alpha_;
    /// The resultant dealias set
    boost::shared_ptr<BasisSet> dealias_;

    // => Parameters <= //
    /// Base for core and diffuse functions (~2)
    double delta_;
    /// Base for even-tempered cap functions (~3.5)
    double beta_;
    /// Number of intercalaters per window (~1)
    int nintercalater_;
    /// Number of core functions per block (~1)
    int ncore_;
    /// Number of diffuse functions per block (~1)
    int ndiffuse_;
    /// Number of cap functions (~1)
    int ncap_;
    /// Number of higher cardinal numbers to cap with
    int nl_; 

    /// Helper functions
    void form_primary_alpha();
    void form_core();
    void form_diffuse();
    void form_intercalater();
    void form_cap();
    void form_basis();

public:
    DealiasBasisSet(boost::shared_ptr<BasisSet> primary_, Options& options);
    virtual ~DealiasBasisSet();

    /// Parameter entry
    void setDelta(double delta) { delta_ = delta; }
    void setBeta(double beta) { beta_ = beta; }
    void setNCore(double n) { ncore_ = n; }
    void setNCap(double n) { ncap_ = n; }
    void setNL(double n) { nl_ = n; }
    void setNIntercalater(double n) { nintercalater_ = n; }
    void setNDiffuse(double n) { ndiffuse_ = n; }
   
    /// Master build routine 
    void buildDealiasBasisSet();

    /// Convenience routine
    boost::shared_ptr<BasisSet> dealiasSet() const { return dealias_; }
};

/*- PSTensorII: Range-Separated Pseudospectral Techniques
 * This class will eventually be refactored to PSTensor
 * -*/
class PSTensorII {

protected:
    /// Debug level
    int debug_;
    /// Print level
    int print_;
    /// Memory available to this algorithm, in doubles
    unsigned long int memory_;

    /// Molecule (for convenience)
    boost::shared_ptr<Molecule> molecule_;
    /// Primary basis set
    boost::shared_ptr<BasisSet> primary_;
    /// Dealias basis set
    boost::shared_ptr<BasisSet> dealias_;
    /// Pseudospectral grid
    boost::shared_ptr<PseudospectralGrid> grid_;
    /// options reference
    Options& options_;

    /// Full C matrix (must provide orthonormal MO basis)
    boost::shared_ptr<Matrix> C_;
    /// Active occupied C Matrix (for convenience)
    boost::shared_ptr<Matrix> Caocc_;
    /// Active virtual C Matrix (for convenience)
    boost::shared_ptr<Matrix> Cavir_;

    /// Dealias or not?
    bool do_dealias_;
    /// Omega or not?
    bool do_omega_;
    /// Renormalize or not? False implies do_dealias_ = false
    bool do_renormalize_;

    /// Omega value to use. 
    double omega_;
    /// Minimum eigenvalue of a function to keep in the dealias basis
    double min_S_dealias_;
    /// Alpha value in rotated fitting. 1.0 is conventional
    double alpha_; 

    /// Build the PseudospectralGrid grid_
    void buildGrid();
    /// Build the BasisSet dealias_
    void buildDealiasSet();

    /// Build using the naive, nondealiased quadrature algorithm 
    void buildQuadrature();
    /// Build using the renormalized, nondealiased algorithm
    void buildRenormalized();
    /// Build using the renormalized, dealiased algorithm
    void buildDealiased();

    /// Print header, including algorithms to use
    void print_header();
    /// Initialize, build R and Q
    void common_init();

    /// Build the primary AO collocation matrix
    void form_Rpao();
    /// Build the primary MO collocation matrix
    void form_Rpmo();
    /// Build the dealias AO collocation matrix
    void form_Rdao();
    /// Build the dealias MO collocation matrix
    void form_Rdmo();

    /// Build the Primary AO x Dealias AO overlap matrix 
    void form_Spdao();
    /// Build the Primary MO x Dealias AO overlap matrix 
    void form_Spdmo();
    /// Build the Dealias AO x Dealias AO overlap matrix 
    void form_Sddao();
    /// Build the Dealias OO x Dealias OO overlap matrix 
    void form_Sddoo();
    /// Build the Dealias orthonormalization matrix
    void form_Cdd();

    /// Build the Q R pair from quadrature
    void form_Q_quadrature();
    /// Build the Q R pair from renormalization 
    void form_Q_renormalized();
    /// Build the Q R pair from dealiased renormalization 
    void form_Q_dealiased();

    /// Primary AO functions
    int nso_;
    /// Primary MO functions
    int nmo_;
    /// Dealias AO functions
    int ndso_;
    /// Dealias MO functions
    int ndmo_;

    /// Number of grid points
    int naux_;
    
    /// Number of frozen occupieds
    int nfocc_; 
    /// Total number of occupieds
    int nocc_; 
    /// Number of active occupieds
    int naocc_; 
    /// Number of frozen virtuals
    int nfvir_; 
    /// Total number of virtuals
    int nvir_; 
    /// Number of active virtuals
    int navir_; 

    /// Grid weights (hopefully SPD) 
    boost::shared_ptr<Vector> w_;
    /// Primary AO collocation matrix
    boost::shared_ptr<Matrix> Rpao_;
    /// Primary MO collocation matrix
    boost::shared_ptr<Matrix> Rpmo_;
    /// Dealias AO collocation matrix
    boost::shared_ptr<Matrix> Rdao_;
    /// Dealias MO collocation matrix
    boost::shared_ptr<Matrix> Rdmo_;
   
    /// Primary AO x Dealias MO overlap matrix
    boost::shared_ptr<Matrix> Spdao_;
    /// Dealias AO x Dealias AO overlap matrix
    boost::shared_ptr<Matrix> Sddao_;
    /// Dealias OO x Dealias OO overlap matrix
    boost::shared_ptr<Matrix> Sddoo_;
 
    /// Primary MO x Dealias AO overlap matrix
    boost::shared_ptr<Matrix> Spdmo_;
    /// Dealias Orthogonalization matrix
    boost::shared_ptr<Matrix> Cdd_;

    /// Target Q tensor (nmo x naux)
    boost::shared_ptr<Matrix> Qmo_;
    /// Target R tensor (nmo x naux)
    boost::shared_ptr<Matrix> Rmo_;

    /// Possible target O tensor (naux x naux)
    boost::shared_ptr<Matrix> O_;

public:
    PSTensorII(boost::shared_ptr<BasisSet> primary, 
             boost::shared_ptr<Matrix> C,
             int nocc,
             int nvir,
             int naocc,
             int navir,
             unsigned long int memory,
             Options& options);
    ~PSTensorII();

    boost::shared_ptr<Matrix> Q();
    boost::shared_ptr<Matrix> Qocc();
    boost::shared_ptr<Matrix> Qvir();
    boost::shared_ptr<Matrix> Qaocc();
    boost::shared_ptr<Matrix> Qavir();

    boost::shared_ptr<Matrix> R();
    boost::shared_ptr<Matrix> Rocc();
    boost::shared_ptr<Matrix> Rvir();
    boost::shared_ptr<Matrix> Raocc();
    boost::shared_ptr<Matrix> Ravir();

    boost::shared_ptr<Matrix> O();

    // Begin naive routines for debugging    
    boost::shared_ptr<Matrix> Aso();
    boost::shared_ptr<Matrix> Amo();
    boost::shared_ptr<Matrix> Aoo();
    boost::shared_ptr<Matrix> Aov();
    boost::shared_ptr<Matrix> Avv();

    boost::shared_ptr<Matrix> Imo();
    boost::shared_ptr<Matrix> Ipsmo();
    boost::shared_ptr<Matrix> Idpsmo();
    // End naive routines for debugging
};

class PSTensor {

protected:

    /// Debug level
    int debug_;
    /// Print level
    int print_;

    /// Molecule (fo convenience)
    boost::shared_ptr<Molecule> molecule_;
    /// Primary basis set
    boost::shared_ptr<BasisSet> primary_;
    /// Dealias basis set
    boost::shared_ptr<BasisSet> dealias_;
    /// Pseudospectral grid
    boost::shared_ptr<PseudospectralGrid> grid_;
    /// options reference
    Options& options_;

    /// Full C matrix (must provide orthonormal MO basis)
    boost::shared_ptr<Matrix> C_;
    /// Active occupied C Matrix (for convenience)
    boost::shared_ptr<Matrix> Caocc_;
    /// Active virtual C Matrix (for convenience)
    boost::shared_ptr<Matrix> Cavir_;

    /// Minimum eigenvalue to keep in the primary basis
    double min_S_primary_;
    /// Minimum eigenvalue to keep in the dealias basis
    double min_S_dealias_;

    /// Omega to use for low-pass smoothing
    double omega_;
    /// Use omega integrals or not?
    bool use_omega_;

    /// Number of AO primary functions
    int nso_;
    /// Number of MO primary functions
    int nmo_;
    /// Number of AO dealias functions
    int ndealias_;
    /// nmo_ + ndealias_
    int naug_;

    /// Number of primary functions felt by the grid
    int nmo2_;
    /// Number of dealias functions felt by the grid
    int ndealias2_;

    /// Number of grid points
    int naux_;

    /// Number of frozen occupieds
    int nfocc_; 
    /// Total number of occupieds
    int nocc_; 
    /// Number of active occupieds
    int naocc_; 
    /// Number of frozen virtuals
    int nfvir_; 
    /// Total number of virtuals
    int nvir_; 
    /// Number of active virtuals
    int navir_; 

    void common_init();
    void print_header();
    void buildDealiasSet();
    void buildGrid();
    void buildR();

    void form_Spdao(); 
    void form_Spdmo();

    void form_Rpao();
    void form_Rdao();
    void form_Rpmo();
    void form_Rdmo();
    void form_Ra();

    void buildQ();
    void buildQ_renormalized();
    void buildQ_canonical();
    void buildQ_quadrature();
    void form_Cpp();
    void form_U();
    void form_Cpd();
    void form_V();
    void form_Cdd();
    void form_W();
    void form_X();
    void form_Q();
    void validate_X();

    /// Grid weights (hopefully SPD) 
    boost::shared_ptr<Vector> w_;

    /// Target Q tensor (nmo x naux)
    boost::shared_ptr<Matrix> Qmo_;
    /// Target R tensor (nmo x naux)
    boost::shared_ptr<Matrix> Rmo_;

    /// AO-basis overlap (nso x dealias)
    boost::shared_ptr<Matrix> Spdao_;
    /// Overlap matrix (nmo x dealias)
    boost::shared_ptr<Matrix> Spdmo_;

    /// AO primary collocation matrix (primary x naux)
    boost::shared_ptr<Matrix> Rpao_;
    /// AO dealias collocation matrix (dealias x naux)
    boost::shared_ptr<Matrix> Rdao_;
    /// MO primary collocation matrix (primary' x naux)
    boost::shared_ptr<Matrix> Rpmo_;
    /// MO dealias collocation matrix (dealias x naux)
    boost::shared_ptr<Matrix> Rdmo_;
    /// Finished augmented collocation matrix (naug x naux)
    boost::shared_ptr<Matrix> Ra_;

    boost::shared_ptr<Matrix> Cpp_;
    boost::shared_ptr<Matrix> Cpd_;
    boost::shared_ptr<Matrix> Cdd_;

    boost::shared_ptr<Matrix> U_;
    boost::shared_ptr<Matrix> V_;
    boost::shared_ptr<Matrix> W_;
    boost::shared_ptr<Matrix> X_;

public:

    PSTensor(boost::shared_ptr<BasisSet> primary, 
             boost::shared_ptr<Matrix> C,
             int nocc,
             int nvir,
             int naocc,
             int navir,
             Options& options,
             double omega = -1.0);
    ~PSTensor();

    boost::shared_ptr<Matrix> Q();
    boost::shared_ptr<Matrix> Qocc();
    boost::shared_ptr<Matrix> Qvir();
    boost::shared_ptr<Matrix> Qaocc();
    boost::shared_ptr<Matrix> Qavir();

    boost::shared_ptr<Matrix> R();
    boost::shared_ptr<Matrix> Rocc();
    boost::shared_ptr<Matrix> Rvir();
    boost::shared_ptr<Matrix> Raocc();
    boost::shared_ptr<Matrix> Ravir();
    
    boost::shared_ptr<Matrix> Aso();
    boost::shared_ptr<Matrix> Amo();
    boost::shared_ptr<Matrix> Aoo();
    boost::shared_ptr<Matrix> Aov();
    boost::shared_ptr<Matrix> Avv();

    boost::shared_ptr<Matrix> Imo();
    boost::shared_ptr<Matrix> Ipsmo();
};


// Class to demonstrate pseudospectral techniques
// given the current molecule and options.
//
// "C'mon you apes, you want to live forever?!"
//
class PseudoTrial {

protected:

    // => Starter stuff <= //
    
    // Debug flag
    int debug_;
    // Print flag
    int print_;
    // options 
    Options& options_;
    // Molecule
    boost::shared_ptr<Molecule> molecule_;

    // => Bases/Grids <= // 

    // Dealias or not? 
    bool do_dealias_;
    // Minimum eigenvalue for a primary basis function
    double min_S_primary_;
    // Minimum eigenvalue for a dealias basis function
    double min_S_dealias_;
    // Primary basis set
    boost::shared_ptr<BasisSet> primary_;
    // Dealias basis set
    boost::shared_ptr<BasisSet> dealias_;
    // Pseudospectral grid
    boost::shared_ptr<PseudospectralGrid> grid_;
    // Number of primary basis functions
    int nso_;
    // Number of orthogonalized primary basis functions
    int nmo_;
    // Number of dealias basis functions
    int ndealias_;
    // Number of orthogonalized dealias basis functions
    int ndealias2_;
    // Number of primary + dealias basis functions (raw)
    int naug_;
    // Number of primary + dealias basis functions (finished)
    int naug2_;
    // Number of grid points
    int naux_;

    // => Temps <= // 

    // Raw S matrices
    
    // Overlap matrix (primary x primary)
    boost::shared_ptr<Matrix> Spp_;
    // Overlap matrix (primary x dealias)
    boost::shared_ptr<Matrix> Spd_;
    // Overlap matrix (dealias x dealias)
    boost::shared_ptr<Matrix> Sdd_;
    // Augmented Overlap matrix (aug x aug)
    boost::shared_ptr<Matrix> Sa_;

    // Orthonormal primary
    
    // Overlap matrix (primary x dealias)
    boost::shared_ptr<Matrix> Spd3_;
    // Augmented Overlap matrix (aug x aug)
    boost::shared_ptr<Matrix> Sa3_;

    // Orthogonal primary - dealias

    // Overlap matrix (dealias x dealias)
    boost::shared_ptr<Matrix> Sdd4_;
    // Augmented Overlap matrix (aug x aug) (finished)
    boost::shared_ptr<Matrix> Sa4_;

    // Orthonormal primary - dealias
    
    // Augmented Overlap matrix (aug' x aug') (finished)
    boost::shared_ptr<Matrix> Sa2_;    

    // X matrix, primary (primary x primary')
    boost::shared_ptr<Matrix> Xpp_;
    // X matrix, dealias (dealias x dealias')
    boost::shared_ptr<Matrix> Xdd_;
    // Orthogonalization matrix (dealias x primary')
    boost::shared_ptr<Matrix> Cdp_;

    // Collocation matrix (primary)
    boost::shared_ptr<Matrix> Rp_;
    // Collocation matrix (dealias)
    boost::shared_ptr<Matrix> Rd_;
    // Collocation matrix (primary')
    boost::shared_ptr<Matrix> Rp2_;
    // Collocation matrix (dealias')
    boost::shared_ptr<Matrix> Rd2_;
    // Collocation matrix (augmented')
    boost::shared_ptr<Matrix> Ra_;

    // Weight Vector 
    boost::shared_ptr<Vector> w_;
    // C matrix
    boost::shared_ptr<Matrix> C_;
    // Cinv matrix
    boost::shared_ptr<Matrix> Cinv_;
    // Full Q matrix
    boost::shared_ptr<Matrix> Qfull_;
    // Full Q matrix
    boost::shared_ptr<Matrix> Qmo_;
    // Projector matrix 
    boost::shared_ptr<Matrix> P_;
    // Transformer matrix
    boost::shared_ptr<Matrix> SX_;

    // => Targets <= //

    // Q_m^P
    boost::shared_ptr<Matrix> Q_;
    // R_n^P
    boost::shared_ptr<Matrix> R_;
    // A_ls^P
    boost::shared_ptr<Matrix> A_;
    // QR_mn^P 
    boost::shared_ptr<Matrix> T_;

    // => Final Targets <= //  

    // AO basis (mn|ls) tensor (exact)
    boost::shared_ptr<Matrix> I_; 
    // AO basis (mn|ls) tensor (PS)
    boost::shared_ptr<Matrix> Ips_; 

    // => Helpers <= //
    
    // Build everything 
    void common_init();
    void print_header();

    void form_molecule();
    void form_bases();
    void form_grid();

    void form_Spp();
    void form_Spd();
    void form_Sdd();

    void form_Sa();
    void form_Sa3();
    void form_Sa4();
    void form_Sa2();

    void form_Xpp();
    void form_Spd3();
    void form_Cdp();
    void form_Sdd4();
    void form_Xdd();

    void form_Rp();
    void form_Rd();
    void form_Rp2();
    void form_Rd2();
    void form_Ra();
    
    void form_Q();
    void form_P();
    void form_SX();
    void form_A(); 

    void form_I();
    void form_Ips();
    void verify(); 

public:

    PseudoTrial();
    ~PseudoTrial();
    
    boost::shared_ptr<Matrix> getI() const;
    boost::shared_ptr<Matrix> getIPS() const;

    boost::shared_ptr<Matrix> getQ() const;
    boost::shared_ptr<Matrix> getR() const;
    boost::shared_ptr<Matrix> getA() const;

};

}
#endif
