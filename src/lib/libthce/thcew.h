#ifndef THCEW_H
#define THCEW_H

#include <libmints/wavefunction.h>
#include <libmints/typedefs.h>
#include <map>

namespace psi {

class THCE;

class THCEW : public Wavefunction {

protected:

    // Energy map
    std::map<std::string, double> energies_;
    // THCE object
    boost::shared_ptr<THCE> thce_;

    void common_init();

public:
    THCEW();
    virtual ~THCEW();

}; 

/**
 * Starter class for RHF-based THCE correlated methods
 *
 * Dimensions in the THCE:
 * -nso
 * -nmo
 * -nocc
 * -nfocc
 * -naocc
 * -navir
 * -nfvir  
 * -nvir
 * -nact
 *
 * Optional Dimensions in the THCE:
 * -nw
 * -naux
 * -ngrid
 * -namp
 *
 * Tensors in the THCE:
 * -Cmo   [nmo   x nso]
 * -Cocc  [nocc  x nso]
 * -Cfocc [nfocc x nso]
 * -Caocc [naocc x nso]
 * -Cavir [navir x nso]
 * -Cfvir [nfvir x nso]
 * -Cvir  [nvir  x nso]
 *
 * -eps_mo   [nmo  ]
 * -eps_occ  [nocc ]
 * -eps_focc [nfocc]
 * -eps_aocc [naocc]
 * -eps_avir [navir]
 * -eps_fvir [nfvir]
 * -eps_vir  [nvir ]
 * 
 * Optional Tensors in the THCE:
 * -pi_i     [nw x naocc]
 * -pi_a     [nw x navir]
 * -Bii      [naocc x naocc x naux] [Disk]
 * -Bia      [naocc x navir x naux] [Disk]
 * -Bai      [navir x naocc x naux] [Disk]
 * -Baa      [navir x navir x naux] [Disk]
 * -Bpp      [nact  x nact  x naux] [Disk]
 * -Xi       [naocc x ngrid]
 * -Xa       [navir x ngrid]
 * -Ziiii    [ngrid x ngrid] [Swapped] 
 * -Ziiia    [ngrid x ngrid] [Swapped] 
 * -Ziiaa    [ngrid x ngrid] [Swapped] 
 * -Ziaia    [ngrid x ngrid] [Swapped] 
 * -Ziaaa    [ngrid x ngrid] [Swapped] 
 * -Zaaaa    [ngrid x ngrid] [Swapped] 
 * -Zpppp    [ngrid x ngrid] [Swapped] 
 * -Sii      [ngrid x ngrid] [Swapped]
 * -Sia      [ngrid x ngrid] [Swapped]
 * -Saa      [ngrid x ngrid] [Swapped]
 * -Spp      [ngrid x ngrid] [Swapped]
 * -Lii      [ngrid x naux]  [Swapped]
 * -Lia      [ngrid x naux]  [Swapped]
 * -Laa      [ngrid x naux]  [Swapped]
 * -Lpp      [ngrid x naux]  [Swapped]
 * -Ti       [naocc x namp]
 * -Ta       [navir x namp]
 * -STia     [namp x namp] [Swapped]
 **/
class RTHCEW : public THCEW {

protected:

    void common_init();    

    // => Laplace <= //
   
    // pi_i and pi_a [nw]
    void build_laplace(double delta, double omega = 0.0); 

    // => DF <= //

    // Bia [naux]
    void build_df_ia(boost::shared_ptr<BasisSet> auxiliary);
    // Bii, Bia, Bai, and Baa [naux]
    void build_df_act(boost::shared_ptr<BasisSet> auxiliary);
    // Bpp [naux]
    void build_df_pp(boost::shared_ptr<BasisSet> auxiliary);

    // => LS-THC <= // 

    // Xi, Xa, Ziaia [swapped], Sia [swapped], Lia [swapped] [naux,ngrid]
    void build_lsthc_ia(boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<Matrix> X);
    // Xi, Xa, Ziiii -> Zaaaa [swapped], Sii -> Saa [swapped], Lii -> Laa [swapped] [naux,ngrid]
    void build_lsthc_act(boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<Matrix> X);
    // Xi, Xa, Zpppp [swapped], Spp [swapped], Lpp [swapped] [naux,ngrid]
    void build_lsthc_pp(boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<Matrix> X);

    // Ti, Ta, STia [swapped] [namp]
    void build_meth_ia(boost::shared_ptr<Matrix> X);

public:
    RTHCEW();
    virtual ~RTHCEW();
};



} // End namespace

#endif

