  #
  # @BEGIN LICENSE
  #
  # Psi4: an open-source quantum chemistry software package
  #
  # Copyright (c) 2007-2016 The Psi4 Developers.
  #
  # The copyrights for code used from other parties are included in
  # the corresponding files.
  #
  # This program is free software; you can redistribute it and/or modify
  # it under the terms of the GNU General Public License as published by
  # the Free Software Foundation; either version 2 of the License, or
  # (at your option) any later version.
  #
  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.
  #
  # You should have received a copy of the GNU General Public License along
  # with this program; if not, write to the Free Software Foundation, Inc.,
  # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
  #
  # @END LICENSE
  #

"""
List of SuperFunctionals that LibXC natively provides.
"""

from psi4 import core

xc_func_list = [
    "XC_LDA_XC_TETER93",             #   20  /*Teter 93 parametrization                                     */
    "XC_LDA_XC_ZLP",                 #   43  /*Zhao, Levy & Parr, Eq. (20)                                  */
    "XC_LDA_XC_KSDT",                #  259  /*Karasiev et al. parametrization                              */
    "XC_GGA_XC_OPBE_D",              #   65  /*oPBE_D functional of Goerigk and Grimme                      */
    "XC_GGA_XC_OPWLYP_D",            #   66  /*oPWLYP-D functional of Goerigk and Grimme                    */
    "XC_GGA_XC_OBLYP_D",             #   67  /*oBLYP-D functional of Goerigk and Grimme                     */
    "XC_GGA_XC_HCTH_407P",           #   93  /*HCTH/407+                                                    */
    "XC_GGA_XC_HCTH_P76",            #   94  /*HCTH p=7/6                                                   */
    "XC_GGA_XC_HCTH_P14",            #   95  /*HCTH p=1/4                                                   */
    "XC_GGA_XC_B97_GGA1",            #   96  /*Becke 97 GGA-1                                               */
    "XC_GGA_XC_KT2",                 #  146  /*Keal and Tozer version 2                                     */
    "XC_GGA_XC_TH1",                 #  154  /*Tozer and Handy v. 1                                         */
    "XC_GGA_XC_TH2",                 #  155  /*Tozer and Handy v. 2                                         */
    "XC_GGA_XC_TH3",                 #  156  /*Tozer and Handy v. 3                                         */
    "XC_GGA_XC_TH4",                 #  157  /*Tozer and Handy v. 4                                         */
    "XC_GGA_XC_HCTH_93",             #  161  /*HCTH functional fitted to  93 molecules                      */
    "XC_GGA_XC_HCTH_120",            #  162  /*HCTH functional fitted to 120 molecules                      */
    "XC_GGA_XC_HCTH_147",            #  163  /*HCTH functional fitted to 147 molecules                      */
    "XC_GGA_XC_HCTH_407",            #  164  /*HCTH functional fitted to 407 molecules                      */
    "XC_GGA_XC_EDF1",                #  165  /*Empirical functionals from Adamson, Gill, and Pople          */
    "XC_GGA_XC_XLYP",                #  166  /*XLYP functional                                              */
    "XC_GGA_XC_B97_D",               #  170  /*Grimme functional to be used with C6 vdW term                */
    "XC_GGA_XC_PBE1W",               #  173  /*Functionals fitted for water                                 */
    "XC_GGA_XC_MPWLYP1W",            #  174  /*Functionals fitted for water                                 */
    "XC_GGA_XC_PBELYP1W",            #  175  /*Functionals fitted for water                                 */
    "XC_GGA_XC_MOHLYP",              #  194  /*Functional for organometallic chemistry                      */
    "XC_GGA_XC_MOHLYP2",             #  195  /*Functional for barrier heights                               */
    "XC_GGA_XC_TH_FL",               #  196  /*Tozer and Handy v. FL                                        */
    "XC_GGA_XC_TH_FC",               #  197  /*Tozer and Handy v. FC                                        */
    "XC_GGA_XC_TH_FCFO",             #  198  /*Tozer and Handy v. FCFO                                      */
    "XC_GGA_XC_TH_FCO",              #  199  /*Tozer and Handy v. FCO                                       */
    "XC_GGA_XC_VV10",                #  255  /*Vydrov and Van Voorhis                                       */
    "XC_HYB_GGA_XC_B97_1p",          #  266  /*version of B97 by Cohen and Handy                            */
    "XC_HYB_GGA_XC_B3PW91",          #  401  /*The original (ACM) hybrid of Becke                           */
    "XC_HYB_GGA_XC_B3LYP",           #  402  /*The (in)famous B3LYP                                         */
    "XC_HYB_GGA_XC_B3P86",           #  403  /*Perdew 86 hybrid similar to B3PW91                           */
    "XC_HYB_GGA_XC_O3LYP",           #  404  /*hybrid using the optx functional                             */
    "XC_HYB_GGA_XC_mPW1K",           #  405  /*mixture of mPW91 and PW91 optimized for kinetics             */
    "XC_HYB_GGA_XC_PBEH",            #  406  /*aka PBE0 or PBE1PBE                                          */
    "XC_HYB_GGA_XC_B97",             #  407  /*Becke 97                                                     */
    "XC_HYB_GGA_XC_B97_1",           #  408  /*Becke 97-1                                                   */
    "XC_HYB_GGA_XC_B97_2",           #  410  /*Becke 97-2                                                   */
    "XC_HYB_GGA_XC_X3LYP",           #  411  /*hybrid by Xu and Goddard                                     */
    "XC_HYB_GGA_XC_B1WC",            #  412  /*Becke 1-parameter mixture of WC and PBE                      */
    "XC_HYB_GGA_XC_B97_K",           #  413  /*Boese-Martin for Kinetics                                    */
    "XC_HYB_GGA_XC_B97_3",           #  414  /*Becke 97-3                                                   */
    "XC_HYB_GGA_XC_MPW3PW",          #  415  /*mixture with the mPW functional                              */
    "XC_HYB_GGA_XC_B1LYP",           #  416  /*Becke 1-parameter mixture of B88 and LYP                     */
    "XC_HYB_GGA_XC_B1PW91",          #  417  /*Becke 1-parameter mixture of B88 and PW91                    */
    "XC_HYB_GGA_XC_mPW1PW",          #  418  /*Becke 1-parameter mixture of mPW91 and PW91                  */
    "XC_HYB_GGA_XC_MPW3LYP",         #  419  /*mixture of mPW and LYP                                       */
    "XC_HYB_GGA_XC_SB98_1a",         #  420  /*Schmider-Becke 98 parameterization 1a                        */
    "XC_HYB_GGA_XC_SB98_1b",         #  421  /*Schmider-Becke 98 parameterization 1b                        */
    "XC_HYB_GGA_XC_SB98_1c",         #  422  /*Schmider-Becke 98 parameterization 1c                        */
    "XC_HYB_GGA_XC_SB98_2a",         #  423  /*Schmider-Becke 98 parameterization 2a                        */
    "XC_HYB_GGA_XC_SB98_2b",         #  424  /*Schmider-Becke 98 parameterization 2b                        */
    "XC_HYB_GGA_XC_SB98_2c",         #  425  /*Schmider-Becke 98 parameterization 2c                        */
    "XC_HYB_GGA_XC_HSE03",           #  427  /*the 2003 version of the screened hybrid HSE                  */
    "XC_HYB_GGA_XC_HSE06",           #  428  /*the 2006 version of the screened hybrid HSE                  */
    "XC_HYB_GGA_XC_HJS_PBE",         #  429  /*HJS hybrid screened exchange PBE version                     */
    "XC_HYB_GGA_XC_HJS_PBE_SOL",     #  430  /*HJS hybrid screened exchange PBE_SOL version                 */
    "XC_HYB_GGA_XC_HJS_B88",         #  431  /*HJS hybrid screened exchange B88 version                     */
    "XC_HYB_GGA_XC_HJS_B97X",        #  432  /*HJS hybrid screened exchange B97x version                    */
    "XC_HYB_GGA_XC_CAM_B3LYP",       #  433  /*CAM version of B3LYP                                         */
    "XC_HYB_GGA_XC_TUNED_CAM_B3LYP", #  434  /*CAM version of B3LYP tuned for excitations                   */
    "XC_HYB_GGA_XC_BHANDH",          #  435  /*Becke half-and-half                                          */
    "XC_HYB_GGA_XC_BHANDHLYP",       #  436  /*Becke half-and-half with B88 exchange                        */
    "XC_HYB_GGA_XC_MB3LYP_RC04",     #  437  /*B3LYP with RC04 LDA                                          */
    "XC_HYB_GGA_XC_MPWLYP1M",        #  453  /*MPW with 1 par. for metals/LYP                               */
    "XC_HYB_GGA_XC_REVB3LYP",        #  454  /*Revised B3LYP                                                */
    "XC_HYB_GGA_XC_CAMY_BLYP",       #  455  /*BLYP with yukawa screening                                   */
    "XC_HYB_GGA_XC_PBE0_13",         #  456  /*PBE0-1/3                                                     */
    "XC_HYB_GGA_XC_B3LYPs",          #  459  /*B3LYP* functional                                            */
    "XC_HYB_GGA_XC_WB97",            #  463  /*Chai and Head-Gordon                                         */
    "XC_HYB_GGA_XC_WB97X",           #  464  /*Chai and Head-Gordon                                         */
    "XC_HYB_GGA_XC_LRC_WPBEH",       #  465  /*Long-range corrected functional by Rorhdanz et al            */
    "XC_HYB_GGA_XC_WB97X_V",         #  466  /*Mardirossian and Head-Gordon                                 */
    "XC_HYB_GGA_XC_LCY_PBE",         #  467  /*PBE with yukawa screening                                    */
    "XC_HYB_GGA_XC_LCY_BLYP",        #  468  /*BLYP with yukawa screening                                   */
    "XC_HYB_GGA_XC_LC_VV10",         #  469  /*Vydrov and Van Voorhis                                       */
    "XC_HYB_GGA_XC_CAMY_B3LYP",      #  470  /*B3LYP with Yukawa screening                                  */
    "XC_HYB_GGA_XC_WB97X_D",         #  471  /*Chai and Head-Gordon                                         */
    "XC_HYB_GGA_XC_HPBEINT",         #  472  /*hPBEint                                                      */
    "XC_HYB_GGA_XC_LRC_WPBE",        #  473  /*Long-range corrected functional by Rorhdanz et al            */
    "XC_HYB_GGA_XC_B3LYP5",          #  475  /*B3LYP with VWN functional 5 instead of RPA                   */
    "XC_HYB_GGA_XC_EDF2",            #  476  /*Empirical functional from Lin, George and Gill               */
    "XC_HYB_GGA_XC_CAP0",            #  477  /*Correct Asymptotic Potential hybrid                          */
    "XC_MGGA_XC_ZLP",                #   42  /*Zhao, Levy & Parr, Eq. (21)                                  */
    "XC_MGGA_XC_OTPSS_D",            #   64  /*oTPSS_D functional of Goerigk and Grimme                     */
    "XC_MGGA_XC_TPSSLYP1W",          #  242  /*Functionals fitted for water                                 */
    "XC_MGGA_XC_B97M_V",             #  254  /*Mardirossian and Head-Gordon                                 */
    "XC_HYB_MGGA_XC_M05",            #  438  /*M05 functional from Minnesota                                */
    "XC_HYB_MGGA_XC_M05_2X",         #  439  /*M05-2X functional from Minnesota                             */
    "XC_HYB_MGGA_XC_B88B95",         #  440  /*Mixture of B88 with BC95 (B1B95)                             */
    "XC_HYB_MGGA_XC_B86B95",         #  441  /*Mixture of B86 with BC95                                     */
    "XC_HYB_MGGA_XC_PW86B95",        #  442  /*Mixture of PW86 with BC95                                    */
    "XC_HYB_MGGA_XC_BB1K",           #  443  /*Mixture of B88 with BC95 from Zhao and Truhlar               */
    "XC_HYB_MGGA_XC_M06_HF",         #  444  /*M06-HF functional from Minnesota                             */
    "XC_HYB_MGGA_XC_MPW1B95",        #  445  /*Mixture of mPW91 with BC95 from Zhao and Truhlar             */
    "XC_HYB_MGGA_XC_MPWB1K",         #  446  /*Mixture of mPW91 with BC95 for kinetics                      */
    "XC_HYB_MGGA_XC_X1B95",          #  447  /*Mixture of X with BC95                                       */
    "XC_HYB_MGGA_XC_XB1K",           #  448  /*Mixture of X with BC95 for kinetics                          */
    "XC_HYB_MGGA_XC_M06",            #  449  /*M06 functional from Minnesota                                */
    "XC_HYB_MGGA_XC_M06_2X",         #  450  /*M06-2X functional from Minnesota                             */
    "XC_HYB_MGGA_XC_PW6B95",         #  451  /*Mixture of PW91 with BC95 from Zhao and Truhlar              */
    "XC_HYB_MGGA_XC_PWB6K",          #  452  /*Mixture of PW91 with BC95 from Zhao and Truhlar for kinetics */
    "XC_HYB_MGGA_XC_TPSSH",          #  457  /*TPSS hybrid                                                  */
    "XC_HYB_MGGA_XC_REVTPSSH",       #  458  /*revTPSS hybrid                                               */
    "XC_HYB_MGGA_XC_M08_HX",         #  460  /*M08-HX functional from Minnesota                             */
    "XC_HYB_MGGA_XC_M08_SO",         #  461  /*M08-SO functional from Minnesota                             */
    "XC_HYB_MGGA_XC_M11",            #  462  /*M11    functional from Minnesota                             */
    "XC_HYB_MGGA_XC_WB97M_V",        #  531  /*Mardirossian and Head-Gordon                                 */
]


# filter out a few -D for now as we are missing them
xc_func_list.remove("XC_GGA_XC_OPWLYP_D")
xc_func_list.remove("XC_GGA_XC_OBLYP_D")

# filter out -V for now
xc_func_list = [x for x in xc_func_list if "_V" not in x]

# Deal with xc mix upper/lower case
lower_to_xc_dict = {x.lower() : x for x in xc_func_list}
# Translate to something a user would want to input
psi_to_xc_translate = {}
for x in xc_func_list:
    key = x.split('_XC_')[-1].replace('_', '-').lower()
    key = key.replace("hcth-", "hcth")

    psi_to_xc_translate[key] = x

# Add extra translations
psi_to_xc_translate["b97-0"] = "XC_HYB_GGA_XC_B97"
psi_to_xc_translate["hcth"] = "XC_GGA_XC_HCTH_93"


def find_xc_func_name(name):
    """
    XC names unfortunately come in a mix of upper and lower case. This will search through
    the translation dictionaries to find a reasonable match.
    """

    if name in xc_func_list:
        return name

    name = name.lower()
    if name in list(psi_to_xc_translate):
        return psi_to_xc_translate[name]

    if name in list(lower_to_xc_dict):
        return lower_to_xc_dict[name]

    raise KeyError("LibXC keyname %s was not found!")



def build_libxc_xc_func(name, npoints, deriv):



    xc_name = find_xc_func_name(name)
    sup = core.SuperFunctional.XC_build(xc_name)

    descr = "    " + name.upper() + " "
    if sup.is_gga():
        if sup.x_alpha() > 0:
            descr += "Hyb-GGA "
        else:
            descr += "GGA "
    descr += "Exchange-Correlation Functional\n"

    sup.set_description(descr)

    sup.set_max_points(npoints)
    sup.set_deriv(deriv)
    sup.set_name(name.upper())
    sup.allocate()

    # Figure out dispersion
    if "_D" in xc_name:
        if xc_name == "XC_HYB_GGA_XC_WB97X_D":
            disp = ('wB97', '-CHG')
        elif xc_name == "XC_GGA_XC_B97_D":
            disp = ('B97-D', '-d3zero')
        elif xc_name == "XC_GGA_XC_OPBE_D":
            disp = ('OPBE', '-d3zero')
        elif xc_name == "XC_MGGA_XC_OTPSS_D":
            disp = ('OTPSS', '-d3zero')
        else:
            raise KeyError("Dispersion for functional %s not found" % name)
    else:
        disp = False

    return (sup, False)


# Translation layer
libxc_xc_functional_list = {x : build_libxc_xc_func for x in xc_func_list}
libxc_xc_functional_list.update({x : build_libxc_xc_func for x in list(psi_to_xc_translate)})





