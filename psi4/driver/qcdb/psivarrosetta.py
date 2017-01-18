#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
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

useme2psivar = {

    # <<<  DFT  >>>

    'DFT.usemeraw': 'DFT FUNCTIONAL TOTAL ENERGY',  # for herding. plays well with other uses?
    #'-nobas.DFTdX.usemedash': 'DISPERSION CORRECTION ENERGY',  # for herding. plays well with other uses?
    #'DHDFT.usemeraw': 'DOUBLE-HYBRID CORRECTION ENERGY',  # for herding. plays well with other uses?  # violation of conventions to get plain dhdft E!

    'blyp.usemeraw': 'BLYP FUNCTIONAL TOTAL ENERGY',
    'blypd2.usemedash': 'BLYP-D2 DISPERSION CORRECTION ENERGY',
    'blypd3.usemedash': 'BLYP-D3 DISPERSION CORRECTION ENERGY',
    'blypd3bj.usemedash': 'BLYP-D3(BJ) DISPERSION CORRECTION ENERGY',
    'blypd3m.usemedash': 'BLYP-D3M DISPERSION CORRECTION ENERGY',
    'blypd3mbj.usemedash': 'BLYP-D3M(BJ) DISPERSION CORRECTION ENERGY',

    'b3lyp.usemeraw': 'B3LYP FUNCTIONAL TOTAL ENERGY',
    'b3lypd2.usemedash': 'B3LYP-D2 DISPERSION CORRECTION ENERGY',
    'b3lypd3.usemedash': 'B3LYP-D3 DISPERSION CORRECTION ENERGY',
    'b3lypd3bj.usemedash': 'B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY',
    'b3lypxdm.usemedash': 'B3LYP-XDM DISPERSION CORRECTION ENERGY',
    'b3lypd3m.usemedash': 'B3LYP-D3M DISPERSION CORRECTION ENERGY',
    'b3lypd3mbj.usemedash': 'B3LYP-D3M(BJ) DISPERSION CORRECTION ENERGY',

    'b2plyp.usemeraw': 'B2PLYP TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def
    'b2plypd2.usemedash': 'B2PLYP-D2 DISPERSION CORRECTION ENERGY',
    'b2plypd3.usemedash': 'B2PLYP-D3 DISPERSION CORRECTION ENERGY',
    'b2plypd3bj.usemedash': 'B2PLYP-D3(BJ) DISPERSION CORRECTION ENERGY',
    'b2plypd3m.usemedash': 'B2PLYP-D3M DISPERSION CORRECTION ENERGY',
    'b2plypd3mbj.usemedash': 'B2PLYP-D3M(BJ) DISPERSION CORRECTION ENERGY',

    'b970.usemeraw': 'B970 FUNCTIONAL TOTAL ENERGY',
    'b970d2.usemedash': 'B970-D2 DISPERSION CORRECTION ENERGY',

    'b97.usemeraw': 'B97 FUNCTIONAL TOTAL ENERGY',
    'b97d2.usemedash': 'B97-D2 DISPERSION CORRECTION ENERGY',
    'b97d3.usemedash': 'B97-D3 DISPERSION CORRECTION ENERGY',
    'b97d3bj.usemedash': 'B97-D3(BJ) DISPERSION CORRECTION ENERGY',
    'b97d3m.usemedash': 'B97-D3M DISPERSION CORRECTION ENERGY',
    'b97d3mbj.usemedash': 'B97-D3M(BJ) DISPERSION CORRECTION ENERGY',

    'bp86.usemeraw': 'BP86 FUNCTIONAL TOTAL ENERGY',
    'bp86d2.usemedash': 'BP86-D2 DISPERSION CORRECTION ENERGY',
    'bp86d3.usemedash': 'BP86-D3 DISPERSION CORRECTION ENERGY',
    'bp86d3bj.usemedash': 'BP86-D3(BJ) DISPERSION CORRECTION ENERGY',
    'bp86d3m.usemedash': 'BP86-D3M DISPERSION CORRECTION ENERGY',
    'bp86d3mbj.usemedash': 'BP86-D3M(BJ) DISPERSION CORRECTION ENERGY',

    'wb97x.usemeraw': 'WB97X FUNCTIONAL TOTAL ENERGY',
    'wb97xd.usemeraw': 'WB97X-D TOTAL ENERGY',
    'wb97xd.usemedash': 'WB97X-D DISPERSION CORRECTION ENERGY',

    'wb97x2.usemeraw': 'WB97X-2 TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def

    'wb97xv.usemeraw': 'WB97X-V TOTAL ENERGY',

    'm052x.usemeraw': 'M05-2X FUNCTIONAL TOTAL ENERGY',
    'm052xd3.usemedash': 'M05-2X-D3 DISPERSION CORRECTION ENERGY',

    'm062x.usemeraw': 'M06-2X FUNCTIONAL TOTAL ENERGY',
    'm062xd3.usemedash': 'M06-2X-D3 DISPERSION CORRECTION ENERGY',

    'pbe.usemeraw': 'PBE FUNCTIONAL TOTAL ENERGY',
    'pbed2.usemedash': 'PBE-D2 DISPERSION CORRECTION ENERGY',
    'pbed3.usemedash': 'PBE-D3 DISPERSION CORRECTION ENERGY',
    'pbed3bj.usemedash': 'PBE-D3(BJ) DISPERSION CORRECTION ENERGY',
    'pbed3m.usemedash': 'PBE-D3M DISPERSION CORRECTION ENERGY',
    'pbed3mbj.usemedash': 'PBE-D3M(BJ) DISPERSION CORRECTION ENERGY',

    'pbe0.usemeraw': 'PBE0 FUNCTIONAL TOTAL ENERGY',
    'pbe0d2.usemedash': 'PBE0-D2 DISPERSION CORRECTION ENERGY',
    'pbe0d3.usemedash': 'PBE0-D3 DISPERSION CORRECTION ENERGY',
    'pbe0d3bj.usemedash': 'PBE0-D3(BJ) DISPERSION CORRECTION ENERGY',
    'pbe0d3m.usemedash': 'PBE0-D3M DISPERSION CORRECTION ENERGY',
    'pbe0d3mbj.usemedash': 'PBE0-D3M(BJ) DISPERSION CORRECTION ENERGY',

    'wpbe.usemeraw': 'WPBE FUNCTIONAL TOTAL ENERGY',
    'wpbed3.usemedash': 'WPBE-D3 DISPERSION CORRECTION ENERGY',
    'wpbed3bj.usemedash': 'WPBE-D3(BJ) DISPERSION CORRECTION ENERGY',
    'wpbed3m.usemedash': 'WPBE-D3M DISPERSION CORRECTION ENERGY',
    'wpbed3mbj.usemedash': 'WPBE-D3M(BJ) DISPERSION CORRECTION ENERGY',

    'xyg3.usemeraw': 'XYG3 TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def

    'vv10.usemeraw': 'VV10 FUNCTIONAL TOTAL ENERGY',

    'lcvv10.usemeraw': 'LC-VV10 FUNCTIONAL TOTAL ENERGY',

    'dsdpbep86.usemeraw': 'DSD-PBEP86 TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def  # also DSD technically implies -D
    'dsdpbep86d2.usemedash': 'DSD-PBEP86-D2 DISPERSION CORRECTION ENERGY',
    'dsdpbep86d3.usemedash': 'DSD-PBEP86-D3 DISPERSION CORRECTION ENERGY',
    'dsdpbep86d3bj.usemedash': 'DSD-PBEP86-D3(BJ) DISPERSION CORRECTION ENERGY',

    'm08hx.usemeraw': 'M08-HX FUNCTIONAL TOTAL ENERGY',
    'm08so.usemeraw': 'M08-SO FUNCTIONAL TOTAL ENERGY',
    'm11.usemeraw': 'M11 FUNCTIONAL TOTAL ENERGY',
    'm11l.usemeraw': 'M11L FUNCTIONAL TOTAL ENERGY',

    'pbe02.usemeraw': 'PBE0-2 TOTAL ENERGY',  # no psivar for fctl + dh, which would be the more restrictive def

    'dldf.usemeraw': 'DLDF FUNCTIONAL TOTAL ENERGY',
    'dldfd.usemedash': 'DLDF+D DISPERSION CORRECTION ENERGY',

    # <<<  WFN  >>>

    #'usemeraw': 'HF TOTAL ENERGY',
    'usemeraw': 'SCF TOTAL ENERGY',

    'mp2.usemecorl': 'MP2 CORRELATION ENERGY',
    'mp3.usemecorl': 'MP3 CORRELATION ENERGY',
    'mp4.usemecorl': 'MP4 CORRELATION ENERGY',
    'ccsd.usemecorl': 'CCSD CORRELATION ENERGY',
    'ccsdt.usemecorl': 'CCSD(T) CORRELATION ENERGY',
    'ccsdfullt.usemecorl': 'CCSDT CORRELATION ENERGY',
    'ccsdtq.usemecorl': 'CCSDT(Q) CORRELATION ENERGY',

    'fno.usemecrct': 'FNO CORRECTION ENERGY',
    'fnomp3.usemecorl': 'MP3 FNO CORRELATION ENERGY',
    'fnoccsd.usemecorl': 'CCSD FNO CORRELATION ENERGY',
    'fnoccsdt.usemecorl': 'CCSD(T) FNO CORRELATION ENERGY',

    'ccsdt.usemecrct': '(T) CORRECTION ENERGY',
    'ccsdtq.usemecrct': '(Q) CORRECTION ENERGY',

    'mp2.usemetrip': 'MP2 SAME-SPIN CORRELATION ENERGY',
    'mp3.usemetrip': 'MP3 SAME-SPIN CORRELATION ENERGY',
    'ccsd.usemetrip': 'CCSD SAME-SPIN CORRELATION ENERGY',

    # <<<  F12  >>>

    'f12.usemeraw': 'HF-CABS TOTAL ENERGY',
    'mp2f12.usemecorl': 'MP2-F12 CORRELATION ENERGY',
    'ccsdaf12.usemecorl': 'CCSD-F12A CORRELATION ENERGY',
    'ccsdbf12.usemecorl': 'CCSD-F12B CORRELATION ENERGY',
    'ccsdcf12.usemecorl': 'CCSD-F12C CORRELATION ENERGY',
    'ccsdnstaf12.usemecorl': 'CCSD(T)-F12A CORRELATION ENERGY',
    'ccsdstaf12.usemecorl': 'CCSD(T*)-F12A CORRELATION ENERGY',
    'ccsdtaf12.usemecorl': 'CCSD(T**)-F12A CORRELATION ENERGY',
    'ccsdnstbf12.usemecorl': 'CCSD(T)-F12B CORRELATION ENERGY',
    'ccsdstbf12.usemecorl': 'CCSD(T*)-F12B CORRELATION ENERGY',
    'ccsdtbf12.usemecorl': 'CCSD(T**)-F12B CORRELATION ENERGY',
    'ccsdnstcf12.usemecorl': 'CCSD(T)-F12C CORRELATION ENERGY',
    'ccsdstcf12.usemecorl': 'CCSD(T*)-F12C CORRELATION ENERGY',
    'ccsdtcf12.usemecorl': 'CCSD(T**)-F12C CORRELATION ENERGY',

    'ccsdnstabf12.usemecrct': '(T)-F12AB CORRECTION ENERGY',
    'ccsdstabf12.usemecrct': '(T*)-F12AB CORRECTION ENERGY',
    'ccsdtabf12.usemecrct': '(T**)-F12AB CORRECTION ENERGY',
    'ccsdnstcf12.usemecrct': '(T)-F12C CORRECTION ENERGY',
    'ccsdstcf12.usemecrct': '(T*)-F12C CORRECTION ENERGY',
    'ccsdtcf12.usemecrct': '(T**)-F12C CORRECTION ENERGY',

    'mp2f12.usemetrip': 'MP2-F12 SAME-SPIN CORRELATION ENERGY',
    'ccsdaf12.usemetrip': 'CCSD-F12A SAME-SPIN CORRELATION ENERGY',
    'ccsdbf12.usemetrip': 'CCSD-F12B SAME-SPIN CORRELATION ENERGY',
    'ccsdcf12.usemetrip': 'CCSD-F12C SAME-SPIN CORRELATION ENERGY',

    # <<<  SAPT  >>>

    'usemesapt': None,
    'usemedftsapt': None,
    'usemempsapt': None,
    #'usemempsapt': 'MP2C DISP20 ENERGY',

    'mp2cDisp20': 'MP2C DISP20 ENERGY',

    'E1pol': 'DFT-SAPT ELST10,R ENERGY',
    'E1exch': 'DFT-SAPT EXCH10 ENERGY',
    'E1exch(S2)': 'DFT-SAPT EXCH10(S^2) ENERGY',  # ne'er used
    'E2ind': 'DFT-SAPT IND20,R ENERGY',
    'E2ind-exch': 'DFT-SAPT EXCH-IND20,R ENERGY',
    'E2disp': 'DFT-SAPT DISP20 ENERGY',
    'E2disp-exch': 'DFT-SAPT EXCH-DISP20 ENERGY',

    'Elst10,r': 'SAPT ELST10,R ENERGY',
    'Elst12,r': 'SAPT ELST12,R ENERGY',
    'Elst13,r': 'SAPT ELST13,R ENERGY',

    'Exch10': 'SAPT EXCH10 ENERGY',
    'Exch10(S^2)': 'SAPT EXCH10(S^2) ENERGY',
    'Exch11(S^2)': 'SAPT EXCH11(S^2) ENERGY',
    'Exch12(S^2)': 'SAPT EXCH12(S^2) ENERGY',

    'Ind20,r': 'SAPT IND20,R ENERGY',
    'Exch-Ind20,r': 'SAPT EXCH-IND20,R ENERGY',
    'Ind22': 'SAPT IND22 ENERGY',
    'Exch-Ind22': 'SAPT EXCH-IND22 ENERGY',
    'Ind30,r': 'SAPT IND30,R ENERGY',
    'Exch-Ind30,r': 'SAPT EXCH-IND30,R ENERGY',
    'Ind-Disp30': 'SAPT IND-DISP30 ENERGY',
    'Exch-Ind-Disp30': 'SAPT EXCH-IND-DISP30 ENERGY',

    'Disp20': 'SAPT DISP20 ENERGY',
    'Exch-Disp20': 'SAPT EXCH-DISP20 ENERGY',
    #'Disp20(OS)': 'SAPT DISP20(OS) ENERGY',
    #'Exch-Disp20(OS)': 'SAPT EXCH-DISP20(OS) ENERGY',
    'Disp20(SS)': 'SAPT SAME-SPIN DISP20 ENERGY',
    'Exch-Disp20(SS)': 'SAPT SAME-SPIN EXCH-DISP20 ENERGY',
    'Disp21': 'SAPT DISP21 ENERGY',
    'Disp22(SDQ)': 'SAPT DISP22(SDQ) ENERGY',  # added for modern parsing, may confuse old usemesapt parsing
    #'Disp22(T)': 'SAPT DISP22(T) ENERGY',  # ditto  # ne'er used
    'Disp22(SDQ).1': 'SAPT DISP22(SDQ) ENERGY',
    #'Disp22(T).1': 'SAPT DISP22(T) ENERGY',  # ne'er used  # edited to remove est
    'Est.Disp22(T)': 'SAPT EST.DISP22(T) ENERGY',
    'Disp2(CCD)': 'SAPT DISP2(CCD) ENERGY',
    'Disp22(S)(CCD)': 'SAPT DISP22(S)(CCD) ENERGY',
    #'Disp22(T)(CCD)': 'SAPT DISP22(T)(CCD) ENERGY',  # ne'er used
    'Est.Disp22(T)(CCD)': 'SAPT EST.DISP22(T)(CCD) ENERGY',
    'Disp30': 'SAPT DISP30 ENERGY',
    'Exch-Disp30': 'SAPT EXCH-DISP30 ENERGY',

    'TotalHF': 'SAPT HF TOTAL ENERGY',
    #'deltaHF,r(2)': None,  # ne'er used
    #'deltaHF,r(3)': None,  # ne'er used
    }

psivar2useme = dict((v, k) for k, v in useme2psivar.items())


optclue2psivar = {
    'full': ['CCSD CORRELATION ENERGY', 'CCSD TOTAL ENERGY',
               'CCSD(T) TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY', '(T) CORRECTION ENERGY',
               'CCSDT TOTAL ENERGY', 'CCSDT CORRELATION ENERGY',
               'CCSDT(Q) TOTAL ENERGY', 'CCSDT(Q) CORRELATION ENERGY', '(Q) CORRECTION ENERGY'],
    'fno1e3': ['CCSD(T) FNO CORRELATION ENERGY', 'CCSD FNO CORRELATION ENERGY', 'MP3 FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY',
               'CCSD CORRELATION ENERGY', 'CCSD TOTAL ENERGY',
               'CCSD(T) TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY', '(T) CORRECTION ENERGY',
               'CCSDT TOTAL ENERGY', 'CCSDT CORRELATION ENERGY',
               'CCSDT(Q) TOTAL ENERGY', 'CCSDT(Q) CORRELATION ENERGY', '(Q) CORRECTION ENERGY'],
    'fno1e4': ['CCSD(T) FNO CORRELATION ENERGY', 'CCSD FNO CORRELATION ENERGY', 'MP3 FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY',
               'CCSD CORRELATION ENERGY', 'CCSD TOTAL ENERGY',
               'CCSD(T) TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY', '(T) CORRECTION ENERGY',
               'CCSDT TOTAL ENERGY', 'CCSDT CORRELATION ENERGY',
               'CCSDT(Q) TOTAL ENERGY', 'CCSDT(Q) CORRELATION ENERGY', '(Q) CORRECTION ENERGY'],
    'fno1e5': ['CCSD(T) FNO CORRELATION ENERGY', 'CCSD FNO CORRELATION ENERGY', 'MP3 FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY',
               'CCSD CORRELATION ENERGY', 'CCSD TOTAL ENERGY',
               'CCSD(T) TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY', '(T) CORRECTION ENERGY',
               'CCSDT TOTAL ENERGY', 'CCSDT CORRELATION ENERGY',
               'CCSDT(Q) TOTAL ENERGY', 'CCSDT(Q) CORRELATION ENERGY', '(Q) CORRECTION ENERGY'],
    'fno5e5': ['CCSD(T) FNO CORRELATION ENERGY', 'CCSD FNO CORRELATION ENERGY', 'MP3 FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY',
               'CCSD CORRELATION ENERGY', 'CCSD TOTAL ENERGY',
               'CCSD(T) TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY', '(T) CORRECTION ENERGY',
               'CCSDT TOTAL ENERGY', 'CCSDT CORRELATION ENERGY',
               'CCSDT(Q) TOTAL ENERGY', 'CCSDT(Q) CORRELATION ENERGY', '(Q) CORRECTION ENERGY'],
    'fno1e6': ['CCSD(T) FNO CORRELATION ENERGY', 'CCSD FNO CORRELATION ENERGY', 'MP3 FNO CORRELATION ENERGY', 'FNO CORRECTION ENERGY',
               'CCSD CORRELATION ENERGY', 'CCSD TOTAL ENERGY',
               'CCSD(T) TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY', '(T) CORRECTION ENERGY',
               'CCSDT TOTAL ENERGY', 'CCSDT CORRELATION ENERGY',
               'CCSDT(Q) TOTAL ENERGY', 'CCSDT(Q) CORRELATION ENERGY', '(Q) CORRECTION ENERGY'],
    'dsrgs0p1':['MP2 CORRELATION ENERGY', 'MP2 TOTAL ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY'],
    'dsrgs0p5':['MP2 CORRELATION ENERGY', 'MP2 TOTAL ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY'],
    'dsrgs1p0':['MP2 CORRELATION ENERGY', 'MP2 TOTAL ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY'],
    'mrcc': ['CCSD CORRELATION ENERGY',
               'CCSD(T) CORRELATION ENERGY', '(T) CORRECTION ENERGY',
               'CCSDT CORRELATION ENERGY',
               'CCSDT(Q) CORRELATION ENERGY', '(Q) CORRECTION ENERGY'],
    'nfc': ['B2PLYP TOTAL ENERGY', 'B2PLYP-D2 TOTAL ENERGY', 'B2PLYP-D3 TOTAL ENERGY', 'B2PLYP-D3(BJ) TOTAL ENERGY',
                'B2PLYP-D3M TOTAL ENERGY', 'B2PLYP-D3M(BJ) TOTAL ENERGY',
            'DSD-PBEP86 TOTAL ENERGY', 'DSD-PBEP86-D2 TOTAL ENERGY', 'DSD-PBEP86-D3 TOTAL ENERGY', 'DSD-PBEP86-D3(BJ) TOTAL ENERGY',
            'WB97X-2 TOTAL ENERGY'],
    'fc': ['B2PLYP TOTAL ENERGY', 'B2PLYP-D2 TOTAL ENERGY', 'B2PLYP-D3 TOTAL ENERGY', 'B2PLYP-D3(BJ) TOTAL ENERGY',
                'B2PLYP-D3M TOTAL ENERGY', 'B2PLYP-D3M(BJ) TOTAL ENERGY',
            'DSD-PBEP86 TOTAL ENERGY', 'DSD-PBEP86-D2 TOTAL ENERGY', 'DSD-PBEP86-D3 TOTAL ENERGY', 'DSD-PBEP86-D3(BJ) TOTAL ENERGY',
            'WB97X-2 TOTAL ENERGY'],
    'dfhf': ['HF-CABS TOTAL ENERGY', 'MP2-F12 TOTAL ENERGY', 'SCS-MP2-F12 TOTAL ENERGY', 'SCS(N)-MP2-F12 TOTAL ENERGY',
             'SCS(MI)-MP2-F12 TOTAL ENERGY', 'DW-MP2-F12 TOTAL ENERGY', 'MP2C-F12 TOTAL ENERGY',
             'SCF TOTAL ENERGY', 'HF TOTAL ENERGY', 'MP2 TOTAL ENERGY', 'SCS-MP2 TOTAL ENERGY', 'SCS(N)-MP2 TOTAL ENERGY',
             'SCS(MI)-MP2 TOTAL ENERGY', 'DW-MP2 TOTAL ENERGY', 'MP2C TOTAL ENERGY',
             'B3LYP FUNCTIONAL TOTAL ENERGY', 'B3LYP TOTAL ENERGY', 'B3LYP-D2 TOTAL ENERGY', 'B3LYP-D3 TOTAL ENERGY', 'B3LYP-D3(BJ) TOTAL ENERGY', 'B3LYP-XDM TOTAL ENERGY',
             'BLYP FUNCTIONAL TOTAL ENERGY', 'BLYP TOTAL ENERGY', 'BLYP-D2 TOTAL ENERGY', 'BLYP-D3 TOTAL ENERGY', 'BLYP-D3(BJ) TOTAL ENERGY',
             'BP86 FUNCTIONAL TOTAL ENERGY', 'BP86 TOTAL ENERGY', 'BP86-D2 TOTAL ENERGY', 'BP86-D3 TOTAL ENERGY', 'BP86-D3(BJ) TOTAL ENERGY',
             'PBE FUNCTIONAL TOTAL ENERGY', 'PBE TOTAL ENERGY', 'PBE-D2 TOTAL ENERGY', 'PBE-D3 TOTAL ENERGY', 'PBE-D3(BJ) TOTAL ENERGY',
             'PBE0 FUNCTIONAL TOTAL ENERGY', 'PBE0 TOTAL ENERGY', 'PBE0-D2 TOTAL ENERGY', 'PBE0-D3 TOTAL ENERGY', 'PBE0-D3(BJ) TOTAL ENERGY',
             'B97 FUNCTIONAL TOTAL ENERGY', 'B97 TOTAL ENERGY', 'B97-D2 TOTAL ENERGY', 'B97-D3 TOTAL ENERGY', 'B97-D3(BJ) TOTAL ENERGY',
             'B2PLYP TOTAL ENERGY', 'B2PLYP-D2 TOTAL ENERGY', 'B2PLYP-D3 TOTAL ENERGY', 'B2PLYP-D3(BJ) TOTAL ENERGY',
             'WPBE FUNCTIONAL TOTAL ENERGY', 'WPBE TOTAL ENERGY', 'WPBE-D3 TOTAL ENERGY', 'WPBE-D3(BJ) TOTAL ENERGY',
             'M05-2X FUNCTIONAL TOTAL ENERGY', 'M05-2X TOTAL ENERGY',
             'WB97X FUNCTIONAL TOTAL ENERGY', 'WB97X-D TOTAL ENERGY',
             'B3LYP-D3M TOTAL ENERGY', 'BLYP-D3M TOTAL ENERGY', 'BP86-D3M TOTAL ENERGY', 'PBE-D3M TOTAL ENERGY',
                'PBE0-D3M TOTAL ENERGY', 'B97-D3M TOTAL ENERGY', 'B2PLYP-D3M TOTAL ENERGY', 'WPBE-D3M TOTAL ENERGY',
             'B3LYP-D3M(BJ) TOTAL ENERGY', 'BLYP-D3M(BJ) TOTAL ENERGY', 'BP86-D3M(BJ) TOTAL ENERGY', 'PBE-D3M(BJ) TOTAL ENERGY',
                'PBE0-D3M(BJ) TOTAL ENERGY', 'B97-D3M(BJ) TOTAL ENERGY', 'B2PLYP-D3M(BJ) TOTAL ENERGY', 'WPBE-D3M(BJ) TOTAL ENERGY',
            ],
    'dfmp': ['MP2-F12 CORRELATION ENERGY', 'MP2-F12 TOTAL ENERGY', 'MP2-F12 SAME-SPIN CORRELATION ENERGY',
             'SCS-MP2-F12 CORRELATION ENERGY', 'SCS-MP2-F12 TOTAL ENERGY',
             'SCS(N)-MP2-F12 CORRELATION ENERGY', 'SCS(N)-MP2-F12 TOTAL ENERGY',
             'SCS(MI)-MP2-F12 CORRELATION ENERGY', 'SCS(MI)-MP2-F12 TOTAL ENERGY',
             'DW-MP2-F12 CORRELATION ENERGY', 'DW-MP2-F12 TOTAL ENERGY',
             'MP2C-F12 CORRELATION ENERGY', 'MP2C-F12 TOTAL ENERGY',
             'MP2 CORRELATION ENERGY', 'MP2 TOTAL ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY',
             'SCS-MP2 CORRELATION ENERGY', 'SCS-MP2 TOTAL ENERGY',
             'SCS(N)-MP2 CORRELATION ENERGY', 'SCS(N)-MP2 TOTAL ENERGY',
             'SCS(MI)-MP2 CORRELATION ENERGY', 'SCS(MI)-MP2 TOTAL ENERGY',
             'DW-MP2 CORRELATION ENERGY', 'DW-MP2 TOTAL ENERGY',
             'MP2C CORRELATION ENERGY', 'MP2C TOTAL ENERGY',
             'SAPT2+DMP2 TOTAL ENERGY', 'SAPT2+(CCD)DMP2 TOTAL ENERGY',
             'SAPT2+(3)DMP2 TOTAL ENERGY', 'SAPT2+(3)(CCD)DMP2 TOTAL ENERGY',
             'SAPT2+3DMP2 TOTAL ENERGY', 'SAPT2+3(CCD)DMP2 TOTAL ENERGY',
             'B2PLYP TOTAL ENERGY', 'B2PLYP-D2 TOTAL ENERGY', 'B2PLYP-D3 TOTAL ENERGY', 'B2PLYP-D3(BJ) TOTAL ENERGY',
                'B2PLYP-D3M TOTAL ENERGY', 'B2PLYP-D3M(BJ) TOTAL ENERGY',
            ],
    }
