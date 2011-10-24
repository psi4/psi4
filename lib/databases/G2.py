import re
import input

# <<< G2 Database Module >>>
# The G2/97 database from Curtiss et al. in JCP 106 1063 (1997) and JCP 109 42 (1998).
# Geometries from MP2 links at http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/G2-97.htm . Specifically, from:
#    http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g2-97_cart_neut.txt
#    http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g2-97small.txt
#    http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g2-97anion.txt
#    http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g2-97aux.txt
#    http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/g2-97extra.txt
# Reference energy sources detailed at BIND dictionary below.
dbse = 'G2'
isOS = 'true'

# <<< Database Members >>>
HRXN = range(1, 303)
HRXN_SM = []
HRXN_LG = []
G21IE = range(1, 39)     # G2-1 Ionization Energies (38)
G21EA = range(39, 64)    # G2-1 Electron Affinities (25)
G21PA = range(64, 72)    # G2-1 Proton Affinities (8)
G21EF = range(72, 127)   # G2-1 Enthalpies of Formation (55)
G22IE = range(127, 177)  # G2-2 Ionization Energies (50)
G22EA = range(177, 210)  # G2-2 Electron Affinities (33)
G22EF = range(210, 303)  # G2-2 Enthalpies of Formation (93)

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
# [  1 -  38] G2-1 Ionization Energies (38)
ACTV['%s-%s'            % (dbse, '1'                     )] = ['%s-%s-reagent'      % (dbse, 'Li_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Li_cation') ]
RXNM['%s-%s'            % (dbse, '1'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '1')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '2'                     )] = ['%s-%s-reagent'      % (dbse, 'Be'),
                                                               '%s-%s-reagent'      % (dbse, 'Be_cation') ]
RXNM['%s-%s'            % (dbse, '2'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '2')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '3'                     )] = ['%s-%s-reagent'      % (dbse, 'B_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'B_cation') ]
RXNM['%s-%s'            % (dbse, '3'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '3')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '4'                     )] = ['%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'C_cation') ]
RXNM['%s-%s'            % (dbse, '4'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '4')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '5'                     )] = ['%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'N_cation') ]
RXNM['%s-%s'            % (dbse, '5'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '5')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '6'                     )] = ['%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'O_cation') ]
RXNM['%s-%s'            % (dbse, '6'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '6')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '7'                     )] = ['%s-%s-reagent'      % (dbse, 'F_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'F_cation') ]
RXNM['%s-%s'            % (dbse, '7'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '7')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '8'                     )] = ['%s-%s-reagent'      % (dbse, 'Na_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Na_cation') ]
RXNM['%s-%s'            % (dbse, '8'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '8')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '9'                     )] = ['%s-%s-reagent'      % (dbse, 'Mg'),
                                                               '%s-%s-reagent'      % (dbse, 'Mg_cation') ]
RXNM['%s-%s'            % (dbse, '9'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '9')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '10'                    )] = ['%s-%s-reagent'      % (dbse, 'Al_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Al_cation') ]
RXNM['%s-%s'            % (dbse, '10'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '10')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '11'                    )] = ['%s-%s-reagent'      % (dbse, 'Si'),
                                                               '%s-%s-reagent'      % (dbse, 'Si_cation') ]
RXNM['%s-%s'            % (dbse, '11'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '11')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '12'                    )] = ['%s-%s-reagent'      % (dbse, 'P_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'P_cation') ]
RXNM['%s-%s'            % (dbse, '12'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '12')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '13'                    )] = ['%s-%s-reagent'      % (dbse, 'S'),
                                                               '%s-%s-reagent'      % (dbse, 'S_cation') ]
RXNM['%s-%s'            % (dbse, '13'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '13')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '14'                    )] = ['%s-%s-reagent'      % (dbse, 'Cl_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_cation') ]
RXNM['%s-%s'            % (dbse, '14'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '14')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '15'                    )] = ['%s-%s-reagent'      % (dbse, 'Methane'),
                                                               '%s-%s-reagent'      % (dbse, 'Methane_cation') ]
RXNM['%s-%s'            % (dbse, '15'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '15')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '16'                    )] = ['%s-%s-reagent'      % (dbse, 'NH3'),
                                                               '%s-%s-reagent'      % (dbse, 'NH3_cation') ]
RXNM['%s-%s'            % (dbse, '16'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '16')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '17'                    )] = ['%s-%s-reagent'      % (dbse, 'OH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'OH_cation') ]
RXNM['%s-%s'            % (dbse, '17'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '17')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '18'                    )] = ['%s-%s-reagent'      % (dbse, 'H2O'),
                                                               '%s-%s-reagent'      % (dbse, 'H2O_cation') ]
RXNM['%s-%s'            % (dbse, '18'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '18')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '19'                    )] = ['%s-%s-reagent'      % (dbse, 'HF'),
                                                               '%s-%s-reagent'      % (dbse, 'HF_cation') ]
RXNM['%s-%s'            % (dbse, '19'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '19')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '20'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH4'),
                                                               '%s-%s-reagent'      % (dbse, 'SiH4_cation') ]
RXNM['%s-%s'            % (dbse, '20'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '20')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '21'                    )] = ['%s-%s-reagent'      % (dbse, 'PH'),
                                                               '%s-%s-reagent'      % (dbse, 'PH_cation') ]
RXNM['%s-%s'            % (dbse, '21'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '21')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '22'                    )] = ['%s-%s-reagent'      % (dbse, 'PH2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'PH2_cation') ]
RXNM['%s-%s'            % (dbse, '22'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '22')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '23'                    )] = ['%s-%s-reagent'      % (dbse, 'PH3'),
                                                               '%s-%s-reagent'      % (dbse, 'PH3_cation') ]
RXNM['%s-%s'            % (dbse, '23'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '23')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '24'                    )] = ['%s-%s-reagent'      % (dbse, 'SH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'SH_cation') ]
RXNM['%s-%s'            % (dbse, '24'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '24')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '25'                    )] = ['%s-%s-reagent'      % (dbse, 'SH2'),
                                                               '%s-%s-reagent'      % (dbse, 'SH2_cation_2B1') ]
RXNM['%s-%s'            % (dbse, '25'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '25')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '26'                    )] = ['%s-%s-reagent'      % (dbse, 'SH2'),
                                                               '%s-%s-reagent'      % (dbse, 'SH2_cation_2A1') ]
RXNM['%s-%s'            % (dbse, '26'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '26')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '27'                    )] = ['%s-%s-reagent'      % (dbse, 'HCl'),
                                                               '%s-%s-reagent'      % (dbse, 'HCl_cation') ]
RXNM['%s-%s'            % (dbse, '27'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '27')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '28'                    )] = ['%s-%s-reagent'      % (dbse, 'Ethyne'),
                                                               '%s-%s-reagent'      % (dbse, 'Ethyne_cation') ]
RXNM['%s-%s'            % (dbse, '28'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '28')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '29'                    )] = ['%s-%s-reagent'      % (dbse, 'Ethene'),
                                                               '%s-%s-reagent'      % (dbse, 'Ethene_cation') ]
RXNM['%s-%s'            % (dbse, '29'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '29')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '30'                    )] = ['%s-%s-reagent'      % (dbse, 'CO'),
                                                               '%s-%s-reagent'      % (dbse, 'CO_cation') ]
RXNM['%s-%s'            % (dbse, '30'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '30')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '31'                    )] = ['%s-%s-reagent'      % (dbse, 'N2'),
                                                               '%s-%s-reagent'      % (dbse, 'N2_cation_2SIGMAg') ]
RXNM['%s-%s'            % (dbse, '31'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '31')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '32'                    )] = ['%s-%s-reagent'      % (dbse, 'N2'),
                                                               '%s-%s-reagent'      % (dbse, 'N2_cation_2PIu') ]
RXNM['%s-%s'            % (dbse, '32'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '32')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '33'                    )] = ['%s-%s-reagent'      % (dbse, 'O2'),
                                                               '%s-%s-reagent'      % (dbse, 'O2_cation') ]
RXNM['%s-%s'            % (dbse, '33'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '33')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '34'                    )] = ['%s-%s-reagent'      % (dbse, 'P2'),
                                                               '%s-%s-reagent'      % (dbse, 'P2_cation') ]
RXNM['%s-%s'            % (dbse, '34'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '34')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '35'                    )] = ['%s-%s-reagent'      % (dbse, 'S2'),
                                                               '%s-%s-reagent'      % (dbse, 'S2_cation') ]
RXNM['%s-%s'            % (dbse, '35'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '35')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '36'                    )] = ['%s-%s-reagent'      % (dbse, 'Cl2'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl2_cation') ]
RXNM['%s-%s'            % (dbse, '36'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '36')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '37'                    )] = ['%s-%s-reagent'      % (dbse, 'ClF'),
                                                               '%s-%s-reagent'      % (dbse, 'ClF_cation') ]
RXNM['%s-%s'            % (dbse, '37'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '37')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '38'                    )] = ['%s-%s-reagent'      % (dbse, 'CS'),
                                                               '%s-%s-reagent'      % (dbse, 'CS_cation') ]
RXNM['%s-%s'            % (dbse, '38'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '38')], [-1, +1]))

# [ 39 -  63] G2-1 Electron Affinities (25)
ACTV['%s-%s'            % (dbse, '39'                    )] = ['%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'C_anion') ]
RXNM['%s-%s'            % (dbse, '39'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '39')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '40'                    )] = ['%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'O_anion') ]
RXNM['%s-%s'            % (dbse, '40'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '40')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '41'                    )] = ['%s-%s-reagent'      % (dbse, 'F_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'F_anion') ]
RXNM['%s-%s'            % (dbse, '41'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '41')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '42'                    )] = ['%s-%s-reagent'      % (dbse, 'Si'),
                                                               '%s-%s-reagent'      % (dbse, 'Si_anion') ]
RXNM['%s-%s'            % (dbse, '42'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '42')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '43'                    )] = ['%s-%s-reagent'      % (dbse, 'P_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'P_anion') ]
RXNM['%s-%s'            % (dbse, '43'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '43')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '44'                    )] = ['%s-%s-reagent'      % (dbse, 'S'),
                                                               '%s-%s-reagent'      % (dbse, 'S_anion') ]
RXNM['%s-%s'            % (dbse, '44'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '44')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '45'                    )] = ['%s-%s-reagent'      % (dbse, 'Cl_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_anion') ]
RXNM['%s-%s'            % (dbse, '45'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '45')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '46'                    )] = ['%s-%s-reagent'      % (dbse, 'CH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CH_anion') ]
RXNM['%s-%s'            % (dbse, '46'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '46')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '47'                    )] = ['%s-%s-reagent'      % (dbse, 'H2C_singlet'),  # ACK 1/3 nc
                                                               '%s-%s-reagent'      % (dbse, 'H2C_anion') ]
RXNM['%s-%s'            % (dbse, '47'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '47')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '48'                    )] = ['%s-%s-reagent'      % (dbse, 'Methyl_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Methyl_anion') ]
RXNM['%s-%s'            % (dbse, '48'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '48')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '49'                    )] = ['%s-%s-reagent'      % (dbse, 'NH'),
                                                               '%s-%s-reagent'      % (dbse, 'NH_anion') ]
RXNM['%s-%s'            % (dbse, '49'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '49')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '50'                    )] = ['%s-%s-reagent'      % (dbse, 'NH2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'NH2_anion') ]
RXNM['%s-%s'            % (dbse, '50'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '50')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '51'                    )] = ['%s-%s-reagent'      % (dbse, 'OH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'OH_anion') ]
RXNM['%s-%s'            % (dbse, '51'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '51')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '52'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'SiH_anion') ]
RXNM['%s-%s'            % (dbse, '52'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '52')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '53'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH2_singlet'),  # ACK 1/3 nc
                                                               '%s-%s-reagent'      % (dbse, 'SiH2_anion') ]
RXNM['%s-%s'            % (dbse, '53'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '53')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '54'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH3_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'SiH3_anion') ]
RXNM['%s-%s'            % (dbse, '54'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '54')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '55'                    )] = ['%s-%s-reagent'      % (dbse, 'PH'),
                                                               '%s-%s-reagent'      % (dbse, 'PH_anion') ]
RXNM['%s-%s'            % (dbse, '55'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '55')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '56'                    )] = ['%s-%s-reagent'      % (dbse, 'PH2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'PH2_anion') ]
RXNM['%s-%s'            % (dbse, '56'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '56')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '57'                    )] = ['%s-%s-reagent'      % (dbse, 'SH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'SH_anion') ]
RXNM['%s-%s'            % (dbse, '57'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '57')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '58'                    )] = ['%s-%s-reagent'      % (dbse, 'O2'),
                                                               '%s-%s-reagent'      % (dbse, 'O2_anion') ]
RXNM['%s-%s'            % (dbse, '58'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '58')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '59'                    )] = ['%s-%s-reagent'      % (dbse, 'NO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'NO_anion') ]
RXNM['%s-%s'            % (dbse, '59'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '59')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '60'                    )] = ['%s-%s-reagent'      % (dbse, 'CN_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CN_anion') ]
RXNM['%s-%s'            % (dbse, '60'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '60')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '61'                    )] = ['%s-%s-reagent'      % (dbse, 'PO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'PO_anion') ]
RXNM['%s-%s'            % (dbse, '61'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '61')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '62'                    )] = ['%s-%s-reagent'      % (dbse, 'S2'),
                                                               '%s-%s-reagent'      % (dbse, 'S2_anion') ]
RXNM['%s-%s'            % (dbse, '62'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '62')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '63'                    )] = ['%s-%s-reagent'      % (dbse, 'Cl2'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl2_anion') ]
RXNM['%s-%s'            % (dbse, '63'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '63')], [+1, -1]))

# [ 64 -  71] G2-1 Proton Affinities (8)
ACTV['%s-%s'            % (dbse, '64'                    )] = ['%s-%s-reagent'      % (dbse, 'NH3'),
                                                               '%s-%s-reagent'      % (dbse, 'NH4_cation') ]
RXNM['%s-%s'            % (dbse, '64'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '64')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '65'                    )] = ['%s-%s-reagent'      % (dbse, 'H2O'),
                                                               '%s-%s-reagent'      % (dbse, 'H3O_cation') ]
RXNM['%s-%s'            % (dbse, '65'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '65')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '66'                    )] = ['%s-%s-reagent'      % (dbse, 'Ethyne'),
                                                               '%s-%s-reagent'      % (dbse, 'Vinyl_ncl_cation') ]
RXNM['%s-%s'            % (dbse, '66'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '66')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '67'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH4'),
                                                               '%s-%s-reagent'      % (dbse, 'SiH5_cation') ]
RXNM['%s-%s'            % (dbse, '67'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '67')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '68'                    )] = ['%s-%s-reagent'      % (dbse, 'PH3'),
                                                               '%s-%s-reagent'      % (dbse, 'PH4_cation') ]
RXNM['%s-%s'            % (dbse, '68'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '68')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '69'                    )] = ['%s-%s-reagent'      % (dbse, 'SH2'),
                                                               '%s-%s-reagent'      % (dbse, 'SH3_cation') ]
RXNM['%s-%s'            % (dbse, '69'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '69')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '70'                    )] = ['%s-%s-reagent'      % (dbse, 'HCl'),
                                                               '%s-%s-reagent'      % (dbse, 'H2Cl_cation') ]
RXNM['%s-%s'            % (dbse, '70'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '70')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '71'                    )] = ['%s-%s-reagent'      % (dbse, 'H2'),
                                                               '%s-%s-reagent'      % (dbse, 'H3_cation') ]
RXNM['%s-%s'            % (dbse, '71'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '71')], [+1, -1]))

# [ 72 - 126] G2-1 Enthalpies of Formation (55)
ACTV['%s-%s'            % (dbse, '72'                    )] = ['%s-%s-reagent'      % (dbse, 'LiH'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Li_radical') ]
RXNM['%s-%s'            % (dbse, '72'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '72')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '73'                    )] = ['%s-%s-reagent'      % (dbse, 'BeH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Be') ]
RXNM['%s-%s'            % (dbse, '73'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '73')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '74'                    )] = ['%s-%s-reagent'      % (dbse, 'CH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '74'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '74')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '75'                    )] = ['%s-%s-reagent'      % (dbse, 'H2C_triplet'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '75'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '75')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '76'                    )] = ['%s-%s-reagent'      % (dbse, 'H2C_singlet'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '76'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '76')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '77'                    )] = ['%s-%s-reagent'      % (dbse, 'Methyl_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '77'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '77')], [-1, +3, +1]))

ACTV['%s-%s'            % (dbse, '78'                    )] = ['%s-%s-reagent'      % (dbse, 'Methane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '78'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '78')], [-1, +4, +1]))

ACTV['%s-%s'            % (dbse, '79'                    )] = ['%s-%s-reagent'      % (dbse, 'NH'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '79'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '79')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '80'                    )] = ['%s-%s-reagent'      % (dbse, 'NH2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '80'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '80')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '81'                    )] = ['%s-%s-reagent'      % (dbse, 'NH3'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '81'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '81')], [-1, +3, +1]))

ACTV['%s-%s'            % (dbse, '82'                    )] = ['%s-%s-reagent'      % (dbse, 'OH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '82'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '82')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '83'                    )] = ['%s-%s-reagent'      % (dbse, 'H2O'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '83'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '83')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '84'                    )] = ['%s-%s-reagent'      % (dbse, 'HF'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '84'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '84')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '85'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH2_singlet'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '85'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '85')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '86'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH2_triplet'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '86'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '86')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '87'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH3_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '87'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '87')], [-1, +3, +1]))

ACTV['%s-%s'            % (dbse, '88'                    )] = ['%s-%s-reagent'      % (dbse, 'SiH4'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '88'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '88')], [-1, +4, +1]))

ACTV['%s-%s'            % (dbse, '89'                    )] = ['%s-%s-reagent'      % (dbse, 'PH2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'P_radical') ]
RXNM['%s-%s'            % (dbse, '89'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '89')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '90'                    )] = ['%s-%s-reagent'      % (dbse, 'PH3'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'P_radical') ]
RXNM['%s-%s'            % (dbse, '90'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '90')], [-1, +3, +1]))

ACTV['%s-%s'            % (dbse, '91'                    )] = ['%s-%s-reagent'      % (dbse, 'SH2'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '91'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '91')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '92'                    )] = ['%s-%s-reagent'      % (dbse, 'HCl'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '92'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '92')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '93'                    )] = ['%s-%s-reagent'      % (dbse, 'Li2'),
                                                               '%s-%s-reagent'      % (dbse, 'Li_radical') ]
RXNM['%s-%s'            % (dbse, '93'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '93')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '94'                    )] = ['%s-%s-reagent'      % (dbse, 'LiF'),
                                                               '%s-%s-reagent'      % (dbse, 'Li_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '94'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '94')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '95'                    )] = ['%s-%s-reagent'      % (dbse, 'Ethyne'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '95'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '95')], [-1, +2, +2]))

ACTV['%s-%s'            % (dbse, '96'                    )] = ['%s-%s-reagent'      % (dbse, 'Ethene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '96'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '96')], [-1, +4, +2]))

ACTV['%s-%s'            % (dbse, '97'                    )] = ['%s-%s-reagent'      % (dbse, 'Ethane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '97'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '97')], [-1, +6, +2]))

ACTV['%s-%s'            % (dbse, '98'                    )] = ['%s-%s-reagent'      % (dbse, 'CN_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '98'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '98')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '99'                    )] = ['%s-%s-reagent'      % (dbse, 'HCN'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '99'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '99')], [-1, +1, +1, +1]))

ACTV['%s-%s'            % (dbse, '100'                   )] = ['%s-%s-reagent'      % (dbse, 'CO'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '100'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '100')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '101'                   )] = ['%s-%s-reagent'      % (dbse, 'HCO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '101'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '101')], [-1, +1, +1, +1]))

ACTV['%s-%s'            % (dbse, '102'                   )] = ['%s-%s-reagent'      % (dbse, 'Formaldehyde'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '102'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '102')], [-1, +2, +1, +1]))

ACTV['%s-%s'            % (dbse, '103'                   )] = ['%s-%s-reagent'      % (dbse, 'MeOH'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '103'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '103')], [-1, +4, +1, +1]))

ACTV['%s-%s'            % (dbse, '104'                   )] = ['%s-%s-reagent'      % (dbse, 'N2'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '104'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '104')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '105'                   )] = ['%s-%s-reagent'      % (dbse, 'Hydrazine'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '105'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '105')], [-1, +4, +2]))

ACTV['%s-%s'            % (dbse, '106'                   )] = ['%s-%s-reagent'      % (dbse, 'NO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '106'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '106')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '107'                   )] = ['%s-%s-reagent'      % (dbse, 'O2'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '107'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '107')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '108'                   )] = ['%s-%s-reagent'      % (dbse, 'HOOH'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '108'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '108')], [-1, +2, +2]))

ACTV['%s-%s'            % (dbse, '109'                   )] = ['%s-%s-reagent'      % (dbse, 'F2'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '109'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '109')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '110'                   )] = ['%s-%s-reagent'      % (dbse, 'CO2'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '110'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '110')], [-1, +1, +2]))

ACTV['%s-%s'            % (dbse, '111'                   )] = ['%s-%s-reagent'      % (dbse, 'Na2'),
                                                               '%s-%s-reagent'      % (dbse, 'Na_radical') ]
RXNM['%s-%s'            % (dbse, '111'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '111')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '112'                   )] = ['%s-%s-reagent'      % (dbse, 'Si2'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '112'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '112')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '113'                   )] = ['%s-%s-reagent'      % (dbse, 'P2'),
                                                               '%s-%s-reagent'      % (dbse, 'P_radical') ]
RXNM['%s-%s'            % (dbse, '113'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '113')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '114'                   )] = ['%s-%s-reagent'      % (dbse, 'S2'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '114'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '114')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '115'                   )] = ['%s-%s-reagent'      % (dbse, 'Cl2'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '115'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '115')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '116'                   )] = ['%s-%s-reagent'      % (dbse, 'NaCl'),
                                                               '%s-%s-reagent'      % (dbse, 'Na_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '116'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '116')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '117'                   )] = ['%s-%s-reagent'      % (dbse, 'SiO'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '117'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '117')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '118'                   )] = ['%s-%s-reagent'      % (dbse, 'CS'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '118'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '118')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '119'                   )] = ['%s-%s-reagent'      % (dbse, 'SO'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '119'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '119')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '120'                   )] = ['%s-%s-reagent'      % (dbse, 'ClO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '120'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '120')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '121'                   )] = ['%s-%s-reagent'      % (dbse, 'ClF'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '121'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '121')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '122'                   )] = ['%s-%s-reagent'      % (dbse, 'Disilane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '122'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '122')], [-1, +6, +2]))

ACTV['%s-%s'            % (dbse, '123'                   )] = ['%s-%s-reagent'      % (dbse, 'MeCl'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '123'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '123')], [-1, +3, +1, +1]))

ACTV['%s-%s'            % (dbse, '124'                   )] = ['%s-%s-reagent'      % (dbse, 'MeSH'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '124'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '124')], [-1, +3, +1, +1]))

ACTV['%s-%s'            % (dbse, '125'                   )] = ['%s-%s-reagent'      % (dbse, 'HOCl'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '125'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '125')], [-1, +1, +1, +1]))

ACTV['%s-%s'            % (dbse, '126'                   )] = ['%s-%s-reagent'      % (dbse, 'SO2'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '126'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '126')], [-1, +2, +1]))

# [127 - 176] G2-2 Ionization Energies (50)
ACTV['%s-%s'            % (dbse, '127'                   )] = ['%s-%s-reagent'      % (dbse, 'H_radical') ]
RXNM['%s-%s'            % (dbse, '127'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '127')], [-1]))

ACTV['%s-%s'            % (dbse, '128'                   )] = ['%s-%s-reagent'      % (dbse, 'He'),
                                                               '%s-%s-reagent'      % (dbse, 'He_cation') ]
RXNM['%s-%s'            % (dbse, '128'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '128')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '129'                   )] = ['%s-%s-reagent'      % (dbse, 'Ne'),
                                                               '%s-%s-reagent'      % (dbse, 'Ne_cation') ]
RXNM['%s-%s'            % (dbse, '129'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '129')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '130'                   )] = ['%s-%s-reagent'      % (dbse, 'Ar'),
                                                               '%s-%s-reagent'      % (dbse, 'Ar_cation') ]
RXNM['%s-%s'            % (dbse, '130'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '130')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '131'                   )] = ['%s-%s-reagent'      % (dbse, 'BF3'),
                                                               '%s-%s-reagent'      % (dbse, 'BF3_cation') ]
RXNM['%s-%s'            % (dbse, '131'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '131')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '132'                   )] = ['%s-%s-reagent'      % (dbse, 'BCl3'),
                                                               '%s-%s-reagent'      % (dbse, 'BCl3_cation') ]
RXNM['%s-%s'            % (dbse, '132'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '132')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '133'                   )] = ['%s-%s-reagent'      % (dbse, 'B2F4'),
                                                               '%s-%s-reagent'      % (dbse, 'B2F4_cation') ]
RXNM['%s-%s'            % (dbse, '133'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '133')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '134'                   )] = ['%s-%s-reagent'      % (dbse, 'CO2'),
                                                               '%s-%s-reagent'      % (dbse, 'CO2_cation') ]
RXNM['%s-%s'            % (dbse, '134'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '134')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '135'                   )] = ['%s-%s-reagent'      % (dbse, 'CF2'),
                                                               '%s-%s-reagent'      % (dbse, 'CF2_cation') ]
RXNM['%s-%s'            % (dbse, '135'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '135')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '136'                   )] = ['%s-%s-reagent'      % (dbse, 'OCS'),
                                                               '%s-%s-reagent'      % (dbse, 'OCS_cation') ]
RXNM['%s-%s'            % (dbse, '136'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '136')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '137'                   )] = ['%s-%s-reagent'      % (dbse, 'CS2'),
                                                               '%s-%s-reagent'      % (dbse, 'CS2_cation') ]
RXNM['%s-%s'            % (dbse, '137'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '137')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '138'                   )] = ['%s-%s-reagent'      % (dbse, 'H2C_triplet'),  # ACK- form not fully specified, now think 3 bc p.43
                                                               '%s-%s-reagent'      % (dbse, 'H2C_cation') ]
RXNM['%s-%s'            % (dbse, '138'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '138')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '139'                   )] = ['%s-%s-reagent'      % (dbse, 'Methyl_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Methyl_cation') ]
RXNM['%s-%s'            % (dbse, '139'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '139')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '140'                   )] = ['%s-%s-reagent'      % (dbse, 'C2H5_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C2H5_ncl_cation') ]
RXNM['%s-%s'            % (dbse, '140'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '140')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '141'                   )] = ['%s-%s-reagent'      % (dbse, 'Cyclopropene'),
                                                               '%s-%s-reagent'      % (dbse, 'Cyclopropenyl_cation') ]
RXNM['%s-%s'            % (dbse, '141'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '141')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '142'                   )] = ['%s-%s-reagent'      % (dbse, 'Allene'),
                                                               '%s-%s-reagent'      % (dbse, 'Allene_cation') ]
RXNM['%s-%s'            % (dbse, '142'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '142')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '143'                   )] = ['%s-%s-reagent'      % (dbse, 'Me2CH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Me2CH_cation') ]
RXNM['%s-%s'            % (dbse, '143'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '143')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '144'                   )] = ['%s-%s-reagent'      % (dbse, 'Benzene'),
                                                               '%s-%s-reagent'      % (dbse, 'Benzene_cation') ]
RXNM['%s-%s'            % (dbse, '144'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '144')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '145'                   )] = ['%s-%s-reagent'      % (dbse, 'Toluene'),
                                                               '%s-%s-reagent'      % (dbse, 'Toluene_cation') ]
RXNM['%s-%s'            % (dbse, '145'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '145')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '146'                   )] = ['%s-%s-reagent'      % (dbse, 'CN_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CN_cation') ]
RXNM['%s-%s'            % (dbse, '146'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '146')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '147'                   )] = ['%s-%s-reagent'      % (dbse, 'HCO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'HCO_cation') ]
RXNM['%s-%s'            % (dbse, '147'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '147')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '148'                   )] = ['%s-%s-reagent'      % (dbse, 'H2COH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H2COH_cation') ]
RXNM['%s-%s'            % (dbse, '148'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '148')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '149'                   )] = ['%s-%s-reagent'      % (dbse, 'CH3O_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CH3O_cation') ]
RXNM['%s-%s'            % (dbse, '149'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '149')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '150'                   )] = ['%s-%s-reagent'      % (dbse, 'MeOH'),
                                                               '%s-%s-reagent'      % (dbse, 'MeOH_cation') ]
RXNM['%s-%s'            % (dbse, '150'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '150')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '151'                   )] = ['%s-%s-reagent'      % (dbse, 'MeF'),
                                                               '%s-%s-reagent'      % (dbse, 'MeF_cation') ]
RXNM['%s-%s'            % (dbse, '151'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '151')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '152'                   )] = ['%s-%s-reagent'      % (dbse, 'Thioformaldehyde'),
                                                               '%s-%s-reagent'      % (dbse, 'Thioformaldehyde_cation') ]
RXNM['%s-%s'            % (dbse, '152'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '152')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '153'                   )] = ['%s-%s-reagent'      % (dbse, 'H2CSH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H2CSH_cation') ]
RXNM['%s-%s'            % (dbse, '153'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '153')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '154'                   )] = ['%s-%s-reagent'      % (dbse, 'MeSH'),
                                                               '%s-%s-reagent'      % (dbse, 'MeSH_cation') ]
RXNM['%s-%s'            % (dbse, '154'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '154')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '155'                   )] = ['%s-%s-reagent'      % (dbse, 'MeCl'),
                                                               '%s-%s-reagent'      % (dbse, 'MeCl_cation') ]
RXNM['%s-%s'            % (dbse, '155'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '155')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '156'                   )] = ['%s-%s-reagent'      % (dbse, 'Ethanol'),
                                                               '%s-%s-reagent'      % (dbse, 'Ethanol_cation') ]
RXNM['%s-%s'            % (dbse, '156'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '156')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '157'                   )] = ['%s-%s-reagent'      % (dbse, 'Acetaldehyde'),
                                                               '%s-%s-reagent'      % (dbse, 'Acetaldehyde_cation') ]
RXNM['%s-%s'            % (dbse, '157'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '157')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '158'                   )] = ['%s-%s-reagent'      % (dbse, 'MeOF'),
                                                               '%s-%s-reagent'      % (dbse, 'MeOF_cation') ]
RXNM['%s-%s'            % (dbse, '158'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '158')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '159'                   )] = ['%s-%s-reagent'      % (dbse, 'Thiooxirane'),
                                                               '%s-%s-reagent'      % (dbse, 'Thiooxirane_cation') ]
RXNM['%s-%s'            % (dbse, '159'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '159')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '160'                   )] = ['%s-%s-reagent'      % (dbse, 'Cyanogen'),
                                                               '%s-%s-reagent'      % (dbse, 'Cyanogen_cation') ]
RXNM['%s-%s'            % (dbse, '160'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '160')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '161'                   )] = ['%s-%s-reagent'      % (dbse, 'Furan'),
                                                               '%s-%s-reagent'      % (dbse, 'Furan_cation') ]
RXNM['%s-%s'            % (dbse, '161'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '161')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '162'                   )] = ['%s-%s-reagent'      % (dbse, 'Pyrrole'),
                                                               '%s-%s-reagent'      % (dbse, 'Pyrrole_cation') ]
RXNM['%s-%s'            % (dbse, '162'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '162')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '163'                   )] = ['%s-%s-reagent'      % (dbse, 'Phenol'),
                                                               '%s-%s-reagent'      % (dbse, 'Phenol_cation') ]
RXNM['%s-%s'            % (dbse, '163'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '163')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '164'                   )] = ['%s-%s-reagent'      % (dbse, 'Aniline'),
                                                               '%s-%s-reagent'      % (dbse, 'Aniline_cation') ]
RXNM['%s-%s'            % (dbse, '164'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '164')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '165'                   )] = ['%s-%s-reagent'      % (dbse, 'B2H4'),
                                                               '%s-%s-reagent'      % (dbse, 'B2H4_cation') ]
RXNM['%s-%s'            % (dbse, '165'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '165')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '166'                   )] = ['%s-%s-reagent'      % (dbse, 'NH'),
                                                               '%s-%s-reagent'      % (dbse, 'NH_cation') ]
RXNM['%s-%s'            % (dbse, '166'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '166')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '167'                   )] = ['%s-%s-reagent'      % (dbse, 'NH2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'NH2_cation') ]
RXNM['%s-%s'            % (dbse, '167'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '167')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '168'                   )] = ['%s-%s-reagent'      % (dbse, 'N2H2'),
                                                               '%s-%s-reagent'      % (dbse, 'N2H2_cation') ]
RXNM['%s-%s'            % (dbse, '168'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '168')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '169'                   )] = ['%s-%s-reagent'      % (dbse, 'N2H3_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'N2H3_cation') ]
RXNM['%s-%s'            % (dbse, '169'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '169')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '170'                   )] = ['%s-%s-reagent'      % (dbse, 'HOF'),
                                                               '%s-%s-reagent'      % (dbse, 'HOF_cation') ]
RXNM['%s-%s'            % (dbse, '170'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '170')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '171'                   )] = ['%s-%s-reagent'      % (dbse, 'SiH2_singlet'),
                                                               '%s-%s-reagent'      % (dbse, 'SiH2_cation') ]
RXNM['%s-%s'            % (dbse, '171'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '171')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '172'                   )] = ['%s-%s-reagent'      % (dbse, 'SiH3_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'SiH3_cation') ]
RXNM['%s-%s'            % (dbse, '172'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '172')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '173'                   )] = ['%s-%s-reagent'      % (dbse, 'Si2H2'),
                                                               '%s-%s-reagent'      % (dbse, 'Si2H2_cation') ]
RXNM['%s-%s'            % (dbse, '173'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '173')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '174'                   )] = ['%s-%s-reagent'      % (dbse, 'Si2H4'),
                                                               '%s-%s-reagent'      % (dbse, 'Si2H4_cation') ]
RXNM['%s-%s'            % (dbse, '174'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '174')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '175'                   )] = ['%s-%s-reagent'      % (dbse, 'Si2H5'),
                                                               '%s-%s-reagent'      % (dbse, 'Si2H5_cation') ]
RXNM['%s-%s'            % (dbse, '175'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '175')], [-1, +1]))

ACTV['%s-%s'            % (dbse, '176'                   )] = ['%s-%s-reagent'      % (dbse, 'Disilane'),
                                                               '%s-%s-reagent'      % (dbse, 'Disilane_cation') ]
RXNM['%s-%s'            % (dbse, '176'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '176')], [-1, +1]))

# [177 - 209] G2-2 Electron Affinities (33)
ACTV['%s-%s'            % (dbse, '177'                   )] = ['%s-%s-reagent'      % (dbse, 'Li_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Li_anion') ]
RXNM['%s-%s'            % (dbse, '177'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '177')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '178'                   )] = ['%s-%s-reagent'      % (dbse, 'B_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'B_anion') ]
RXNM['%s-%s'            % (dbse, '178'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '178')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '179'                   )] = ['%s-%s-reagent'      % (dbse, 'Na_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Na_anion') ]
RXNM['%s-%s'            % (dbse, '179'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '179')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '180'                   )] = ['%s-%s-reagent'      % (dbse, 'Al_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Al_anion') ]
RXNM['%s-%s'            % (dbse, '180'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '180')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '181'                   )] = ['%s-%s-reagent'      % (dbse, 'CC'),
                                                               '%s-%s-reagent'      % (dbse, 'CC_anion') ]
RXNM['%s-%s'            % (dbse, '181'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '181')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '182'                   )] = ['%s-%s-reagent'      % (dbse, 'CCO'),
                                                               '%s-%s-reagent'      % (dbse, 'CCO_anion') ]
RXNM['%s-%s'            % (dbse, '182'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '182')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '183'                   )] = ['%s-%s-reagent'      % (dbse, 'CF2'),
                                                               '%s-%s-reagent'      % (dbse, 'CF2_anion') ]
RXNM['%s-%s'            % (dbse, '183'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '183')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '184'                   )] = ['%s-%s-reagent'      % (dbse, 'NCO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'NCO_anion') ]
RXNM['%s-%s'            % (dbse, '184'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '184')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '185'                   )] = ['%s-%s-reagent'      % (dbse, 'NO2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'NO2_anion') ]
RXNM['%s-%s'            % (dbse, '185'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '185')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '186'                   )] = ['%s-%s-reagent'      % (dbse, 'O3'),
                                                               '%s-%s-reagent'      % (dbse, 'O3_anion') ]
RXNM['%s-%s'            % (dbse, '186'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '186')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '187'                   )] = ['%s-%s-reagent'      % (dbse, 'OF_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'OF_anion') ]
RXNM['%s-%s'            % (dbse, '187'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '187')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '188'                   )] = ['%s-%s-reagent'      % (dbse, 'SO2'),
                                                               '%s-%s-reagent'      % (dbse, 'SO2_anion') ]
RXNM['%s-%s'            % (dbse, '188'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '188')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '189'                   )] = ['%s-%s-reagent'      % (dbse, 'S2O'),
                                                               '%s-%s-reagent'      % (dbse, 'S2O_anion') ]
RXNM['%s-%s'            % (dbse, '189'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '189')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '190'                   )] = ['%s-%s-reagent'      % (dbse, 'CCH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CCH_anion') ]
RXNM['%s-%s'            % (dbse, '190'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '190')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '191'                   )] = ['%s-%s-reagent'      % (dbse, 'H2CCH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H2CCH_anion') ]
RXNM['%s-%s'            % (dbse, '191'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '191')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '192'                   )] = ['%s-%s-reagent'      % (dbse, 'H2CCC'),
                                                               '%s-%s-reagent'      % (dbse, 'H2CCC_anion') ]
RXNM['%s-%s'            % (dbse, '192'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '192')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '193'                   )] = ['%s-%s-reagent'      % (dbse, 'H2CCCH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H2CCCH_anion') ]
RXNM['%s-%s'            % (dbse, '193'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '193')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '194'                   )] = ['%s-%s-reagent'      % (dbse, 'CH2CHCH2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CH2CHCH2_anion') ]
RXNM['%s-%s'            % (dbse, '194'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '194')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '195'                   )] = ['%s-%s-reagent'      % (dbse, 'HCO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'HCO_anion') ]
RXNM['%s-%s'            % (dbse, '195'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '195')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '196'                   )] = ['%s-%s-reagent'      % (dbse, 'HCF'),
                                                               '%s-%s-reagent'      % (dbse, 'HCF_anion') ]
RXNM['%s-%s'            % (dbse, '196'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '196')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '197'                   )] = ['%s-%s-reagent'      % (dbse, 'CH3O_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CH3O_anion') ]
RXNM['%s-%s'            % (dbse, '197'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '197')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '198'                   )] = ['%s-%s-reagent'      % (dbse, 'CH3S_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CH3S_anion') ]
RXNM['%s-%s'            % (dbse, '198'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '198')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '199'                   )] = ['%s-%s-reagent'      % (dbse, 'Thioformaldehyde'),
                                                               '%s-%s-reagent'      % (dbse, 'Thioformaldehyde_anion') ]
RXNM['%s-%s'            % (dbse, '199'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '199')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '200'                   )] = ['%s-%s-reagent'      % (dbse, 'H2CCN_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H2CCN_anion') ]
RXNM['%s-%s'            % (dbse, '200'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '200')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '201'                   )] = ['%s-%s-reagent'      % (dbse, 'H2CNC_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H2CNC_anion') ]
RXNM['%s-%s'            % (dbse, '201'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '201')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '202'                   )] = ['%s-%s-reagent'      % (dbse, 'HCCO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'HCCO_anion') ]
RXNM['%s-%s'            % (dbse, '202'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '202')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '203'                   )] = ['%s-%s-reagent'      % (dbse, 'H2CCHO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H2CCHO_anion') ]
RXNM['%s-%s'            % (dbse, '203'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '203')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '204'                   )] = ['%s-%s-reagent'      % (dbse, 'cAcetyl_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'cAcetyl_anion') ]
RXNM['%s-%s'            % (dbse, '204'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '204')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '205'                   )] = ['%s-%s-reagent'      % (dbse, 'CH3CH2O_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CH3CH2O_anion') ]
RXNM['%s-%s'            % (dbse, '205'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '205')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '206'                   )] = ['%s-%s-reagent'      % (dbse, 'CH3CH2S_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'CH3CH2S_anion') ]
RXNM['%s-%s'            % (dbse, '206'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '206')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '207'                   )] = ['%s-%s-reagent'      % (dbse, 'LiH'),
                                                               '%s-%s-reagent'      % (dbse, 'LiH_anion') ]
RXNM['%s-%s'            % (dbse, '207'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '207')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '208'                   )] = ['%s-%s-reagent'      % (dbse, 'HNO'),
                                                               '%s-%s-reagent'      % (dbse, 'HNO_anion') ]
RXNM['%s-%s'            % (dbse, '208'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '208')], [+1, -1]))

ACTV['%s-%s'            % (dbse, '209'                   )] = ['%s-%s-reagent'      % (dbse, 'HOO_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'HOO_anion') ]
RXNM['%s-%s'            % (dbse, '209'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '209')], [+1, -1]))

# [210 - 302] G2-2 Enthalpies of Formation (93)
ACTV['%s-%s'            % (dbse, '210'                   )] = ['%s-%s-reagent'      % (dbse, 'BF3'),
                                                               '%s-%s-reagent'      % (dbse, 'B_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '210'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '210')], [-1, +1, +3]))

ACTV['%s-%s'            % (dbse, '211'                   )] = ['%s-%s-reagent'      % (dbse, 'BCl3'),
                                                               '%s-%s-reagent'      % (dbse, 'B_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '211'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '211')], [-1, +1, +3]))

ACTV['%s-%s'            % (dbse, '212'                   )] = ['%s-%s-reagent'      % (dbse, 'AlF3'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Al_radical') ]
RXNM['%s-%s'            % (dbse, '212'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '212')], [-1, +3, +1]))

ACTV['%s-%s'            % (dbse, '213'                   )] = ['%s-%s-reagent'      % (dbse, 'AlCl3'),
                                                               '%s-%s-reagent'      % (dbse, 'Al_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '213'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '213')], [-1, +1, +3]))

ACTV['%s-%s'            % (dbse, '214'                   )] = ['%s-%s-reagent'      % (dbse, 'CF4'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '214'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '214')], [-1, +1, +4]))

ACTV['%s-%s'            % (dbse, '215'                   )] = ['%s-%s-reagent'      % (dbse, 'CCl4'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '215'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '215')], [-1, +1, +4]))

ACTV['%s-%s'            % (dbse, '216'                   )] = ['%s-%s-reagent'      % (dbse, 'OCS'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '216'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '216')], [-1, +1, +1, +1]))

ACTV['%s-%s'            % (dbse, '217'                   )] = ['%s-%s-reagent'      % (dbse, 'CS2'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '217'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '217')], [-1, +1, +2]))

ACTV['%s-%s'            % (dbse, '218'                   )] = ['%s-%s-reagent'      % (dbse, 'COF2'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '218'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '218')], [-1, +1, +1, +2]))

ACTV['%s-%s'            % (dbse, '219'                   )] = ['%s-%s-reagent'      % (dbse, 'SiF4'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '219'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '219')], [-1, +4, +1]))

ACTV['%s-%s'            % (dbse, '220'                   )] = ['%s-%s-reagent'      % (dbse, 'SiCl4'),
                                                               '%s-%s-reagent'      % (dbse, 'Si'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '220'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '220')], [-1, +1, +4]))

ACTV['%s-%s'            % (dbse, '221'                   )] = ['%s-%s-reagent'      % (dbse, 'N2O'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '221'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '221')], [-1, +2, +1]))

ACTV['%s-%s'            % (dbse, '222'                   )] = ['%s-%s-reagent'      % (dbse, 'ClNO'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '222'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '222')], [-1, +1, +1, +1]))

ACTV['%s-%s'            % (dbse, '223'                   )] = ['%s-%s-reagent'      % (dbse, 'NF3'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '223'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '223')], [-1, +1, +3]))

ACTV['%s-%s'            % (dbse, '224'                   )] = ['%s-%s-reagent'      % (dbse, 'PF3'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'P_radical') ]
RXNM['%s-%s'            % (dbse, '224'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '224')], [-1, +3, +1]))

ACTV['%s-%s'            % (dbse, '225'                   )] = ['%s-%s-reagent'      % (dbse, 'O3'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '225'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '225')], [-1, +3]))

ACTV['%s-%s'            % (dbse, '226'                   )] = ['%s-%s-reagent'      % (dbse, 'F2O'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '226'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '226')], [-1, +1, +2]))

ACTV['%s-%s'            % (dbse, '227'                   )] = ['%s-%s-reagent'      % (dbse, 'ClF3'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '227'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '227')], [-1, +3, +1]))

ACTV['%s-%s'            % (dbse, '228'                   )] = ['%s-%s-reagent'      % (dbse, 'Tetrafluoroethene'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '228'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '228')], [-1, +2, +4]))

ACTV['%s-%s'            % (dbse, '229'                   )] = ['%s-%s-reagent'      % (dbse, 'Tetrachloroethene'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '229'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '229')], [-1, +2, +4]))

ACTV['%s-%s'            % (dbse, '230'                   )] = ['%s-%s-reagent'      % (dbse, 'CF3CN'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '230'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '230')], [-1, +2, +1, +3]))

ACTV['%s-%s'            % (dbse, '231'                   )] = ['%s-%s-reagent'      % (dbse, 'Propyne'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '231'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '231')], [-1, +4, +3]))

ACTV['%s-%s'            % (dbse, '232'                   )] = ['%s-%s-reagent'      % (dbse, 'Allene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '232'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '232')], [-1, +4, +3]))

ACTV['%s-%s'            % (dbse, '233'                   )] = ['%s-%s-reagent'      % (dbse, 'Cyclopropene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '233'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '233')], [-1, +4, +3]))

ACTV['%s-%s'            % (dbse, '234'                   )] = ['%s-%s-reagent'      % (dbse, 'Propene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '234'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '234')], [-1, +6, +3]))

ACTV['%s-%s'            % (dbse, '235'                   )] = ['%s-%s-reagent'      % (dbse, 'Cyclopropane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '235'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '235')], [-1, +6, +3]))

ACTV['%s-%s'            % (dbse, '236'                   )] = ['%s-%s-reagent'      % (dbse, 'Propane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '236'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '236')], [-1, +8, +3]))

ACTV['%s-%s'            % (dbse, '237'                   )] = ['%s-%s-reagent'      % (dbse, 't13Butadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '237'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '237')], [-1, +6, +4]))

ACTV['%s-%s'            % (dbse, '238'                   )] = ['%s-%s-reagent'      % (dbse, '2Butyne'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '238'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '238')], [-1, +6, +4]))

ACTV['%s-%s'            % (dbse, '239'                   )] = ['%s-%s-reagent'      % (dbse, 'Methylenecyclopropane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '239'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '239')], [-1, +6, +4]))

ACTV['%s-%s'            % (dbse, '240'                   )] = ['%s-%s-reagent'      % (dbse, 'Bicyclo110butane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '240'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '240')], [-1, +6, +4]))

ACTV['%s-%s'            % (dbse, '241'                   )] = ['%s-%s-reagent'      % (dbse, 'Cyclobutene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '241'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '241')], [-1, +6, +4]))

ACTV['%s-%s'            % (dbse, '242'                   )] = ['%s-%s-reagent'      % (dbse, 'Cyclobutane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '242'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '242')], [-1, +8, +4]))

ACTV['%s-%s'            % (dbse, '243'                   )] = ['%s-%s-reagent'      % (dbse, 'Isobutene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '243'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '243')], [-1, +8, +4]))

ACTV['%s-%s'            % (dbse, '244'                   )] = ['%s-%s-reagent'      % (dbse, 'tButane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '244'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '244')], [-1, +10, +4]))

ACTV['%s-%s'            % (dbse, '245'                   )] = ['%s-%s-reagent'      % (dbse, 'Isobutane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '245'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '245')], [-1, +10, +4]))

ACTV['%s-%s'            % (dbse, '246'                   )] = ['%s-%s-reagent'      % (dbse, 'Spiropentane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '246'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '246')], [-1, +8, +5]))

ACTV['%s-%s'            % (dbse, '247'                   )] = ['%s-%s-reagent'      % (dbse, 'Benzene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '247'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '247')], [-1, +6, +6]))

ACTV['%s-%s'            % (dbse, '248'                   )] = ['%s-%s-reagent'      % (dbse, 'CH2F2'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '248'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '248')], [-1, +2, +1, +2]))

ACTV['%s-%s'            % (dbse, '249'                   )] = ['%s-%s-reagent'      % (dbse, 'CHF3'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '249'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '249')], [-1, +1, +1, +3]))

ACTV['%s-%s'            % (dbse, '250'                   )] = ['%s-%s-reagent'      % (dbse, 'CH2Cl2'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '250'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '250')], [-1, +2, +1, +2]))

ACTV['%s-%s'            % (dbse, '251'                   )] = ['%s-%s-reagent'      % (dbse, 'CHCl3'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '251'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '251')], [-1, +1, +1, +3]))

ACTV['%s-%s'            % (dbse, '252'                   )] = ['%s-%s-reagent'      % (dbse, 'Methylamine'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '252'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '252')], [-1, +5, +1, +1]))

ACTV['%s-%s'            % (dbse, '253'                   )] = ['%s-%s-reagent'      % (dbse, 'MeCN'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '253'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '253')], [-1, +3, +2, +1]))

ACTV['%s-%s'            % (dbse, '254'                   )] = ['%s-%s-reagent'      % (dbse, 'Nitromethane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '254'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '254')], [-1, +3, +1, +1, +2]))

ACTV['%s-%s'            % (dbse, '255'                   )] = ['%s-%s-reagent'      % (dbse, 'Methylnitrite'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '255'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '255')], [-1, +3, +1, +1, +2]))

ACTV['%s-%s'            % (dbse, '256'                   )] = ['%s-%s-reagent'      % (dbse, 'Methylsilane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Si') ]
RXNM['%s-%s'            % (dbse, '256'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '256')], [-1, +6, +1, +1]))

ACTV['%s-%s'            % (dbse, '257'                   )] = ['%s-%s-reagent'      % (dbse, 'FormicAcid'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '257'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '257')], [-1, +2, +1, +2]))

ACTV['%s-%s'            % (dbse, '258'                   )] = ['%s-%s-reagent'      % (dbse, 'Methylformate'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '258'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '258')], [-1, +4, +2, +2]))

ACTV['%s-%s'            % (dbse, '259'                   )] = ['%s-%s-reagent'      % (dbse, 'Acetamide'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '259'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '259')], [-1, +5, +2, +1, +1]))

ACTV['%s-%s'            % (dbse, '260'                   )] = ['%s-%s-reagent'      % (dbse, 'Aziridine'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '260'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '260')], [-1, +5, +2, +1]))

ACTV['%s-%s'            % (dbse, '261'                   )] = ['%s-%s-reagent'      % (dbse, 'Cyanogen'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '261'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '261')], [-1, +2, +2]))

ACTV['%s-%s'            % (dbse, '262'                   )] = ['%s-%s-reagent'      % (dbse, 'Dimethylamine'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '262'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '262')], [-1, +7, +2, +1]))

ACTV['%s-%s'            % (dbse, '263'                   )] = ['%s-%s-reagent'      % (dbse, 'tEthylamine'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '263'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '263')], [-1, +7, +2, +1]))

ACTV['%s-%s'            % (dbse, '264'                   )] = ['%s-%s-reagent'      % (dbse, 'Ketene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '264'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '264')], [-1, +2, +2, +1]))

ACTV['%s-%s'            % (dbse, '265'                   )] = ['%s-%s-reagent'      % (dbse, 'Oxirane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '265'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '265')], [-1, +4, +2, +1]))

ACTV['%s-%s'            % (dbse, '266'                   )] = ['%s-%s-reagent'      % (dbse, 'Acetaldehyde'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '266'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '266')], [-1, +4, +2, +1]))

ACTV['%s-%s'            % (dbse, '267'                   )] = ['%s-%s-reagent'      % (dbse, 'Glyoxal'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '267'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '267')], [-1, +2, +2, +2]))

ACTV['%s-%s'            % (dbse, '268'                   )] = ['%s-%s-reagent'      % (dbse, 'Ethanol'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '268'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '268')], [-1, +6, +2, +1]))

ACTV['%s-%s'            % (dbse, '269'                   )] = ['%s-%s-reagent'      % (dbse, 'Dimethylether'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '269'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '269')], [-1, +6, +2, +1]))

ACTV['%s-%s'            % (dbse, '270'                   )] = ['%s-%s-reagent'      % (dbse, 'Thiooxirane'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '270'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '270')], [-1, +4, +2, +1]))

ACTV['%s-%s'            % (dbse, '271'                   )] = ['%s-%s-reagent'      % (dbse, 'Dimethylsulfoxide'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '271'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '271')], [-1, +6, +2, +1, +1]))

ACTV['%s-%s'            % (dbse, '272'                   )] = ['%s-%s-reagent'      % (dbse, 'Thioethanol'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '272'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '272')], [-1, +6, +2, +1]))

ACTV['%s-%s'            % (dbse, '273'                   )] = ['%s-%s-reagent'      % (dbse, 'Dimethylthioether'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '273'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '273')], [-1, +6, +2, +1]))

ACTV['%s-%s'            % (dbse, '274'                   )] = ['%s-%s-reagent'      % (dbse, 'Vinylfluoride'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '274'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '274')], [-1, +3, +2, +1]))

ACTV['%s-%s'            % (dbse, '275'                   )] = ['%s-%s-reagent'      % (dbse, 'Ethylchloride'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '275'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '275')], [-1, +5, +2, +1]))

ACTV['%s-%s'            % (dbse, '276'                   )] = ['%s-%s-reagent'      % (dbse, 'Vinylchloride'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '276'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '276')], [-1, +3, +2, +1]))

ACTV['%s-%s'            % (dbse, '277'                   )] = ['%s-%s-reagent'      % (dbse, 'Vinylcyanide'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '277'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '277')], [-1, +3, +3, +1]))

ACTV['%s-%s'            % (dbse, '278'                   )] = ['%s-%s-reagent'      % (dbse, 'Acetone'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '278'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '278')], [-1, +6, +3, +1]))

ACTV['%s-%s'            % (dbse, '279'                   )] = ['%s-%s-reagent'      % (dbse, 'AceticAcid'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '279'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '279')], [-1, +4, +2, +2]))

ACTV['%s-%s'            % (dbse, '280'                   )] = ['%s-%s-reagent'      % (dbse, 'Acetylfluoride'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'F_radical') ]
RXNM['%s-%s'            % (dbse, '280'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '280')], [-1, +3, +2, +1, +1]))

ACTV['%s-%s'            % (dbse, '281'                   )] = ['%s-%s-reagent'      % (dbse, 'Acetylchloride'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '281'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '281')], [-1, +3, +2, +1, +1]))

ACTV['%s-%s'            % (dbse, '282'                   )] = ['%s-%s-reagent'      % (dbse, 'Propylchloride'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'Cl_radical') ]
RXNM['%s-%s'            % (dbse, '282'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '282')], [-1, +7, +3, +1]))

ACTV['%s-%s'            % (dbse, '283'                   )] = ['%s-%s-reagent'      % (dbse, 'Isopropanol'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '283'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '283')], [-1, +8, +3, +1]))

ACTV['%s-%s'            % (dbse, '284'                   )] = ['%s-%s-reagent'      % (dbse, 'Methylethylether'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '284'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '284')], [-1, +8, +3, +1]))

ACTV['%s-%s'            % (dbse, '285'                   )] = ['%s-%s-reagent'      % (dbse, 'Trimethylamine'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '285'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '285')], [-1, +9, +3, +1]))

ACTV['%s-%s'            % (dbse, '286'                   )] = ['%s-%s-reagent'      % (dbse, 'Furan'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '286'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '286')], [-1, +4, +4, +1]))

ACTV['%s-%s'            % (dbse, '287'                   )] = ['%s-%s-reagent'      % (dbse, 'Thiophene'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '287'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '287')], [-1, +4, +4, +1]))

ACTV['%s-%s'            % (dbse, '288'                   )] = ['%s-%s-reagent'      % (dbse, 'Pyrrole'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '288'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '288')], [-1, +5, +4, +1]))

ACTV['%s-%s'            % (dbse, '289'                   )] = ['%s-%s-reagent'      % (dbse, 'Pyridine'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical') ]
RXNM['%s-%s'            % (dbse, '289'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '289')], [-1, +5, +5, +1]))

ACTV['%s-%s'            % (dbse, '290'                   )] = ['%s-%s-reagent'      % (dbse, 'H2'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical') ]
RXNM['%s-%s'            % (dbse, '290'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '290')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '291'                   )] = ['%s-%s-reagent'      % (dbse, 'SH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '291'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '291')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '292'                   )] = ['%s-%s-reagent'      % (dbse, 'CCH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '292'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '292')], [-1, +1, +2]))

ACTV['%s-%s'            % (dbse, '293'                   )] = ['%s-%s-reagent'      % (dbse, 'H2CCH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '293'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '293')], [-1, +3, +2]))

ACTV['%s-%s'            % (dbse, '294'                   )] = ['%s-%s-reagent'      % (dbse, 'cAcetyl_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '294'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '294')], [-1, +3, +2, +1]))

ACTV['%s-%s'            % (dbse, '295'                   )] = ['%s-%s-reagent'      % (dbse, 'H2COH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '295'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '295')], [-1, +3, +1, +1]))

ACTV['%s-%s'            % (dbse, '296'                   )] = ['%s-%s-reagent'      % (dbse, 'CH3O_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '296'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '296')], [-1, +3, +1, +1]))

ACTV['%s-%s'            % (dbse, '297'                   )] = ['%s-%s-reagent'      % (dbse, 'CH3CH2O_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '297'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '297')], [-1, +5, +2, +1]))

ACTV['%s-%s'            % (dbse, '298'                   )] = ['%s-%s-reagent'      % (dbse, 'CH3S_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C'),
                                                               '%s-%s-reagent'      % (dbse, 'S') ]
RXNM['%s-%s'            % (dbse, '298'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '298')], [-1, +3, +1, +1]))

ACTV['%s-%s'            % (dbse, '299'                   )] = ['%s-%s-reagent'      % (dbse, 'C2H5_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '299'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '299')], [-1, +5, +2]))

ACTV['%s-%s'            % (dbse, '300'                   )] = ['%s-%s-reagent'      % (dbse, 'Me2CH_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '300'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '300')], [-1, +7, +3]))

ACTV['%s-%s'            % (dbse, '301'                   )] = ['%s-%s-reagent'      % (dbse, 'Me3C_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'H_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'C') ]
RXNM['%s-%s'            % (dbse, '301'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '301')], [-1, +9, +4]))

ACTV['%s-%s'            % (dbse, '302'                   )] = ['%s-%s-reagent'      % (dbse, 'NO2_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'N_radical'),
                                                               '%s-%s-reagent'      % (dbse, 'O') ]
RXNM['%s-%s'            % (dbse, '302'                   )] = dict(zip(ACTV['%s-%s' % (dbse, '302')], [-1, +1, +2]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
# [  1 -  38] G2-1 Ionization Energies (38)
# Ionization Energy for species M is the energy of the reaction M --> M+ + e-  [E(cation) - E(neutral)]
# Reference ionization energies [eV] from Table II of Curtiss et al. JCP 109 42 (1998).
BIND['%s-%s'            % (dbse, '1'                     )] =    5.39
BIND['%s-%s'            % (dbse, '2'                     )] =    9.32
BIND['%s-%s'            % (dbse, '3'                     )] =    8.30
BIND['%s-%s'            % (dbse, '4'                     )] =   11.26
BIND['%s-%s'            % (dbse, '5'                     )] =   14.54
BIND['%s-%s'            % (dbse, '6'                     )] =   13.61
BIND['%s-%s'            % (dbse, '7'                     )] =   17.42
BIND['%s-%s'            % (dbse, '8'                     )] =    5.14
BIND['%s-%s'            % (dbse, '9'                     )] =    7.65
BIND['%s-%s'            % (dbse, '10'                    )] =    5.98
BIND['%s-%s'            % (dbse, '11'                    )] =    8.15
BIND['%s-%s'            % (dbse, '12'                    )] =   10.49
BIND['%s-%s'            % (dbse, '13'                    )] =   10.36
BIND['%s-%s'            % (dbse, '14'                    )] =   12.97
BIND['%s-%s'            % (dbse, '15'                    )] =   12.62
BIND['%s-%s'            % (dbse, '16'                    )] =   10.18
BIND['%s-%s'            % (dbse, '17'                    )] =   13.01
BIND['%s-%s'            % (dbse, '18'                    )] =   12.62
BIND['%s-%s'            % (dbse, '19'                    )] =   16.04
BIND['%s-%s'            % (dbse, '20'                    )] =   11.00
BIND['%s-%s'            % (dbse, '21'                    )] =   10.15
BIND['%s-%s'            % (dbse, '22'                    )] =    9.82
BIND['%s-%s'            % (dbse, '23'                    )] =    9.87
BIND['%s-%s'            % (dbse, '24'                    )] =   10.37
BIND['%s-%s'            % (dbse, '25'                    )] =   10.47
BIND['%s-%s'            % (dbse, '26'                    )] =   12.78
BIND['%s-%s'            % (dbse, '27'                    )] =   12.75
BIND['%s-%s'            % (dbse, '28'                    )] =   11.40
BIND['%s-%s'            % (dbse, '29'                    )] =   10.51
BIND['%s-%s'            % (dbse, '30'                    )] =   14.01
BIND['%s-%s'            % (dbse, '31'                    )] =   15.58
BIND['%s-%s'            % (dbse, '32'                    )] =   16.70
BIND['%s-%s'            % (dbse, '33'                    )] =   12.07
BIND['%s-%s'            % (dbse, '34'                    )] =   10.53
BIND['%s-%s'            % (dbse, '35'                    )] =    9.36
BIND['%s-%s'            % (dbse, '36'                    )] =   11.50
BIND['%s-%s'            % (dbse, '37'                    )] =   12.66
BIND['%s-%s'            % (dbse, '38'                    )] =   11.33
for ii in range(1, 39):
   BIND['%s-%s'            % (dbse, str(ii)                 )] *=  23.06  # convert from [eV] to [kcal/mol]
# [ 39 -  63] G2-1 Electron Affinities (25)
# Electron Affinity for species M is the negative of the reaction M + e- --> M-  [E(neutral) - E(anion)]
# Reference electron affinities [eV] from Table III of Curtiss et al. JCP 109 42 (1998).
BIND['%s-%s'            % (dbse, '39'                    )] =    1.26
BIND['%s-%s'            % (dbse, '40'                    )] =    1.46
BIND['%s-%s'            % (dbse, '41'                    )] =    3.40
BIND['%s-%s'            % (dbse, '42'                    )] =    1.39
BIND['%s-%s'            % (dbse, '43'                    )] =    0.75
BIND['%s-%s'            % (dbse, '44'                    )] =    2.08
BIND['%s-%s'            % (dbse, '45'                    )] =    3.62
BIND['%s-%s'            % (dbse, '46'                    )] =    1.24
BIND['%s-%s'            % (dbse, '47'                    )] =    0.65
BIND['%s-%s'            % (dbse, '48'                    )] =    0.08
BIND['%s-%s'            % (dbse, '49'                    )] =    0.38
BIND['%s-%s'            % (dbse, '50'                    )] =    0.77
BIND['%s-%s'            % (dbse, '51'                    )] =    1.83
BIND['%s-%s'            % (dbse, '52'                    )] =    1.28
BIND['%s-%s'            % (dbse, '53'                    )] =    1.12
BIND['%s-%s'            % (dbse, '54'                    )] =    1.41
BIND['%s-%s'            % (dbse, '55'                    )] =    1.03
BIND['%s-%s'            % (dbse, '56'                    )] =    1.27
BIND['%s-%s'            % (dbse, '57'                    )] =    2.36
BIND['%s-%s'            % (dbse, '58'                    )] =    0.44
BIND['%s-%s'            % (dbse, '59'                    )] =    0.02
BIND['%s-%s'            % (dbse, '60'                    )] =    3.86
BIND['%s-%s'            % (dbse, '61'                    )] =    1.09
BIND['%s-%s'            % (dbse, '62'                    )] =    1.66
BIND['%s-%s'            % (dbse, '63'                    )] =    2.39
for ii in range(39, 64):
   BIND['%s-%s'            % (dbse, str(ii)                 )] *=  23.06  # convert from [eV] to [kcal/mol]
# [ 64 -  71] G2-1 Proton Affinities (8) in kcal/mol from
# J. Chem. Phys. 109, 7764 (1998) doi: 10.1063/1.477422
# Proton affinity for species M is the negative enthalpy of the reaction M + H+ --> MH+
BIND['%s-%s'            % (dbse, '64'                    )] =    202.5
BIND['%s-%s'            % (dbse, '65'                    )] =    165.1
BIND['%s-%s'            % (dbse, '66'                    )] =    152.3
BIND['%s-%s'            % (dbse, '67'                    )] =    154.0
BIND['%s-%s'            % (dbse, '68'                    )] =    187.1
BIND['%s-%s'            % (dbse, '69'                    )] =    168.8
BIND['%s-%s'            % (dbse, '70'                    )] =    133.6
BIND['%s-%s'            % (dbse, '71'                    )] =    100.8
# [ 72 - 126] G2-1 Enthalpies of Formation (55)
BIND['%s-%s'            % (dbse, '72'                    )] =    32.7
BIND['%s-%s'            % (dbse, '73'                    )] =    82.6
BIND['%s-%s'            % (dbse, '74'                    )] =    141.1
BIND['%s-%s'            % (dbse, '75'                    )] =    94.6
BIND['%s-%s'            % (dbse, '76'                    )] =    101.3
BIND['%s-%s'            % (dbse, '77'                    )] =    35.7
BIND['%s-%s'            % (dbse, '78'                    )] =   -16.7
BIND['%s-%s'            % (dbse, '79'                    )] =    86.2
BIND['%s-%s'            % (dbse, '80'                    )] =    45.7
BIND['%s-%s'            % (dbse, '81'                    )] =   -9.1 
BIND['%s-%s'            % (dbse, '82'                    )] =    9.0
BIND['%s-%s'            % (dbse, '83'                    )] =   -57.4 
BIND['%s-%s'            % (dbse, '84'                    )] =   -66.2
BIND['%s-%s'            % (dbse, '85'                    )] =    62.7
BIND['%s-%s'            % (dbse, '86'                    )] =    86.1
BIND['%s-%s'            % (dbse, '87'                    )] =    48.0
BIND['%s-%s'            % (dbse, '88'                    )] =    8.3
BIND['%s-%s'            % (dbse, '89'                    )] =    33.8
BIND['%s-%s'            % (dbse, '90'                    )] =    3.9
BIND['%s-%s'            % (dbse, '91'                    )] =   -4.1 
BIND['%s-%s'            % (dbse, '92'                    )] =   -22.4 
BIND['%s-%s'            % (dbse, '93'                    )] =    49.5
BIND['%s-%s'            % (dbse, '94'                    )] =   -81.4
BIND['%s-%s'            % (dbse, '95'                    )] =    56.0
BIND['%s-%s'            % (dbse, '96'                    )] =    14.8
BIND['%s-%s'            % (dbse, '97'                    )] =   -16.8 
BIND['%s-%s'            % (dbse, '98'                    )] =    106.5
BIND['%s-%s'            % (dbse, '99'                    )] =    31.3
BIND['%s-%s'            % (dbse, '100'                   )] =   -29.0 
BIND['%s-%s'            % (dbse, '101'                   )] =    9.2
BIND['%s-%s'            % (dbse, '102'                   )] =   -27.0
BIND['%s-%s'            % (dbse, '103'                   )] =   -46.8
BIND['%s-%s'            % (dbse, '104'                   )] =    1.3
BIND['%s-%s'            % (dbse, '105'                   )] =    27.2
BIND['%s-%s'            % (dbse, '106'                   )] =    21.0
BIND['%s-%s'            % (dbse, '107'                   )] =    2.4
BIND['%s-%s'            % (dbse, '108'                   )] =   -30.8
BIND['%s-%s'            % (dbse, '109'                   )] =    0.3
BIND['%s-%s'            % (dbse, '110'                   )] =   -96.7
BIND['%s-%s'            % (dbse, '111'                   )] =    32.2
BIND['%s-%s'            % (dbse, '112'                   )] =    139.6
BIND['%s-%s'            % (dbse, '113'                   )] =    36.1
BIND['%s-%s'            % (dbse, '114'                   )] =    33.9
BIND['%s-%s'            % (dbse, '115'                   )] =    1.4
BIND['%s-%s'            % (dbse, '116'                   )] =   -44.5
BIND['%s-%s'            % (dbse, '117'                   )] =   -23.2
BIND['%s-%s'            % (dbse, '118'                   )] =    65.1
BIND['%s-%s'            % (dbse, '119'                   )] =    3.8
BIND['%s-%s'            % (dbse, '120'                   )] =    26.4
BIND['%s-%s'            % (dbse, '121'                   )] =   -13.9
BIND['%s-%s'            % (dbse, '122'                   )] =    20.0
BIND['%s-%s'            % (dbse, '123'                   )] =   -18.6
BIND['%s-%s'            % (dbse, '124'                   )] =   -2.9
BIND['%s-%s'            % (dbse, '125'                   )] =   -17.6
BIND['%s-%s'            % (dbse, '126'                   )] =   -65.3
# [127 - 176] G2-2 Ionization Energies (50)
# Ionization Energy for species M is the energy of the reaction M --> M+ + e-  [E(cation) - E(neutral)]
# Reference ionization energies [eV] from Table II of Curtiss et al. JCP 109 42 (1998).
BIND['%s-%s'            % (dbse, '127'                   )] =   13.60  
BIND['%s-%s'            % (dbse, '128'                   )] =   24.59 
BIND['%s-%s'            % (dbse, '129'                   )] =   21.56 
BIND['%s-%s'            % (dbse, '130'                   )] =   15.76 
BIND['%s-%s'            % (dbse, '131'                   )] =   15.56
BIND['%s-%s'            % (dbse, '132'                   )] =   11.60
BIND['%s-%s'            % (dbse, '133'                   )] =   12.07
BIND['%s-%s'            % (dbse, '134'                   )] =   13.77
BIND['%s-%s'            % (dbse, '135'                   )] =   11.42
BIND['%s-%s'            % (dbse, '136'                   )] =   11.17
BIND['%s-%s'            % (dbse, '137'                   )] =   10.07
BIND['%s-%s'            % (dbse, '138'                   )] =   10.40
BIND['%s-%s'            % (dbse, '139'                   )] =    9.84
BIND['%s-%s'            % (dbse, '140'                   )] =    8.12
BIND['%s-%s'            % (dbse, '141'                   )] =    9.67
BIND['%s-%s'            % (dbse, '142'                   )] =    9.69
BIND['%s-%s'            % (dbse, '143'                   )] =    7.37
BIND['%s-%s'            % (dbse, '144'                   )] =    9.25
BIND['%s-%s'            % (dbse, '145'                   )] =    8.83
BIND['%s-%s'            % (dbse, '146'                   )] =   13.60
BIND['%s-%s'            % (dbse, '147'                   )] =    8.14
BIND['%s-%s'            % (dbse, '148'                   )] =    7.55
BIND['%s-%s'            % (dbse, '149'                   )] =   10.73
BIND['%s-%s'            % (dbse, '150'                   )] =   10.85
BIND['%s-%s'            % (dbse, '151'                   )] =   12.47
BIND['%s-%s'            % (dbse, '152'                   )] =    9.38
BIND['%s-%s'            % (dbse, '153'                   )] =    7.54
BIND['%s-%s'            % (dbse, '154'                   )] =    9.44
BIND['%s-%s'            % (dbse, '155'                   )] =   11.22
BIND['%s-%s'            % (dbse, '156'                   )] =   10.47
BIND['%s-%s'            % (dbse, '157'                   )] =   10.23
BIND['%s-%s'            % (dbse, '158'                   )] =   11.34
BIND['%s-%s'            % (dbse, '159'                   )] =    9.05
BIND['%s-%s'            % (dbse, '160'                   )] =   13.37
BIND['%s-%s'            % (dbse, '161'                   )] =    8.83
BIND['%s-%s'            % (dbse, '162'                   )] =    8.21
BIND['%s-%s'            % (dbse, '163'                   )] =    8.51
BIND['%s-%s'            % (dbse, '164'                   )] =    7.72
BIND['%s-%s'            % (dbse, '165'                   )] =    9.70
BIND['%s-%s'            % (dbse, '166'                   )] =   13.49
BIND['%s-%s'            % (dbse, '167'                   )] =   11.14
BIND['%s-%s'            % (dbse, '168'                   )] =    9.59
BIND['%s-%s'            % (dbse, '169'                   )] =    7.61
BIND['%s-%s'            % (dbse, '170'                   )] =   12.71
BIND['%s-%s'            % (dbse, '171'                   )] =    9.15
BIND['%s-%s'            % (dbse, '172'                   )] =    8.14
BIND['%s-%s'            % (dbse, '173'                   )] =    8.20
BIND['%s-%s'            % (dbse, '174'                   )] =    8.09
BIND['%s-%s'            % (dbse, '175'                   )] =    7.60
BIND['%s-%s'            % (dbse, '176'                   )] =    9.74
for ii in range(127, 177):
   BIND['%s-%s'            % (dbse, str(ii)                 )] *=  23.06  # convert from [eV] to [kcal/mol]
# [177 - 209] G2-2 Electron Affinities (33)
# Electron Affinity for species M is the negative of the reaction M + e- --> M-  [E(neutral) - E(anion)]
# Reference electron affinities [eV] from Table III of Curtiss et al. JCP 109 42 (1998).
BIND['%s-%s'            % (dbse, '177'                   )] =    0.62
BIND['%s-%s'            % (dbse, '178'                   )] =    0.28
BIND['%s-%s'            % (dbse, '179'                   )] =    0.55
BIND['%s-%s'            % (dbse, '180'                   )] =    0.44
BIND['%s-%s'            % (dbse, '181'                   )] =    3.27
BIND['%s-%s'            % (dbse, '182'                   )] =    2.29
BIND['%s-%s'            % (dbse, '183'                   )] =    0.18
BIND['%s-%s'            % (dbse, '184'                   )] =    3.61
BIND['%s-%s'            % (dbse, '185'                   )] =    2.27
BIND['%s-%s'            % (dbse, '186'                   )] =    2.10
BIND['%s-%s'            % (dbse, '187'                   )] =    2.27
BIND['%s-%s'            % (dbse, '188'                   )] =    1.11
BIND['%s-%s'            % (dbse, '189'                   )] =    1.88
BIND['%s-%s'            % (dbse, '190'                   )] =    2.97
BIND['%s-%s'            % (dbse, '191'                   )] =    0.67
BIND['%s-%s'            % (dbse, '192'                   )] =    1.79
BIND['%s-%s'            % (dbse, '193'                   )] =    0.89
BIND['%s-%s'            % (dbse, '194'                   )] =    0.47
BIND['%s-%s'            % (dbse, '195'                   )] =    0.31
BIND['%s-%s'            % (dbse, '196'                   )] =    0.54
BIND['%s-%s'            % (dbse, '197'                   )] =    1.57
BIND['%s-%s'            % (dbse, '198'                   )] =    1.87
BIND['%s-%s'            % (dbse, '199'                   )] =    0.47
BIND['%s-%s'            % (dbse, '200'                   )] =    1.54
BIND['%s-%s'            % (dbse, '201'                   )] =    1.06
BIND['%s-%s'            % (dbse, '202'                   )] =    2.35
BIND['%s-%s'            % (dbse, '203'                   )] =    1.82
BIND['%s-%s'            % (dbse, '204'                   )] =    0.42
BIND['%s-%s'            % (dbse, '205'                   )] =    1.71
BIND['%s-%s'            % (dbse, '206'                   )] =    1.95
BIND['%s-%s'            % (dbse, '207'                   )] =    0.34
BIND['%s-%s'            % (dbse, '208'                   )] =    0.34
BIND['%s-%s'            % (dbse, '209'                   )] =    1.08
for ii in range(177, 210):
   BIND['%s-%s'            % (dbse, str(ii)                 )] *=  23.06  # convert from [eV] to [kcal/mol]
# [210 - 302] G2-2 Enthalpies of Formation (93)
BIND['%s-%s'            % (dbse, '210'                   )] =   -270.8
BIND['%s-%s'            % (dbse, '211'                   )] =   -98.2
BIND['%s-%s'            % (dbse, '212'                   )] =   -286.8
BIND['%s-%s'            % (dbse, '213'                   )] =   -142.1
BIND['%s-%s'            % (dbse, '214'                   )] =   -227.2
BIND['%s-%s'            % (dbse, '215'                   )] =   -25.2
BIND['%s-%s'            % (dbse, '216'                   )] =   -35.9
BIND['%s-%s'            % (dbse, '217'                   )] =    25.6
BIND['%s-%s'            % (dbse, '218'                   )] =   -147.8
BIND['%s-%s'            % (dbse, '219'                   )] =   -377.7
BIND['%s-%s'            % (dbse, '220'                   )] =   -161.8
BIND['%s-%s'            % (dbse, '221'                   )] =    21.0
BIND['%s-%s'            % (dbse, '222'                   )] =    12.0
BIND['%s-%s'            % (dbse, '223'                   )] =   -33.8
BIND['%s-%s'            % (dbse, '224'                   )] =   -222.4
BIND['%s-%s'            % (dbse, '225'                   )] =    33.7
BIND['%s-%s'            % (dbse, '226'                   )] =    5.9
BIND['%s-%s'            % (dbse, '227'                   )] =   -37.4
BIND['%s-%s'            % (dbse, '228'                   )] =   -164.8
BIND['%s-%s'            % (dbse, '229'                   )] =   -7.3
BIND['%s-%s'            % (dbse, '230'                   )] =   -122.3
BIND['%s-%s'            % (dbse, '231'                   )] =    47.4
BIND['%s-%s'            % (dbse, '232'                   )] =    48.2
BIND['%s-%s'            % (dbse, '233'                   )] =    71.2
BIND['%s-%s'            % (dbse, '234'                   )] =    8.9
BIND['%s-%s'            % (dbse, '235'                   )] =    17.7
BIND['%s-%s'            % (dbse, '236'                   )] =   -20.0 
BIND['%s-%s'            % (dbse, '237'                   )] =    31.5
BIND['%s-%s'            % (dbse, '238'                   )] =    39.7
BIND['%s-%s'            % (dbse, '239'                   )] =    51.3
BIND['%s-%s'            % (dbse, '240'                   )] =    59.0
BIND['%s-%s'            % (dbse, '241'                   )] =    44.3
BIND['%s-%s'            % (dbse, '242'                   )] =    12.8
BIND['%s-%s'            % (dbse, '243'                   )] =    1.7
BIND['%s-%s'            % (dbse, '244'                   )] =   -23.6 
BIND['%s-%s'            % (dbse, '245'                   )] =   -25.5
BIND['%s-%s'            % (dbse, '246'                   )] =    51.3
BIND['%s-%s'            % (dbse, '247'                   )] =    27.8
BIND['%s-%s'            % (dbse, '248'                   )] =   -109.0
BIND['%s-%s'            % (dbse, '249'                   )] =   -169.2
BIND['%s-%s'            % (dbse, '250'                   )] =   -21.8
BIND['%s-%s'            % (dbse, '251'                   )] =   -24.5
BIND['%s-%s'            % (dbse, '252'                   )] =   -1.9 
BIND['%s-%s'            % (dbse, '253'                   )] =    19.8
BIND['%s-%s'            % (dbse, '254'                   )] =   -17.4
BIND['%s-%s'            % (dbse, '255'                   )] =   -15.5
BIND['%s-%s'            % (dbse, '256'                   )] =   -3.6
BIND['%s-%s'            % (dbse, '257'                   )] =   -90.8
BIND['%s-%s'            % (dbse, '258'                   )] =   -85.6
BIND['%s-%s'            % (dbse, '259'                   )] =   -53.5
BIND['%s-%s'            % (dbse, '260'                   )] =    34.5
BIND['%s-%s'            % (dbse, '261'                   )] =    74.4
BIND['%s-%s'            % (dbse, '262'                   )] =    0.5
BIND['%s-%s'            % (dbse, '263'                   )] =   -6.8
BIND['%s-%s'            % (dbse, '264'                   )] =   -11.4
BIND['%s-%s'            % (dbse, '265'                   )] =   -10.9
BIND['%s-%s'            % (dbse, '266'                   )] =   -38.5
BIND['%s-%s'            % (dbse, '267'                   )] =   -52.2
BIND['%s-%s'            % (dbse, '268'                   )] =   -52.9
BIND['%s-%s'            % (dbse, '269'                   )] =   -41.8
BIND['%s-%s'            % (dbse, '270'                   )] =    21.7
BIND['%s-%s'            % (dbse, '271'                   )] =   -30.2
BIND['%s-%s'            % (dbse, '272'                   )] =   -6.7
BIND['%s-%s'            % (dbse, '273'                   )] =   -5.1
BIND['%s-%s'            % (dbse, '274'                   )] =   -33.0
BIND['%s-%s'            % (dbse, '275'                   )] =   -24.1
BIND['%s-%s'            % (dbse, '276'                   )] =    7.0
BIND['%s-%s'            % (dbse, '277'                   )] =    47.6
BIND['%s-%s'            % (dbse, '278'                   )] =   -49.2
BIND['%s-%s'            % (dbse, '279'                   )] =   -101.8
BIND['%s-%s'            % (dbse, '280'                   )] =   -105.4
BIND['%s-%s'            % (dbse, '281'                   )] =   -57.7
BIND['%s-%s'            % (dbse, '282'                   )] =   -27.7
BIND['%s-%s'            % (dbse, '283'                   )] =   -60.6
BIND['%s-%s'            % (dbse, '284'                   )] =   -48.3
BIND['%s-%s'            % (dbse, '285'                   )] =   -0.3
BIND['%s-%s'            % (dbse, '286'                   )] =   -4.2
BIND['%s-%s'            % (dbse, '287'                   )] =    32.8
BIND['%s-%s'            % (dbse, '288'                   )] =    32.1
BIND['%s-%s'            % (dbse, '289'                   )] =    39.8
BIND['%s-%s'            % (dbse, '290'                   )] =   -1.1
BIND['%s-%s'            % (dbse, '291'                   )] =    34.4
BIND['%s-%s'            % (dbse, '292'                   )] =    137.8
BIND['%s-%s'            % (dbse, '293'                   )] =    73.7
BIND['%s-%s'            % (dbse, '294'                   )] =   -1.3
BIND['%s-%s'            % (dbse, '295'                   )] =   -2.2
BIND['%s-%s'            % (dbse, '296'                   )] =    6.6
BIND['%s-%s'            % (dbse, '297'                   )] =    1.1
BIND['%s-%s'            % (dbse, '298'                   )] =    31.6
BIND['%s-%s'            % (dbse, '299'                   )] =    32.4
BIND['%s-%s'            % (dbse, '300'                   )] =    26.8
BIND['%s-%s'            % (dbse, '301'                   )] =    19.7
BIND['%s-%s'            % (dbse, '302'                   )] =    7.9
#G21PA = range(64, 72)    # G2-1 Proton Affinities (8)
#G21EF = range(72, 127)   # G2-1 Enthalpies of Formation (55)
#G22EF = range(210, 303)  # G2-2 Enthalpies of Formation (93)

# <<< Comment Lines >>>
TAGL = {}
# [  1 -  38] G2-1 Ionization Energies (38)
TAGL['%s-%s'            % (dbse, '1'                     )] = """Ionization Energy of Li (lithium atom)"""
TAGL['%s-%s'            % (dbse, '2'                     )] = """Ionization Energy of Be (berylliun atom)"""
TAGL['%s-%s'            % (dbse, '3'                     )] = """Ionization Energy of B (boron atom)"""
TAGL['%s-%s'            % (dbse, '4'                     )] = """Ionization Energy of C (carbon atom)"""
TAGL['%s-%s'            % (dbse, '5'                     )] = """Ionization Energy of N (nitrogen atom)"""
TAGL['%s-%s'            % (dbse, '6'                     )] = """Ionization Energy of O (oxygen atom)"""
TAGL['%s-%s'            % (dbse, '7'                     )] = """Ionization Energy of F (fluorine atom)"""
TAGL['%s-%s'            % (dbse, '8'                     )] = """Ionization Energy of Na (sodium atom)"""
TAGL['%s-%s'            % (dbse, '9'                     )] = """Ionization Energy of Mg (magnesium atom)"""
TAGL['%s-%s'            % (dbse, '10'                    )] = """Ionization Energy of Al (aluminium atom)"""
TAGL['%s-%s'            % (dbse, '11'                    )] = """Ionization Energy of Si (silicon atom)"""
TAGL['%s-%s'            % (dbse, '12'                    )] = """Ionization Energy of P (phosphorus atom)"""
TAGL['%s-%s'            % (dbse, '13'                    )] = """Ionization Energy of S (sulfur atom)"""
TAGL['%s-%s'            % (dbse, '14'                    )] = """Ionization Energy of Cl (chlorine atom)"""
TAGL['%s-%s'            % (dbse, '15'                    )] = """Ionization Energy of CH4 (methane)"""
TAGL['%s-%s'            % (dbse, '16'                    )] = """Ionization Energy of NH3 (ammonia)"""
TAGL['%s-%s'            % (dbse, '17'                    )] = """Ionization Energy of OH (hydroxyl radical)"""
TAGL['%s-%s'            % (dbse, '18'                    )] = """Ionization Energy of H2O (water)"""
TAGL['%s-%s'            % (dbse, '19'                    )] = """Ionization Energy of HF (hydrogen fluoride)"""
TAGL['%s-%s'            % (dbse, '20'                    )] = """Ionization Energy of SiH4 (silane)"""
TAGL['%s-%s'            % (dbse, '21'                    )] = """Ionization Energy of PH (phosphorus monohydride)"""
TAGL['%s-%s'            % (dbse, '22'                    )] = """Ionization Energy of PH2 (phosphino radical)"""
TAGL['%s-%s'            % (dbse, '23'                    )] = """Ionization Energy of PH3 (phoshine)"""
TAGL['%s-%s'            % (dbse, '24'                    )] = """Ionization Energy of SH (mercapato radical)"""
TAGL['%s-%s'            % (dbse, '25'                    )] = """Ionization Energy of H2S (hydrogen sulfide with 2B1 cation)"""
TAGL['%s-%s'            % (dbse, '26'                    )] = """Ionization Energy of H2S (hydrogen sulfide with 2A1 cation)"""
TAGL['%s-%s'            % (dbse, '27'                    )] = """Ionization Energy of HCl (hydrogen chloride)"""
TAGL['%s-%s'            % (dbse, '28'                    )] = """Ionization Energy of HCCH (ethyne)"""
TAGL['%s-%s'            % (dbse, '29'                    )] = """Ionization Energy of H2C=CH2 (ethene)"""
TAGL['%s-%s'            % (dbse, '30'                    )] = """Ionization Energy of CO (carbon monoxide)"""
TAGL['%s-%s'            % (dbse, '31'                    )] = """Ionization Energy of N2 (nitrogen diatomic with 2Sigma cation)"""
TAGL['%s-%s'            % (dbse, '32'                    )] = """Ionization Energy of N2 (nitrogen diatomic with 2Pi cation)"""
TAGL['%s-%s'            % (dbse, '33'                    )] = """Ionization Energy of O2 (oxygen diatomic)"""
TAGL['%s-%s'            % (dbse, '34'                    )] = """Ionization Energy of P2 (phosphorus diatomic)"""
TAGL['%s-%s'            % (dbse, '35'                    )] = """Ionization Energy of S2 (sulfur diatomic)"""
TAGL['%s-%s'            % (dbse, '36'                    )] = """Ionization Energy of Cl2 (chlorine diatomic)"""
TAGL['%s-%s'            % (dbse, '37'                    )] = """Ionization Energy of ClF (chlorine monofluoride)"""
TAGL['%s-%s'            % (dbse, '38'                    )] = """Ionization Energy of SC (carbon monosulfide)"""
# [ 39 -  63] G2-1 Electron Affinities (25)
TAGL['%s-%s'            % (dbse, '39'                    )] = """Electron Affinity of C (carbon atom)"""
TAGL['%s-%s'            % (dbse, '40'                    )] = """Electron Affinity of O (oxygen atom)"""
TAGL['%s-%s'            % (dbse, '41'                    )] = """Electron Affinity of F (fluorine atom)"""
TAGL['%s-%s'            % (dbse, '42'                    )] = """Electron Affinity of Si (silicon atom)"""
TAGL['%s-%s'            % (dbse, '43'                    )] = """Electron Affinity of P (phosphorus atom)"""
TAGL['%s-%s'            % (dbse, '44'                    )] = """Electron Affinity of S (sulfur atom)"""
TAGL['%s-%s'            % (dbse, '45'                    )] = """Electron Affinity of Cl (chlorine atom)"""
TAGL['%s-%s'            % (dbse, '46'                    )] = """Electron Affinity of CH (methylidyne)"""
TAGL['%s-%s'            % (dbse, '47'                    )] = """Electron Affinity of CH2 (methylene)"""
TAGL['%s-%s'            % (dbse, '48'                    )] = """Electron Affinity of CH3 (methyl radical)"""
TAGL['%s-%s'            % (dbse, '49'                    )] = """Electron Affinity of NH (imidogen)"""
TAGL['%s-%s'            % (dbse, '50'                    )] = """Electron Affinity of NH2 (amino radical)"""
TAGL['%s-%s'            % (dbse, '51'                    )] = """Electron Affinity of OH (hydroxyl radical)"""
TAGL['%s-%s'            % (dbse, '52'                    )] = """Electron Affinity of SiH (silylidyne)"""
TAGL['%s-%s'            % (dbse, '53'                    )] = """Electron Affinity of SiH2 (silicon dihydride)"""
TAGL['%s-%s'            % (dbse, '54'                    )] = """Electron Affinity of SiH3 (silyl)"""
TAGL['%s-%s'            % (dbse, '55'                    )] = """Electron Affinity of PH (phosphorus monohydride)"""
TAGL['%s-%s'            % (dbse, '56'                    )] = """Electron Affinity of PH2 (phosphino radical)"""
TAGL['%s-%s'            % (dbse, '57'                    )] = """Electron Affinity of SH (mercapato radical)"""
TAGL['%s-%s'            % (dbse, '58'                    )] = """Electron Affinity of O2 (oxygen diatomic)"""
TAGL['%s-%s'            % (dbse, '59'                    )] = """Electron Affinity of NO (nitric oxide)"""
TAGL['%s-%s'            % (dbse, '60'                    )] = """Electron Affinity of CN (cyano radical)"""
TAGL['%s-%s'            % (dbse, '61'                    )] = """Electron Affinity of PO (phosphorus monoxide)"""
TAGL['%s-%s'            % (dbse, '62'                    )] = """Electron Affinity of S2 (sulfur diatomic)"""
TAGL['%s-%s'            % (dbse, '63'                    )] = """Electron Affinity of Cl2 (chlorine diatomic)"""
# [ 64 -  71] G2-1 Proton Affinities (8)
TAGL['%s-%s'            % (dbse, '64'                    )] = """Proton Affinity of NH3 (ammonia)"""
TAGL['%s-%s'            % (dbse, '65'                    )] = """Proton Affinity of H2O (water)"""
TAGL['%s-%s'            % (dbse, '66'                    )] = """Proton Affinity of HCCH (ethyne)"""
TAGL['%s-%s'            % (dbse, '67'                    )] = """Proton Affinity of SiH4 (silane)"""
TAGL['%s-%s'            % (dbse, '68'                    )] = """Proton Affinity of PH3 (phosphine)"""
TAGL['%s-%s'            % (dbse, '69'                    )] = """Proton Affinity of H2S (hydrogen sulfide)"""
TAGL['%s-%s'            % (dbse, '70'                    )] = """Proton Affinity of HCl (hydrogen chloride)"""
TAGL['%s-%s'            % (dbse, '71'                    )] = """Proton Affinity of H2 (hydrogen diatomic)"""
# [ 72 - 126] G2-1 Enthalpies of Formation (55)
TAGL['%s-%s'            % (dbse, '72'                    )] = """Atomization Energy of LiH (lithium hydride)"""
TAGL['%s-%s'            % (dbse, '73'                    )] = """Atomization Energy of BeH (beryllium monohydride)"""
TAGL['%s-%s'            % (dbse, '74'                    )] = """Atomization Energy of CH (methylidyne)"""
TAGL['%s-%s'            % (dbse, '75'                    )] = """Atomization Energy of CH2 (triplet methylene)"""
TAGL['%s-%s'            % (dbse, '76'                    )] = """Atomization Energy of CH2 (singlet methylene)"""
TAGL['%s-%s'            % (dbse, '77'                    )] = """Atomization Energy of CH3 (methyl radical)"""
TAGL['%s-%s'            % (dbse, '78'                    )] = """Atomization Energy of CH4 (methane)"""
TAGL['%s-%s'            % (dbse, '79'                    )] = """Atomization Energy of NH (imidogen)"""
TAGL['%s-%s'            % (dbse, '80'                    )] = """Atomization Energy of NH2 (amino radical)"""
TAGL['%s-%s'            % (dbse, '81'                    )] = """Atomization Energy of NH3 (ammonia)"""
TAGL['%s-%s'            % (dbse, '82'                    )] = """Atomization Energy of OH (hydroxyl radical)"""
TAGL['%s-%s'            % (dbse, '83'                    )] = """Atomization Energy of H2O (water)"""
TAGL['%s-%s'            % (dbse, '84'                    )] = """Atomization Energy of HF (hydrogen fluoride)"""
TAGL['%s-%s'            % (dbse, '85'                    )] = """Atomization Energy of SiH2 (singlet silicon dihydride)"""
TAGL['%s-%s'            % (dbse, '86'                    )] = """Atomization Energy of SiH2 (triplet silicon dihydride)"""
TAGL['%s-%s'            % (dbse, '87'                    )] = """Atomization Energy of SiH3 (silyl)"""
TAGL['%s-%s'            % (dbse, '88'                    )] = """Atomization Energy of SiH4 (silane)"""
TAGL['%s-%s'            % (dbse, '89'                    )] = """Atomization Energy of PH2 (phosphino radical)"""
TAGL['%s-%s'            % (dbse, '90'                    )] = """Atomization Energy of PH3 (phosphine)"""
TAGL['%s-%s'            % (dbse, '91'                    )] = """Atomization Energy of H2S (hydrogen sulfide)"""
TAGL['%s-%s'            % (dbse, '92'                    )] = """Atomization Energy of HCl (hydrogen chloride)"""
TAGL['%s-%s'            % (dbse, '93'                    )] = """Atomization Energy of Li2 (lithium diatonic)"""
TAGL['%s-%s'            % (dbse, '94'                    )] = """Atomization Energy of LiF (lithium fluoride)"""
TAGL['%s-%s'            % (dbse, '95'                    )] = """Atomization Energy of HCCH (ethyne)"""
TAGL['%s-%s'            % (dbse, '96'                    )] = """Atomization Energy of H2C=CH2 (ethene)"""
TAGL['%s-%s'            % (dbse, '97'                    )] = """Atomization Energy of H3C-CH3 (ethane)"""
TAGL['%s-%s'            % (dbse, '98'                    )] = """Atomization Energy of CN (cyano radical)"""
TAGL['%s-%s'            % (dbse, '99'                    )] = """Atomization Energy of HCN (hydrogen cyanide)"""
TAGL['%s-%s'            % (dbse, '100'                   )] = """Atomization Energy of CO (carbon monoxide)"""
TAGL['%s-%s'            % (dbse, '101'                   )] = """Atomization Energy of HCO (formyl radical)"""
TAGL['%s-%s'            % (dbse, '102'                   )] = """Atomization Energy of H2C=O (formaldehyde)"""
TAGL['%s-%s'            % (dbse, '103'                   )] = """Atomization Energy of H3COH (methanol)"""
TAGL['%s-%s'            % (dbse, '104'                   )] = """Atomization Energy of N2 (nitrogen diatomic)"""
TAGL['%s-%s'            % (dbse, '105'                   )] = """Atomization Energy of H2N-NH2 (hydrazine)"""
TAGL['%s-%s'            % (dbse, '106'                   )] = """Atomization Energy of NO (nitric oxide)"""
TAGL['%s-%s'            % (dbse, '107'                   )] = """Atomization Energy of O2 (oxygen diatomic)"""
TAGL['%s-%s'            % (dbse, '108'                   )] = """Atomization Energy of HOOH (hydrogen peroxide)"""
TAGL['%s-%s'            % (dbse, '109'                   )] = """Atomization Energy of F2 (fluorine diatomic)"""
TAGL['%s-%s'            % (dbse, '110'                   )] = """Atomization Energy of CO2 (carbon dioxide)"""
TAGL['%s-%s'            % (dbse, '111'                   )] = """Atomization Energy of Na2 (sodium diatomic)"""
TAGL['%s-%s'            % (dbse, '112'                   )] = """Atomization Energy of Si2 (silicon diatomic)"""
TAGL['%s-%s'            % (dbse, '113'                   )] = """Atomization Energy of P2 (phosphorus diatomic)"""
TAGL['%s-%s'            % (dbse, '114'                   )] = """Atomization Energy of S2 (sulfur diatomic)"""
TAGL['%s-%s'            % (dbse, '115'                   )] = """Atomization Energy of Cl2 (chlorine diatomic)"""
TAGL['%s-%s'            % (dbse, '116'                   )] = """Atomization Energy of NaCl (sodium chloride)"""
TAGL['%s-%s'            % (dbse, '117'                   )] = """Atomization Energy of SiO (silicon monoxide)"""
TAGL['%s-%s'            % (dbse, '118'                   )] = """Atomization Energy of SC (carbon monosulfide)"""
TAGL['%s-%s'            % (dbse, '119'                   )] = """Atomization Energy of SO (sulfur monoxide)"""
TAGL['%s-%s'            % (dbse, '120'                   )] = """Atomization Energy of ClO (monochlorine monoxide)"""
TAGL['%s-%s'            % (dbse, '121'                   )] = """Atomization Energy of FCl (chlorine monofluoride)"""
TAGL['%s-%s'            % (dbse, '122'                   )] = """Atomization Energy of H3Si-SiH3 (disilane)"""
TAGL['%s-%s'            % (dbse, '123'                   )] = """Atomization Energy of H3CCl (methyl chloride)"""
TAGL['%s-%s'            % (dbse, '124'                   )] = """Atomization Energy of H3CSH (methanethiol)"""
TAGL['%s-%s'            % (dbse, '125'                   )] = """Atomization Energy of HOCl (hypochlorous acid)"""
TAGL['%s-%s'            % (dbse, '126'                   )] = """Atomization Energy of SO2 (sulfur dioxide)"""
# [127 - 176] G2-2 Ionization Energies (50)
TAGL['%s-%s'            % (dbse, '127'                   )] = """Ionization Energy of H (hydrogen atom)"""
TAGL['%s-%s'            % (dbse, '128'                   )] = """Ionization Energy of He (helium atom)"""
TAGL['%s-%s'            % (dbse, '129'                   )] = """Ionization Energy of Ne (neon atom)"""
TAGL['%s-%s'            % (dbse, '130'                   )] = """Ionization Energy of Ar (argon atom)"""
TAGL['%s-%s'            % (dbse, '131'                   )] = """Ionization Energy of BF3 (trifluoroborane)"""
TAGL['%s-%s'            % (dbse, '132'                   )] = """Ionization Energy of BCl3 (trichloroborane)"""
TAGL['%s-%s'            % (dbse, '133'                   )] = """Ionization Energy of B2F4 (diboron tetrafluoride)"""
TAGL['%s-%s'            % (dbse, '134'                   )] = """Ionization Energy of CO2 (carbon dioxide)"""
TAGL['%s-%s'            % (dbse, '135'                   )] = """Ionization Energy of CF2 (difluoromethylene)"""
TAGL['%s-%s'            % (dbse, '136'                   )] = """Ionization Energy of OCS (carbonyl sulfide)"""
TAGL['%s-%s'            % (dbse, '137'                   )] = """Ionization Energy of CS2 (carbon disulfide)"""
TAGL['%s-%s'            % (dbse, '138'                   )] = """Ionization Energy of CH2 (methylene)"""
TAGL['%s-%s'            % (dbse, '139'                   )] = """Ionization Energy of CH3 (methyl radical)"""
TAGL['%s-%s'            % (dbse, '140'                   )] = """Ionization Energy of C2H5 (ethyl radical)"""
TAGL['%s-%s'            % (dbse, '141'                   )] = """Ionization Energy of C3H4 (cyclopropene)"""
TAGL['%s-%s'            % (dbse, '142'                   )] = """Ionization Energy of H2C=C=CH2 (allene)"""
TAGL['%s-%s'            % (dbse, '143'                   )] = """Ionization Energy of C3H7 (isopropyl radical)"""
TAGL['%s-%s'            % (dbse, '144'                   )] = """Ionization Energy of C6H6 (benzene)"""
TAGL['%s-%s'            % (dbse, '145'                   )] = """Ionization Energy of C6H5-CH3 (toluene)"""
TAGL['%s-%s'            % (dbse, '146'                   )] = """Ionization Energy of CN (cyano radical)"""
TAGL['%s-%s'            % (dbse, '147'                   )] = """Ionization Energy of HCO (formyl radical)"""
TAGL['%s-%s'            % (dbse, '148'                   )] = """Ionization Energy of CH2OH (hydroxymethyl radical)"""
TAGL['%s-%s'            % (dbse, '149'                   )] = """Ionization Energy of CH3O (methoxy radical)"""
TAGL['%s-%s'            % (dbse, '150'                   )] = """Ionization Energy of H3COH (methanol)"""
TAGL['%s-%s'            % (dbse, '151'                   )] = """Ionization Energy of H3CF (methyl fluoride)"""
TAGL['%s-%s'            % (dbse, '152'                   )] = """Ionization Energy of H2C=S (thioformaldehyde)"""
TAGL['%s-%s'            % (dbse, '153'                   )] = """Ionization Energy of CH2SH (mercapatomethyl radical*)"""
TAGL['%s-%s'            % (dbse, '154'                   )] = """Ionization Energy of H3CSH (methanethiol)"""
TAGL['%s-%s'            % (dbse, '155'                   )] = """Ionization Energy of H3CCl (chloromethane)"""
TAGL['%s-%s'            % (dbse, '156'                   )] = """Ionization Energy of C2H2OH (????)"""
TAGL['%s-%s'            % (dbse, '157'                   )] = """Ionization Energy of CH3CHO (acetaldehyde)"""
TAGL['%s-%s'            % (dbse, '158'                   )] = """Ionization Energy of CH3OF (methyl hypofluorite)"""
TAGL['%s-%s'            % (dbse, '159'                   )] = """Ionization Energy of C2H4S (thiirane)"""
TAGL['%s-%s'            % (dbse, '160'                   )] = """Ionization Energy of NCCN (cyanogen)"""
TAGL['%s-%s'            % (dbse, '161'                   )] = """Ionization Energy of C4H4O (furan)"""
TAGL['%s-%s'            % (dbse, '162'                   )] = """Ionization Energy of C4H5N (pyrrole)"""
TAGL['%s-%s'            % (dbse, '163'                   )] = """Ionization Energy of C6H5-OH (phenol)"""
TAGL['%s-%s'            % (dbse, '164'                   )] = """Ionization Energy of C6H5-NH2 (aniline)"""
TAGL['%s-%s'            % (dbse, '165'                   )] = """Ionization Energy of B2H4 (diborane)"""
TAGL['%s-%s'            % (dbse, '166'                   )] = """Ionization Energy of NH (imidogen)"""
TAGL['%s-%s'            % (dbse, '167'                   )] = """Ionization Energy of NH2 (amino radical)"""
TAGL['%s-%s'            % (dbse, '168'                   )] = """Ionization Energy of N2H2 (trans diazine)"""
TAGL['%s-%s'            % (dbse, '169'                   )] = """Ionization Energy of N2H3 (hydrazinyl radical)"""
TAGL['%s-%s'            % (dbse, '170'                   )] = """Ionization Energy of HOF (hypofluorous acid)"""
TAGL['%s-%s'            % (dbse, '171'                   )] = """Ionization Energy of SiH2 (silicon dihydride)"""
TAGL['%s-%s'            % (dbse, '172'                   )] = """Ionization Energy of SiH3 (silyl)"""
TAGL['%s-%s'            % (dbse, '173'                   )] = """Ionization Energy of Si2H2 (disilyne)"""
TAGL['%s-%s'            % (dbse, '174'                   )] = """Ionization Energy of Si2H4 (disilene)"""
TAGL['%s-%s'            % (dbse, '175'                   )] = """Ionization Energy of Si2H5 (disilanyl radical)"""
TAGL['%s-%s'            % (dbse, '176'                   )] = """Ionization Energy of Si2H6 (disilane)"""
# [177 - 209] G2-2 Electron Affinities (33)
TAGL['%s-%s'            % (dbse, '177'                   )] = """Electron Affinity Li (lithium atom)"""
TAGL['%s-%s'            % (dbse, '178'                   )] = """Electron Affinity B (boron atom)"""
TAGL['%s-%s'            % (dbse, '179'                   )] = """Electron Affinity Na (sodium atom)"""
TAGL['%s-%s'            % (dbse, '180'                   )] = """Electron Affinity Al (aluminium atom)"""
TAGL['%s-%s'            % (dbse, '181'                   )] = """Electron Affinity C2 (dicarbon)"""
TAGL['%s-%s'            % (dbse, '182'                   )] = """Electron Affinity C2O (dicarbon monoxide)"""
TAGL['%s-%s'            % (dbse, '183'                   )] = """Electron Affinity CF2 (difluoromethylene)"""
TAGL['%s-%s'            % (dbse, '184'                   )] = """Electron Affinity NCO (isocyanato radical)"""
TAGL['%s-%s'            % (dbse, '185'                   )] = """Electron Affinity NO2 (nitrogen dioxide)"""
TAGL['%s-%s'            % (dbse, '186'                   )] = """Electron Affinity O3 (ozone)"""
TAGL['%s-%s'            % (dbse, '187'                   )] = """Electron Affinity OF (oxygen monofluoride)"""
TAGL['%s-%s'            % (dbse, '188'                   )] = """Electron Affinity SO2 (sulfur dioxide)"""
TAGL['%s-%s'            % (dbse, '189'                   )] = """Electron Affinity S2O (disulfur monoxide)"""
TAGL['%s-%s'            % (dbse, '190'                   )] = """Electron Affinity CCH (ethynyl radical)"""
TAGL['%s-%s'            % (dbse, '191'                   )] = """Electron Affinity H2CCH (vinyl radical)"""
TAGL['%s-%s'            % (dbse, '192'                   )] = """Electron Affinity H2C=C=C (propadienylidene)"""
TAGL['%s-%s'            % (dbse, '193'                   )] = """Electron Affinity H2C=C=CH (????)"""
TAGL['%s-%s'            % (dbse, '194'                   )] = """Electron Affinity H2C=CH-CH2 (allyl radical)"""
TAGL['%s-%s'            % (dbse, '195'                   )] = """Electron Affinity HCO (formyl radical)"""
TAGL['%s-%s'            % (dbse, '196'                   )] = """Electron Affinity HCF (fluoromethylene)"""
TAGL['%s-%s'            % (dbse, '197'                   )] = """Electron Affinity CH3O (methoxy radical)"""
TAGL['%s-%s'            % (dbse, '198'                   )] = """Electron Affinity CH3S (thiomethoxy radical)"""
TAGL['%s-%s'            % (dbse, '199'                   )] = """Electron Affinity H2C=S (thioformaldehyde)"""
TAGL['%s-%s'            % (dbse, '200'                   )] = """Electron Affinity H2CCN (cyanomethyl radical)"""
TAGL['%s-%s'            % (dbse, '201'                   )] = """Electron Affinity H2CNC (????)"""
TAGL['%s-%s'            % (dbse, '202'                   )] = """Electron Affinity HC=C=O (ketenyl radical)"""
TAGL['%s-%s'            % (dbse, '203'                   )] = """Electron Affinity H2CCHO (????)"""
TAGL['%s-%s'            % (dbse, '204'                   )] = """Electron Affinity H3CCO (acetyl radical)"""
TAGL['%s-%s'            % (dbse, '205'                   )] = """Electron Affinity H3C-CH2O (ethoxy radical)"""
TAGL['%s-%s'            % (dbse, '206'                   )] = """Electron Affinity H3C-CH2S (methylthiomethyl radical)"""
TAGL['%s-%s'            % (dbse, '207'                   )] = """Electron Affinity LiH (lithium hydride)"""
TAGL['%s-%s'            % (dbse, '208'                   )] = """Electron Affinity HNO (nitrosyl hydride)"""
TAGL['%s-%s'            % (dbse, '209'                   )] = """Electron Affinity HOO (hydroperoxy radical)"""
# [210 - 302] G2-2 Enthalpies of Formation (93)
TAGL['%s-%s'            % (dbse, '210'                   )] = """Atomization Energy of BF3 (trifluoroborane)"""
TAGL['%s-%s'            % (dbse, '211'                   )] = """Atomization Energy of BCl3 (trichloroborane)"""
TAGL['%s-%s'            % (dbse, '212'                   )] = """Atomization Energy of AlF3 (aluminum trifluoride)"""
TAGL['%s-%s'            % (dbse, '213'                   )] = """Atomization Energy of AlCl3 (aluminum trichloride)"""
TAGL['%s-%s'            % (dbse, '214'                   )] = """Atomization Energy of CF4 (carbon tetrafluoride)"""
TAGL['%s-%s'            % (dbse, '215'                   )] = """Atomization Energy of CCl4 (carbon tetrachloride)"""
TAGL['%s-%s'            % (dbse, '216'                   )] = """Atomization Energy of O=C=S (carbonyl sulfide)"""
TAGL['%s-%s'            % (dbse, '217'                   )] = """Atomization Energy of CS2 (carbon disulfide)"""
TAGL['%s-%s'            % (dbse, '218'                   )] = """Atomization Energy of CF2O (carbonic difluoride)"""
TAGL['%s-%s'            % (dbse, '219'                   )] = """Atomization Energy of SiF4 (tetrafluorosilane)"""
TAGL['%s-%s'            % (dbse, '220'                   )] = """Atomization Energy of SiCl4 (tetrachlorosilane)"""
TAGL['%s-%s'            % (dbse, '221'                   )] = """Atomization Energy of N2O (nitrous oxide)"""
TAGL['%s-%s'            % (dbse, '222'                   )] = """Atomization Energy of ClNO (nitrosyl chloride)"""
TAGL['%s-%s'            % (dbse, '223'                   )] = """Atomization Energy of NF3 (nitrogen trifluoride)"""
TAGL['%s-%s'            % (dbse, '224'                   )] = """Atomization Energy of PF3 (phosphorus trifluoride)"""
TAGL['%s-%s'            % (dbse, '225'                   )] = """Atomization Energy of O3 (ozone)"""
TAGL['%s-%s'            % (dbse, '226'                   )] = """Atomization Energy of F2O (difluorine monoxide)"""
TAGL['%s-%s'            % (dbse, '227'                   )] = """Atomization Energy of ClF3 (chlorine trifluoride)"""
TAGL['%s-%s'            % (dbse, '228'                   )] = """Atomization Energy of F2C=CF2 (tetrafluoroethene)"""
TAGL['%s-%s'            % (dbse, '229'                   )] = """Atomization Energy of Cl2C=CCl2 (tetrachloroethene)"""
TAGL['%s-%s'            % (dbse, '230'                   )] = """Atomization Energy of F3CCN (trifluoroacetonitrile)"""
TAGL['%s-%s'            % (dbse, '231'                   )] = """Atomization Energy of H3CCCH (propyne)"""
TAGL['%s-%s'            % (dbse, '232'                   )] = """Atomization Energy of H2C=C=CH2 (allene)"""
TAGL['%s-%s'            % (dbse, '233'                   )] = """Atomization Energy of C3H4 (cyclopropene)"""
TAGL['%s-%s'            % (dbse, '234'                   )] = """Atomization Energy of H3C-CH=CH2 (propene)"""
TAGL['%s-%s'            % (dbse, '235'                   )] = """Atomization Energy of C3H6 (cyclopropane)"""
TAGL['%s-%s'            % (dbse, '236'                   )] = """Atomization Energy of C3H8 (propane)"""
TAGL['%s-%s'            % (dbse, '237'                   )] = """Atomization Energy of H2C=CH-CH=CH2 (butadiene)"""
TAGL['%s-%s'            % (dbse, '238'                   )] = """Atomization Energy of C4H6 (2-butyne)"""
TAGL['%s-%s'            % (dbse, '239'                   )] = """Atomization Energy of C4H6 (methylene cyclopropane)"""
TAGL['%s-%s'            % (dbse, '240'                   )] = """Atomization Energy of C4H6 (bicyclobutane)"""
TAGL['%s-%s'            % (dbse, '241'                   )] = """Atomization Energy of C4H6 (cyclobutene)"""
TAGL['%s-%s'            % (dbse, '242'                   )] = """Atomization Energy of C4H8 (cyclobutane)"""
TAGL['%s-%s'            % (dbse, '243'                   )] = """Atomization Energy of C4H8 (isobutene)"""
TAGL['%s-%s'            % (dbse, '244'                   )] = """Atomization Energy of C4H10 (trans butane)"""
TAGL['%s-%s'            % (dbse, '245'                   )] = """Atomization Energy of C4H10 (isobutane)"""
TAGL['%s-%s'            % (dbse, '246'                   )] = """Atomization Energy of C5H8 (spiropentane)"""
TAGL['%s-%s'            % (dbse, '247'                   )] = """Atomization Energy of C6H6 (benzene)"""
TAGL['%s-%s'            % (dbse, '248'                   )] = """Atomization Energy of H2CF2 (methylene fluoride)"""
TAGL['%s-%s'            % (dbse, '249'                   )] = """Atomization Energy of HCF3 (trifluoromethane)"""
TAGL['%s-%s'            % (dbse, '250'                   )] = """Atomization Energy of H2CCl2 (methylene chloride)"""
TAGL['%s-%s'            % (dbse, '251'                   )] = """Atomization Energy of HCCl3 (chloroform)"""
TAGL['%s-%s'            % (dbse, '252'                   )] = """Atomization Energy of H3CNH2 (methylamine)"""
TAGL['%s-%s'            % (dbse, '253'                   )] = """Atomization Energy of H3CCN (methyl cyanide)"""
TAGL['%s-%s'            % (dbse, '254'                   )] = """Atomization Energy of H3CNO2 (nitromethane)"""
TAGL['%s-%s'            % (dbse, '255'                   )] = """Atomization Energy of H3CONO (methyl nitrite)"""
TAGL['%s-%s'            % (dbse, '256'                   )] = """Atomization Energy of H3CSiH3 (methyl silane)"""
TAGL['%s-%s'            % (dbse, '257'                   )] = """Atomization Energy of HCOOH (formic acid)"""
TAGL['%s-%s'            % (dbse, '258'                   )] = """Atomization Energy of H3C-O-CHO (methyl formate)"""
TAGL['%s-%s'            % (dbse, '259'                   )] = """Atomization Energy of H3C-CONH2 (acetamide)"""
TAGL['%s-%s'            % (dbse, '260'                   )] = """Atomization Energy of C2H4NH (aziridine)"""
TAGL['%s-%s'            % (dbse, '261'                   )] = """Atomization Energy of NCCN (cyanogen)"""
TAGL['%s-%s'            % (dbse, '262'                   )] = """Atomization Energy of (H3C)2NH (dimethylamine)"""
TAGL['%s-%s'            % (dbse, '263'                   )] = """Atomization Energy of H3C-CH2NH2 (trans ethylamine)"""
TAGL['%s-%s'            % (dbse, '264'                   )] = """Atomization Energy of H2C=C=O (ketene)"""
TAGL['%s-%s'            % (dbse, '265'                   )] = """Atomization Energy of C2H4O (oxirane)"""
TAGL['%s-%s'            % (dbse, '266'                   )] = """Atomization Energy of H3C-CHO (acetaldehyde)"""
TAGL['%s-%s'            % (dbse, '267'                   )] = """Atomization Energy of OHC-CHO (glyoxal)"""
TAGL['%s-%s'            % (dbse, '268'                   )] = """Atomization Energy of H3C-CH2OH (ethanol)"""
TAGL['%s-%s'            % (dbse, '269'                   )] = """Atomization Energy of H3C-O-CH3 (dimethylether)"""
TAGL['%s-%s'            % (dbse, '270'                   )] = """Atomization Energy of C2H4S (thiooxirane)"""
TAGL['%s-%s'            % (dbse, '271'                   )] = """Atomization Energy of (H3C)2SO (dimethyl sulfoxide)"""
TAGL['%s-%s'            % (dbse, '272'                   )] = """Atomization Energy of H3C-CH2SH (ethanethiol)"""
TAGL['%s-%s'            % (dbse, '273'                   )] = """Atomization Energy of H3C-S-CH3 (dimethyl sulfide)"""
TAGL['%s-%s'            % (dbse, '274'                   )] = """Atomization Energy of H2C=CHF (vinyl fluoride)"""
TAGL['%s-%s'            % (dbse, '275'                   )] = """Atomization Energy of H3C-CH2Cl (ethyl chloride)"""
TAGL['%s-%s'            % (dbse, '276'                   )] = """Atomization Energy of H2C-CHCl (vinyl chloride)"""
TAGL['%s-%s'            % (dbse, '277'                   )] = """Atomization Energy of H2C=CHCN (acrylonitrile)"""
TAGL['%s-%s'            % (dbse, '278'                   )] = """Atomization Energy of H3C-COCH3 (acetone)"""
TAGL['%s-%s'            % (dbse, '279'                   )] = """Atomization Energy of H3C-COOH (acetic acid)"""
TAGL['%s-%s'            % (dbse, '280'                   )] = """Atomization Energy of H3C-COF (acetyl fluoride)"""
TAGL['%s-%s'            % (dbse, '281'                   )] = """Atomization Energy of H3C-COCl (acetyl chloride)"""
TAGL['%s-%s'            % (dbse, '282'                   )] = """Atomization Energy of H3C-CH2-CH2Cl (propyl chloride)"""
TAGL['%s-%s'            % (dbse, '283'                   )] = """Atomization Energy of (H3C)2CHOH (isopropanol)"""
TAGL['%s-%s'            % (dbse, '284'                   )] = """Atomization Energy of C2H5-O-CH3 (methyl ethyl ether)"""
TAGL['%s-%s'            % (dbse, '285'                   )] = """Atomization Energy of (H3C)3N (trimethylamine)"""
TAGL['%s-%s'            % (dbse, '286'                   )] = """Atomization Energy of C4H4O (furan)"""
TAGL['%s-%s'            % (dbse, '287'                   )] = """Atomization Energy of C4H4S (thiophene)"""
TAGL['%s-%s'            % (dbse, '288'                   )] = """Atomization Energy of C4H5N (pyrrole)"""
TAGL['%s-%s'            % (dbse, '289'                   )] = """Atomization Energy of C5H5N (pyridine)"""
TAGL['%s-%s'            % (dbse, '290'                   )] = """Atomization Energy of H2 (hydrogen diatomic)"""
TAGL['%s-%s'            % (dbse, '291'                   )] = """Atomization Energy of SH (mercapato radical)"""
TAGL['%s-%s'            % (dbse, '292'                   )] = """Atomization Energy of CCH (ethylnyl radical)"""
TAGL['%s-%s'            % (dbse, '293'                   )] = """Atomization Energy of C2H3 (vinyl radical)"""
TAGL['%s-%s'            % (dbse, '294'                   )] = """Atomization Energy of H3CCO (acetyl radical)"""
TAGL['%s-%s'            % (dbse, '295'                   )] = """Atomization Energy of H2COH (hydroxymethyl radical)"""
TAGL['%s-%s'            % (dbse, '296'                   )] = """Atomization Energy of CH3O (methoxy radical)"""
TAGL['%s-%s'            % (dbse, '297'                   )] = """Atomization Energy of H3C-CH2O (ethoxy radical)"""
TAGL['%s-%s'            % (dbse, '298'                   )] = """Atomization Energy of CH3S (thiomthoxy radical)"""
TAGL['%s-%s'            % (dbse, '299'                   )] = """Atomization Energy of C2H5 (ethyl radical)"""
TAGL['%s-%s'            % (dbse, '300'                   )] = """Atomization Energy of (H3C)2CH (isopropyl radical)"""
TAGL['%s-%s'            % (dbse, '301'                   )] = """Atomization Energy of (H3C)3C (t-butyl radical)"""
TAGL['%s-%s'            % (dbse, '302'                   )] = """Atomization Energy of NO2 (nitrogen dioxide)"""
TAGL['%s-%s-reagent'    % (dbse, '2Butyne'               )] = """2-Butyne (C4H6) eclipsed, D3h """
TAGL['%s-%s-reagent'    % (dbse, 'Acetaldehyde'          )] = """Acetaldehyde (CH3CHO), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Acetaldehyde_cation'   )] = """Acetaldehyde Radical Cation (CH3CHO+) """
TAGL['%s-%s-reagent'    % (dbse, 'Acetamide'             )] = """Acetamide (CH3-CONH2), C1 """
TAGL['%s-%s-reagent'    % (dbse, 'AceticAcid'            )] = """Acetic Acid (CH3-COOH) single bonds trans, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Acetone'               )] = """Acetone (CH3-CO-CH3), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Acetylchloride'        )] = """Acetyl Chloride (H3C-COCl) HCCO cis, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Acetylfluoride'        )] = """Acetyl Fluoride (H3C-COF) HCCO cis, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Al_anion'              )] = """Aluminum Anion (Al-) """
TAGL['%s-%s-reagent'    % (dbse, 'Al_cation'             )] = """Aluminum Atom Cation (Al+) """
TAGL['%s-%s-reagent'    % (dbse, 'Al_radical'            )] = """Aluminum Atom Radical (Al) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'AlCl3'                 )] = """Aluminum Chloride (AlCl3) planar, D3h """
TAGL['%s-%s-reagent'    % (dbse, 'AlF3'                  )] = """Aluminum Fluoride (AlF3) planar, D3h """
TAGL['%s-%s-reagent'    % (dbse, 'Allene'                )] = """Allene (C3H4), D2d """
TAGL['%s-%s-reagent'    % (dbse, 'Allene_cation'         )] = """Allene Radical Cation (C3H4+) twisted """
TAGL['%s-%s-reagent'    % (dbse, 'Aniline'               )] = """Aniline (C6H5-NH2), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Aniline_cation'        )] = """Aniline Radical Cation (C6H5-NH2+), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Ar'                    )] = """Argon Atom (Ar) """
TAGL['%s-%s-reagent'    % (dbse, 'Ar_cation'             )] = """Argon Radical Cation (Ar+) """
TAGL['%s-%s-reagent'    % (dbse, 'Aziridine'             )] = """Aziridine (cyclic CH2-NH-CH2), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'B2F4'                  )] = """B2F4, D2d """
TAGL['%s-%s-reagent'    % (dbse, 'B2F4_cation'           )] = """B2F4+ Radical Cation, D2d """
TAGL['%s-%s-reagent'    % (dbse, 'B2H4'                  )] = """B2H4 nonplanar doubly bridged, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'B2H4_cation'           )] = """B2H4+ Radical Cation nonplanar doubly bridged, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'B_anion'               )] = """Boron Anion (B-) """
TAGL['%s-%s-reagent'    % (dbse, 'B_cation'              )] = """Boron Cation (B+) """
TAGL['%s-%s-reagent'    % (dbse, 'B_radical'             )] = """Boron Atom Radical (B) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'BCl3'                  )] = """Boron Trichloride (BCl3) planar, D3h """
TAGL['%s-%s-reagent'    % (dbse, 'BCl3_cation'           )] = """Boron Trichloride Radical Cation (BCl3+), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Be'                    )] = """Beryllium Atom (Be) """
TAGL['%s-%s-reagent'    % (dbse, 'Be_cation'             )] = """Beryllium Radical Cation (Be+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'BeH_radical'           )] = """Beryllium Hydride Radical (BeH), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'Benzene'               )] = """Benzene (C6H6), D6h """
TAGL['%s-%s-reagent'    % (dbse, 'Benzene_cation'        )] = """Benzene Radical Cation (C6H6+) compressed 2B2g """
TAGL['%s-%s-reagent'    % (dbse, 'BF3'                   )] = """Boron Trifluoride (BF3) planar, D3h """
TAGL['%s-%s-reagent'    % (dbse, 'BF3_cation'            )] = """Boron Trifluoride Radical Cation (BF3+) 2B2, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Bicyclo110butane'      )] = """Bicyclo[1.1.0]butane (C4H6), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'C'                     )] = """Carbon Atom (C) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'C2H5_ncl_cation'       )] = """C2H5+ Cation nonclassical bridged structure, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'C2H5_radical'          )] = """C2H5 Radical staggered, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'C_anion'               )] = """Carbon Atom Radical Anion (C-) quartet """
TAGL['%s-%s-reagent'    % (dbse, 'C_cation'              )] = """Carbon Atom Radical Cation (C+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'cAcetyl_anion'         )] = """Acetyl Anion (CH3CO-) HCCO cis, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'cAcetyl_radical'       )] = """Acetyl Radical (CH3CO) HCCO cis 2A', Cs """
TAGL['%s-%s-reagent'    % (dbse, 'CC'                    )] = """CC """
TAGL['%s-%s-reagent'    % (dbse, 'CC_anion'              )] = """CC- Radical Anion """
TAGL['%s-%s-reagent'    % (dbse, 'CCH_anion'             )] = """CCH- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'CCH_radical'           )] = """CCH Radical, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'CCl4'                  )] = """Carbon Tetrachloride (CCl4), Td """
TAGL['%s-%s-reagent'    % (dbse, 'CCO'                   )] = """CCO linear triplet """
TAGL['%s-%s-reagent'    % (dbse, 'CCO_anion'             )] = """CCO- Radical Anion """
TAGL['%s-%s-reagent'    % (dbse, 'CF2'                   )] = """CF2 """
TAGL['%s-%s-reagent'    % (dbse, 'CF2_anion'             )] = """CF2- Radical Anion 2B1 """
TAGL['%s-%s-reagent'    % (dbse, 'CF2_cation'            )] = """CF2+ Radical Cation """
TAGL['%s-%s-reagent'    % (dbse, 'CF3CN'                 )] = """CF3CN, C3v """
TAGL['%s-%s-reagent'    % (dbse, 'CF4'                   )] = """Carbon Tetrafluoride (CF4), Td """
TAGL['%s-%s-reagent'    % (dbse, 'CH2CHCH2_anion'        )] = """Prop-2-enyl Anion (CH2CHCH2-), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'CH2CHCH2_radical'      )] = """Prop-2-enyl Radical (CH2CHCH2) 2A2, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'CH2Cl2'                )] = """Dichloromethane (H2CCl2), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'CH2F2'                 )] = """Difluoromethane (H2CF2), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'CH3CH2O_anion'         )] = """CH3CH2O- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'CH3CH2O_radical'       )] = """CH3CH2O Radical 2A'', Cs """
TAGL['%s-%s-reagent'    % (dbse, 'CH3CH2S_anion'         )] = """CH3CH2S- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'CH3CH2S_radical'       )] = """CH3CH2S Radical """
TAGL['%s-%s-reagent'    % (dbse, 'CH3O_anion'            )] = """Methoxy Anion (CH3O-) """
TAGL['%s-%s-reagent'    % (dbse, 'CH3O_cation'           )] = """Methoxy Cation (CH3O+) 3A1, C3v """
TAGL['%s-%s-reagent'    % (dbse, 'CH3O_radical'          )] = """Methoxy Radical (CH3O) 2A', Cs """
TAGL['%s-%s-reagent'    % (dbse, 'CH3S_anion'            )] = """Methylthiyl Anion (CH3S-) """
TAGL['%s-%s-reagent'    % (dbse, 'CH3S_radical'          )] = """Methylthiyl Radical (CH3S) 2A', Cs """
TAGL['%s-%s-reagent'    % (dbse, 'CH_anion'              )] = """CH- Anion triplet """
TAGL['%s-%s-reagent'    % (dbse, 'CH_radical'            )] = """CH Radical doublet, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'CHCl3'                 )] = """Chloroform (HCCl3), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'CHF3'                  )] = """Trifluoromethane (HCF3), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'Cl2'                   )] = """Chlorine Molecule (Cl2), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'Cl2_anion'             )] = """Chlorine Molecule Radical Anion (Cl2-) """
TAGL['%s-%s-reagent'    % (dbse, 'Cl2_cation'            )] = """Chlorine Molecule Radical Cation (Cl2+) 2PIg """
TAGL['%s-%s-reagent'    % (dbse, 'Cl_anion'              )] = """Chloride Anion (Cl-) singlet """
TAGL['%s-%s-reagent'    % (dbse, 'Cl_cation'             )] = """Chlorine Atom Cation (Cl+) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'Cl_radical'            )] = """Chlorine Atom Radical (Cl) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'ClF'                   )] = """Chlorine Monofluoride (ClF) 1SG, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'ClF3'                  )] = """Chlorine Trifluoride (ClF3), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'ClF_cation'            )] = """Chlorine Monofluoride Radical Cation (ClF+) """
TAGL['%s-%s-reagent'    % (dbse, 'ClNO'                  )] = """Nitrosyl Chloride (ClNO), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'ClO_radical'           )] = """ClO Radical 2PI, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'CN_anion'              )] = """Cyanide Anion (CN-) """
TAGL['%s-%s-reagent'    % (dbse, 'CN_cation'             )] = """Cyanide Cation (CN+) """
TAGL['%s-%s-reagent'    % (dbse, 'CN_radical'            )] = """Cyano Radical (CN) 2Sigma+, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'CO'                    )] = """Carbon Monoxide (CO), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'CO2'                   )] = """Carbon Dioxide (CO2), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'CO2_cation'            )] = """Carbon Dioxide Radical Cation (CO2+), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'CO_cation'             )] = """Carbon Monoxide Radical Cation (CO+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'COF2'                  )] = """Carbonyl Fluoride (COF2), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'CS'                    )] = """Carbon Monosulfide (CS), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'CS2'                   )] = """Carbon Disulfide (CS2) linear, D*h """
TAGL['%s-%s-reagent'    % (dbse, 'CS2_cation'            )] = """Carbon Disulfide Radical Cation (CS2+), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'CS_cation'             )] = """Carbon Monosulfide Radical Cation (CS+) 2Sg doublet """
TAGL['%s-%s-reagent'    % (dbse, 'Cyanogen'              )] = """Cyanogen (NCCN), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'Cyanogen_cation'       )] = """Cyanogen Radical Cation (NCCN+), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'Cyclobutane'           )] = """Cyclobutane (C4H8), D2d """
TAGL['%s-%s-reagent'    % (dbse, 'Cyclobutene'           )] = """Cyclobutene (C4H6), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Cyclopropane'          )] = """Cyclopropane (C3H6), D3h """
TAGL['%s-%s-reagent'    % (dbse, 'Cyclopropene'          )] = """Cyclopropene (C3H4), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Cyclopropenyl_cation'  )] = """Cyclopropenyl Radical Cation (C3H4+), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Dimethylamine'         )] = """Dimethylamine ((CH3)2NH), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Dimethylether'         )] = """Dimethylether ((CH3)2O), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Dimethylsulfoxide'     )] = """Dimethylsulfoxide ((CH3)2SO), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Dimethylthioether'     )] = """Dimethyl Thioether ((CH3)2S), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Disilane'              )] = """Disilane (H3Si-SiH3), D3d """
TAGL['%s-%s-reagent'    % (dbse, 'Disilane_cation'       )] = """Disilane Radical Cation (H3Si-SiH3+) 2Ag, D3d """
TAGL['%s-%s-reagent'    % (dbse, 'Ethane'                )] = """Ethane (H3C-CH3), D3d """
TAGL['%s-%s-reagent'    % (dbse, 'Ethanol'               )] = """Ethanol (CH3CH2OH) trans, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Ethanol_cation'        )] = """Ethanol Radical Cation (CH3CH2OH+), C1 """
TAGL['%s-%s-reagent'    % (dbse, 'Ethene'                )] = """Ethene (H2C=CH2), D2h """
TAGL['%s-%s-reagent'    % (dbse, 'Ethene_cation'         )] = """Ethene Radical Cation (H2C=CH2+) doublet, D2h """
TAGL['%s-%s-reagent'    % (dbse, 'Ethylchloride'         )] = """Ethyl Chloride (CH3-CH2-Cl), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Ethyne'                )] = """Ethyne (HCCH), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'Ethyne_cation'         )] = """Ethyne Radical Cation (HCCH+) 2PIu linear doublet, D*h """
TAGL['%s-%s-reagent'    % (dbse, 'F2'                    )] = """Fluorine Molecule (F2), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'F2O'                   )] = """Oxygen Difluoride (F2O), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'F_anion'               )] = """Fluoride Anion (F-) singlet """
TAGL['%s-%s-reagent'    % (dbse, 'F_cation'              )] = """Fluorine Atom Cation (F+) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'F_radical'             )] = """Fluorine Atom Radical (F) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'Formaldehyde'          )] = """Formaldehyde (H2C=O), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'FormicAcid'            )] = """Formic Acid (HCOOH) HOCO cis, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Furan'                 )] = """Furan (cyclic C4H4O), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Furan_cation'          )] = """Furan Radical Cation (cyclic C4H4O+) """
TAGL['%s-%s-reagent'    % (dbse, 'Glyoxal'               )] = """Glyoxal (O=CH-CH=O) trans, C2h """
TAGL['%s-%s-reagent'    % (dbse, 'H2'                    )] = """Hydrogen Molecule (H2), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'H2C_anion'             )] = """Methylene Radical Anion (CH2-) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'H2C_cation'            )] = """Methylene Radical Cation (CH2+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'H2C_singlet'           )] = """Methylene (CH2) singlet, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'H2C_triplet'           )] = """Methylene (CH2) triplet, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCC'                 )] = """H2CCC """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCC_anion'           )] = """H2CCC- Radical Anion """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCCH_anion'          )] = """H2CCCH- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCCH_radical'        )] = """H2CCCH Radical """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCH_anion'           )] = """H2CCH- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCH_radical'         )] = """H2CCH Radical 2A', Cs """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCHO_anion'          )] = """H2CCHO- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCHO_radical'        )] = """H2CCHO Radical """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCN_anion'           )] = """H2CCN- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'H2CCN_radical'         )] = """H2CCN Radical """
TAGL['%s-%s-reagent'    % (dbse, 'H2Cl_cation'           )] = """H2Cl+ Cation """
TAGL['%s-%s-reagent'    % (dbse, 'H2CNC_anion'           )] = """H2CNC- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'H2CNC_radical'         )] = """H2CNC Radical """
TAGL['%s-%s-reagent'    % (dbse, 'H2COH_cation'          )] = """H2COH+ Cation """
TAGL['%s-%s-reagent'    % (dbse, 'H2COH_radical'         )] = """H2COH Radical, C1 """
TAGL['%s-%s-reagent'    % (dbse, 'H2CSH_cation'          )] = """H2CSH+ Cation 1A', Cs """
TAGL['%s-%s-reagent'    % (dbse, 'H2CSH_radical'         )] = """H2CSH Radical, C1 """
TAGL['%s-%s-reagent'    % (dbse, 'H2O'                   )] = """Water (H2O), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'H2O_cation'            )] = """Water Radical Cation (H2O+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'H3_cation'             )] = """H3+ Cation, D3h """
TAGL['%s-%s-reagent'    % (dbse, 'H3O_cation'            )] = """Hydronium Cation (H3O+) """
TAGL['%s-%s-reagent'    % (dbse, 'H_radical'             )] = """Hydrogen Atom Radical (H) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'HCCO_anion'            )] = """HCCO- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'HCCO_radical'          )] = """HCCO Radical nonlinear planar """
TAGL['%s-%s-reagent'    % (dbse, 'HCF'                   )] = """HCF 1A' """
TAGL['%s-%s-reagent'    % (dbse, 'HCF_anion'             )] = """HCF- Radical Anion 2A'' """
TAGL['%s-%s-reagent'    % (dbse, 'HCl'                   )] = """Hydrogen Chloride (HCl), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'HCl_cation'            )] = """Hydrogen Chloride Radical Cation (HCl+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'HCN'                   )] = """Hydrogen Cyanide (HCN), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'HCO_anion'             )] = """HCO- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'HCO_cation'            )] = """HCO+ Cation """
TAGL['%s-%s-reagent'    % (dbse, 'HCO_radical'           )] = """HCO Radical bent, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'He'                    )] = """Helium Atom (He) """
TAGL['%s-%s-reagent'    % (dbse, 'He_cation'             )] = """Helium Radical Cation (He+) """
TAGL['%s-%s-reagent'    % (dbse, 'HF'                    )] = """Hydrogen Fluoride (HF), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'HF_cation'             )] = """Hydrogen Fluoride Radical Cation (HF+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'HNO'                   )] = """HNO """
TAGL['%s-%s-reagent'    % (dbse, 'HNO_anion'             )] = """HNO- Radical Anion """
TAGL['%s-%s-reagent'    % (dbse, 'HOCl'                  )] = """Hypochlorous Acid (HOCl), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'HOF'                   )] = """HOF """
TAGL['%s-%s-reagent'    % (dbse, 'HOF_cation'            )] = """HOF+ Radical Cation """
TAGL['%s-%s-reagent'    % (dbse, 'HOO_anion'             )] = """Hydroperoxyl Anion (HOO-) 1A' """
TAGL['%s-%s-reagent'    % (dbse, 'HOO_radical'           )] = """Hydroperoxyl Radical (HOO) 2A'' """
TAGL['%s-%s-reagent'    % (dbse, 'HOOH'                  )] = """Hydrogen Peroxide (HOOH), C2 """
TAGL['%s-%s-reagent'    % (dbse, 'Hydrazine'             )] = """Hydrazine (H2N-NH2), C2 """
TAGL['%s-%s-reagent'    % (dbse, 'Isobutane'             )] = """Isobutane (C4H10), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'Isobutene'             )] = """Isobutene (C4H8) single bonds trans, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Isopropanol'           )] = """Isopropyl Alcohol ((CH3)2CH-OH) gauche isomer, C1 """
TAGL['%s-%s-reagent'    % (dbse, 'Ketene'                )] = """Ketene (H2C=C=O), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Li2'                   )] = """Dilithium (Li2), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'Li_anion'              )] = """Lithium Anion (Li-) """
TAGL['%s-%s-reagent'    % (dbse, 'Li_cation'             )] = """Lithium Cation (Li+) """
TAGL['%s-%s-reagent'    % (dbse, 'Li_radical'            )] = """Lithium Atom Radical (Li) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'LiF'                   )] = """Lithium Fluoride (LiF), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'LiH'                   )] = """Lithium Hydride (LiH), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'LiH_anion'             )] = """Lithium Hydride Radical Anion (LiH-) 2Sigma """
TAGL['%s-%s-reagent'    % (dbse, 'Me2CH_cation'          )] = """(CH3)2CH+ Cation """
TAGL['%s-%s-reagent'    % (dbse, 'Me2CH_radical'         )] = """(CH3)2CH Radical 2A', Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Me3C_radical'          )] = """t-Butyl Radical ((CH3)3C), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'MeCl'                  )] = """Methyl Chloride (CH3Cl), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'MeCl_cation'           )] = """Methyl Chloride Radical Cation (CH3Cl+) """
TAGL['%s-%s-reagent'    % (dbse, 'MeCN'                  )] = """Acetonitrile (H3C-CN), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'MeF'                   )] = """Methyl Fluoride (CH3F) """
TAGL['%s-%s-reagent'    % (dbse, 'MeF_cation'            )] = """Methyl Fluoride Radical Cation (CH3F+), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'MeOF'                  )] = """Methyl Hypochlorite (CH3-OF) """
TAGL['%s-%s-reagent'    % (dbse, 'MeOF_cation'           )] = """Methyl Hypochlorite Radical Cation (CH3-OF+) """
TAGL['%s-%s-reagent'    % (dbse, 'MeOH'                  )] = """Methanol (CH3-OH), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'MeOH_cation'           )] = """Methanol Radical Cation (CH3-OH+) """
TAGL['%s-%s-reagent'    % (dbse, 'MeSH'                  )] = """Methanethiol (CH3-SH) staggered, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'MeSH_cation'           )] = """Methanethiol Radical Cation (CH3-SH+) """
TAGL['%s-%s-reagent'    % (dbse, 'Methane'               )] = """Methane (CH4), Td """
TAGL['%s-%s-reagent'    % (dbse, 'Methane_cation'        )] = """Methane Radical Cation (CH4+) 2B2 doublet """
TAGL['%s-%s-reagent'    % (dbse, 'Methyl_anion'          )] = """Methyl Anion (CH3-), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'Methyl_cation'         )] = """Methyl Cation (CH3+), D3h """
TAGL['%s-%s-reagent'    % (dbse, 'Methyl_radical'        )] = """Methyl Radical (CH3), D3h """
TAGL['%s-%s-reagent'    % (dbse, 'Methylamine'           )] = """Methylamine (CH3-NH2), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Methylenecyclopropane' )] = """Methylenecyclopropane (C4H6), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Methylethylether'      )] = """Methyl Ethyl Ether (CH3-CH2-O-CH3) trans, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Methylformate'         )] = """Methyl Formate (CH3-O-CHO), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Methylnitrite'         )] = """Methylnitrite (CH3-O-N=O) NOCH trans ONOC cis, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Methylsilane'          )] = """Methylsilane (CH3-SiH3), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'Mg'                    )] = """Magnesium Atom (Mg) """
TAGL['%s-%s-reagent'    % (dbse, 'Mg_cation'             )] = """Magnesium Atom Radical Cation (Mg+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'N2'                    )] = """Nitrogen Molecule (N2), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'N2_cation_2PIu'        )] = """Nitrogen Molecule Radical Cation (N2+) 2PIu """
TAGL['%s-%s-reagent'    % (dbse, 'N2_cation_2SIGMAg'     )] = """Nitrogen Molecule Radical Cation (N2+) 2SIGMAg """
TAGL['%s-%s-reagent'    % (dbse, 'N2H2'                  )] = """HNNH """
TAGL['%s-%s-reagent'    % (dbse, 'N2H2_cation'           )] = """HNNH+ Radical Cation trans """
TAGL['%s-%s-reagent'    % (dbse, 'N2H3_cation'           )] = """N2H3+ Cation """
TAGL['%s-%s-reagent'    % (dbse, 'N2H3_radical'          )] = """N2H3 Radical """
TAGL['%s-%s-reagent'    % (dbse, 'N2O'                   )] = """Nitrous Oxide (N2O), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'N_cation'              )] = """Nitrogen Atom Cation (N+) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'N_radical'             )] = """Nitrogen Atom Radical (N) quartet """
TAGL['%s-%s-reagent'    % (dbse, 'Na2'                   )] = """Disodium (Na2), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'Na_anion'              )] = """Sodium Atom Anion (Na-) """
TAGL['%s-%s-reagent'    % (dbse, 'Na_cation'             )] = """Sodium Atom Cation (Na+) """
TAGL['%s-%s-reagent'    % (dbse, 'Na_radical'            )] = """Sodium Atom Radical (Na) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'NaCl'                  )] = """Sodium Chloride (NaCl), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'NCO_anion'             )] = """NCO- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'NCO_radical'           )] = """NCO Radical linear """
TAGL['%s-%s-reagent'    % (dbse, 'Ne'                    )] = """Neon Atom (Ne) """
TAGL['%s-%s-reagent'    % (dbse, 'Ne_cation'             )] = """Neon Atom Radical Cation (Ne+) """
TAGL['%s-%s-reagent'    % (dbse, 'NF3'                   )] = """Nitrogen Trifluoride (NF3), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'NH'                    )] = """NH triplet, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'NH2_anion'             )] = """NH2- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'NH2_cation'            )] = """NH2+ Cation 3B2, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'NH2_radical'           )] = """NH2 Radical 2B1, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'NH3'                   )] = """Ammonia (NH3), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'NH3_cation'            )] = """Ammonia Radical Cation (NH3+) planar doublet """
TAGL['%s-%s-reagent'    % (dbse, 'NH4_cation'            )] = """Ammonium Cation (NH4+) """
TAGL['%s-%s-reagent'    % (dbse, 'NH_anion'              )] = """NH- Radical Anion doublet """
TAGL['%s-%s-reagent'    % (dbse, 'NH_cation'             )] = """NH+ Radical Cation """
TAGL['%s-%s-reagent'    % (dbse, 'Nitromethane'          )] = """Nitromethane (CH3-NO2), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'NO2_anion'             )] = """Nitrite Anion (NO2-) """
TAGL['%s-%s-reagent'    % (dbse, 'NO2_radical'           )] = """Nitrogen Dioxide Radical (NO2) 2A1, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'NO_anion'              )] = """NO- Anion 3Sigma- triplet """
TAGL['%s-%s-reagent'    % (dbse, 'NO_radical'            )] = """NO Radical 2PI, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'O'                     )] = """Oxygen Atom (O) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'O2'                    )] = """Oxygen Molecule (O2) triplet, D*h """
TAGL['%s-%s-reagent'    % (dbse, 'O2_anion'              )] = """Oxygen Molecule Radical Anion (O2-) """
TAGL['%s-%s-reagent'    % (dbse, 'O2_cation'             )] = """Oxygen Molecule Radical Cation (O2+) """
TAGL['%s-%s-reagent'    % (dbse, 'O3'                    )] = """Ozone (O3), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'O3_anion'              )] = """Oxone Radical Anion (O3-), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'O_anion'               )] = """Oxygen Atom Radical Anion (O-) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'O_cation'              )] = """Oxygen Atom Radical Cation (O+) quartet """
TAGL['%s-%s-reagent'    % (dbse, 'OCS'                   )] = """Carbonyl Sulfide (O=C=S) linear, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'OCS_cation'            )] = """Carbonyl Sulfide Radical Cation (O=C=S+) """
TAGL['%s-%s-reagent'    % (dbse, 'OF_anion'              )] = """OF- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'OF_radical'            )] = """OF Radical """
TAGL['%s-%s-reagent'    % (dbse, 'OH_anion'              )] = """Hydroxide Anion (OH-) """
TAGL['%s-%s-reagent'    % (dbse, 'OH_cation'             )] = """Hydroxyl Cation (OH+) 3SIGMA triplet """
TAGL['%s-%s-reagent'    % (dbse, 'OH_radical'            )] = """Hydroxyl Radical (OH), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'Oxirane'               )] = """Oxirane (cyclic CH2-O-CH2), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'P2'                    )] = """Phosphorus Molecule (P2), D*h """
TAGL['%s-%s-reagent'    % (dbse, 'P2_cation'             )] = """Phosphorus Molecule Radical Cation (P2+) 2PIu """
TAGL['%s-%s-reagent'    % (dbse, 'P_anion'               )] = """Phosphorus Atom Anion (P-) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'P_cation'              )] = """Phosphorus Atom Cation (P+) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'P_radical'             )] = """Phosphorus Atom Radical (P) quartet """
TAGL['%s-%s-reagent'    % (dbse, 'PF3'                   )] = """Phosphorus Trifluoride (PF3), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'PH'                    )] = """PH triplet """
TAGL['%s-%s-reagent'    % (dbse, 'PH2_anion'             )] = """PH2- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'PH2_cation'            )] = """PH2+ Cation 1A1 """
TAGL['%s-%s-reagent'    % (dbse, 'PH2_radical'           )] = """PH2 Radical, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'PH3'                   )] = """Phosphine (PH3), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'PH3_cation'            )] = """Phosphine Radical Cation (PH3+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'PH4_cation'            )] = """PH4+ Cation """
TAGL['%s-%s-reagent'    % (dbse, 'PH_anion'              )] = """PH- Radical Anion doublet """
TAGL['%s-%s-reagent'    % (dbse, 'PH_cation'             )] = """PH+ Radical Cation doublet """
TAGL['%s-%s-reagent'    % (dbse, 'Phenol'                )] = """Phenol (C6H5-OH) planar, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Phenol_cation'         )] = """Phenol Radical Cation (C6H5-OH+) 2A'', Cs """
TAGL['%s-%s-reagent'    % (dbse, 'PO_anion'              )] = """PO- Anion 3SIGMA triplet """
TAGL['%s-%s-reagent'    % (dbse, 'PO_radical'            )] = """PO Radical 2PI """
TAGL['%s-%s-reagent'    % (dbse, 'Propane'               )] = """Propane (C3H8), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Propene'               )] = """Propene (C3H6), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Propylchloride'        )] = """Propyl Chloride (CH3CH2CH2Cl), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Propyne'               )] = """Propyne (C3H4), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'Pyridine'              )] = """Pyridine (cyclic C5H5N), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Pyrrole'               )] = """Pyrrole (cyclic C4H4NH) planar, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Pyrrole_cation'        )] = """Pyrrole Radical Cation (cyclic C4H4NH+) planar 2A2, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'S'                     )] = """Sulfur Atom (S) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'S2'                    )] = """Sulfur Molecule (S2) triplet, D*h """
TAGL['%s-%s-reagent'    % (dbse, 'S2_anion'              )] = """Sulfur Molecule Radical Anion (S2-) 2SIGMA doublet """
TAGL['%s-%s-reagent'    % (dbse, 'S2_cation'             )] = """Sulfur Molecule Radical Cation (S2+) 2PIg """
TAGL['%s-%s-reagent'    % (dbse, 'S2O'                   )] = """Disulfur Monoxide (S2O), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'S2O_anion'             )] = """Disulfur Monoxide Radical Anion (S2O-) A'' doublet, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'S_anion'               )] = """Sulfur Atom Radical Anion (S-) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'S_cation'              )] = """Sulfur Atom Radical Cation (S+) quartet """
TAGL['%s-%s-reagent'    % (dbse, 'SH2'                   )] = """Hydrogen Sulfide (H2S), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'SH2_cation_2A1'        )] = """Hydrogen Sulfide Radical Cation (SH2+) 2A1 doublet """
TAGL['%s-%s-reagent'    % (dbse, 'SH2_cation_2B1'        )] = """Hydrogen Sulfide Radical Cation (SH2+) 2B1 doublet """
TAGL['%s-%s-reagent'    % (dbse, 'SH3_cation'            )] = """SH3+ Cation """
TAGL['%s-%s-reagent'    % (dbse, 'SH_anion'              )] = """SH- Anion """
TAGL['%s-%s-reagent'    % (dbse, 'SH_cation'             )] = """SH+ Cation triplet """
TAGL['%s-%s-reagent'    % (dbse, 'SH_radical'            )] = """SH Radical, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'Si'                    )] = """Silicon Atom (Si) triplet """
TAGL['%s-%s-reagent'    % (dbse, 'Si2'                   )] = """Silicon Molecule (Si2) triplet 3SigmaG-, D*h """
TAGL['%s-%s-reagent'    % (dbse, 'Si2H2'                 )] = """Si2H2 """
TAGL['%s-%s-reagent'    % (dbse, 'Si2H2_cation'          )] = """Si2H2+ Radical Cation """
TAGL['%s-%s-reagent'    % (dbse, 'Si2H4'                 )] = """Si2H4 1Ag, C2h """
TAGL['%s-%s-reagent'    % (dbse, 'Si2H4_cation'          )] = """Si2H4+ Radical Cation 2B3u, D2h """
TAGL['%s-%s-reagent'    % (dbse, 'Si2H5'                 )] = """Si2H5 Radical """
TAGL['%s-%s-reagent'    % (dbse, 'Si2H5_cation'          )] = """Si2H5+ Cation """
TAGL['%s-%s-reagent'    % (dbse, 'Si_anion'              )] = """Silicon Atom Radical Anion (Si-) quartet """
TAGL['%s-%s-reagent'    % (dbse, 'Si_cation'             )] = """Silicon Atom Radical Cation (Si+) doublet """
TAGL['%s-%s-reagent'    % (dbse, 'SiCl4'                 )] = """Silicon Tetrachloride (SiCl4), Td """
TAGL['%s-%s-reagent'    % (dbse, 'SiF4'                  )] = """Silicon Tetrafluoride (SiF4), Td """
TAGL['%s-%s-reagent'    % (dbse, 'SiH2_anion'            )] = """Silylene Radical Anion (SiH2-) 2B1 doublet """
TAGL['%s-%s-reagent'    % (dbse, 'SiH2_cation'           )] = """Silylene Radical Cation (SiH2+) """
TAGL['%s-%s-reagent'    % (dbse, 'SiH2_singlet'          )] = """Silylene (SiH2) singlet 1A1, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'SiH2_triplet'          )] = """Silylene (SiH2) triplet 3B1, C2v """
TAGL['%s-%s-reagent'    % (dbse, 'SiH3_anion'            )] = """Silyl Anion (SiH3-), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'SiH3_cation'           )] = """Silyl Cation (SiH3+), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'SiH3_radical'          )] = """Silyl Radical (SiH3), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'SiH4'                  )] = """Silane (SiH4), Td """
TAGL['%s-%s-reagent'    % (dbse, 'SiH4_cation'           )] = """Silane Radical Cation (SiH4+) doublet, Cs """
TAGL['%s-%s-reagent'    % (dbse, 'SiH5_cation'           )] = """SiH5+ Cation, C1 """
TAGL['%s-%s-reagent'    % (dbse, 'SiH_anion'             )] = """SiH- Anion 3SIGMA """
TAGL['%s-%s-reagent'    % (dbse, 'SiH_radical'           )] = """SiH Radical """
TAGL['%s-%s-reagent'    % (dbse, 'SiO'                   )] = """Silicon Monoxide (SiO), C*v """
TAGL['%s-%s-reagent'    % (dbse, 'SO'                    )] = """Sulfur Monoxide (SO) triplet, C*v """
TAGL['%s-%s-reagent'    % (dbse, 'SO2'                   )] = """Sulfur Dioxide (SO2), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'SO2_anion'             )] = """Sulfur Dioxide Radical Anion (SO2-) """
TAGL['%s-%s-reagent'    % (dbse, 'Spiropentane'          )] = """Spiropentane (C5H8), D2d """
TAGL['%s-%s-reagent'    % (dbse, 't13Butadiene'          )] = """Trans-1,3-Butadiene (C4H6), C2h """
TAGL['%s-%s-reagent'    % (dbse, 'tButane'               )] = """Trans-Butane (C4H10), C2h """
TAGL['%s-%s-reagent'    % (dbse, 'tEthylamine'           )] = """Trans-Ethylamine (CH3-CH2-NH2), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Tetrachloroethene'     )] = """Tetrachloroethene (Cl2C=CCl2), D2h """
TAGL['%s-%s-reagent'    % (dbse, 'Tetrafluoroethene'     )] = """Tetrafluoroethene (F2C=CF2), D2h """
TAGL['%s-%s-reagent'    % (dbse, 'Thioethanol'           )] = """Thio Ethanol (CH3-CH2-SH), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Thioformaldehyde'      )] = """Thioformaldehyde (CH2S) """
TAGL['%s-%s-reagent'    % (dbse, 'Thioformaldehyde_anion' )] = """Thioformaldehyde Radical Anion (CH2S-) """
TAGL['%s-%s-reagent'    % (dbse, 'Thioformaldehyde_cation' )] = """Thioformaldehyde Radical Cation (CH2S+) """
TAGL['%s-%s-reagent'    % (dbse, 'Thiooxirane'           )] = """Thiooxirane (cyclic CH2-S-CH2), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Thiooxirane_cation'    )] = """Thiooxirane Radical Cation (cyclic CH2-S-CH2+), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Thiophene'             )] = """Thiophene (cyclic C4H4S), C2v """
TAGL['%s-%s-reagent'    % (dbse, 'Toluene'               )] = """Toluene (C6H5-CH3), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Toluene_cation'        )] = """Toluene Radical Cation (C6H5-CH3+), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Trimethylamine'        )] = """Trimethyl Amine ((CH3)3N), C3v """
TAGL['%s-%s-reagent'    % (dbse, 'Vinyl_ncl_cation'      )] = """Nonclassical Vinyl Cation (C2H3+) """
TAGL['%s-%s-reagent'    % (dbse, 'Vinylchloride'         )] = """Vinyl Chloride (H2C=CHCl), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Vinylcyanide'          )] = """Vinyl Cyanide (H2C=CHCN), Cs """
TAGL['%s-%s-reagent'    % (dbse, 'Vinylfluoride'         )] = """Vinyl Fluoride (H2C=CHF), Cs """

# <<< Molecule Specifications >>>
G2_2Butyne = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     2.07195500
C        0.00000000     0.00000000     0.60997000
C        0.00000000     0.00000000    -0.60997000
C        0.00000000     0.00000000    -2.07195500
H        0.00000000     1.02069600     2.46456200
H       -0.88394900    -0.51034800     2.46456200
H        0.88394900    -0.51034800     2.46456200
H        0.00000000     1.02069600    -2.46456200
H        0.88394900    -0.51034800    -2.46456200
H       -0.88394900    -0.51034800    -2.46456200
units angstrom
}
""")

G2_Acetaldehyde = input.process_input("""
molecule dimer {
0 1
O        1.21805500     0.36124000     0.00000000
C        0.00000000     0.46413300     0.00000000
H       -0.47724100     1.46529500     0.00000000
C       -0.94810200    -0.70013800     0.00000000
H       -0.38594600    -1.63423600     0.00000000
H       -1.59632100    -0.65247500     0.88094600
H       -1.59632100    -0.65247500    -0.88094600
units angstrom
}
""")

G2_Acetaldehyde_cation = input.process_input("""
molecule dimer {
1 2
O        1.20836682     0.32239770     0.00000000
C        0.00000000     0.46429315     0.00000000
H       -0.34038411     1.51908997     0.00000000
C       -0.95327027    -0.68482580     0.00000000
H       -0.43297255    -1.64073184     0.00000000
H       -1.58697814    -0.56717191     0.88745207
H       -1.58697814    -0.56717191    -0.88745207
units angstrom
}
""")

G2_Acetamide = input.process_input("""
molecule dimer {
0 1
O        0.42454600     1.32702400     0.00803400
C        0.07715800     0.14978900    -0.00424900
N        0.98551800    -0.87853700    -0.04891000
C       -1.37147500    -0.28866500    -0.00014400
H        0.70795200    -1.82424900     0.16994200
H       -1.99722900     0.58492200    -0.17547700
H       -1.56084200    -1.03927000    -0.77168600
H       -1.63211300    -0.72300700     0.96981400
H        1.95313300    -0.63157400     0.11186600
units angstrom
}
""")

G2_AceticAcid = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.15456000     0.00000000
O        0.16638400     1.36008400     0.00000000
O       -1.23644900    -0.41503600     0.00000000
H       -1.86764600     0.33358200     0.00000000
C        1.07377600    -0.89274800     0.00000000
H        2.04818900    -0.40813500     0.00000000
H        0.96866100    -1.52835300     0.88174700
H        0.96866100    -1.52835300    -0.88174700
units angstrom
}
""")

G2_Acetone = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000     1.40559100
C        0.00000000     0.00000000     0.17906000
C        0.00000000     1.28549000    -0.61634200
C        0.00000000    -1.28549000    -0.61634200
H        0.00000000     2.13491700     0.06653500
H        0.00000000    -2.13491700     0.06653500
H       -0.88108600     1.33154800    -1.26401300
H        0.88108600     1.33154800    -1.26401300
H        0.88108600    -1.33154800    -1.26401300
H       -0.88108600    -1.33154800    -1.26401300
units angstrom
}
""")

G2_Acetylchloride = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.52387800     0.00000000
C        1.48607500     0.71637700     0.00000000
Cl      -0.45228600    -1.21799900     0.00000000
O       -0.84553900     1.37494000     0.00000000
H        1.70102700     1.78479300     0.00000000
H        1.91784700     0.24006700     0.88267900
H        1.91784700     0.24006700    -0.88267900
units angstrom
}
""")

G2_Acetylfluoride = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.18639600     0.00000000
O        0.12665100     1.37719900     0.00000000
F       -1.24395000    -0.38274500     0.00000000
C        1.04945400    -0.87622400     0.00000000
H        2.03588300    -0.41709900     0.00000000
H        0.92486900    -1.50840700     0.88154900
H        0.92486900    -1.50840700    -0.88154900
units angstrom
}
""")

G2_Al_anion = input.process_input("""
molecule dimer {
-1 3
Al       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Al_cation = input.process_input("""
molecule dimer {
1 1
Al       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Al_radical = input.process_input("""
molecule dimer {
0 2
Al       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_AlCl3 = input.process_input("""
molecule dimer {
0 1
Al       0.00000000     0.00000000     0.00000000
Cl       0.00000000     2.06904100     0.00000000
Cl       1.79184200    -1.03452000     0.00000000
Cl      -1.79184200    -1.03452000     0.00000000
units angstrom
}
""")

G2_AlF3 = input.process_input("""
molecule dimer {
0 1
Al       0.00000000     0.00000000     0.00000000
F        0.00000000     1.64472000     0.00000000
F        1.42436900    -0.82236000     0.00000000
F       -1.42436900    -0.82236000     0.00000000
units angstrom
}
""")

G2_Allene = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.00000000
C        0.00000000     0.00000000     1.31119000
C        0.00000000     0.00000000    -1.31119000
H        0.00000000     0.92677800     1.87664200
H        0.00000000    -0.92677800     1.87664200
H        0.92677800     0.00000000    -1.87664200
H       -0.92677800     0.00000000    -1.87664200
units angstrom
}
""")

G2_Allene_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
C        0.00000000     1.31169643     0.00000000
C        0.00000000    -1.31169643    -0.00000000
H        0.84425956     1.87480477     0.40303304
H       -0.84422501     1.87486427    -0.40301655
H        0.84425956    -1.87480477    -0.40303304
H       -0.84422501    -1.87486427     0.40301655
units angstrom
}
""")

G2_Aniline = input.process_input("""
molecule dimer {
0 1
C       -0.00673100     0.93476600     0.00000000
C        0.00665400     0.22192400     1.20574700
C        0.00665400     0.22192400    -1.20574700
C        0.00665400    -1.17061600     1.20247300
C        0.00665400    -1.17061600    -1.20247300
C        0.01050300    -1.87723200     0.00000000
N        0.06416100     2.33791900     0.00000000
H       -0.34004400     2.75863300    -0.83040200
H       -0.34004400     2.75863300     0.83040200
H        0.01054300    -1.70521400     2.14934400
H        0.01054300    -1.70521400    -2.14934400
H        0.00997700    -2.96354500     0.00000000
H        0.00878800     0.76518900     2.14913900
H        0.00878800     0.76518900    -2.14913900
units angstrom
}
""")

G2_Aniline_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.92274200
C        0.00000000     1.23865500     0.20874000
C        0.00000000    -1.23865500     0.20874000
C        0.00000000     1.22606600    -1.13664100
C        0.00000000    -1.22606600    -1.13664100
C        0.00000000     0.00000000    -1.83111700
N        0.00000000     0.00000000     2.25588200
H        0.00000000    -0.86394700     2.78904100
H        0.00000000     0.86394700     2.78904100
H        0.00000000     2.15538700    -1.69656200
H        0.00000000    -2.15538700    -1.69656200
H        0.00000000     0.00000000    -2.91640500
H        0.00000000     2.17311400     0.76266700
H        0.00000000    -2.17311400     0.76266700
units angstrom
}
""")

G2_Ar = input.process_input("""
molecule dimer {
0 1
Ar       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Ar_cation = input.process_input("""
molecule dimer {
1 2
Ar       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Aziridine = input.process_input("""
molecule dimer {
0 1
C       -0.03845000    -0.39732600     0.73942100
N       -0.03845000     0.87518900     0.00000000
C       -0.03845000    -0.39732600    -0.73942100
H        0.90305200     1.26823900     0.00000000
H       -0.95566100    -0.60492600     1.28004700
H       -0.95566100    -0.60492600    -1.28004700
H        0.86940900    -0.70839900     1.24903300
H        0.86940900    -0.70839900    -1.24903300
units angstrom
}
""")

G2_B2F4 = input.process_input("""
molecule dimer {
0 1
B        0.00000000     0.00000000     0.00000000
B        1.71325435     0.00000000     0.00000000
F        2.40684300     1.13348094     0.00000000
F        2.40684300    -1.13348094     0.00000000
F       -0.69358865     0.00000000     1.13348094
F       -0.69358865     0.00000000    -1.13348094
units angstrom
}
""")

G2_B2F4_cation = input.process_input("""
molecule dimer {
1 2
B        0.00000000     0.00000000     0.00000000
B        2.03812360     0.00000000     0.00000000
F        2.48301612     1.19811530     0.00000000
F        2.48301612    -1.19811530     0.00000000
F       -0.44489252     0.00000000     1.19811530
F       -0.44489252     0.00000000    -1.19811530
units angstrom
}
""")

G2_B2H4 = input.process_input("""
molecule dimer {
0 1
B        0.00000000     0.72934423    -0.11467318
B        0.00000000    -0.72934423    -0.11467318
H       -0.90314031     0.00000000     0.54872914
H        0.90314031     0.00000000     0.54872914
H        0.00000000     1.89702991     0.02463674
H        0.00000000    -1.89702991     0.02463674
units angstrom
}
""")

G2_B2H4_cation = input.process_input("""
molecule dimer {
1 2
B        0.00000000     0.76879700    -0.04735558
B        0.00000000    -0.76879700    -0.04735558
H       -1.00855278     0.00000000     0.31086854
H        1.00855278     0.00000000     0.31086854
H        0.00000000     1.94037629    -0.07409066
H        0.00000000    -1.94037629    -0.07409066
units angstrom
}
""")

G2_B_anion = input.process_input("""
molecule dimer {
-1 3
B        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_B_cation = input.process_input("""
molecule dimer {
1 1
B        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_B_radical = input.process_input("""
molecule dimer {
0 2
B        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_BCl3 = input.process_input("""
molecule dimer {
0 1
B        0.00000000     0.00000000     0.00000000
Cl       0.00000000     1.73535200     0.00000000
Cl       1.50285900    -0.86767600     0.00000000
Cl      -1.50285900    -0.86767600     0.00000000
units angstrom
}
""")

G2_BCl3_cation = input.process_input("""
molecule dimer {
1 2
B        0.00000000     0.00000000     0.18895042
Cl       0.00000000     0.00000000     1.84953710
Cl       0.00000000    -1.33808914    -0.95255538
Cl       0.00000000     1.33808914    -0.95255538
units angstrom
}
""")

G2_Be = input.process_input("""
molecule dimer {
0 1
Be       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Be_cation = input.process_input("""
molecule dimer {
1 2
Be       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_BeH_radical = input.process_input("""
molecule dimer {
0 2
Be       0.00000000     0.00000000     0.26965400
H        0.00000000     0.00000000    -1.07861600
units angstrom
}
""")

G2_Benzene = input.process_input("""
molecule dimer {
0 1
C        0.00000000     1.39524800     0.00000000
C        1.20832000     0.69762400     0.00000000
C        1.20832000    -0.69762400     0.00000000
C        0.00000000    -1.39524800     0.00000000
C       -1.20832000    -0.69762400     0.00000000
C       -1.20832000     0.69762400     0.00000000
H        0.00000000     2.48236000     0.00000000
H        2.14978700     1.24118000     0.00000000
H        2.14978700    -1.24118000     0.00000000
H        0.00000000    -2.48236000     0.00000000
H       -2.14978700    -1.24118000     0.00000000
H       -2.14978700     1.24118000     0.00000000
units angstrom
}
""")

G2_Benzene_cation = input.process_input("""
molecule dimer {
1 2
C       -0.69395173     0.00000000     0.00000000
C        2.06751942     0.00000000    -0.00000000
C        0.00000000     0.00000000    -1.24659620
C        0.00000000     0.00000000     1.24659620
C        1.37356769     0.00000000    -1.24659620
C        1.37356769     0.00000000     1.24659620
H       -1.78116635     0.00000000     0.00000000
H        3.15473404     0.00000000    -0.00000000
H       -0.56406266    -0.00000000    -2.17394095
H       -0.56406266     0.00000000     2.17394095
H        1.93763035     0.00000000    -2.17394095
H        1.93763035    -0.00000000     2.17394095
units angstrom
}
""")

G2_BF3 = input.process_input("""
molecule dimer {
0 1
B        0.00000000     0.00000000     0.00000000
F        0.00000000     1.32176000     0.00000000
F        1.14467800    -0.66088000     0.00000000
F       -1.14467800    -0.66088000     0.00000000
units angstrom
}
""")

G2_BF3_cation = input.process_input("""
molecule dimer {
1 2
B        0.00000000     0.00000000     0.00000000
F        1.72487052     0.00000000     0.00000000
F       -0.33409601     1.21351798     0.00000000
F       -0.33409601    -1.21351798    -0.00000000
units angstrom
}
""")

G2_Bicyclo110butane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     1.13134300     0.31042400
C        0.00000000    -1.13134300     0.31042400
C        0.74795200     0.00000000    -0.31181200
C       -0.74795200     0.00000000    -0.31181200
H        0.00000000     1.23703300     1.39761700
H        0.00000000     2.07737500    -0.22766800
H        0.00000000    -1.23703300     1.39761700
H        0.00000000    -2.07737500    -0.22766800
H        1.41441000     0.00000000    -1.16162600
H       -1.41441000     0.00000000    -1.16162600
units angstrom
}
""")

G2_C = input.process_input("""
molecule dimer {
0 3
C        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_C2H5_ncl_cation = input.process_input("""
molecule dimer {
1 1
C        0.00000000     0.69019560    -0.06247487
C        0.00000000    -0.69019560    -0.06247487
H        0.00000000     0.00000000     1.04469040
H        0.93570131     1.24482654    -0.07374799
H       -0.93570131     1.24482654    -0.07374799
H       -0.93570131    -1.24482654    -0.07374799
H        0.93570131    -1.24482654    -0.07374799
units angstrom
}
""")

G2_C2H5_radical = input.process_input("""
molecule dimer {
0 2
C       -0.01435900    -0.69461700     0.00000000
C       -0.01435900     0.79447300     0.00000000
H        1.00610100    -1.10404200     0.00000000
H       -0.51703700    -1.09361300     0.88483900
H       -0.51703700    -1.09361300    -0.88483900
H        0.10013700     1.34606500     0.92370500
H        0.10013700     1.34606500    -0.92370500
units angstrom
}
""")

G2_C_anion = input.process_input("""
molecule dimer {
-1 4
C        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_C_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_cAcetyl_anion = input.process_input("""
molecule dimer {
-1 1
C       -0.32865919     1.08560372    -0.24213691
C       -0.33348486    -0.54490805    -0.24569218
O        0.62651950    -0.95914372     0.46158300
H        0.49460891     1.50305765     0.36439898
H       -0.24433064     1.46295905    -1.27441047
H       -1.28956995     1.46295905     0.14432204
units angstrom
}
""")

G2_cAcetyl_radical = input.process_input("""
molecule dimer {
0 2
C       -0.97829100    -0.64781400     0.00000000
C        0.00000000     0.50628300     0.00000000
H       -0.45555100    -1.60783700     0.00000000
H       -1.61762600    -0.56327100     0.88106100
H       -1.61762600    -0.56327100    -0.88106100
O        1.19506900     0.44794500     0.00000000
units angstrom
}
""")

G2_CC = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -0.63191401
C        0.00000000     0.00000000     0.63191401
units angstrom
}
""")

G2_CC_anion = input.process_input("""
molecule dimer {
-1 2
C        0.00000000     0.00000000    -0.64146028
C        0.00000000     0.00000000     0.64146028
units angstrom
}
""")

G2_CCH_anion = input.process_input("""
molecule dimer {
-1 1
C        0.00000000     0.00000000    -0.49917112
C        0.00000000     0.00000000     0.76151220
H        0.00000000     0.00000000    -1.57404646
units angstrom
}
""")

G2_CCH_radical = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000    -0.46262800
C        0.00000000     0.00000000     0.71716200
H        0.00000000     0.00000000    -1.52719800
units angstrom
}
""")

G2_CCl4 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.00000000
Cl       1.02134000     1.02134000     1.02134000
Cl      -1.02134000    -1.02134000     1.02134000
Cl      -1.02134000     1.02134000    -1.02134000
Cl       1.02134000    -1.02134000    -1.02134000
units angstrom
}
""")

G2_CCO = input.process_input("""
molecule dimer {
0 3
C        0.05585944     0.00000000     0.00000000
C        1.43205267     0.00000000     0.00000000
O       -1.11593408     0.00000000     0.00000000
units angstrom
}
""")

G2_CCO_anion = input.process_input("""
molecule dimer {
-1 2
C        0.10036789     0.00000000     0.00000000
C        1.41076260     0.00000000     0.00000000
O       -1.13334786     0.00000000     0.00000000
units angstrom
}
""")

G2_CF2 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.60464260
F        0.00000000     1.03621209    -0.20154753
F        0.00000000    -1.03621209    -0.20154753
units angstrom
}
""")

G2_CF2_anion = input.process_input("""
molecule dimer {
-1 2
C       -0.54354195     0.00000000    -0.45608585
F       -0.53586937     0.00000000     1.00657555
F        0.89823067     0.00000000    -0.70251832
units angstrom
}
""")

G2_CF2_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.42685029
F        0.00000000     1.09029674    -0.14228343
F        0.00000000    -1.09029674    -0.14228343
units angstrom
}
""")

G2_CF3CN = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -0.32635000
C        0.00000000     0.00000000     1.15083000
F        0.00000000     1.25757900    -0.78722500
F        1.08909600    -0.62879000    -0.78722500
F       -1.08909600    -0.62879000    -0.78722500
N        0.00000000     0.00000000     2.32974100
units angstrom
}
""")

G2_CF4 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.00000000
F        0.76743600     0.76743600     0.76743600
F       -0.76743600    -0.76743600     0.76743600
F       -0.76743600     0.76743600    -0.76743600
F        0.76743600    -0.76743600    -0.76743600
units angstrom
}
""")

G2_CH2CHCH2_anion = input.process_input("""
molecule dimer {
-1 1
C        0.00000000     0.00000000     0.37533829
H        0.00000000     0.00000000     1.47780278
C        0.00000000    -1.27671168    -0.17698893
C        0.00000000     1.27671168    -0.17698893
H        0.00000000    -2.16314629     0.45117377
H        0.00000000     2.16314629     0.45117377
H        0.00000000    -1.43844444    -1.25415647
H        0.00000000     1.43844444    -1.25415647
units angstrom
}
""")

G2_CH2CHCH2_radical = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000     0.44601039
H        0.00000000     0.00000000     1.53402965
C        0.00000000    -1.21762516    -0.19619829
C        0.00000000     1.21762516    -0.19619829
H        0.00000000    -2.15128955     0.35090337
H        0.00000000     2.15128955     0.35090337
H        0.00000000    -1.27757339    -1.27875965
H        0.00000000     1.27757339    -1.27875965
units angstrom
}
""")

G2_CH2Cl2 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.75994500
Cl       0.00000000     1.47420000    -0.21511500
Cl       0.00000000    -1.47420000    -0.21511500
H       -0.89458500     0.00000000     1.37712700
H        0.89458500     0.00000000     1.37712700
units angstrom
}
""")

G2_CH2F2 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.50290300
F        0.00000000     1.10971600    -0.29060100
F        0.00000000    -1.10971600    -0.29060100
H       -0.90836900     0.00000000     1.10669900
H        0.90836900     0.00000000     1.10669900
units angstrom
}
""")

G2_CH3CH2O_anion = input.process_input("""
molecule dimer {
-1 1
C        0.26874680    -1.13905086     0.21591472
C        0.26136244     0.41854997     0.20998203
O       -0.68820968     0.93883794    -0.55291674
H        0.23368504     0.69143815     1.32379502
H        1.34305261     0.69143815    -0.05702318
H        1.05724526    -1.58144458     0.84940479
H       -0.70800466    -1.49456494     0.56485417
H        0.39904377    -1.49456494    -0.81307743
units angstrom
}
""")

G2_CH3CH2O_radical = input.process_input("""
molecule dimer {
0 2
C        1.00475700    -0.56826300     0.00000000
C        0.00000000     0.58869100     0.00000000
O       -1.26006200     0.00072900     0.00000000
H        0.14695600     1.20468100     0.89652900
H        0.14695600     1.20468100    -0.89652900
H        2.01936300    -0.16410000     0.00000000
H        0.86934000    -1.18683200     0.88807100
H        0.86934000    -1.18683200    -0.88807100
units angstrom
}
""")

G2_CH3CH2S_anion = input.process_input("""
molecule dimer {
-1 1
C        0.07347988    -1.44165111     0.77722101
C        0.07325218     0.08629022     0.77481250
S       -0.08444864     0.79854140    -0.89324124
H        0.16936507    -1.86415764     1.79143037
H       -0.85344950    -1.81352024     0.33111302
H        0.90038154    -1.81352024     0.16530266
H       -0.74597370     0.42335050     1.42846195
H        1.00046254     0.42335050     1.26335071
units angstrom
}
""")

G2_CH3CH2S_radical = input.process_input("""
molecule dimer {
0 2
C       -1.59365494    -0.41513574     0.00415508
C       -0.46652748     0.60724751     0.01166674
S        1.19129188    -0.10375186    -0.01168753
H       -2.56442407     0.08843110     0.04326521
H       -1.56179132    -1.02311363    -0.90261391
H       -1.51471663    -1.08529801     0.86235819
H       -0.55234711     1.30274318    -0.82910813
H       -0.50629650     1.22459645     0.91816812
units angstrom
}
""")

G2_CH3O_anion = input.process_input("""
molecule dimer {
-1 1
C       -0.00020322     0.52996763     0.00000000
O        0.00030352    -0.79323293     0.00000000
H       -1.02182517     1.05496336     0.00000000
H        0.51030818     1.05554717     0.88457559
H        0.51030818     1.05554717    -0.88457559
units angstrom
}
""")

G2_CH3O_cation = input.process_input("""
molecule dimer {
1 3
C        0.00000000     0.00000000    -0.55388766
O        0.00000000     0.00000000     0.75485421
H        0.00000000    -1.06283967    -0.90516924
H       -0.92044615     0.53141983    -0.90516924
H        0.92044615     0.53141983    -0.90516924
units angstrom
}
""")

G2_CH3O_radical = input.process_input("""
molecule dimer {
0 2
C       -0.00861800    -0.58647500     0.00000000
O       -0.00861800     0.79954100     0.00000000
H        1.05536300    -0.86875600     0.00000000
H       -0.46735800    -1.00436300     0.90327900
H       -0.46735800    -1.00436300    -0.90327900
units angstrom
}
""")

G2_CH3S_anion = input.process_input("""
molecule dimer {
-1 1
C       -0.00020053    -1.11779899     0.00000000
S        0.00012933     0.70779900     0.00000000
H        1.01520150    -1.53954878     0.00000000
H       -0.50803380    -1.53922061    -0.87943251
H       -0.50803380    -1.53922061     0.87943251
units angstrom
}
""")

G2_CH3S_radical = input.process_input("""
molecule dimer {
0 2
C       -0.00385600     1.10622200     0.00000000
S       -0.00385600    -0.69257900     0.00000000
H        1.04326900     1.42705700     0.00000000
H       -0.47921700     1.50843700     0.89519700
H       -0.47921700     1.50843700    -0.89519700
units angstrom
}
""")

G2_CH_anion = input.process_input("""
molecule dimer {
-1 3
C        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.16101900
units angstrom
}
""")

G2_CH_radical = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000     0.16007400
H        0.00000000     0.00000000    -0.96044600
units angstrom
}
""")

G2_CHCl3 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.45167900
H        0.00000000     0.00000000     1.53758600
Cl       0.00000000     1.68172300    -0.08328700
Cl       1.45641500    -0.84086200    -0.08328700
Cl      -1.45641500    -0.84086200    -0.08328700
units angstrom
}
""")

G2_CHF3 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.34102300
H        0.00000000     0.00000000     1.42948500
F        0.00000000     1.25820000    -0.12872700
F        1.08963300    -0.62910000    -0.12872700
F       -1.08963300    -0.62910000    -0.12872700
units angstrom
}
""")

G2_Cl2 = input.process_input("""
molecule dimer {
0 1
Cl       0.00000000     0.00000000     1.00754100
Cl       0.00000000     0.00000000    -1.00754100
units angstrom
}
""")

G2_Cl2_anion = input.process_input("""
molecule dimer {
-1 2
Cl       0.00000000     0.00000000     0.00000000
Cl       0.00000000     0.00000000     2.65181100
units angstrom
}
""")

G2_Cl2_cation = input.process_input("""
molecule dimer {
1 2
Cl       0.00000000     0.00000000     0.00000000
Cl       0.00000000     0.00000000     1.92809700
units angstrom
}
""")

G2_Cl_anion = input.process_input("""
molecule dimer {
-1 1
Cl       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Cl_cation = input.process_input("""
molecule dimer {
1 3
Cl       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Cl_radical = input.process_input("""
molecule dimer {
0 2
Cl       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_ClF = input.process_input("""
molecule dimer {
0 1
F        0.00000000     0.00000000    -1.08479400
Cl       0.00000000     0.00000000     0.57430200
units angstrom
}
""")

G2_ClF3 = input.process_input("""
molecule dimer {
0 1
Cl       0.00000000     0.00000000     0.37679600
F        0.00000000     0.00000000    -1.25834600
F        0.00000000     1.71454400     0.27331000
F        0.00000000    -1.71454400     0.27331000
units angstrom
}
""")

G2_ClF_cation = input.process_input("""
molecule dimer {
1 2
Cl       0.00000000     0.00000000     0.00000000
F        0.00000000     0.00000000     1.55480000
units angstrom
}
""")

G2_ClNO = input.process_input("""
molecule dimer {
0 1
Cl      -0.53772400    -0.96129100     0.00000000
N        0.00000000     0.99703700     0.00000000
O        1.14266400     1.17033500     0.00000000
units angstrom
}
""")

G2_ClO_radical = input.process_input("""
molecule dimer {
0 2
Cl       0.00000000     0.00000000     0.51417200
O        0.00000000     0.00000000    -1.09261500
units angstrom
}
""")

G2_CN_anion = input.process_input("""
molecule dimer {
-1 1
C        0.00000000     0.00000000     0.00000000
N        0.00000000     0.00000000     1.20015400
units angstrom
}
""")

G2_CN_cation = input.process_input("""
molecule dimer {
1 1
C        0.00000000     0.00000000    -0.65732686
N        0.00000000     0.00000000     0.56342302
units angstrom
}
""")

G2_CN_radical = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000    -0.61104600
N        0.00000000     0.00000000     0.52375300
units angstrom
}
""")

G2_CO = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000     0.49300300
C        0.00000000     0.00000000    -0.65733700
units angstrom
}
""")

G2_CO2 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.17865800
O        0.00000000     0.00000000    -1.17865800
units angstrom
}
""")

G2_CO2_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.17344582
O        0.00000000     0.00000000    -1.17344582
units angstrom
}
""")

G2_CO_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.10256500
units angstrom
}
""")

G2_COF2 = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000     1.33071500
C        0.00000000     0.00000000     0.14435800
F        0.00000000     1.06949000    -0.63954800
F        0.00000000    -1.06949000    -0.63954800
units angstrom
}
""")

G2_CS = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -1.12338200
S        0.00000000     0.00000000     0.42126800
units angstrom
}
""")

G2_CS2 = input.process_input("""
molecule dimer {
0 1
S        0.00000000     0.00000000     1.56111700
C        0.00000000     0.00000000     0.00000000
S        0.00000000     0.00000000    -1.56111700
units angstrom
}
""")

G2_CS2_cation = input.process_input("""
molecule dimer {
1 2
S        0.00000000     0.00000000     1.54500793
C        0.00000000     0.00000000     0.00000000
S        0.00000000     0.00000000    -1.54500793
units angstrom
}
""")

G2_CS_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
S        0.00000000     0.00000000     1.45878100
units angstrom
}
""")

G2_Cyanogen = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     1.87587500
C        0.00000000     0.00000000     0.69057300
C        0.00000000     0.00000000    -0.69057300
N        0.00000000     0.00000000    -1.87587500
units angstrom
}
""")

G2_Cyanogen_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
C        1.33970000     0.00000000     0.00000000
N       -1.18470000     0.00000000     0.00000000
N        2.52440000     0.00000000    -0.00000000
units angstrom
}
""")

G2_Cyclobutane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     1.07114200     0.14762600
C        0.00000000    -1.07114200     0.14762600
C       -1.07114200     0.00000000    -0.14762600
C        1.07114200     0.00000000    -0.14762600
H        0.00000000     1.98685800    -0.45007700
H        0.00000000     1.34292100     1.20752000
H        0.00000000    -1.98685800    -0.45007700
H        0.00000000    -1.34292100     1.20752000
H       -1.98685800     0.00000000     0.45007700
H       -1.34292100     0.00000000    -1.20752000
H        1.98685800     0.00000000     0.45007700
H        1.34292100     0.00000000    -1.20752000
units angstrom
}
""")

G2_Cyclobutene = input.process_input("""
molecule dimer {
0 1
C        0.00000000    -0.67276200     0.81121700
C        0.00000000     0.67276200     0.81121700
C        0.00000000    -0.78198000    -0.69664800
C        0.00000000     0.78198000    -0.69664800
H        0.00000000    -1.42239300     1.59776300
H        0.00000000     1.42239300     1.59776300
H       -0.88931000    -1.23924200    -1.14259100
H        0.88931000    -1.23924200    -1.14259100
H        0.88931000     1.23924200    -1.14259100
H       -0.88931000     1.23924200    -1.14259100
units angstrom
}
""")

G2_Cyclopropane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.86699800     0.00000000
C        0.75084200    -0.43349900     0.00000000
C       -0.75084200    -0.43349900     0.00000000
H        0.00000000     1.45576200     0.91052600
H        0.00000000     1.45576200    -0.91052600
H        1.26072700    -0.72788100    -0.91052600
H        1.26072700    -0.72788100     0.91052600
H       -1.26072700    -0.72788100     0.91052600
H       -1.26072700    -0.72788100    -0.91052600
units angstrom
}
""")

G2_Cyclopropene = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.85829900
C        0.00000000    -0.65054500    -0.49880200
C        0.00000000     0.65054500    -0.49880200
H        0.91243800     0.00000000     1.45638700
H       -0.91243800     0.00000000     1.45638700
H        0.00000000    -1.58409800    -1.03846900
H        0.00000000     1.58409800    -1.03846900
units angstrom
}
""")

G2_Cyclopropenyl_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.80338719
C        0.00000000    -0.68296602    -0.47254244
C        0.00000000     0.68296602    -0.47254244
H        0.90070718     0.00000000     1.45514987
H       -0.90070718     0.00000000     1.45514987
H        0.00000000    -1.61264009    -1.03005682
H        0.00000000     1.61264009    -1.03005682
units angstrom
}
""")

G2_Dimethylamine = input.process_input("""
molecule dimer {
0 1
C       -0.02753000    -0.22470200     1.20488000
N       -0.02753000     0.59247000     0.00000000
C       -0.02753000    -0.22470200    -1.20488000
H        0.79150100    -0.96274200     1.24850600
H        0.03959800     0.42118200     2.08340500
H       -0.97222000    -0.77298700     1.26175000
H        0.80530300     1.17822000     0.00000000
H        0.79150100    -0.96274200    -1.24850600
H        0.03959800     0.42118200    -2.08340500
H       -0.97222000    -0.77298700    -1.26175000
units angstrom
}
""")

G2_Dimethylether = input.process_input("""
molecule dimer {
0 1
C        0.00000000     1.16572500    -0.19995000
O        0.00000000     0.00000000     0.60011000
C        0.00000000    -1.16572500    -0.19995000
H        0.00000000     2.01776900     0.48020300
H        0.89178400     1.21432000    -0.84047400
H       -0.89178400     1.21432000    -0.84047400
H        0.00000000    -2.01776900     0.48020300
H       -0.89178400    -1.21432000    -0.84047400
H        0.89178400    -1.21432000    -0.84047400
units angstrom
}
""")

G2_Dimethylsulfoxide = input.process_input("""
molecule dimer {
0 1
S        0.00000200     0.23183800    -0.43864300
O        0.00002000     1.50074200     0.37981900
C        1.33952800    -0.80902200     0.18071700
C       -1.33954800    -0.80899200     0.18071800
H        1.25583500    -0.89638500     1.26682500
H       -2.27940400    -0.31392400    -0.06867400
H        1.30440700    -1.79332700    -0.29258900
H        2.27939500    -0.31397400    -0.06867400
H       -1.30444700    -1.79329800    -0.29258700
H       -1.25585700    -0.89635500     1.26682600
units angstrom
}
""")

G2_Dimethylthioether = input.process_input("""
molecule dimer {
0 1
C        0.00000000     1.36666800    -0.51371300
S        0.00000000     0.00000000     0.66427300
C        0.00000000    -1.36666800    -0.51371300
H        0.00000000     2.29668700     0.05728400
H        0.89164400     1.34568000    -1.14459600
H       -0.89164400     1.34568000    -1.14459600
H        0.00000000    -2.29668700     0.05728400
H       -0.89164400    -1.34568000    -1.14459600
H        0.89164400    -1.34568000    -1.14459600
units angstrom
}
""")

G2_Disilane = input.process_input("""
molecule dimer {
0 1
Si       0.00000000     0.00000000     1.16768300
Si       0.00000000     0.00000000    -1.16768300
H        0.00000000     1.39328600     1.68602000
H       -1.20662100    -0.69664300     1.68602000
H        1.20662100    -0.69664300     1.68602000
H        0.00000000    -1.39328600    -1.68602000
H       -1.20662100     0.69664300    -1.68602000
H        1.20662100     0.69664300    -1.68602000
units angstrom
}
""")

G2_Disilane_cation = input.process_input("""
molecule dimer {
1 2
Si       0.00000000     0.00000000     1.32944157
Si       0.00000000     0.00000000    -1.32944157
H        0.00000000     1.45532487     1.54277503
H        0.00000000    -1.45532487    -1.54277503
H       -1.26034831    -0.72766244     1.54277503
H        1.26034831    -0.72766244     1.54277503
H       -1.26034831     0.72766244    -1.54277503
H        1.26034831     0.72766244    -1.54277503
units angstrom
}
""")

G2_Ethane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.76220900
C        0.00000000     0.00000000    -0.76220900
H        0.00000000     1.01895700     1.15722900
H       -0.88244300    -0.50947900     1.15722900
H        0.88244300    -0.50947900     1.15722900
H        0.00000000    -1.01895700    -1.15722900
H       -0.88244300     0.50947900    -1.15722900
H        0.88244300     0.50947900    -1.15722900
units angstrom
}
""")

G2_Ethanol = input.process_input("""
molecule dimer {
0 1
C        1.16818100    -0.40038200     0.00000000
C        0.00000000     0.55946200     0.00000000
O       -1.19008300    -0.22766900     0.00000000
H       -1.94662300     0.38152500     0.00000000
H        0.04255700     1.20750800     0.88693300
H        0.04255700     1.20750800    -0.88693300
H        2.11589100     0.14480000     0.00000000
H        1.12859900    -1.03723400     0.88588100
H        1.12859900    -1.03723400    -0.88588100
units angstrom
}
""")

G2_Ethanol_cation = input.process_input("""
molecule dimer {
1 2
C       -1.27755510    -0.35448093     0.13064827
C        0.17745959     0.62587690    -0.02914754
O        1.16434004    -0.25472312    -0.20495340
H       -2.02803293     0.43385286     0.18345496
H        1.62011973    -0.49472058     0.63732092
H       -1.29845462    -0.95606148    -0.77352467
H       -1.13097427    -0.91279506     1.05147367
H       -0.05751305     1.13505646    -0.96449924
H        0.18070788     1.20407692     0.89639723
units angstrom
}
""")

G2_Ethene = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.66748000
C        0.00000000     0.00000000    -0.66748000
H        0.00000000     0.92283200     1.23769500
H        0.00000000    -0.92283200     1.23769500
H        0.00000000     0.92283200    -1.23769500
H        0.00000000    -0.92283200    -1.23769500
units angstrom
}
""")

G2_Ethene_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
C        1.41632000     0.00000000     0.00000000
H       -0.55421249     0.93421938     0.00000000
H       -0.55421249    -0.93421938    -0.00000000
H        1.97053249    -0.93421938     0.00000000
H        1.97053249     0.93421938    -0.00000000
units angstrom
}
""")

G2_Ethylchloride = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.80763600     0.00000000
C        1.50582700     0.64783200     0.00000000
Cl      -0.82355300    -0.77997000     0.00000000
H       -0.34497900     1.34164900     0.88524800
H       -0.34497900     1.34164900    -0.88524800
H        1.97690300     1.63487700     0.00000000
H        1.83924600     0.10425000     0.88539800
H        1.83924600     0.10425000    -0.88539800
units angstrom
}
""")

G2_Ethyne = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.60808000
C        0.00000000     0.00000000    -0.60808000
H        0.00000000     0.00000000    -1.67399000
H        0.00000000     0.00000000     1.67399000
units angstrom
}
""")

G2_Ethyne_cation = input.process_input("""
molecule dimer {
1 2
H        0.00000000     0.00000000     0.00000000
C        1.08138200     0.00000000     0.00000000
C        2.33919900     0.00000000    -0.00000000
H        3.42058100     0.00000000    -0.00000000
units angstrom
}
""")

G2_F2 = input.process_input("""
molecule dimer {
0 1
F        0.00000000     0.00000000     0.71030400
F        0.00000000     0.00000000    -0.71030400
units angstrom
}
""")

G2_F2O = input.process_input("""
molecule dimer {
0 1
F        0.00000000     1.11057600    -0.27372900
O        0.00000000     0.00000000     0.61589000
F        0.00000000    -1.11057600    -0.27372900
units angstrom
}
""")

G2_F_anion = input.process_input("""
molecule dimer {
-1 1
F        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_F_cation = input.process_input("""
molecule dimer {
1 3
F        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_F_radical = input.process_input("""
molecule dimer {
0 2
F        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Formaldehyde = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000     0.68350100
C        0.00000000     0.00000000    -0.53661400
H        0.00000000     0.93439000    -1.12416400
H        0.00000000    -0.93439000    -1.12416400
units angstrom
}
""")

G2_FormicAcid = input.process_input("""
molecule dimer {
0 1
O       -1.04094500    -0.43643200     0.00000000
C        0.00000000     0.42394900     0.00000000
O        1.16937200     0.10374100     0.00000000
H       -0.64957000    -1.33513400     0.00000000
H       -0.37784700     1.45296700     0.00000000
units angstrom
}
""")

G2_Furan = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000     1.16333900
C        0.00000000     1.09470000     0.34803900
C        0.00000000    -1.09470000     0.34803900
C        0.00000000     0.71320000    -0.96216100
C        0.00000000    -0.71320000    -0.96216100
H        0.00000000     2.04935900     0.85111300
H        0.00000000    -2.04935900     0.85111300
H        0.00000000     1.37082800    -1.81973800
H        0.00000000    -1.37082800    -1.81973800
units angstrom
}
""")

G2_Furan_cation = input.process_input("""
molecule dimer {
1 2
O        0.00000000     0.00000000     0.00000000
C        1.34662106     0.00000000     0.00000000
C       -0.38613584     1.29007263     0.00000000
C        1.84560785     1.32657778     0.00000000
C        0.74165360     2.14849410     0.00000000
units angstrom
}
""")

G2_Glyoxal = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.75643000     0.00000000
C        0.00000000    -0.75643000     0.00000000
O        1.04609000     1.38991600     0.00000000
H       -0.99994000     1.22819100     0.00000000
O       -1.04609000    -1.38991600     0.00000000
H        0.99994000    -1.22819100     0.00000000
units angstrom
}
""")

G2_H2 = input.process_input("""
molecule dimer {
0 1
H        0.00000000     0.00000000     0.36858300
H        0.00000000     0.00000000    -0.36858300
units angstrom
}
""")

G2_H2C_anion = input.process_input("""
molecule dimer {
-1 2
C        0.00000000     0.00000000     0.00000000
H        1.14059000     0.00000000     0.00000000
H       -0.15710864     1.12971785     0.00000000
units angstrom
}
""")

G2_H2C_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.09390812
H        0.00000000     1.02331672    -0.28172436
H        0.00000000    -1.02331672    -0.28172436
units angstrom
}
""")

G2_H2C_singlet = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.17434300
H        0.00000000     0.86223200    -0.52302900
H        0.00000000    -0.86223200    -0.52302900
units angstrom
}
""")

G2_H2C_triplet = input.process_input("""
molecule dimer {
0 3
C        0.00000000     0.00000000     0.11038100
H        0.00000000     0.98262200    -0.33114200
H        0.00000000    -0.98262200    -0.33114200
units angstrom
}
""")

G2_H2CCC = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.20355795
C        0.00000000     0.00000000    -1.13011303
C        0.00000000     0.00000000     1.49496484
H        0.00000000     0.92522527    -1.70522929
H        0.00000000    -0.92522527    -1.70522929
units angstrom
}
""")

G2_H2CCC_anion = input.process_input("""
molecule dimer {
-1 2
C        0.00000000     0.00000000     0.23109742
C        0.00000000     0.00000000    -1.14336514
C        0.00000000     0.00000000     1.48817762
H        0.00000000     0.91958993    -1.72772970
H        0.00000000    -0.91958993    -1.72772970
units angstrom
}
""")

G2_H2CCCH_anion = input.process_input("""
molecule dimer {
-1 1
C       -0.11338103     0.00000000    -0.02227231
C        1.24246914     0.00000000     0.04217717
C       -1.38964290     0.00000000     0.08478279
H        1.81867488    -0.92326226     0.06899318
H        1.81867488     0.92326226     0.06899318
H       -2.07402100     0.00000000    -0.76611232
units angstrom
}
""")

G2_H2CCCH_radical = input.process_input("""
molecule dimer {
0 2
C       -0.13224389     0.00000000     0.00001297
C        1.25863724     0.00000000     0.00053148
C       -1.33048125     0.00000000    -0.00073708
H        1.81004529    -0.93008428     0.00076185
H        1.81004529     0.93008428     0.00076185
H       -2.39556318     0.00000000    -0.00036790
units angstrom
}
""")

G2_H2CCH_anion = input.process_input("""
molecule dimer {
-1 1
C        0.06694586     0.00000000    -0.57047905
C        0.06674714     0.00000000     0.78902847
H        1.01411341     0.00000000    -1.13884045
H       -0.79736875     0.00000000    -1.27040824
H       -1.01890262     0.00000000     1.09795216
units angstrom
}
""")

G2_H2CCH_radical = input.process_input("""
molecule dimer {
0 2
C        0.04979800    -0.57627200     0.00000000
C        0.04979800     0.71098800     0.00000000
H       -0.87675000    -1.15184400     0.00000000
H        0.96918300    -1.15463900     0.00000000
H       -0.69001300     1.49818500     0.00000000
units angstrom
}
""")

G2_H2CCHO_anion = input.process_input("""
molecule dimer {
-1 1
C        0.17757961    -1.16024834     0.22778254
C        0.17773353     0.21976718     0.22797997
O       -0.41448289     1.05501459    -0.53165995
H        0.82405484     0.64165137     1.05702061
H       -0.39362611    -1.72667418    -0.50490682
H        0.75355556    -1.71220695     0.96659071
units angstrom
}
""")

G2_H2CCHO_radical = input.process_input("""
molecule dimer {
0 2
C        0.19873443     0.30761603    -0.22344843
C        0.21674045    -1.14794868    -0.24369362
O       -0.47328792     0.95362861     0.53214453
H       -0.40714611    -1.68582854     0.45777754
H        0.83803177    -1.69372739    -0.94224681
H        0.86256838     0.79252299    -0.96983471
units angstrom
}
""")

G2_H2CCN_anion = input.process_input("""
molecule dimer {
-1 1
C        0.12179259    -1.21645767     0.00000000
C        0.00000000     0.17630466     0.00000000
N       -0.03199487     1.37561272     0.00000000
H       -0.25339573    -1.69418549     0.90657667
H       -0.25339573    -1.69418549    -0.90657667
units angstrom
}
""")

G2_H2CCN_radical = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000    -1.20572174
C        0.00000000     0.00000000     0.20676250
N        0.00000000     0.00000000     1.35385287
H        0.00000000     0.93667732    -1.74160730
H        0.00000000    -0.93667732    -1.74160730
units angstrom
}
""")

G2_H2Cl_cation = input.process_input("""
molecule dimer {
1 1
Cl       0.00000000     0.00000000     0.00000000
H        1.30217600     0.00000000     0.00000000
H       -0.15046684     1.29345354     0.00000000
units angstrom
}
""")

G2_H2CNC_anion = input.process_input("""
molecule dimer {
-1 1
C        0.12302436    -1.17796966     0.00000000
N        0.00000000     0.22822032     0.00000000
C        0.01265719     1.42836755     0.00000000
H       -0.40704464    -1.54996477     0.88896653
H       -0.40704464    -1.54996477    -0.88896653
units angstrom
}
""")

G2_H2CNC_radical = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000    -1.13119329
N        0.00000000     0.00000000     0.22896335
C        0.00000000     0.00000000     1.41241212
H        0.00000000     0.94793542    -1.64502826
H        0.00000000    -0.94793542    -1.64502826
units angstrom
}
""")

G2_H2COH_cation = input.process_input("""
molecule dimer {
1 1
C        0.04202796     0.62939692     0.00000000
O        0.05694409    -0.62660169     0.00000000
H       -0.89083622     1.19308632     0.00000000
H        1.01901996     1.10567005     0.00000000
H       -0.83590423    -1.06232437     0.00000000
units angstrom
}
""")

G2_H2COH_radical = input.process_input("""
molecule dimer {
0 2
C        0.68744800     0.02962600    -0.08201400
O       -0.67209400    -0.12564800     0.03040500
H       -1.09185000     0.74028200    -0.09516700
H        1.12278300     0.97526300     0.22599300
H        1.22113100    -0.88811600     0.11801500
units angstrom
}
""")

G2_H2CSH_cation = input.process_input("""
molecule dimer {
1 1
C        0.04491991     1.06518235     0.00000000
S        0.05432920    -0.55333918     0.00000000
H       -1.27629878    -0.76427991     0.00000000
H       -0.86921199     1.65363261     0.00000000
H        1.00672407     1.57298010     0.00000000
units angstrom
}
""")

G2_H2CSH_radical = input.process_input("""
molecule dimer {
0 2
C        0.15016787     1.13217860    -0.00041860
S       -0.04879566    -0.58303958     0.07121611
H       -0.22950033    -0.79221027    -1.23942783
H       -0.23173829     1.69900065    -0.83641384
H        0.34096202     1.62877129     0.93889555
units angstrom
}
""")

G2_H2O = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000     0.11926200
H        0.00000000     0.76323900    -0.47704700
H        0.00000000    -0.76323900    -0.47704700
units angstrom
}
""")

G2_H2O_cation = input.process_input("""
molecule dimer {
1 2
O        0.00000000     0.00000000     0.00000000
H        1.01033500     0.00000000     0.00000000
H       -0.34485109     0.94966022     0.00000000
units angstrom
}
""")

G2_H3_cation = input.process_input("""
molecule dimer {
1 1
H        0.00000000     0.00000000     0.00000000
H        0.85108753     0.00000000     0.00000000
H        0.42554377     0.73706342     0.00000000
units angstrom
}
""")

G2_H3O_cation = input.process_input("""
molecule dimer {
1 1
O        0.00000000     0.00000000     0.00000000
H        0.99082000     0.00000000     0.00000000
H       -0.36188688     0.92236769     0.00000000
H       -0.36188688    -0.53072856     0.75438012
units angstrom
}
""")

G2_H_radical = input.process_input("""
molecule dimer {
0 2
H        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_HCCO_anion = input.process_input("""
molecule dimer {
-1 1
C        0.00000000     0.00000000    -0.00983257
C        0.00000000     0.00000000    -1.25591173
O        0.00000000     0.00000000     1.23885717
H        0.00000000     0.00000000    -2.31639152
units angstrom
}
""")

G2_HCCO_radical = input.process_input("""
molecule dimer {
0 2
C       -0.00000000    -0.01702509    -0.01311780
C        0.00000004     0.02272529     1.26016544
O       -0.00000001     0.03284574    -1.18842432
H       -0.00000002    -0.11316596     2.31555938
units angstrom
}
""")

G2_HCF = input.process_input("""
molecule dimer {
0 1
C       -0.06642886     0.00000000    -0.72745147
F       -0.07020641     0.00000000     0.59126029
H        1.03043081     0.00000000    -0.95663377
units angstrom
}
""")

G2_HCF_anion = input.process_input("""
molecule dimer {
-1 2
C       -0.07623690     0.00000000    -0.83262070
F       -0.06665839     0.00000000     0.66612563
H        1.05734694     0.00000000    -0.99940644
units angstrom
}
""")

G2_HCl = input.process_input("""
molecule dimer {
0 1
Cl       0.00000000     0.00000000     0.07111000
H        0.00000000     0.00000000    -1.20886800
units angstrom
}
""")

G2_HCl_cation = input.process_input("""
molecule dimer {
1 2
Cl       0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.30980000
units angstrom
}
""")

G2_HCN = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -0.51174700
N        0.00000000     0.00000000     0.66446100
H        0.00000000     0.00000000    -1.58074600
units angstrom
}
""")

G2_HCO_anion = input.process_input("""
molecule dimer {
-1 1
C        0.07897248     0.64323816     0.00000000
O        0.08157896    -0.61218252     0.00000000
H       -1.12646660     1.03803117     0.00000000
units angstrom
}
""")

G2_HCO_cation = input.process_input("""
molecule dimer {
1 1
C        0.00000000     0.00000000    -0.52997659
O        0.00000000     0.00000000     0.60056506
H        0.00000000     0.00000000    -1.62466099
units angstrom
}
""")

G2_HCO_radical = input.process_input("""
molecule dimer {
0 2
C        0.06256000     0.59392600     0.00000000
O        0.06256000    -0.59691400     0.00000000
H       -0.87583500     1.21175500     0.00000000
units angstrom
}
""")

G2_He = input.process_input("""
molecule dimer {
0 1
He       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_He_cation = input.process_input("""
molecule dimer {
1 2
He       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_HF = input.process_input("""
molecule dimer {
0 1
F        0.00000000     0.00000000     0.09338900
H        0.00000000     0.00000000    -0.84050200
units angstrom
}
""")

G2_HF_cation = input.process_input("""
molecule dimer {
1 2
F        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.02964800
units angstrom
}
""")

G2_HNO = input.process_input("""
molecule dimer {
0 1
N        0.06635412     0.59787261     0.00000000
O        0.05962810    -0.63812682     0.00000000
H       -0.94150364     0.91990629     0.00000000
units angstrom
}
""")

G2_HNO_anion = input.process_input("""
molecule dimer {
-1 2
N        0.06647936     0.65447831     0.00000000
O        0.06212066    -0.69010522     0.00000000
H       -0.96232075     0.93949356     0.00000000
units angstrom
}
""")

G2_HOCl = input.process_input("""
molecule dimer {
0 1
O        0.03670200     1.11351700     0.00000000
H       -0.91754800     1.32887900     0.00000000
Cl       0.03670200    -0.60217700     0.00000000
units angstrom
}
""")

G2_HOF = input.process_input("""
molecule dimer {
0 1
O        0.06191044     0.71445274     0.00000000
F        0.04586815    -0.72913728     0.00000000
H       -0.90809687     0.84661356     0.00000000
units angstrom
}
""")

G2_HOF_cation = input.process_input("""
molecule dimer {
1 2
O        0.06732272     0.62836922     0.00000000
F        0.04263136    -0.65754465     0.00000000
H       -0.92226404     0.89094807     0.00000000
units angstrom
}
""")

G2_HOO_anion = input.process_input("""
molecule dimer {
-1 1
O       -0.06588374     0.00000000    -0.70589603
O       -0.04617814     0.00000000     0.80980467
H        0.89649510     0.00000000    -0.83126911
units angstrom
}
""")

G2_HOO_radical = input.process_input("""
molecule dimer {
0 2
O       -0.05863326     0.00000000    -0.60894669
O       -0.05283955     0.00000000     0.71641592
H        0.89178251     0.00000000    -0.85975382
units angstrom
}
""")

G2_HOOH = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.73405800    -0.05275000
O        0.00000000    -0.73405800    -0.05275000
H        0.83954700     0.88075200     0.42200100
H       -0.83954700    -0.88075200     0.42200100
units angstrom
}
""")

G2_Hydrazine = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.71895900    -0.07768700
N        0.00000000    -0.71895900    -0.07768700
H        0.21108200     1.09275200     0.84788700
H       -0.94821400     1.00502600    -0.30407800
H       -0.21108200    -1.09275200     0.84788700
H        0.94821400    -1.00502600    -0.30407800
units angstrom
}
""")

G2_Isobutane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.37694900
H        0.00000000     0.00000000     1.47526900
C        0.00000000     1.45029000    -0.09623400
H        0.00000000     1.49399700    -1.19084700
H       -0.88548200     1.98469500     0.26129700
H        0.88548200     1.98469500     0.26129700
C        1.25598800    -0.72514500    -0.09623400
H        1.29383900    -0.74699800    -1.19084700
H        2.16153700    -0.22549800     0.26129700
H        1.27605500    -1.75919800     0.26129700
C       -1.25598800    -0.72514500    -0.09623400
H       -1.29383900    -0.74699800    -1.19084700
H       -1.27605500    -1.75919800     0.26129700
H       -2.16153700    -0.22549800     0.26129700
units angstrom
}
""")

G2_Isobutene = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     1.45880700
C        0.00000000     0.00000000     0.11958800
H        0.00000000     0.92430200     2.02840900
H        0.00000000    -0.92430200     2.02840900
C        0.00000000     1.27268300    -0.67880300
H        0.00000000     2.15304200    -0.03158800
H        0.88021100     1.32354200    -1.32959200
H       -0.88021100     1.32354200    -1.32959200
C        0.00000000    -1.27268300    -0.67880300
H        0.00000000    -2.15304200    -0.03158800
H       -0.88021100    -1.32354200    -1.32959200
H        0.88021100    -1.32354200    -1.32959200
units angstrom
}
""")

G2_Isopropanol = input.process_input("""
molecule dimer {
0 1
O        0.02719100     1.36369100    -0.16751600
C       -0.00092600     0.03645900     0.37012800
H        0.85946500     1.77564700     0.12130700
H        0.00737100     0.08214500     1.47050600
C       -1.31327500    -0.56351400    -0.08897900
C        1.20072100    -0.76448000    -0.10492000
H       -1.33400500    -0.60725300    -1.18100900
H        1.20284300    -0.80781700    -1.19718900
H       -2.14781200     0.05499300     0.24767600
H        2.13646200    -0.29932400     0.22316400
H       -1.43870900    -1.57427500     0.30834000
H        1.17773600    -1.78443600     0.28996700
units angstrom
}
""")

G2_Ketene = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -1.21934000
C        0.00000000     0.00000000     0.09892000
H        0.00000000     0.93884700    -1.75322400
H        0.00000000    -0.93884700    -1.75322400
O        0.00000000     0.00000000     1.27862000
units angstrom
}
""")

G2_Li2 = input.process_input("""
molecule dimer {
0 1
Li       0.00000000     0.00000000     1.38653000
Li       0.00000000     0.00000000    -1.38653000
units angstrom
}
""")

G2_Li_anion = input.process_input("""
molecule dimer {
-1 1
Li       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Li_cation = input.process_input("""
molecule dimer {
1 1
Li       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Li_radical = input.process_input("""
molecule dimer {
0 2
Li       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_LiF = input.process_input("""
molecule dimer {
0 1
Li       0.00000000     0.00000000    -1.17496500
F        0.00000000     0.00000000     0.39165500
units angstrom
}
""")

G2_LiH = input.process_input("""
molecule dimer {
0 1
Li       0.00000000     0.00000000     0.41000000
H        0.00000000     0.00000000    -1.23000000
units angstrom
}
""")

G2_LiH_anion = input.process_input("""
molecule dimer {
-1 2
Li       0.00000000     0.00000000    -0.44012991
H        0.00000000     0.00000000     1.32038973
units angstrom
}
""")

G2_Me2CH_cation = input.process_input("""
molecule dimer {
1 1
C        0.00000000     0.00000000     0.46238641
H        0.00000000     0.00000000     1.55508714
C        0.01808188    -1.27724332    -0.19629942
C       -0.01808188     1.27724332    -0.19629942
H       -1.04610355    -1.59473319    -0.07247287
H        1.04610355     1.59473319    -0.07247287
H        0.23561147    -1.24162496    -1.26312400
H       -0.23561147     1.24162496    -1.26312400
H        0.58617938    -2.03701192     0.34869060
H       -0.58617938     2.03701192     0.34869060
units angstrom
}
""")

G2_Me2CH_radical = input.process_input("""
molecule dimer {
0 2
C        0.01422300     0.54385000     0.00000000
C        0.01422300    -0.19974200     1.29157200
C        0.01422300    -0.19974200    -1.29157200
H       -0.32289000     1.57532900     0.00000000
H        0.22141700     0.45917400     2.13847700
H        0.22141700     0.45917400    -2.13847700
H       -0.95515700    -0.68462900     1.48463300
H        0.76718100    -0.99530800     1.28623900
H        0.76718100    -0.99530800    -1.28623900
H       -0.95515700    -0.68462900    -1.48463300
units angstrom
}
""")

G2_Me3C_radical = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000     0.19192900
C        0.00000000     1.47818700    -0.02086600
C        1.28014700    -0.73909300    -0.02086600
C       -1.28014700    -0.73909300    -0.02086600
H        0.00000000     1.73149600    -1.09379200
H       -0.88704300     1.94576900     0.41756500
H        0.88704300     1.94576900     0.41756500
H        1.49952000    -0.86574800    -1.09379200
H        2.12860700    -0.20468300     0.41756500
H        1.24156400    -1.74108600     0.41756500
H       -1.49952000    -0.86574800    -1.09379200
H       -1.24156400    -1.74108600     0.41756500
H       -2.12860700    -0.20468300     0.41756500
units angstrom
}
""")

G2_MeCl = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -1.12138900
Cl       0.00000000     0.00000000     0.65595100
H        0.00000000     1.02931800    -1.47428000
H        0.89141500    -0.51465900    -1.47428000
H       -0.89141500    -0.51465900    -1.47428000
units angstrom
}
""")

G2_MeCl_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
Cl       1.76344369     0.00000000     0.00000000
H       -0.26136894     1.07451322     0.00000000
H       -0.34617271    -0.46625400    -0.92240155
H       -0.34617271    -0.46625400     0.92240155
units angstrom
}
""")

G2_MeCN = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -1.18693000
C        0.00000000     0.00000000     0.27387400
N        0.00000000     0.00000000     1.45220600
H        0.00000000     1.02498600    -1.56237000
H        0.88766400    -0.51249300    -1.56237000
H       -0.88766400    -0.51249300    -1.56237000
units angstrom
}
""")

G2_MeF = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -0.63546475
F        0.00000000     0.00000000     0.75474977
H        0.00000000     1.03138672    -0.99331981
H        0.89320710    -0.51569336    -0.99331981
H       -0.89320710    -0.51569336    -0.99331981
units angstrom
}
""")

G2_MeF_cation = input.process_input("""
molecule dimer {
1 2
C        0.01766066     0.59617828     0.00000000
F        0.05696446    -0.71700089     0.00000000
H       -1.17801636     0.82006046     0.00000000
H        0.27968612     1.02793896     0.97712866
H        0.27968612     1.02793896    -0.97712866
units angstrom
}
""")

G2_MeOF = input.process_input("""
molecule dimer {
0 1
C        0.30473301    -1.01308477     0.34830197
O        0.32535782     0.40597033     0.37187559
F       -0.59964837     0.76654085    -0.68538261
H        1.00155834    -1.28466084     1.14475534
H       -0.69510341    -1.39173044     0.57184777
H        0.65911976    -1.39173044    -0.61297615
units angstrom
}
""")

G2_MeOF_cation = input.process_input("""
molecule dimer {
1 2
C        0.28601767    -1.07137862     0.30353160
O        0.29579027     0.37257114     0.31390261
F       -0.54709347     0.83567608    -0.58059405
H        1.01264262    -1.32486519     1.07465051
H       -0.74522057    -1.37425842     0.53068889
H        0.57399100    -1.37425842    -0.71240349
units angstrom
}
""")

G2_MeOH = input.process_input("""
molecule dimer {
0 1
C       -0.04713100     0.66438900     0.00000000
O       -0.04713100    -0.75855100     0.00000000
H       -1.09299500     0.96978500     0.00000000
H        0.87853400    -1.04845800     0.00000000
H        0.43714500     1.08037600     0.89177200
H        0.43714500     1.08037600    -0.89177200
units angstrom
}
""")

G2_MeOH_cation = input.process_input("""
molecule dimer {
1 2
C       -0.65846995     0.04286789    -0.01889908
O        0.71345467    -0.12908510    -0.01148190
H       -1.08431476    -0.77325594    -0.61536758
H        1.23301665     0.72250885     0.03164603
H       -0.92473360    -0.22904599     1.04201453
H       -0.98078593     1.05526652    -0.25304326
units angstrom
}
""")

G2_MeSH = input.process_input("""
molecule dimer {
0 1
C       -0.04795300     1.14951900     0.00000000
S       -0.04795300    -0.66485600     0.00000000
H        1.28307600    -0.82324900     0.00000000
H       -1.09260100     1.46142800     0.00000000
H        0.43224900     1.55120700     0.89225900
H        0.43224900     1.55120700    -0.89225900
units angstrom
}
""")

G2_MeSH_cation = input.process_input("""
molecule dimer {
1 2
C       -0.04651116     1.13532092     0.00000000
S       -0.05508984    -0.65086427     0.00000000
H       -1.08039657     1.47742687     0.00000000
H        1.27782168    -0.85826440     0.00000000
H        0.48153961     1.49137019     0.89158515
H        0.48153961     1.49137019    -0.89158515
units angstrom
}
""")

G2_Methane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.00000000
H        0.62911800     0.62911800     0.62911800
H       -0.62911800    -0.62911800     0.62911800
H        0.62911800    -0.62911800    -0.62911800
H       -0.62911800     0.62911800    -0.62911800
units angstrom
}
""")

G2_Methane_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000     0.00000000
H        1.02338607     0.57384238     0.00000000
H        1.02338607    -0.57384238    -0.00000000
H       -0.51334464     0.00000000    -0.95458697
H       -0.51334464     0.00000000     0.95458697
units angstrom
}
""")

G2_Methyl_anion = input.process_input("""
molecule dimer {
-1 1
C        0.00000000     0.00000000     0.00000000
H       -0.50094623     1.00271770     0.00000000
H       -0.50094623    -0.50135885    -0.86837900
H       -0.50094623    -0.50135885     0.86837900
units angstrom
}
""")

G2_Methyl_cation = input.process_input("""
molecule dimer {
1 1
C        0.00000000     0.00000000     0.00000000
H        0.00000000     1.08827618     0.00000000
H        0.94247482    -0.54413809     0.00000000
H       -0.94247482    -0.54413809     0.00000000
units angstrom
}
""")

G2_Methyl_radical = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000     0.00000000
H        0.00000000     1.07841000     0.00000000
H        0.93393000    -0.53920500     0.00000000
H       -0.93393000    -0.53920500     0.00000000
units angstrom
}
""")

G2_Methylamine = input.process_input("""
molecule dimer {
0 1
C        0.05173600     0.70442200     0.00000000
N        0.05173600    -0.75961600     0.00000000
H       -0.94173500     1.17619200     0.00000000
H       -0.45818100    -1.09943300     0.81237000
H       -0.45818100    -1.09943300    -0.81237000
H        0.59276300     1.05672700     0.88067000
H        0.59276300     1.05672700    -0.88067000
units angstrom
}
""")

G2_Methylenecyclopropane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.31502600
C        0.00000000    -0.76792000    -0.93203200
C        0.00000000     0.76792000    -0.93203200
C        0.00000000     0.00000000     1.64002700
H       -0.91279400    -1.27178900    -1.23930300
H        0.91279400    -1.27178900    -1.23930300
H        0.91279400     1.27178900    -1.23930300
H       -0.91279400     1.27178900    -1.23930300
H        0.00000000    -0.92690800     2.20564000
H        0.00000000     0.92690800     2.20564000
units angstrom
}
""")

G2_Methylethylether = input.process_input("""
molecule dimer {
0 1
O        0.00642900    -0.71274100     0.00000000
C        0.00000000     0.70584500     0.00000000
C        1.32451800    -1.22602900     0.00000000
C       -1.44216900     1.16032500     0.00000000
H        0.53096200     1.08648400     0.88688100
H        0.53096200     1.08648400    -0.88688100
H        1.24164800    -2.31332500     0.00000000
H        1.88132900    -0.90592500    -0.89171000
H        1.88132900    -0.90592500     0.89171000
H       -1.95486300     0.78060500    -0.88585500
H       -1.95486300     0.78060500     0.88585500
H       -1.50202500     2.25208300     0.00000000
units angstrom
}
""")

G2_Methylformate = input.process_input("""
molecule dimer {
0 1
C       -0.93120900    -0.08386600     0.00000000
O       -0.71101900    -1.27820900     0.00000000
O        0.00000000     0.88684100     0.00000000
H       -1.92836000     0.37459800     0.00000000
C        1.35689900     0.39728700     0.00000000
H        1.98013400     1.28816400     0.00000000
H        1.54112100    -0.20617200     0.88939700
H        1.54112100    -0.20617200    -0.88939700
units angstrom
}
""")

G2_Methylnitrite = input.process_input("""
molecule dimer {
0 1
C       -1.31620800     0.30924700     0.00000000
O        0.00000000     0.89685200     0.00000000
H       -1.98553800     1.16601300     0.00000000
H       -1.46433600    -0.30463700     0.89067200
H       -1.46433600    -0.30463700    -0.89067200
N        1.04533400    -0.02281500     0.00000000
O        0.68676400    -1.17841600     0.00000000
units angstrom
}
""")

G2_Methylsilane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -1.24446600
Si       0.00000000     0.00000000     0.63570300
H        0.00000000    -1.01976200    -1.63636300
H       -0.88314000     0.50988100    -1.63636300
H        0.88314000     0.50988100    -1.63636300
H        0.00000000     1.39123400     1.15868200
H       -1.20484400    -0.69561700     1.15868200
H        1.20484400    -0.69561700     1.15868200
units angstrom
}
""")

G2_Mg = input.process_input("""
molecule dimer {
0 1
Mg       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Mg_cation = input.process_input("""
molecule dimer {
1 2
Mg       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_N2 = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     0.56499000
N        0.00000000     0.00000000    -0.56499000
units angstrom
}
""")

G2_N2_cation_2PIu = input.process_input("""
molecule dimer {
1 2
N        0.00000000     0.00000000     0.00000000
N        0.00000000     0.00000000     1.20284500
units angstrom
}
""")

G2_N2_cation_2SIGMAg = input.process_input("""
molecule dimer {
1 2
N        0.00000000     0.00000000     0.00000000
N        0.00000000     0.00000000     1.14657200
units angstrom
}
""")

G2_N2H2 = input.process_input("""
molecule dimer {
0 1
N       -0.00975033    -0.63297235     0.00000000
N        0.00975033     0.63297235     0.00000000
H        0.98467144    -0.92398353     0.00000000
H       -0.98467144     0.92398353     0.00000000
units angstrom
}
""")

G2_N2H2_cation = input.process_input("""
molecule dimer {
1 2
N       -0.00004969     0.00000000    -0.59290896
N        0.00004969     0.00000000     0.59290896
H        0.87096664     0.00000000    -1.16977136
H       -0.87096664     0.00000000     1.16977136
units angstrom
}
""")

G2_N2H3_cation = input.process_input("""
molecule dimer {
1 1
N        0.04639114    -0.54270550     0.00000000
N        0.06062446     0.69678269     0.00000000
H       -0.80622572    -1.13004319     0.00000000
H        0.96248895    -1.01715357     0.00000000
H       -0.90537243     1.06865640     0.00000000
units angstrom
}
""")

G2_N2H3_radical = input.process_input("""
molecule dimer {
0 2
N        0.09044714     0.59626812    -0.01441033
N       -0.00881478    -0.74612571     0.11418998
H       -0.19584965    -1.09429590    -0.83230609
H        0.05664631     1.07092832     0.87883499
H       -0.43222319     1.07237066    -0.74498645
units angstrom
}
""")

G2_N2O = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000    -1.23196900
N        0.00000000     0.00000000    -0.06085100
O        0.00000000     0.00000000     1.13121800
units angstrom
}
""")

G2_N_cation = input.process_input("""
molecule dimer {
1 3
N        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_N_radical = input.process_input("""
molecule dimer {
0 4
N        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Na2 = input.process_input("""
molecule dimer {
0 1
Na       0.00000000     0.00000000     1.57626200
Na       0.00000000     0.00000000    -1.57626200
units angstrom
}
""")

G2_Na_anion = input.process_input("""
molecule dimer {
-1 1
Na       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Na_cation = input.process_input("""
molecule dimer {
1 1
Na       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Na_radical = input.process_input("""
molecule dimer {
0 2
Na       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_NaCl = input.process_input("""
molecule dimer {
0 1
Na       0.00000000     0.00000000    -1.45166000
Cl       0.00000000     0.00000000     0.93931000
units angstrom
}
""")

G2_NCO_anion = input.process_input("""
molecule dimer {
-1 1
C        0.06956384     0.00000000     0.00000000
N        1.27930575     0.00000000     0.00000000
O       -1.17156541     0.00000000     0.00000000
units angstrom
}
""")

G2_NCO_radical = input.process_input("""
molecule dimer {
0 2
C        0.02649011     0.00000000     0.00000000
N        1.27950532     0.00000000     0.00000000
O       -1.13943473     0.00000000     0.00000000
units angstrom
}
""")

G2_Ne = input.process_input("""
molecule dimer {
0 1
Ne       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Ne_cation = input.process_input("""
molecule dimer {
1 2
Ne       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_NF3 = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     0.48967200
F        0.00000000     1.23821800    -0.12695200
F        1.07232800    -0.61910900    -0.12695200
F       -1.07232800    -0.61910900    -0.12695200
units angstrom
}
""")

G2_NH = input.process_input("""
molecule dimer {
0 3
N        0.00000000     0.00000000     0.12992900
H        0.00000000     0.00000000    -0.90950100
units angstrom
}
""")

G2_NH2_anion = input.process_input("""
molecule dimer {
-1 1
N        0.00000000     0.00000000     0.00000000
H        1.04263200     0.00000000     0.00000000
H       -0.15289436     1.03136066     0.00000000
units angstrom
}
""")

G2_NH2_cation = input.process_input("""
molecule dimer {
1 3
N       -0.05533633     0.00000000    -0.01422153
H       -0.05544990     0.00000000     1.01913509
H        0.44280424     0.00000000    -0.91958438
units angstrom
}
""")

G2_NH2_radical = input.process_input("""
molecule dimer {
0 2
N        0.00000000     0.00000000     0.14169000
H        0.00000000     0.80644200    -0.49591300
H        0.00000000    -0.80644200    -0.49591300
units angstrom
}
""")

G2_NH3 = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     0.11648900
H        0.00000000     0.93973100    -0.27180800
H        0.81383100    -0.46986500    -0.27180800
H       -0.81383100    -0.46986500    -0.27180800
units angstrom
}
""")

G2_NH3_cation = input.process_input("""
molecule dimer {
1 2
N        0.00000000     0.00000000     0.00000000
H        1.02598500     0.00000000     0.00000000
H       -0.51299250     0.88852907     0.00000000
H       -0.51299250    -0.88852907    -0.00000000
units angstrom
}
""")

G2_NH4_cation = input.process_input("""
molecule dimer {
1 1
N        0.00000000     0.00000000     0.00000000
H        1.02853000     0.00000000     0.00000000
H       -0.34284332     0.96970739     0.00000000
H       -0.34284332    -0.48485367     0.83979124
H       -0.34284332    -0.48485367    -0.83979124
units angstrom
}
""")

G2_NH_anion = input.process_input("""
molecule dimer {
-1 2
N        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.05745100
units angstrom
}
""")

G2_NH_cation = input.process_input("""
molecule dimer {
1 2
N        0.00000000     0.00000000     0.13295316
H        0.00000000     0.00000000    -0.93067213
units angstrom
}
""")

G2_Nitromethane = input.process_input("""
molecule dimer {
0 1
C       -0.11428200    -1.31456500     0.00000000
N        0.00000000     0.16648000     0.00000000
H        0.89956500    -1.71525600     0.00000000
H       -0.64092100    -1.60721200     0.90495600
H       -0.64092100    -1.60721200    -0.90495600
O        0.06674800     0.72823200    -1.10377500
O        0.06674800     0.72823200     1.10377500
units angstrom
}
""")

G2_NO2_anion = input.process_input("""
molecule dimer {
-1 1
N       -0.43818164     0.00000000    -0.18718380
O       -0.23350772     0.00000000     1.07727914
O        0.61691665     0.00000000    -0.91349332
units angstrom
}
""")

G2_NO2_radical = input.process_input("""
molecule dimer {
0 2
N        0.00000000     0.00000000     0.33227300
O        0.00000000     1.11812200    -0.14537000
O        0.00000000    -1.11812200    -0.14537000
units angstrom
}
""")

G2_NO_anion = input.process_input("""
molecule dimer {
-1 3
N        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.28013000
units angstrom
}
""")

G2_NO_radical = input.process_input("""
molecule dimer {
0 2
N        0.00000000     0.00000000    -0.60944200
O        0.00000000     0.00000000     0.53326100
units angstrom
}
""")

G2_O = input.process_input("""
molecule dimer {
0 3
O        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_O2 = input.process_input("""
molecule dimer {
0 3
O        0.00000000     0.00000000     0.62297800
O        0.00000000     0.00000000    -0.62297800
units angstrom
}
""")

G2_O2_anion = input.process_input("""
molecule dimer {
-1 2
O        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.37974500
units angstrom
}
""")

G2_O2_cation = input.process_input("""
molecule dimer {
1 2
O        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.18689800
units angstrom
}
""")

G2_O3 = input.process_input("""
molecule dimer {
0 1
O        0.00000000     1.10381000    -0.22854200
O        0.00000000     0.00000000     0.45708400
O        0.00000000    -1.10381000    -0.22854200
units angstrom
}
""")

G2_O3_anion = input.process_input("""
molecule dimer {
-1 2
O        0.00000000     0.00000000     0.00000000
O        1.33689965     0.00000000     0.00000000
O       -0.55967183     1.21411207     0.00000000
units angstrom
}
""")

G2_O_anion = input.process_input("""
molecule dimer {
-1 2
O        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_O_cation = input.process_input("""
molecule dimer {
1 4
O        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_OCS = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000    -1.69924300
C        0.00000000     0.00000000    -0.52049200
S        0.00000000     0.00000000     1.04480600
units angstrom
}
""")

G2_OCS_cation = input.process_input("""
molecule dimer {
1 2
O        0.00000000     0.00000000    -1.72239961
C        0.00000000     0.00000000    -0.59652704
S        0.00000000     0.00000000     1.08489745
units angstrom
}
""")

G2_OF_anion = input.process_input("""
molecule dimer {
-1 1
O        0.00000000     0.00000000    -0.79316756
F        0.00000000     0.00000000     0.70503784
units angstrom
}
""")

G2_OF_radical = input.process_input("""
molecule dimer {
0 2
O        0.00000000     0.00000000    -0.71131831
F        0.00000000     0.00000000     0.63228295
units angstrom
}
""")

G2_OH_anion = input.process_input("""
molecule dimer {
-1 1
O        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     0.98034800
units angstrom
}
""")

G2_OH_cation = input.process_input("""
molecule dimer {
1 3
H        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.03458100
units angstrom
}
""")

G2_OH_radical = input.process_input("""
molecule dimer {
0 2
O        0.00000000     0.00000000     0.10878600
H        0.00000000     0.00000000    -0.87028400
units angstrom
}
""")

G2_Oxirane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.73158000    -0.37567400
O        0.00000000     0.00000000     0.86095000
C        0.00000000    -0.73158000    -0.37567400
H        0.91956800     1.26882100    -0.59487800
H       -0.91956800     1.26882100    -0.59487800
H       -0.91956800    -1.26882100    -0.59487800
H        0.91956800    -1.26882100    -0.59487800
units angstrom
}
""")

G2_P2 = input.process_input("""
molecule dimer {
0 1
P        0.00000000     0.00000000     0.96614400
P        0.00000000     0.00000000    -0.96614400
units angstrom
}
""")

G2_P2_cation = input.process_input("""
molecule dimer {
1 2
P        0.00000000     0.00000000     0.00000000
P        0.00000000     0.00000000     2.00985700
units angstrom
}
""")

G2_P_anion = input.process_input("""
molecule dimer {
-1 3
P        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_P_cation = input.process_input("""
molecule dimer {
1 3
P        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_P_radical = input.process_input("""
molecule dimer {
0 4
P        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_PF3 = input.process_input("""
molecule dimer {
0 1
P        0.00000000     0.00000000     0.50676700
F        0.00000000     1.38386100    -0.28153700
F        1.19845900    -0.69193100    -0.28153700
F       -1.19845900    -0.69193100    -0.28153700
units angstrom
}
""")

G2_PH = input.process_input("""
molecule dimer {
0 3
P        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.42563500
units angstrom
}
""")

G2_PH2_anion = input.process_input("""
molecule dimer {
-1 1
P        0.00000000     0.00000000     0.00000000
H        0.99985462     1.03186776     0.00000000
H        0.99985462    -1.03186776    -0.00000000
units angstrom
}
""")

G2_PH2_cation = input.process_input("""
molecule dimer {
1 1
P        0.00000000     0.00000000     0.00000000
H        1.41543800     0.00000000     0.00000000
H       -0.08420138     1.41293130     0.00000000
units angstrom
}
""")

G2_PH2_radical = input.process_input("""
molecule dimer {
0 2
P        0.00000000     0.00000000     0.11539600
H        0.00000000     1.02564200    -0.86546800
H        0.00000000    -1.02564200    -0.86546800
units angstrom
}
""")

G2_PH3 = input.process_input("""
molecule dimer {
0 1
P        0.00000000     0.00000000     0.12461900
H        0.00000000     1.20064700    -0.62309500
H        1.03979100    -0.60032300    -0.62309500
H       -1.03979100    -0.60032300    -0.62309500
units angstrom
}
""")

G2_PH3_cation = input.process_input("""
molecule dimer {
1 2
P        0.00000000     0.00000000     0.00000000
H       -0.36803084     1.34567787     0.00000000
H       -0.36803084    -0.67283894    -1.16539122
H       -0.36803084    -0.67283894     1.16539122
units angstrom
}
""")

G2_PH4_cation = input.process_input("""
molecule dimer {
1 1
P        0.00000000     0.00000000     0.00000000
H        1.39292000     0.00000000     0.00000000
H       -0.46430668     1.31325757     0.00000000
H       -0.46430668    -0.65662878    -1.13731441
H       -0.46430668    -0.65662878     1.13731441
units angstrom
}
""")

G2_PH_anion = input.process_input("""
molecule dimer {
-1 2
P        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.44571400
units angstrom
}
""")

G2_PH_cation = input.process_input("""
molecule dimer {
1 2
P        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.42000000
units angstrom
}
""")

G2_Phenol = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.93714600     0.00000000
C       -1.20585300     0.23463800     0.00000000
C       -1.19204000    -1.16037200     0.00000000
C        0.01613400    -1.85456200     0.00000000
C        1.21705800    -1.14180000     0.00000000
C        1.21522200     0.24983700     0.00000000
O        0.06273500     2.30930300     0.00000000
H       -0.84570300     2.65841400     0.00000000
H       -2.15261800     0.77304100     0.00000000
H       -2.13442700    -1.70220300     0.00000000
H        0.02324400    -2.94075300     0.00000000
H        2.16546100    -1.67304600     0.00000000
H        2.13903000     0.82081000     0.00000000
units angstrom
}
""")

G2_Phenol_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.92495500     0.00000000
C        1.24598600     0.26736000     0.00000000
C        1.25209800    -1.07968300     0.00000000
C        0.01559700    -1.81397000     0.00000000
C       -1.21847400    -1.15824500     0.00000000
C       -1.23163500     0.19287300     0.00000000
O       -0.15124700     2.24466400     0.00000000
H        0.71017100     2.71546400     0.00000000
H        2.17500200     0.83102200     0.00000000
H        2.18542600    -1.63493100     0.00000000
H        0.05432600    -2.89963200     0.00000000
H       -2.13940500    -1.73096100     0.00000000
H       -2.15697900     0.76198400     0.00000000
units angstrom
}
""")

G2_PO_anion = input.process_input("""
molecule dimer {
-1 3
P        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.58066300
units angstrom
}
""")

G2_PO_radical = input.process_input("""
molecule dimer {
0 2
P        0.00000000     0.00000000     0.00000000
O        0.00000000     0.00000000     1.53771300
units angstrom
}
""")

G2_Propane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.58771600
C        0.00000000     1.26685700    -0.26018600
C        0.00000000    -1.26685700    -0.26018600
H       -0.87689800     0.00000000     1.24471300
H        0.87689800     0.00000000     1.24471300
H        0.00000000     2.16615000     0.36206600
H        0.00000000    -2.16615000     0.36206600
H        0.88361900     1.30423400    -0.90440500
H       -0.88361900     1.30423400    -0.90440500
H       -0.88361900    -1.30423400    -0.90440500
H        0.88361900    -1.30423400    -0.90440500
units angstrom
}
""")

G2_Propene = input.process_input("""
molecule dimer {
0 1
C        1.29129000     0.13368200     0.00000000
C        0.00000000     0.47915900     0.00000000
H        1.60116000    -0.90742000     0.00000000
H        2.08080000     0.87733700     0.00000000
H       -0.26322100     1.53609800     0.00000000
C       -1.13975700    -0.49234100     0.00000000
H       -0.77685900    -1.52329100     0.00000000
H       -1.77554000    -0.35286100     0.88042000
H       -1.77554000    -0.35286100    -0.88042000
units angstrom
}
""")

G2_Propylchloride = input.process_input("""
molecule dimer {
0 1
C        0.89262900    -0.64234400     0.00000000
C        2.36558700    -0.24516800     0.00000000
C        0.00000000     0.58292100     0.00000000
H        0.66373100    -1.25211700     0.87920100
H        0.66373100    -1.25211700    -0.87920100
H        3.00547600    -1.13092400     0.00000000
Cl      -1.73281000     0.13974300     0.00000000
H        2.61488200     0.34770400    -0.88473000
H        2.61488200     0.34770400     0.88473000
H        0.17288100     1.19583600     0.88646000
H        0.17288100     1.19583600    -0.88646000
units angstrom
}
""")

G2_Propyne = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.21494700
C        0.00000000     0.00000000     1.43313000
C        0.00000000     0.00000000    -1.24647600
H        0.00000000     0.00000000     2.49888700
H        0.00000000     1.02114500    -1.63616700
H        0.88433700    -0.51057200    -1.63616700
H       -0.88433700    -0.51057200    -1.63616700
units angstrom
}
""")

G2_Pyridine = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     1.42467200
C        0.00000000     0.00000000    -1.38617800
C        0.00000000     1.14427700     0.72030600
C        0.00000000    -1.14427700     0.72030600
C        0.00000000    -1.19640400    -0.67291700
C        0.00000000     1.19640400    -0.67291700
H        0.00000000     0.00000000    -2.47305200
H        0.00000000     2.06072300     1.30747700
H        0.00000000    -2.06072300     1.30747700
H        0.00000000    -2.15529300    -1.18310300
H        0.00000000     2.15529300    -1.18310300
units angstrom
}
""")

G2_Pyrrole = input.process_input("""
molecule dimer {
0 1
H        0.00000000     0.00000000     2.12929600
N        0.00000000     0.00000000     1.11868400
C        0.00000000     1.12451600     0.33356500
C        0.00000000    -1.12451600     0.33356500
C        0.00000000     0.70840700    -0.98380700
C        0.00000000    -0.70840700    -0.98380700
H        0.00000000     2.11287200     0.77049600
H        0.00000000    -2.11287200     0.77049600
H        0.00000000     1.35725200    -1.84908500
H        0.00000000    -1.35725200    -1.84908500
units angstrom
}
""")

G2_Pyrrole_cation = input.process_input("""
molecule dimer {
1 2
H        0.00000000     0.00000000     0.00000000
N        1.01807805     0.00000000     0.00000000
C        1.80467644     1.10768558    -0.00000000
C        1.80467644    -1.10768558     0.00000000
C        3.17132832     0.68721929    -0.00000000
C        3.17132832    -0.68721929    -0.00000000
H        1.38237969     2.10515126    -0.00000000
H        1.38237969    -2.10515126     0.00000000
H        4.02343369     1.35368728    -0.00000000
H        4.02343369    -1.35368728    -0.00000000
units angstrom
}
""")

G2_S = input.process_input("""
molecule dimer {
0 3
S        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_S2 = input.process_input("""
molecule dimer {
0 3
S        0.00000000     0.00000000     0.96011300
S        0.00000000     0.00000000    -0.96011300
units angstrom
}
""")

G2_S2_anion = input.process_input("""
molecule dimer {
-1 2
S        0.00000000     0.00000000     0.00000000
S        0.00000000     0.00000000     2.02918900
units angstrom
}
""")

G2_S2_cation = input.process_input("""
molecule dimer {
1 2
S        0.00000000     0.00000000     0.00000000
S        0.00000000     0.00000000     1.86912000
units angstrom
}
""")

G2_S2O = input.process_input("""
molecule dimer {
0 1
S        0.00000000     0.00000000     0.00000000
S        1.91668572     0.00000000     0.00000000
O       -0.70383898     1.32955542     0.00000000
units angstrom
}
""")

G2_S2O_anion = input.process_input("""
molecule dimer {
-1 2
S        0.00000000     0.00000000     0.00000000
S        1.97970430     0.00000000     0.00000000
O       -0.69690029     1.36154184     0.00000000
units angstrom
}
""")

G2_S_anion = input.process_input("""
molecule dimer {
-1 2
S        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_S_cation = input.process_input("""
molecule dimer {
1 4
S        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_SH2 = input.process_input("""
molecule dimer {
0 1
S        0.00000000     0.00000000     0.10213500
H        0.00000000     0.97426900    -0.81708300
H        0.00000000    -0.97426900    -0.81708300
units angstrom
}
""")

G2_SH2_cation_2A1 = input.process_input("""
molecule dimer {
1 2
S        0.00000000     0.00000000     0.00000000
H        1.35315800     0.00000000     0.00000000
H       -0.80688145     1.08626834     0.00000000
units angstrom
}
""")

G2_SH2_cation_2B1 = input.process_input("""
molecule dimer {
1 2
S        0.00000000     0.00000000     0.00000000
H        1.35079300     0.00000000     0.00000000
H       -0.11124606     1.34620431     0.00000000
units angstrom
}
""")

G2_SH3_cation = input.process_input("""
molecule dimer {
1 1
S        0.00000000     0.00000000     0.00000000
H       -0.69077431     1.15635360     0.00000000
H       -0.69077431    -0.57817680    -1.00143159
H       -0.69077431    -0.57817680     1.00143159
units angstrom
}
""")

G2_SH_anion = input.process_input("""
molecule dimer {
-1 1
S        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.35334200
units angstrom
}
""")

G2_SH_cation = input.process_input("""
molecule dimer {
1 3
S        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.35666400
units angstrom
}
""")

G2_SH_radical = input.process_input("""
molecule dimer {
0 2
S        0.00000000     0.00000000     0.07908300
H        0.00000000     0.00000000    -1.26533000
units angstrom
}
""")

G2_Si = input.process_input("""
molecule dimer {
0 3
Si       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Si2 = input.process_input("""
molecule dimer {
0 3
Si       0.00000000     0.00000000     1.13005400
Si       0.00000000     0.00000000    -1.13005400
units angstrom
}
""")

G2_Si2H2 = input.process_input("""
molecule dimer {
0 1
Si       0.00000000    -1.10090263    -0.05106826
Si       0.00000000     1.10090263    -0.05106826
H        0.99022426     0.00000000     0.71495558
H       -0.99022426     0.00000000     0.71495558
units angstrom
}
""")

G2_Si2H2_cation = input.process_input("""
molecule dimer {
1 2
Si       0.00000000    -1.10734754    -0.05303623
Si       0.00000000     1.10734754    -0.05303623
H        0.98436532     0.00000000     0.74250727
H       -0.98436532     0.00000000     0.74250727
units angstrom
}
""")

G2_Si2H4 = input.process_input("""
molecule dimer {
0 1
Si      -0.11087467    -1.07619651     0.00000000
Si       0.11087467     1.07619651     0.00000000
H        0.20995958    -1.82970509    -1.23530661
H        0.20995958    -1.82970509     1.23530661
H       -0.20995958     1.82970509    -1.23530661
H       -0.20995958     1.82970509     1.23530661
units angstrom
}
""")

G2_Si2H4_cation = input.process_input("""
molecule dimer {
1 2
Si       0.00000000     0.00000000     1.11634468
Si       0.00000000     0.00000000    -1.11634468
H        0.00000000     1.27334707     1.84910344
H        0.00000000    -1.27334707     1.84910344
H        0.00000000     1.27334707    -1.84910344
H        0.00000000    -1.27334707    -1.84910344
units angstrom
}
""")

G2_Si2H5 = input.process_input("""
molecule dimer {
0 2
Si      -0.03757885     0.00000000    -1.11751590
Si      -0.03224950     0.00000000     1.20811067
H        1.33764501     0.00000000    -1.68519426
H       -0.74955969    -1.20867191    -1.60635859
H       -0.74955969     1.20867191    -1.60635859
H        0.56953558     1.21651442     1.81479233
H        0.56953558    -1.21651442     1.81479233
units angstrom
}
""")

G2_Si2H5_cation = input.process_input("""
molecule dimer {
1 1
Si      -0.00520479     0.00000000    -1.16733535
Si      -0.00139822     0.00000000     1.20169631
H        1.44354959     0.00000000    -1.45956642
H       -0.69922946    -1.24610706    -1.53641389
H       -0.69922946     1.24610706    -1.53641389
H        0.02367574     1.22338859     2.02567038
H        0.02367574    -1.22338859     2.02567038
units angstrom
}
""")

G2_Si_anion = input.process_input("""
molecule dimer {
-1 4
Si       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_Si_cation = input.process_input("""
molecule dimer {
1 2
Si       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

G2_SiCl4 = input.process_input("""
molecule dimer {
0 1
Si       0.00000000     0.00000000     0.00000000
Cl       1.16934900     1.16934900     1.16934900
Cl      -1.16934900    -1.16934900     1.16934900
Cl       1.16934900    -1.16934900    -1.16934900
Cl      -1.16934900     1.16934900    -1.16934900
units angstrom
}
""")

G2_SiF4 = input.process_input("""
molecule dimer {
0 1
Si       0.00000000     0.00000000     0.00000000
F        0.91280600     0.91280600     0.91280600
F       -0.91280600    -0.91280600     0.91280600
F       -0.91280600     0.91280600    -0.91280600
F        0.91280600    -0.91280600    -0.91280600
units angstrom
}
""")

G2_SiH2_anion = input.process_input("""
molecule dimer {
-1 2
Si       0.00000000     0.00000000     0.00000000
H        1.07327087     1.11873975     0.00000000
H        1.07327087    -1.11873975    -0.00000000
units angstrom
}
""")

G2_SiH2_cation = input.process_input("""
molecule dimer {
1 2
Si       0.00000000     0.00000000     0.09232471
H        0.00000000     1.28084868    -0.64627299
H        0.00000000    -1.28084868    -0.64627299
units angstrom
}
""")

G2_SiH2_singlet = input.process_input("""
molecule dimer {
0 1
Si       0.00000000     0.00000000     0.13127200
H        0.00000000     1.09693800    -0.91890500
H        0.00000000    -1.09693800    -0.91890500
units angstrom
}
""")

G2_SiH2_triplet = input.process_input("""
molecule dimer {
0 3
Si       0.00000000     0.00000000     0.09486900
H        0.00000000     1.27186200    -0.66408300
H        0.00000000    -1.27186200    -0.66408300
units angstrom
}
""")

G2_SiH3_anion = input.process_input("""
molecule dimer {
-1 1
Si       0.00000000     0.00000000     0.00000000
H       -0.81185801     1.31138715     0.00000000
H       -0.81185801    -0.65569358    -1.13569459
H       -0.81185801    -0.65569358     1.13569459
units angstrom
}
""")

G2_SiH3_cation = input.process_input("""
molecule dimer {
1 1
Si       0.00000000     0.00000000     0.00000128
H        1.46484612     0.00000000    -0.00000599
H       -0.73242306    -1.26859395    -0.00000599
H       -0.73242306     1.26859395    -0.00000599
units angstrom
}
""")

G2_SiH3_radical = input.process_input("""
molecule dimer {
0 2
Si       0.00000000     0.00000000     0.07929900
H        0.00000000     1.41328000    -0.37006100
H        1.22393700    -0.70664000    -0.37006100
H       -1.22393700    -0.70664000    -0.37006100
units angstrom
}
""")

G2_SiH4 = input.process_input("""
molecule dimer {
0 1
Si       0.00000000     0.00000000     0.00000000
H        0.85613500     0.85613500     0.85613500
H       -0.85613500    -0.85613500     0.85613500
H       -0.85613500     0.85613500    -0.85613500
H        0.85613500    -0.85613500    -0.85613500
units angstrom
}
""")

G2_SiH4_cation = input.process_input("""
molecule dimer {
1 2
Si       0.00000000     0.00000000     0.00000000
H        0.72095241     1.28122568     0.00000000
H        0.72095241    -1.28122568    -0.00000000
H        0.04972279     0.00000000    -1.88980298
H       -0.69387845     0.00000000    -1.69890303
units angstrom
}
""")

G2_SiH5_cation = input.process_input("""
molecule dimer {
1 1
Si       0.00000000     0.00000000     0.00000000
H        1.83135683     0.38053894     0.00000000
H        1.83047766    -0.38038045     0.00000000
H       -0.16190003     0.02577135    -1.45374758
H       -0.10782357    -1.27405826     0.71377205
H       -0.10956073     1.24796593     0.75809888
units angstrom
}
""")

G2_SiH_anion = input.process_input("""
molecule dimer {
-1 3
Si       0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.56083600
units angstrom
}
""")

G2_SiH_radical = input.process_input("""
molecule dimer {
0 2
Si       0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.52625900
units angstrom
}
""")

G2_SiO = input.process_input("""
molecule dimer {
0 1
Si       0.00000000     0.00000000     0.56084600
O        0.00000000     0.00000000    -0.98148000
units angstrom
}
""")

G2_SO = input.process_input("""
molecule dimer {
0 3
O        0.00000000     0.00000000    -1.01599200
S        0.00000000     0.00000000     0.50799600
units angstrom
}
""")

G2_SO2 = input.process_input("""
molecule dimer {
0 1
S        0.00000000     0.00000000     0.37026800
O        0.00000000     1.27761700    -0.37026800
O        0.00000000    -1.27761700    -0.37026800
units angstrom
}
""")

G2_SO2_anion = input.process_input("""
molecule dimer {
-1 2
S        0.18146836     0.34051234     0.11429205
O       -1.12219364     0.35969221    -0.70677779
O        0.75925692    -1.04071688     0.47819370
units angstrom
}
""")

G2_Spiropentane = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.00000000
C        0.00000000     0.76201400     1.26575200
C        0.00000000    -0.76201400     1.26575200
C        0.76201400     0.00000000    -1.26575200
C       -0.76201400     0.00000000    -1.26575200
H       -0.91402300     1.26507500     1.56809000
H        0.91402300     1.26507500     1.56809000
H       -0.91402300    -1.26507500     1.56809000
H        0.91402300    -1.26507500     1.56809000
H        1.26507500    -0.91402300    -1.56809000
H        1.26507500     0.91402300    -1.56809000
H       -1.26507500    -0.91402300    -1.56809000
H       -1.26507500     0.91402300    -1.56809000
units angstrom
}
""")

G2_t13Butadiene = input.process_input("""
molecule dimer {
0 1
C        0.60571100     1.74655000     0.00000000
C        0.60571100     0.40408300     0.00000000
C       -0.60571100    -0.40408300     0.00000000
C       -0.60571100    -1.74655000     0.00000000
H        1.52761700     2.31744300     0.00000000
H       -0.32113200     2.31311600     0.00000000
H        1.55350300    -0.13364000     0.00000000
H       -1.55350300     0.13364000     0.00000000
H        0.32113200    -2.31311600     0.00000000
H       -1.52761700    -2.31744300     0.00000000
units angstrom
}
""")

G2_tButane = input.process_input("""
molecule dimer {
0 1
C        0.70258100     1.82087300     0.00000000
C        0.70258100     0.29632500     0.00000000
C       -0.70258100    -0.29632500     0.00000000
C       -0.70258100    -1.82087300     0.00000000
H        1.71980900     2.22234000     0.00000000
H       -1.71980900    -2.22234000     0.00000000
H        0.18815400     2.21036200     0.88361400
H        0.18815400     2.21036200    -0.88361400
H       -0.18815400    -2.21036200     0.88361400
H       -0.18815400    -2.21036200    -0.88361400
H        1.24770700    -0.07266000    -0.87756900
H        1.24770700    -0.07266000     0.87756900
H       -1.24770700     0.07266000    -0.87756900
H       -1.24770700     0.07266000     0.87756900
units angstrom
}
""")

G2_tEthylamine = input.process_input("""
molecule dimer {
0 1
C        1.21001400    -0.35359800     0.00000000
C        0.00000000     0.57595100     0.00000000
N       -1.30535100    -0.08747800     0.00000000
H        2.14931000     0.20849800     0.00000000
H        1.20179600    -0.99776000     0.88490900
H        1.20179600    -0.99776000    -0.88490900
H        0.03456100     1.23096300    -0.87647800
H        0.03456100     1.23096300     0.87647800
H       -1.37232600    -0.69834000     0.81313200
H       -1.37232600    -0.69834000    -0.81313200
units angstrom
}
""")

G2_Tetrachloroethene = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.67540200
C        0.00000000     0.00000000    -0.67540200
Cl       0.00000000     1.44893900     1.58970100
Cl       0.00000000    -1.44893900     1.58970100
Cl       0.00000000    -1.44893900    -1.58970100
Cl       0.00000000     1.44893900    -1.58970100
units angstrom
}
""")

G2_Tetrafluoroethene = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.66323000
C        0.00000000     0.00000000    -0.66323000
F        0.00000000     1.11266500     1.38565200
F        0.00000000    -1.11266500     1.38565200
F        0.00000000     1.11266500    -1.38565200
F        0.00000000    -1.11266500    -1.38565200
units angstrom
}
""")

G2_Thioethanol = input.process_input("""
molecule dimer {
0 1
C        1.51434300     0.67941200     0.00000000
C        0.00000000     0.82641200     0.00000000
S       -0.75606800    -0.83128400     0.00000000
H       -2.03534600    -0.42773800     0.00000000
H       -0.32497000     1.37648200     0.88579300
H       -0.32497000     1.37648200    -0.88579300
H        1.98650300     1.66508200     0.00000000
H        1.85490400     0.13764500     0.88549400
H        1.85490400     0.13764500    -0.88549400
units angstrom
}
""")

G2_Thioformaldehyde = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -1.02886345
S        0.00000000     0.00000000     0.58670648
H        0.00000000    -0.92404263    -1.60706147
H        0.00000000     0.92404263    -1.60706147
units angstrom
}
""")

G2_Thioformaldehyde_anion = input.process_input("""
molecule dimer {
-1 2
C       -1.09430652     0.00000000     0.00132556
S        0.61830485     0.00000000    -0.01659267
H       -1.66351921    -0.92115781     0.12876470
H       -1.66351921     0.92115781     0.12876470
units angstrom
}
""")

G2_Thioformaldehyde_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000     0.00000000    -1.00623692
S        0.00000000     0.00000000     0.57366489
H        0.00000000    -0.93845146    -1.57060835
H        0.00000000     0.93845146    -1.57060835
units angstrom
}
""")

G2_Thiooxirane = input.process_input("""
molecule dimer {
0 1
C        0.00000000    -0.73971900    -0.79233400
S        0.00000000     0.00000000     0.86347400
C        0.00000000     0.73971900    -0.79233400
H       -0.91394000    -1.25014200    -1.07689400
H        0.91394000    -1.25014200    -1.07689400
H        0.91394000     1.25014200    -1.07689400
H       -0.91394000     1.25014200    -1.07689400
units angstrom
}
""")

G2_Thiooxirane_cation = input.process_input("""
molecule dimer {
1 2
C        0.00000000    -0.73376849    -0.80927133
S        0.00000000     0.00000000     0.86360677
C        0.00000000     0.73376849    -0.80927133
H       -0.92104658    -1.27251972    -1.02661309
H        0.92104658    -1.27251972    -1.02661309
H        0.92104658     1.27251972    -1.02661309
H       -0.92104658     1.27251972    -1.02661309
units angstrom
}
""")

G2_Thiophene = input.process_input("""
molecule dimer {
0 1
S        0.00000000     0.00000000     1.18975300
C        0.00000000     1.23387600    -0.00147400
C        0.00000000    -1.23387600    -0.00147400
C        0.00000000     0.70917300    -1.27232200
C        0.00000000    -0.70917300    -1.27232200
H        0.00000000     2.27534300     0.29198400
H        0.00000000    -2.27534300     0.29198400
H        0.00000000     1.32193400    -2.16723100
H        0.00000000    -1.32193400    -2.16723100
units angstrom
}
""")

G2_Toluene = input.process_input("""
molecule dimer {
0 1
C       -0.00968800     0.91276000     0.00000000
C       -0.00686100     0.19598300     1.20155000
C       -0.00686100     0.19598300    -1.20155000
C       -0.00686100    -1.19832800     1.20495300
C       -0.00686100    -1.19832800    -1.20495300
C       -0.00476500    -1.90064300     0.00000000
C        0.02954600     2.41769400     0.00000000
H        1.06081300     2.78643800     0.00000000
H       -0.46716400     2.82493900    -0.88487700
H       -0.46716400     2.82493900     0.88487700
H       -0.00974000    -1.73669700     2.14965400
H       -0.00974000    -1.73669700    -2.14965400
H       -0.00748700    -2.98761700     0.00000000
H       -0.01271200     0.73698500     2.14630200
H       -0.01271200     0.73698500    -2.14630200
units angstrom
}
""")

G2_Toluene_cation = input.process_input("""
molecule dimer {
1 2
C       -0.02936000    -1.87603900     0.00000000
C        0.00000000     0.91035800     0.00000000
C        1.21110500     0.13620700     0.00000000
C       -1.24271400     0.23517700     0.00000000
C       -1.27610600    -1.13910000     0.00000000
C        1.20517500    -1.22216900     0.00000000
H        2.15871800     0.67115600     0.00000000
H       -2.16665300     0.80572400     0.00000000
H       -2.21234300    -1.68886800     0.00000000
H        2.12659900    -1.79436700     0.00000000
H       -0.06611600    -2.96290000     0.00000000
C        0.09108900     2.39292500     0.00000000
H       -0.89244200     2.86321400     0.00000000
H        0.64855000     2.74094200     0.87894000
H        0.64855000     2.74094200    -0.87894000
units angstrom
}
""")

G2_Trimethylamine = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     0.39584600
C        0.00000000     1.37802100    -0.06517500
C        1.19340100    -0.68901100    -0.06517500
C       -1.19340100    -0.68901100    -0.06517500
H        0.00000000     1.46114200    -1.16789900
H        0.88615600     1.89105200     0.31765500
H       -0.88615600     1.89105200     0.31765500
H        1.26538600    -0.73057100    -1.16789900
H        1.19462100    -1.71296000     0.31765500
H        2.08077700    -0.17809200     0.31765500
H       -1.26538600    -0.73057100    -1.16789900
H       -2.08077700    -0.17809200     0.31765500
H       -1.19462100    -1.71296000     0.31765500
units angstrom
}
""")

G2_Vinyl_ncl_cation = input.process_input("""
molecule dimer {
1 1
C        0.61597000     0.00000000     0.00000000
H        0.00000000     1.11709000     0.00000000
C       -0.61597000     0.00000000     0.00000000
H        1.69719828     0.00503444     0.00000000
H       -1.69719828     0.00503444     0.00000000
units angstrom
}
""")

G2_Vinylchloride = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.75601600     0.00000000
C        1.30322300     1.02850700     0.00000000
Cl      -0.63155500    -0.85498000     0.00000000
H       -0.77109800     1.51696300     0.00000000
H        2.05609500     0.24942700     0.00000000
H        1.63209600     2.06112500     0.00000000
units angstrom
}
""")

G2_Vinylcyanide = input.process_input("""
molecule dimer {
0 1
C       -0.16159400    -1.63862500     0.00000000
C        0.58495700    -0.52496100     0.00000000
C        0.00000000     0.78225300     0.00000000
H       -1.24520300    -1.59816900     0.00000000
H        0.30597300    -2.61640500     0.00000000
H        1.66986300    -0.57210700     0.00000000
N       -0.46725900     1.86781100     0.00000000
units angstrom
}
""")

G2_Vinylfluoride = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.43771400     0.00000000
C        1.19192300    -0.14508700     0.00000000
F       -1.14892900    -0.27833200     0.00000000
H       -0.18644500     1.50577800     0.00000000
H        1.29134800    -1.22283300     0.00000000
H        2.08392400     0.46627900     0.00000000
units angstrom
}
""")

# <<< Geometry Specification Strings >>>
rxnpattern = re.compile(r'^(.+)-(.+)-(.+)$')
GEOS = {}
for rxn in HRXN:
   for rgt in ACTV['%s-%s' % (dbse, rxn)]:

            molname = rxnpattern.match(rgt)
            GEOS['%s' % (rgt)] = eval('%s_%s' % (dbse, molname.group(2)))

