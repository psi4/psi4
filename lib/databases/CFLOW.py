import re
import input

# <<< CFLOW Database Module >>>
# Geometries and Reference energies from.
dbse = 'CFLOW'

# <<< Database Members >>>
HRXN = ['2Ae2Ae-3.8', '4Ae4Ae-3.8', '6Ae6Ae-3.8', '8Ae8Ae-3.8', '10Ae10Ae-3.8', '12Ae12Ae-3.8', 
        'BzBz_S-3.9', '2BzBz_S-3.5', '2Bz2Bz_S-3.8', '3Bz2Bz_S-3.5', '3Bz3Bz_S-3.7', '4Bz3Bz_S-3.9', 
        'BkybowlBkybowl-3.54', 'BkybowlBkybowl-3.64', 'BkybowlBkybowl-3.74', 'BkybowlBkybowl-3.84', 'BkybowlBkybowl-3.94', 
        'BkybowlBkybowl-3.63', 'C60Bkybowl', 'C60Bkycatch', ]
HRXN_SM = []
HRXN_LG = []
Alkenes = ['2Ae2Ae-3.8', '4Ae4Ae-3.8', '6Ae6Ae-3.8', '8Ae8Ae-3.8', '10Ae10Ae-3.8', '12Ae12Ae-3.8',]
Arenes = ['BzBz_S-3.9', '2BzBz_S-3.5', '2Bz2Bz_S-3.8', '3Bz2Bz_S-3.5', '3Bz3Bz_S-3.7', '4Bz3Bz_S-3.9',]
Pulay = ['BkybowlBkybowl-3.54', 'BkybowlBkybowl-3.64', 'BkybowlBkybowl-3.74', 'BkybowlBkybowl-3.84', 'BkybowlBkybowl-3.94',]
Grimme = ['BkybowlBkybowl-3.63', 'C60Bkybowl', 'C60Bkycatch', ]
Dimers = ['2Ae2Ae-3.8', '4Ae4Ae-3.8', '6Ae6Ae-3.8', '8Ae8Ae-3.8', '10Ae10Ae-3.8', '12Ae12Ae-3.8',
          'BzBz_S-3.9', '2Bz2Bz_S-3.8', '3Bz3Bz_S-3.7', 
          'BkybowlBkybowl-3.54', 'BkybowlBkybowl-3.64', 'BkybowlBkybowl-3.74', 'BkybowlBkybowl-3.84', 'BkybowlBkybowl-3.94',
          'BkybowlBkybowl-3.63', ]

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supermolecular calculations
for rxn in HRXN:

   if (rxn in Dimers):
      RXNM[      '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                           '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                           '%s-%s-monoA-unCP' % (dbse, rxn) : -2 }

      ACTV_SA[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

      ACTV_CP[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-CP'   % (dbse, rxn) ]

      ACTV[      '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-unCP' % (dbse, rxn) ]

   else:
      RXNM[      '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                           '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                           '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                           '%s-%s-monoA-unCP' % (dbse, rxn) : -1,
                                           '%s-%s-monoB-unCP' % (dbse, rxn) : -1 }

      ACTV_SA[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

      ACTV_CP[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-CP'   % (dbse, rxn),
                                           '%s-%s-monoB-CP'   % (dbse, rxn) ]

      ACTV[      '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                           '%s-%s-monoA-unCP' % (dbse, rxn),
                                           '%s-%s-monoB-unCP' % (dbse, rxn) ]

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, '2Ae2Ae-3.8'            )] =    0.000
BIND['%s-%s'            % (dbse, '4Ae4Ae-3.8'            )] =    0.000
BIND['%s-%s'            % (dbse, '6Ae6Ae-3.8'            )] =    0.000
BIND['%s-%s'            % (dbse, '8Ae8Ae-3.8'            )] =    0.000
BIND['%s-%s'            % (dbse, '10Ae10Ae-3.8'          )] =    0.000
BIND['%s-%s'            % (dbse, '12Ae12Ae-3.8'          )] =    0.000
BIND['%s-%s'            % (dbse, 'BzBz_S-3.9'            )] =    0.000
BIND['%s-%s'            % (dbse, '2BzBz_S-3.5'           )] =    0.000
BIND['%s-%s'            % (dbse, '2Bz2Bz_S-3.8'          )] =    0.000
BIND['%s-%s'            % (dbse, '3Bz2Bz_S-3.5'          )] =    0.000
BIND['%s-%s'            % (dbse, '3Bz3Bz_S-3.7'          )] =    0.000
BIND['%s-%s'            % (dbse, '4Bz3Bz_S-3.9'          )] =    0.000
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.54'   )] =    0.000
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.64'   )] =    0.000
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.74'   )] =    0.000
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.84'   )] =    0.000
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.94'   )] =    0.000
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.63'   )] =    0.000
BIND['%s-%s'            % (dbse, 'C60Bkybowl'            )] =    0.000
BIND['%s-%s'            % (dbse, 'C60Bkycatch'           )] =    0.000

# <<< Comment Lines >>>
TAGL = {}
rxnpattern = re.compile(r'^(.+)-(.+)$')

TAGL['%s-%s'            % (dbse, '10Ae10Ae-3.8'          )] = """Decene (C10H12) Dimer, stacked at 3.8 A """
TAGL['%s-%s-dimer'      % (dbse, '10Ae10Ae-3.8'          )] = """Dimer from Decene (C10H12) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-CP'   % (dbse, '10Ae10Ae-3.8'          )] = """Monomer A from Decene (C10H12) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-CP'   % (dbse, '10Ae10Ae-3.8'          )] = """Monomer B from Decene (C10H12) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-unCP' % (dbse, '10Ae10Ae-3.8'          )] = """Monomer A from Decene (C10H12) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-unCP' % (dbse, '10Ae10Ae-3.8'          )] = """Monomer B from Decene (C10H12) Dimer, stacked at 3.8 A """
TAGL['%s-%s'            % (dbse, '12Ae12Ae-3.8'          )] = """Dodecene (C12H14) Dimer, stacked at 3.8 A """
TAGL['%s-%s-dimer'      % (dbse, '12Ae12Ae-3.8'          )] = """Dimer from Dodecene (C12H14) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-CP'   % (dbse, '12Ae12Ae-3.8'          )] = """Monomer A from Dodecene (C12H14) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-CP'   % (dbse, '12Ae12Ae-3.8'          )] = """Monomer B from Dodecene (C12H14) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-unCP' % (dbse, '12Ae12Ae-3.8'          )] = """Monomer A from Dodecene (C12H14) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-unCP' % (dbse, '12Ae12Ae-3.8'          )] = """Monomer B from Dodecene (C12H14) Dimer, stacked at 3.8 A """
TAGL['%s-%s'            % (dbse, '2Ae2Ae-3.8'            )] = """Ethene (C2H4) Dimer, stacked at 3.8 A """
TAGL['%s-%s-dimer'      % (dbse, '2Ae2Ae-3.8'            )] = """Dimer from Ethene (C2H4) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-CP'   % (dbse, '2Ae2Ae-3.8'            )] = """Monomer A from Ethene (C2H4) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-CP'   % (dbse, '2Ae2Ae-3.8'            )] = """Monomer B from Ethene (C2H4) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-unCP' % (dbse, '2Ae2Ae-3.8'            )] = """Monomer A from Ethene (C2H4) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-unCP' % (dbse, '2Ae2Ae-3.8'            )] = """Monomer B from Ethene (C2H4) Dimer, stacked at 3.8 A """
TAGL['%s-%s'            % (dbse, '2Bz2Bz_S-3.8'          )] = """Napthalene Dimer, stacked at 3.8 A """
TAGL['%s-%s-dimer'      % (dbse, '2Bz2Bz_S-3.8'          )] = """Dimer from Napthalene Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-CP'   % (dbse, '2Bz2Bz_S-3.8'          )] = """Monomer A from Napthalene Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-CP'   % (dbse, '2Bz2Bz_S-3.8'          )] = """Monomer B from Napthalene Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-unCP' % (dbse, '2Bz2Bz_S-3.8'          )] = """Monomer A from Napthalene Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-unCP' % (dbse, '2Bz2Bz_S-3.8'          )] = """Monomer B from Napthalene Dimer, stacked at 3.8 A """
TAGL['%s-%s'            % (dbse, '2BzBz_S-3.5'           )] = """Napthalene-Benzene Complex, stacked at 3.5 A, NOT YET COMPLETELY CERTIAN THIS IS CORRECT SCANNING MIN """
TAGL['%s-%s-dimer'      % (dbse, '2BzBz_S-3.5'           )] = """Dimer from Napthalene-Benzene Complex, stacked at 3.5 A, NOT YET COMPLETELY CERTIAN THIS IS CORRECT SCANNING MIN """
TAGL['%s-%s-monoA-CP'   % (dbse, '2BzBz_S-3.5'           )] = """Monomer A from Napthalene-Benzene Complex, stacked at 3.5 A, NOT YET COMPLETELY CERTIAN THIS IS CORRECT SCANNING MIN """
TAGL['%s-%s-monoB-CP'   % (dbse, '2BzBz_S-3.5'           )] = """Monomer B from Napthalene-Benzene Complex, stacked at 3.5 A, NOT YET COMPLETELY CERTIAN THIS IS CORRECT SCANNING MIN """
TAGL['%s-%s-monoA-unCP' % (dbse, '2BzBz_S-3.5'           )] = """Monomer A from Napthalene-Benzene Complex, stacked at 3.5 A, NOT YET COMPLETELY CERTIAN THIS IS CORRECT SCANNING MIN """
TAGL['%s-%s-monoB-unCP' % (dbse, '2BzBz_S-3.5'           )] = """Monomer B from Napthalene-Benzene Complex, stacked at 3.5 A, NOT YET COMPLETELY CERTIAN THIS IS CORRECT SCANNING MIN """
TAGL['%s-%s'            % (dbse, '3Bz2Bz_S-3.5'          )] = """Anthracene-Napthalene Complex, stacked at 3.5 A """
TAGL['%s-%s-dimer'      % (dbse, '3Bz2Bz_S-3.5'          )] = """Dimer from Anthracene-Napthalene Complex, stacked at 3.5 A """
TAGL['%s-%s-monoA-CP'   % (dbse, '3Bz2Bz_S-3.5'          )] = """Monomer A from Anthracene-Napthalene Complex, stacked at 3.5 A """
TAGL['%s-%s-monoB-CP'   % (dbse, '3Bz2Bz_S-3.5'          )] = """Monomer B from Anthracene-Napthalene Complex, stacked at 3.5 A """
TAGL['%s-%s-monoA-unCP' % (dbse, '3Bz2Bz_S-3.5'          )] = """Monomer A from Anthracene-Napthalene Complex, stacked at 3.5 A """
TAGL['%s-%s-monoB-unCP' % (dbse, '3Bz2Bz_S-3.5'          )] = """Monomer B from Anthracene-Napthalene Complex, stacked at 3.5 A """
TAGL['%s-%s'            % (dbse, '3Bz3Bz_S-3.7'          )] = """Anthracene Dimer, stacked at 3.7 A, NOT YET ABSOLUTETLY CERTAIN THIS IS SCANNING MIN """
TAGL['%s-%s-dimer'      % (dbse, '3Bz3Bz_S-3.7'          )] = """Dimer from Anthracene Dimer, stacked at 3.7 A, NOT YET ABSOLUTETLY CERTAIN THIS IS SCANNING MIN """
TAGL['%s-%s-monoA-CP'   % (dbse, '3Bz3Bz_S-3.7'          )] = """Monomer A from Anthracene Dimer, stacked at 3.7 A, NOT YET ABSOLUTETLY CERTAIN THIS IS SCANNING MIN """
TAGL['%s-%s-monoB-CP'   % (dbse, '3Bz3Bz_S-3.7'          )] = """Monomer B from Anthracene Dimer, stacked at 3.7 A, NOT YET ABSOLUTETLY CERTAIN THIS IS SCANNING MIN """
TAGL['%s-%s-monoA-unCP' % (dbse, '3Bz3Bz_S-3.7'          )] = """Monomer A from Anthracene Dimer, stacked at 3.7 A, NOT YET ABSOLUTETLY CERTAIN THIS IS SCANNING MIN """
TAGL['%s-%s-monoB-unCP' % (dbse, '3Bz3Bz_S-3.7'          )] = """Monomer B from Anthracene Dimer, stacked at 3.7 A, NOT YET ABSOLUTETLY CERTAIN THIS IS SCANNING MIN """
TAGL['%s-%s'            % (dbse, '4Ae4Ae-3.8'            )] = """Butene (C4H6) Dimer, stacked at 3.8 A """
TAGL['%s-%s-dimer'      % (dbse, '4Ae4Ae-3.8'            )] = """Dimer from Butene (C4H6) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-CP'   % (dbse, '4Ae4Ae-3.8'            )] = """Monomer A from Butene (C4H6) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-CP'   % (dbse, '4Ae4Ae-3.8'            )] = """Monomer B from Butene (C4H6) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-unCP' % (dbse, '4Ae4Ae-3.8'            )] = """Monomer A from Butene (C4H6) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-unCP' % (dbse, '4Ae4Ae-3.8'            )] = """Monomer B from Butene (C4H6) Dimer, stacked at 3.8 A """
TAGL['%s-%s'            % (dbse, '4Bz3Bz_S-3.9'          )] = """Naphthacene-Anthracene Complex, stacked at 3.9 A, PLACEHOLDER GEOMETRY """
TAGL['%s-%s-dimer'      % (dbse, '4Bz3Bz_S-3.9'          )] = """Dimer from Naphthacene-Anthracene Complex, stacked at 3.9 A, PLACEHOLDER GEOMETRY """
TAGL['%s-%s-monoA-CP'   % (dbse, '4Bz3Bz_S-3.9'          )] = """Monomer A from Naphthacene-Anthracene Complex, stacked at 3.9 A, PLACEHOLDER GEOMETRY """
TAGL['%s-%s-monoB-CP'   % (dbse, '4Bz3Bz_S-3.9'          )] = """Monomer B from Naphthacene-Anthracene Complex, stacked at 3.9 A, PLACEHOLDER GEOMETRY """
TAGL['%s-%s-monoA-unCP' % (dbse, '4Bz3Bz_S-3.9'          )] = """Monomer A from Naphthacene-Anthracene Complex, stacked at 3.9 A, PLACEHOLDER GEOMETRY """
TAGL['%s-%s-monoB-unCP' % (dbse, '4Bz3Bz_S-3.9'          )] = """Monomer B from Naphthacene-Anthracene Complex, stacked at 3.9 A, PLACEHOLDER GEOMETRY """
TAGL['%s-%s'            % (dbse, '6Ae6Ae-3.8'            )] = """Hexene (C6H8) Dimer, stacked at 3.8 A """
TAGL['%s-%s-dimer'      % (dbse, '6Ae6Ae-3.8'            )] = """Dimer from Hexene (C6H8) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-CP'   % (dbse, '6Ae6Ae-3.8'            )] = """Monomer A from Hexene (C6H8) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-CP'   % (dbse, '6Ae6Ae-3.8'            )] = """Monomer B from Hexene (C6H8) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-unCP' % (dbse, '6Ae6Ae-3.8'            )] = """Monomer A from Hexene (C6H8) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-unCP' % (dbse, '6Ae6Ae-3.8'            )] = """Monomer B from Hexene (C6H8) Dimer, stacked at 3.8 A """
TAGL['%s-%s'            % (dbse, '8Ae8Ae-3.8'            )] = """Octene (C8H10) Dimer, stacked at 3.8 A """
TAGL['%s-%s-dimer'      % (dbse, '8Ae8Ae-3.8'            )] = """Dimer from Octene (C8H10) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-CP'   % (dbse, '8Ae8Ae-3.8'            )] = """Monomer A from Octene (C8H10) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-CP'   % (dbse, '8Ae8Ae-3.8'            )] = """Monomer B from Octene (C8H10) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoA-unCP' % (dbse, '8Ae8Ae-3.8'            )] = """Monomer A from Octene (C8H10) Dimer, stacked at 3.8 A """
TAGL['%s-%s-monoB-unCP' % (dbse, '8Ae8Ae-3.8'            )] = """Monomer B from Octene (C8H10) Dimer, stacked at 3.8 A """
TAGL['%s-%s'            % (dbse, 'BkybowlBkybowl-3.54'   )] = """Corannulene Dimer, stacked at 3.54 A, Pulay geometry """
TAGL['%s-%s-dimer'      % (dbse, 'BkybowlBkybowl-3.54'   )] = """Dimer from Corannulene Dimer, stacked at 3.54 A, Pulay geometry """
TAGL['%s-%s-monoA-CP'   % (dbse, 'BkybowlBkybowl-3.54'   )] = """Monomer A from Corannulene Dimer, stacked at 3.54 A, Pulay geometry """
TAGL['%s-%s-monoB-CP'   % (dbse, 'BkybowlBkybowl-3.54'   )] = """Monomer B from Corannulene Dimer, stacked at 3.54 A, Pulay geometry """
TAGL['%s-%s-monoA-unCP' % (dbse, 'BkybowlBkybowl-3.54'   )] = """Monomer A from Corannulene Dimer, stacked at 3.54 A, Pulay geometry """
TAGL['%s-%s-monoB-unCP' % (dbse, 'BkybowlBkybowl-3.54'   )] = """Monomer B from Corannulene Dimer, stacked at 3.54 A, Pulay geometry """
TAGL['%s-%s'            % (dbse, 'BkybowlBkybowl-3.63'   )] = """Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-dimer'      % (dbse, 'BkybowlBkybowl-3.63'   )] = """Dimer from Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-monoA-CP'   % (dbse, 'BkybowlBkybowl-3.63'   )] = """Monomer A from Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-monoB-CP'   % (dbse, 'BkybowlBkybowl-3.63'   )] = """Monomer B from Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-monoA-unCP' % (dbse, 'BkybowlBkybowl-3.63'   )] = """Monomer A from Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-monoB-unCP' % (dbse, 'BkybowlBkybowl-3.63'   )] = """Monomer B from Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s'            % (dbse, 'BkybowlBkybowl-3.64'   )] = """Corannulene Dimer, stacked at 3.64 A, Pulay geometry """
TAGL['%s-%s-dimer'      % (dbse, 'BkybowlBkybowl-3.64'   )] = """Dimer from Corannulene Dimer, stacked at 3.64 A, Pulay geometry """
TAGL['%s-%s-monoA-CP'   % (dbse, 'BkybowlBkybowl-3.64'   )] = """Monomer A from Corannulene Dimer, stacked at 3.64 A, Pulay geometry """
TAGL['%s-%s-monoB-CP'   % (dbse, 'BkybowlBkybowl-3.64'   )] = """Monomer B from Corannulene Dimer, stacked at 3.64 A, Pulay geometry """
TAGL['%s-%s-monoA-unCP' % (dbse, 'BkybowlBkybowl-3.64'   )] = """Monomer A from Corannulene Dimer, stacked at 3.64 A, Pulay geometry """
TAGL['%s-%s-monoB-unCP' % (dbse, 'BkybowlBkybowl-3.64'   )] = """Monomer B from Corannulene Dimer, stacked at 3.64 A, Pulay geometry """
TAGL['%s-%s'            % (dbse, 'BkybowlBkybowl-3.74'   )] = """Corannulene Dimer, stacked at 3.74 A, Pulay geometry """
TAGL['%s-%s-dimer'      % (dbse, 'BkybowlBkybowl-3.74'   )] = """Dimer from Corannulene Dimer, stacked at 3.74 A, Pulay geometry """
TAGL['%s-%s-monoA-CP'   % (dbse, 'BkybowlBkybowl-3.74'   )] = """Monomer A from Corannulene Dimer, stacked at 3.74 A, Pulay geometry """
TAGL['%s-%s-monoB-CP'   % (dbse, 'BkybowlBkybowl-3.74'   )] = """Monomer B from Corannulene Dimer, stacked at 3.74 A, Pulay geometry """
TAGL['%s-%s-monoA-unCP' % (dbse, 'BkybowlBkybowl-3.74'   )] = """Monomer A from Corannulene Dimer, stacked at 3.74 A, Pulay geometry """
TAGL['%s-%s-monoB-unCP' % (dbse, 'BkybowlBkybowl-3.74'   )] = """Monomer B from Corannulene Dimer, stacked at 3.74 A, Pulay geometry """
TAGL['%s-%s'            % (dbse, 'BkybowlBkybowl-3.84'   )] = """Corannulene Dimer, stacked at 3.84 A, Pulay geometry """
TAGL['%s-%s-dimer'      % (dbse, 'BkybowlBkybowl-3.84'   )] = """Dimer from Corannulene Dimer, stacked at 3.84 A, Pulay geometry """
TAGL['%s-%s-monoA-CP'   % (dbse, 'BkybowlBkybowl-3.84'   )] = """Monomer A from Corannulene Dimer, stacked at 3.84 A, Pulay geometry """
TAGL['%s-%s-monoB-CP'   % (dbse, 'BkybowlBkybowl-3.84'   )] = """Monomer B from Corannulene Dimer, stacked at 3.84 A, Pulay geometry """
TAGL['%s-%s-monoA-unCP' % (dbse, 'BkybowlBkybowl-3.84'   )] = """Monomer A from Corannulene Dimer, stacked at 3.84 A, Pulay geometry """
TAGL['%s-%s-monoB-unCP' % (dbse, 'BkybowlBkybowl-3.84'   )] = """Monomer B from Corannulene Dimer, stacked at 3.84 A, Pulay geometry """
TAGL['%s-%s'            % (dbse, 'BkybowlBkybowl-3.94'   )] = """Corannulene Dimer, stacked at 3.94 A, Pulay geometry """
TAGL['%s-%s-dimer'      % (dbse, 'BkybowlBkybowl-3.94'   )] = """Dimer from Corannulene Dimer, stacked at 3.94 A, Pulay geometry """
TAGL['%s-%s-monoA-CP'   % (dbse, 'BkybowlBkybowl-3.94'   )] = """Monomer A from Corannulene Dimer, stacked at 3.94 A, Pulay geometry """
TAGL['%s-%s-monoB-CP'   % (dbse, 'BkybowlBkybowl-3.94'   )] = """Monomer B from Corannulene Dimer, stacked at 3.94 A, Pulay geometry """
TAGL['%s-%s-monoA-unCP' % (dbse, 'BkybowlBkybowl-3.94'   )] = """Monomer A from Corannulene Dimer, stacked at 3.94 A, Pulay geometry """
TAGL['%s-%s-monoB-unCP' % (dbse, 'BkybowlBkybowl-3.94'   )] = """Monomer B from Corannulene Dimer, stacked at 3.94 A, Pulay geometry """
TAGL['%s-%s'            % (dbse, 'BzBz_S-3.9'            )] = """Benzene Dimer, stacked at 3.9 A """
TAGL['%s-%s-dimer'      % (dbse, 'BzBz_S-3.9'            )] = """Dimer from Benzene Dimer, stacked at 3.9 A """
TAGL['%s-%s-monoA-CP'   % (dbse, 'BzBz_S-3.9'            )] = """Monomer A from Benzene Dimer, stacked at 3.9 A """
TAGL['%s-%s-monoB-CP'   % (dbse, 'BzBz_S-3.9'            )] = """Monomer B from Benzene Dimer, stacked at 3.9 A """
TAGL['%s-%s-monoA-unCP' % (dbse, 'BzBz_S-3.9'            )] = """Monomer A from Benzene Dimer, stacked at 3.9 A """
TAGL['%s-%s-monoB-unCP' % (dbse, 'BzBz_S-3.9'            )] = """Monomer B from Benzene Dimer, stacked at 3.9 A """
TAGL['%s-%s'            % (dbse, 'C60Bkybowl'            )] = """C60 @ Corannulene Buckybowl """
TAGL['%s-%s-dimer'      % (dbse, 'C60Bkybowl'            )] = """Dimer from C60 @ Corannulene Buckybowl """
TAGL['%s-%s-monoA-CP'   % (dbse, 'C60Bkybowl'            )] = """Monomer A from C60 @ Corannulene Buckybowl """
TAGL['%s-%s-monoB-CP'   % (dbse, 'C60Bkybowl'            )] = """Monomer B from C60 @ Corannulene Buckybowl """
TAGL['%s-%s-monoA-unCP' % (dbse, 'C60Bkybowl'            )] = """Monomer A from C60 @ Corannulene Buckybowl """
TAGL['%s-%s-monoB-unCP' % (dbse, 'C60Bkybowl'            )] = """Monomer B from C60 @ Corannulene Buckybowl """
TAGL['%s-%s'            % (dbse, 'C60Bkycatch'           )] = """C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-dimer'      % (dbse, 'C60Bkycatch'           )] = """Dimer from C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-monoA-CP'   % (dbse, 'C60Bkycatch'           )] = """Monomer A from C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-monoB-CP'   % (dbse, 'C60Bkycatch'           )] = """Monomer B from C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-monoA-unCP' % (dbse, 'C60Bkycatch'           )] = """Monomer A from C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-monoB-unCP' % (dbse, 'C60Bkycatch'           )] = """Monomer B from C60 @ C60H28 Buckycatcher """

# <<< Molecule Specifications >>>
monoA_unCP = 'monoA = dimer.extract_subsets(1)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_unCP = 'monoB = dimer.extract_subsets(2)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

CFLOW_10Ae10Ae_3p8 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     1.08000000
H        0.00000000    -0.93530700     1.62000000
C        0.00000000     0.00000000    -0.25197000
H        0.00000000     0.94916400    -0.79997000
C        0.00000000    -1.32588500    -1.01747000
H        0.00000000    -2.27504800    -0.46947000
C        0.00000000    -1.32588500    -2.34944000
H        0.00000000    -0.37672100    -2.89744000
C        0.00000000     1.32588500     1.84550000
H        0.00000000     2.27504800     1.29750000
C        0.00000000     1.32588500     3.17747000
H        0.00000000     0.37672100     3.72547000
C        0.00000000     2.65176900     3.94297100
H        0.00000000     3.60093300     3.39497000
C        0.00000000     2.65176900     5.27494100
H        0.00000000     3.60093300     5.82294100
H        0.00000000     1.70260600     5.82294100
C        0.00000000    -2.65176900    -3.11494100
H        0.00000000    -3.60093300    -2.56694000
C        0.00000000    -2.65176900    -4.44691100
H        0.00000000    -3.60093300    -4.99491100
H        0.00000000    -1.70260600    -4.99491100
--
0 1
C        3.80000000     0.00000000     1.08000000
H        3.80000000    -0.93530700     1.62000000
C        3.80000000     0.00000000    -0.25197000
H        3.80000000     0.94916400    -0.79997000
C        3.80000000    -1.32588500    -1.01747000
H        3.80000000    -2.27504800    -0.46947000
C        3.80000000    -1.32588500    -2.34944000
H        3.80000000    -0.37672100    -2.89744000
C        3.80000000     1.32588500     1.84550000
H        3.80000000     2.27504800     1.29750000
C        3.80000000     1.32588500     3.17747000
H        3.80000000     0.37672100     3.72547000
C        3.80000000     2.65176900     3.94297100
H        3.80000000     3.60093300     3.39497000
C        3.80000000     2.65176900     5.27494100
H        3.80000000     3.60093300     5.82294100
H        3.80000000     1.70260600     5.82294100
C        3.80000000    -2.65176900    -3.11494100
H        3.80000000    -3.60093300    -2.56694000
C        3.80000000    -2.65176900    -4.44691100
H        3.80000000    -3.60093300    -4.99491100
H        3.80000000    -1.70260600    -4.99491100
units angstrom
}
""")

CFLOW_12Ae12Ae_3p8 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     1.08000000
H        0.00000000    -0.93530700     1.62000000
C        0.00000000     0.00000000    -0.25197000
H        0.00000000     0.94916400    -0.79997000
C        0.00000000    -1.32588500    -1.01747000
H        0.00000000    -2.27504800    -0.46947000
C        0.00000000    -1.32588500    -2.34944000
H        0.00000000    -0.37672100    -2.89744000
C        0.00000000     1.32588500     1.84550000
H        0.00000000     2.27504800     1.29750000
C        0.00000000     1.32588500     3.17747000
H        0.00000000     0.37672100     3.72547000
C        0.00000000     2.65176900     3.94297100
H        0.00000000     3.60093300     3.39497000
C        0.00000000     2.65176900     5.27494100
H        0.00000000     1.70260600     5.82294100
C        0.00000000    -2.65176900    -3.11494100
H        0.00000000    -3.60093300    -2.56694000
C        0.00000000    -2.65176900    -4.44691100
H        0.00000000    -3.60093300    -4.99491100
H        0.00000000    -1.70260600    -4.99491100
C        0.00000000     3.97765400     6.04044100
H        0.00000000     4.92681800     5.49244100
C        0.00000000     3.97765400     7.37241100
H        0.00000000     4.92681800     7.92041100
H        0.00000000     3.02849000     7.92041100
--
0 1
C        3.80000000     0.00000000     1.08000000
H        3.80000000    -0.93530700     1.62000000
C        3.80000000     0.00000000    -0.25197000
H        3.80000000     0.94916400    -0.79997000
C        3.80000000    -1.32588500    -1.01747000
H        3.80000000    -2.27504800    -0.46947000
C        3.80000000    -1.32588500    -2.34944000
H        3.80000000    -0.37672100    -2.89744000
C        3.80000000     1.32588500     1.84550000
H        3.80000000     2.27504800     1.29750000
C        3.80000000     1.32588500     3.17747000
H        3.80000000     0.37672100     3.72547000
C        3.80000000     2.65176900     3.94297100
H        3.80000000     3.60093300     3.39497000
C        3.80000000     2.65176900     5.27494100
H        3.80000000     1.70260600     5.82294100
C        3.80000000    -2.65176900    -3.11494100
H        3.80000000    -3.60093300    -2.56694000
C        3.80000000    -2.65176900    -4.44691100
H        3.80000000    -3.60093300    -4.99491100
H        3.80000000    -1.70260600    -4.99491100
C        3.80000000     3.97765400     6.04044100
H        3.80000000     4.92681800     5.49244100
C        3.80000000     3.97765400     7.37241100
H        3.80000000     4.92681800     7.92041100
H        3.80000000     3.02849000     7.92041100
units angstrom
}
""")

CFLOW_2Ae2Ae_3p8 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     1.08000000
H        0.00000000     0.93530700     1.62000000
H        0.00000000    -0.93530700     1.62000000
C        0.00000000     0.00000000    -0.25197000
H        0.00000000     0.94916400    -0.79997000
H        0.00000000    -0.94916400    -0.79997000
--
0 1
C        3.80000000     0.00000000     1.08000000
H        3.80000000     0.93530700     1.62000000
H        3.80000000    -0.93530700     1.62000000
C        3.80000000     0.00000000    -0.25197000
H        3.80000000     0.94916400    -0.79997000
H        3.80000000    -0.94916400    -0.79997000
units angstrom
}
""")

CFLOW_2Bz2Bz_S_3p8 = input.process_input("""
molecule dimer {
0 1
C       -2.41591300    -0.70276500     0.00000000
C       -1.23100600    -1.39300200     0.00000000
C        0.00000000    -0.71089400     0.00000000
C        0.00000000     0.71089400     0.00000000
C       -1.23100600     1.39300200     0.00000000
C       -2.41591300     0.70276500     0.00000000
C        1.23100600    -1.39300200     0.00000000
C        1.23100600     1.39300200     0.00000000
C        2.41591300     0.70276500     0.00000000
C        2.41591300    -0.70276500     0.00000000
H        1.22889800    -2.47183700     0.00000000
H       -3.35016300    -1.23782500     0.00000000
H       -1.22889800    -2.47183700     0.00000000
H       -1.22889800     2.47183700     0.00000000
H       -3.35016300     1.23782500     0.00000000
H        1.22889800     2.47183700     0.00000000
H        3.35016300     1.23782500     0.00000000
H        3.35016300    -1.23782500     0.00000000
--
0 1
C       -2.41591300    -0.70276500     3.80000000
C       -1.23100600    -1.39300200     3.80000000
C        0.00000000    -0.71089400     3.80000000
C        0.00000000     0.71089400     3.80000000
C       -1.23100600     1.39300200     3.80000000
C       -2.41591300     0.70276500     3.80000000
C        1.23100600    -1.39300200     3.80000000
C        1.23100600     1.39300200     3.80000000
C        2.41591300     0.70276500     3.80000000
C        2.41591300    -0.70276500     3.80000000
H        1.22889800    -2.47183700     3.80000000
H       -3.35016300    -1.23782500     3.80000000
H       -1.22889800    -2.47183700     3.80000000
H       -1.22889800     2.47183700     3.80000000
H       -3.35016300     1.23782500     3.80000000
H        1.22889800     2.47183700     3.80000000
H        3.35016300     1.23782500     3.80000000
H        3.35016300    -1.23782500     3.80000000
units angstrom
}
""")

CFLOW_2BzBz_S_3p5 = input.process_input("""
molecule dimer {
0 1
C       -2.41591300    -0.70276500     0.00000000
C       -1.23100600    -1.39300200     0.00000000
C        0.00000000    -0.71089400     0.00000000
C        0.00000000     0.71089400     0.00000000
C       -1.23100600     1.39300200     0.00000000
C       -2.41591300     0.70276500     0.00000000
C        1.23100600    -1.39300200     0.00000000
C        1.23100600     1.39300200     0.00000000
C        2.41591300     0.70276500     0.00000000
C        2.41591300    -0.70276500     0.00000000
H        1.22889800    -2.47183700     0.00000000
H       -3.35016300    -1.23782500     0.00000000
H       -1.22889800    -2.47183700     0.00000000
H       -1.22889800     2.47183700     0.00000000
H       -3.35016300     1.23782500     0.00000000
H        1.22889800     2.47183700     0.00000000
H        3.35016300     1.23782500     0.00000000
H        3.35016300    -1.23782500     0.00000000
--
0 1
C        0.00000000    -1.39150000     3.50000000
C        1.20507400    -0.69575000     3.50000000
C        1.20507400     0.69575000     3.50000000
C        0.00000000     1.39150000     3.50000000
C       -1.20507400     0.69575000     3.50000000
C       -1.20507400    -0.69575000     3.50000000
H        0.00000000    -2.47150000     3.50000000
H        2.14038200    -1.23575000     3.50000000
H        2.14038200     1.23575000     3.50000000
H        0.00000000     2.47150000     3.50000000
H       -2.14038200     1.23575000     3.50000000
H       -2.14038200    -1.23575000     3.50000000
units angstrom
}
""")

CFLOW_3Bz2Bz_S_3p5 = input.process_input("""
molecule dimer {
0 1
C       -3.63206100    -0.70629600     0.00000000
C       -2.45292600    -1.39698600     0.00000000
C       -1.21398000    -0.71604800     0.00000000
C       -1.21398000     0.71604800     0.00000000
C       -2.45292600     1.39698600     0.00000000
C       -3.63206100     0.70629600     0.00000000
C        0.00000000    -1.39478500     0.00000000
C        0.00000000     1.39478500     0.00000000
C        1.21398000     0.71604800     0.00000000
C        1.21398000    -0.71604800     0.00000000
C        2.45292600    -1.39698600     0.00000000
H        2.45137100    -2.47594200     0.00000000
C        3.63206100    -0.70629600     0.00000000
C        3.63206100     0.70629600     0.00000000
C        2.45292600     1.39698600     0.00000000
H        0.00000000    -2.47562800     0.00000000
H       -4.56737900    -1.23960400     0.00000000
H       -2.45137100    -2.47594200     0.00000000
H       -2.45137100     2.47594200     0.00000000
H       -4.56737900     1.23960400     0.00000000
H        0.00000000     2.47562800     0.00000000
H        4.56737900    -1.23960400     0.00000000
H        4.56737900     1.23960400     0.00000000
H        2.45137100     2.47594200     0.00000000
--
0 1
C       -2.41591300    -0.70276500     3.50000000
C       -1.23100600    -1.39300200     3.50000000
C        0.00000000    -0.71089400     3.50000000
C        0.00000000     0.71089400     3.50000000
C       -1.23100600     1.39300200     3.50000000
C       -2.41591300     0.70276500     3.50000000
C        1.23100600    -1.39300200     3.50000000
C        1.23100600     1.39300200     3.50000000
C        2.41591300     0.70276500     3.50000000
C        2.41591300    -0.70276500     3.50000000
H        1.22889800    -2.47183700     3.50000000
H       -3.35016300    -1.23782500     3.50000000
H       -1.22889800    -2.47183700     3.50000000
H       -1.22889800     2.47183700     3.50000000
H       -3.35016300     1.23782500     3.50000000
H        1.22889800     2.47183700     3.50000000
H        3.35016300     1.23782500     3.50000000
H        3.35016300    -1.23782500     3.50000000
units angstrom
}
""")

CFLOW_3Bz3Bz_S_3p7 = input.process_input("""
molecule dimer {
0 1
C       -3.63206100    -0.70629600     0.00000000
C       -2.45292600    -1.39698600     0.00000000
C       -1.21398000    -0.71604800     0.00000000
C       -1.21398000     0.71604800     0.00000000
C       -2.45292600     1.39698600     0.00000000
C       -3.63206100     0.70629600     0.00000000
C        0.00000000    -1.39478500     0.00000000
C        0.00000000     1.39478500     0.00000000
C        1.21398000     0.71604800     0.00000000
C        1.21398000    -0.71604800     0.00000000
C        2.45292600    -1.39698600     0.00000000
H        2.45137100    -2.47594200     0.00000000
C        3.63206100    -0.70629600     0.00000000
C        3.63206100     0.70629600     0.00000000
C        2.45292600     1.39698600     0.00000000
H        0.00000000    -2.47562800     0.00000000
H       -4.56737900    -1.23960400     0.00000000
H       -2.45137100    -2.47594200     0.00000000
H       -2.45137100     2.47594200     0.00000000
H       -4.56737900     1.23960400     0.00000000
H        0.00000000     2.47562800     0.00000000
H        4.56737900    -1.23960400     0.00000000
H        4.56737900     1.23960400     0.00000000
H        2.45137100     2.47594200     0.00000000
--
0 1
C       -3.63206100    -0.70629600     3.70000000
C       -2.45292600    -1.39698600     3.70000000
C       -1.21398000    -0.71604800     3.70000000
C       -1.21398000     0.71604800     3.70000000
C       -2.45292600     1.39698600     3.70000000
C       -3.63206100     0.70629600     3.70000000
C        0.00000000    -1.39478500     3.70000000
C        0.00000000     1.39478500     3.70000000
C        1.21398000     0.71604800     3.70000000
C        1.21398000    -0.71604800     3.70000000
C        2.45292600    -1.39698600     3.70000000
H        2.45137100    -2.47594200     3.70000000
C        3.63206100    -0.70629600     3.70000000
C        3.63206100     0.70629600     3.70000000
C        2.45292600     1.39698600     3.70000000
H        0.00000000    -2.47562800     3.70000000
H       -4.56737900    -1.23960400     3.70000000
H       -2.45137100    -2.47594200     3.70000000
H       -2.45137100     2.47594200     3.70000000
H       -4.56737900     1.23960400     3.70000000
H        0.00000000     2.47562800     3.70000000
H        4.56737900    -1.23960400     3.70000000
H        4.56737900     1.23960400     3.70000000
H        2.45137100     2.47594200     3.70000000
units angstrom
}
""")

CFLOW_4Ae4Ae_3p8 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     1.08000000
H        0.00000000     0.93530700     1.62000000
H        0.00000000    -0.93530700     1.62000000
C        0.00000000     0.00000000    -0.25197000
H        0.00000000     0.94916400    -0.79997000
C        0.00000000    -1.32588500    -1.01747000
H        0.00000000    -2.27504800    -0.46947000
C        0.00000000    -1.32588500    -2.34944000
H        0.00000000    -2.27504800    -2.89744000
H        0.00000000    -0.37672100    -2.89744000
--
0 1
C        3.80000000     0.00000000     1.08000000
H        3.80000000     0.93530700     1.62000000
H        3.80000000    -0.93530700     1.62000000
C        3.80000000     0.00000000    -0.25197000
H        3.80000000     0.94916400    -0.79997000
C        3.80000000    -1.32588500    -1.01747000
H        3.80000000    -2.27504800    -0.46947000
C        3.80000000    -1.32588500    -2.34944000
H        3.80000000    -2.27504800    -2.89744000
H        3.80000000    -0.37672100    -2.89744000
units angstrom
}
""")

CFLOW_4Bz3Bz_S_3p9 = input.process_input("""
molecule dimer {
0 1
C        0.69575000    -1.20507400     0.00000000
C       -0.69575000    -1.20507400     0.00000000
C        1.39150000     0.00000000     0.00000000
C       -1.39150000     0.00000000     0.00000000
C       -0.69575000     1.20507400     0.00000000
C        0.69575000     1.20507400     0.00000000
C        1.39150000     2.41014800     0.00000000
C       -1.39150000     2.41014800     0.00000000
C       -0.69575000     3.61522200     0.00000000
C        0.69575000     3.61522200     0.00000000
C        1.39150000     4.82029600     0.00000000
C       -1.39150000     4.82029600     0.00000000
C       -0.69575000     6.02537000     0.00000000
C        0.69575000     6.02537000     0.00000000
C        1.39150000     7.23044400     0.00000000
C       -1.39150000     7.23044400     0.00000000
C       -0.69575000     8.43551800     0.00000000
C        0.69575000     8.43551800     0.00000000
H        1.23575000    -2.14038200     0.00000000
H       -1.23575000    -2.14038200     0.00000000
H        2.47150000     0.00000000     0.00000000
H       -2.47150000     0.00000000     0.00000000
H        2.47150000     2.41014800     0.00000000
H       -2.47150000     2.41014800     0.00000000
H        2.47150000     4.82029600     0.00000000
H       -2.47150000     4.82029600     0.00000000
H        2.47150000     7.23044400     0.00000000
H       -2.47150000     7.23044400     0.00000000
H       -1.23575000     9.37082600     0.00000000
H        1.23575000     9.37082600     0.00000000
--
0 1
C       -3.63206100    -0.70629600     3.90000000
C       -2.45292600    -1.39698600     3.90000000
C       -1.21398000    -0.71604800     3.90000000
C       -1.21398000     0.71604800     3.90000000
C       -2.45292600     1.39698600     3.90000000
C       -3.63206100     0.70629600     3.90000000
C        0.00000000    -1.39478500     3.90000000
C        0.00000000     1.39478500     3.90000000
C        1.21398000     0.71604800     3.90000000
C        1.21398000    -0.71604800     3.90000000
C        2.45292600    -1.39698600     3.90000000
H        2.45137100    -2.47594200     3.90000000
C        3.63206100    -0.70629600     3.90000000
C        3.63206100     0.70629600     3.90000000
C        2.45292600     1.39698600     3.90000000
H        0.00000000    -2.47562800     3.90000000
H       -4.56737900    -1.23960400     3.90000000
H       -2.45137100    -2.47594200     3.90000000
H       -2.45137100     2.47594200     3.90000000
H       -4.56737900     1.23960400     3.90000000
H        0.00000000     2.47562800     3.90000000
H        4.56737900    -1.23960400     3.90000000
H        4.56737900     1.23960400     3.90000000
H        2.45137100     2.47594200     3.90000000
units angstrom
}
""")

CFLOW_6Ae6Ae_3p8 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     1.08000000
H        0.00000000    -0.93530700     1.62000000
C        0.00000000     0.00000000    -0.25197000
H        0.00000000     0.94916400    -0.79997000
C        0.00000000    -1.32588500    -1.01747000
H        0.00000000    -2.27504800    -0.46947000
C        0.00000000    -1.32588500    -2.34944000
H        0.00000000    -2.27504800    -2.89744000
H        0.00000000    -0.37672100    -2.89744000
C        0.00000000     1.32588500     1.84550000
H        0.00000000     2.27504800     1.29750000
C        0.00000000     1.32588500     3.17747000
H        0.00000000     2.27504800     3.72547000
H        0.00000000     0.37672100     3.72547000
--
0 1
C        3.80000000     0.00000000     1.08000000
H        3.80000000    -0.93530700     1.62000000
C        3.80000000     0.00000000    -0.25197000
H        3.80000000     0.94916400    -0.79997000
C        3.80000000    -1.32588500    -1.01747000
H        3.80000000    -2.27504800    -0.46947000
C        3.80000000    -1.32588500    -2.34944000
H        3.80000000    -2.27504800    -2.89744000
H        3.80000000    -0.37672100    -2.89744000
C        3.80000000     1.32588500     1.84550000
H        3.80000000     2.27504800     1.29750000
C        3.80000000     1.32588500     3.17747000
H        3.80000000     2.27504800     3.72547000
H        3.80000000     0.37672100     3.72547000
units angstrom
}
""")

CFLOW_8Ae8Ae_3p8 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     1.08000000
H        0.00000000    -0.93530700     1.62000000
C        0.00000000     0.00000000    -0.25197000
H        0.00000000     0.94916400    -0.79997000
C        0.00000000    -1.32588500    -1.01747000
H        0.00000000    -2.27504800    -0.46947000
C        0.00000000    -1.32588500    -2.34944000
H        0.00000000    -2.27504800    -2.89744000
H        0.00000000    -0.37672100    -2.89744000
C        0.00000000     1.32588500     1.84550000
H        0.00000000     2.27504800     1.29750000
C        0.00000000     1.32588500     3.17747000
H        0.00000000     0.37672100     3.72547000
C        0.00000000     2.65176900     3.94297100
H        0.00000000     3.60093300     3.39497000
C        0.00000000     2.65176900     5.27494100
H        0.00000000     3.60093300     5.82294100
H        0.00000000     1.70260600     5.82294100
--
0 1
C        3.80000000     0.00000000     1.08000000
H        3.80000000    -0.93530700     1.62000000
C        3.80000000     0.00000000    -0.25197000
H        3.80000000     0.94916400    -0.79997000
C        3.80000000    -1.32588500    -1.01747000
H        3.80000000    -2.27504800    -0.46947000
C        3.80000000    -1.32588500    -2.34944000
H        3.80000000    -2.27504800    -2.89744000
H        3.80000000    -0.37672100    -2.89744000
C        3.80000000     1.32588500     1.84550000
H        3.80000000     2.27504800     1.29750000
C        3.80000000     1.32588500     3.17747000
H        3.80000000     0.37672100     3.72547000
C        3.80000000     2.65176900     3.94297100
H        3.80000000     3.60093300     3.39497000
C        3.80000000     2.65176900     5.27494100
H        3.80000000     3.60093300     5.82294100
H        3.80000000     1.70260600     5.82294100
units angstrom
}
""")

CFLOW_BkybowlBkybowl_3p54 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97212500     2.41862700
C       -0.70622800     0.97212500     2.41862700
C       -1.14280400    -0.37137200     2.41849100
C        0.00000000    -1.20165400     2.41827400
C        1.14280400    -0.37137200     2.41849100
C        1.45779000     2.00650700     1.89581300
C       -1.45779000     2.00650700     1.89581300
C       -2.35873800    -0.76639200     1.89565100
C        0.00000000    -2.48003500     1.89534800
C        2.35873800    -0.76639200     1.89565100
C        0.69261800     3.17924500     1.54846400
C       -0.69261800     3.17924500     1.54846400
C       -2.80958100     1.64120300     1.54875100
C       -3.23765700     0.32374300     1.54864100
C       -2.42918200    -2.16498400     1.54865300
C       -1.30841500    -2.97916300     1.54840200
C        1.30841500    -2.97916300     1.54840200
C        2.42918200    -2.16498400     1.54865300
C        3.23765700     0.32374300     1.54864100
C        2.80958100     1.64120300     1.54875100
H        1.20851300     4.06642600     1.18749100
H       -1.20851300     4.06642600     1.18749100
H       -3.49401500     2.40602700     1.18800700
H       -4.24094400     0.10730100     1.18793900
H       -3.36816400    -2.57958300     1.18817300
H       -1.41248600    -4.00023700     1.18769900
H        1.41248600    -4.00023700     1.18769900
H        3.36816400    -2.57958300     1.18817300
H        4.24094400     0.10730100     1.18793900
H        3.49401500     2.40602700     1.18800700
--
0 1
C        0.70622800     0.97212500     5.95862700
C       -0.70622800     0.97212500     5.95862700
C       -1.14280400    -0.37137200     5.95849100
C        0.00000000    -1.20165400     5.95827400
C        1.14280400    -0.37137200     5.95849100
C        1.45779000     2.00650700     5.43581300
C       -1.45779000     2.00650700     5.43581300
C       -2.35873800    -0.76639200     5.43565100
C        0.00000000    -2.48003500     5.43534800
C        2.35873800    -0.76639200     5.43565100
C        0.69261800     3.17924500     5.08846400
C       -0.69261800     3.17924500     5.08846400
C       -2.80958100     1.64120300     5.08875100
C       -3.23765700     0.32374300     5.08864100
C       -2.42918200    -2.16498400     5.08865300
C       -1.30841500    -2.97916300     5.08840200
C        1.30841500    -2.97916300     5.08840200
C        2.42918200    -2.16498400     5.08865300
C        3.23765700     0.32374300     5.08864100
C        2.80958100     1.64120300     5.08875100
H        1.20851300     4.06642600     4.72749100
H       -1.20851300     4.06642600     4.72749100
H       -3.49401500     2.40602700     4.72800700
H       -4.24094400     0.10730100     4.72793900
H       -3.36816400    -2.57958300     4.72817300
H       -1.41248600    -4.00023700     4.72769900
H        1.41248600    -4.00023700     4.72769900
H        3.36816400    -2.57958300     4.72817300
H        4.24094400     0.10730100     4.72793900
H        3.49401500     2.40602700     4.72800700
units angstrom
}
""")

CFLOW_BkybowlBkybowl_3p63 = input.process_input("""
molecule dimer {
0 1
C        0.70893000     0.97575900     2.47883800
C       -0.70893000     0.97575900     2.47883800
C       -1.14707300    -0.37270700     2.47883800
C        0.00000000    -1.20610400     2.47883800
C        1.14707300    -0.37270700     2.47883800
C        1.45861200     2.00760700     1.92583300
C       -1.45861200     2.00760700     1.92583300
C       -2.36008300    -0.76683800     1.92583300
C        0.00000000    -2.48153800     1.92583300
C        2.36008300    -0.76683800     1.92583300
C        0.69534300     3.17631000     1.54839900
C       -0.69534300     3.17631000     1.54839900
C       -2.80597800     1.64284500     1.54839900
C       -3.23572300     0.32022300     1.54839900
C       -2.42953300    -2.16097600     1.54839900
C       -1.30444400    -2.97840200     1.54839900
C        1.30444400    -2.97840200     1.54839900
C        2.42953300    -2.16097600     1.54839900
C        3.23572300     0.32022300     1.54839900
C        2.80597800     1.64284500     1.54839900
H        1.21381900     4.04388400     1.14278800
H       -1.21381900     4.04388400     1.14278800
H       -3.47087200     2.40403900     1.14278800
H       -4.22105300     0.09521800     1.14278800
H       -3.35893500    -2.55810600     1.14278800
H       -1.39493500    -3.98503600     1.14278800
H        1.39493500    -3.98503600     1.14278800
H        3.35893500    -2.55810600     1.14278800
H        4.22105300     0.09521800     1.14278800
H        3.47087200     2.40403900     1.14278800
--
0 1
C        0.70881600     0.97560100    -1.15499800
C       -0.70881600     0.97560100    -1.15499800
C       -1.14688800    -0.37264600    -1.15499800
C        0.00000000    -1.20590900    -1.15499800
C        1.14688800    -0.37264600    -1.15499800
C        1.45684600     2.00517600    -1.71512200
C       -1.45684600     2.00517600    -1.71512200
C       -2.35722600    -0.76590900    -1.71512200
C        0.00000000    -2.47853400    -1.71512200
C        2.35722600    -0.76590900    -1.71512200
C        0.69580300     3.17487300    -2.09209200
C       -0.69580300     3.17487300    -2.09209200
C       -2.80446800     1.64283800    -2.09209200
C       -3.23449900     0.31934100    -2.09209200
C       -2.42906000    -2.15954300    -2.09209200
C       -1.30322700    -2.97750900    -2.09209200
C        1.30322700    -2.97750900    -2.09209200
C        2.42906000    -2.15954300    -2.09209200
C        3.23449900     0.31934100    -2.09209200
C        2.80446800     1.64283800    -2.09209200
H        1.21657800     4.04807100    -2.48428500
H       -1.21657800     4.04807100    -2.48428500
H       -3.47400100     2.40795700    -2.48428500
H       -4.22588700     0.09388800    -2.48428500
H       -3.36362900    -2.55987100    -2.48428500
H       -1.39516400    -3.99004500    -2.48428500
H        1.39516400    -3.99004500    -2.48428500
H        3.36362900    -2.55987100    -2.48428500
H        4.22588700     0.09388800    -2.48428500
H        3.47400100     2.40795700    -2.48428500
units angstrom
}
""")

CFLOW_BkybowlBkybowl_3p64 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97212500     2.41862700
C       -0.70622800     0.97212500     2.41862700
C       -1.14280400    -0.37137200     2.41849100
C        0.00000000    -1.20165400     2.41827400
C        1.14280400    -0.37137200     2.41849100
C        1.45779000     2.00650700     1.89581300
C       -1.45779000     2.00650700     1.89581300
C       -2.35873800    -0.76639200     1.89565100
C        0.00000000    -2.48003500     1.89534800
C        2.35873800    -0.76639200     1.89565100
C        0.69261800     3.17924500     1.54846400
C       -0.69261800     3.17924500     1.54846400
C       -2.80958100     1.64120300     1.54875100
C       -3.23765700     0.32374300     1.54864100
C       -2.42918200    -2.16498400     1.54865300
C       -1.30841500    -2.97916300     1.54840200
C        1.30841500    -2.97916300     1.54840200
C        2.42918200    -2.16498400     1.54865300
C        3.23765700     0.32374300     1.54864100
C        2.80958100     1.64120300     1.54875100
H        1.20851300     4.06642600     1.18749100
H       -1.20851300     4.06642600     1.18749100
H       -3.49401500     2.40602700     1.18800700
H       -4.24094400     0.10730100     1.18793900
H       -3.36816400    -2.57958300     1.18817300
H       -1.41248600    -4.00023700     1.18769900
H        1.41248600    -4.00023700     1.18769900
H        3.36816400    -2.57958300     1.18817300
H        4.24094400     0.10730100     1.18793900
H        3.49401500     2.40602700     1.18800700
--
0 1
C        0.70622800     0.97212500     6.05862700
C       -0.70622800     0.97212500     6.05862700
C       -1.14280400    -0.37137200     6.05849100
C        0.00000000    -1.20165400     6.05827400
C        1.14280400    -0.37137200     6.05849100
C        1.45779000     2.00650700     5.53581300
C       -1.45779000     2.00650700     5.53581300
C       -2.35873800    -0.76639200     5.53565100
C        0.00000000    -2.48003500     5.53534800
C        2.35873800    -0.76639200     5.53565100
C        0.69261800     3.17924500     5.18846400
C       -0.69261800     3.17924500     5.18846400
C       -2.80958100     1.64120300     5.18875100
C       -3.23765700     0.32374300     5.18864100
C       -2.42918200    -2.16498400     5.18865300
C       -1.30841500    -2.97916300     5.18840200
C        1.30841500    -2.97916300     5.18840200
C        2.42918200    -2.16498400     5.18865300
C        3.23765700     0.32374300     5.18864100
C        2.80958100     1.64120300     5.18875100
H        1.20851300     4.06642600     4.82749100
H       -1.20851300     4.06642600     4.82749100
H       -3.49401500     2.40602700     4.82800700
H       -4.24094400     0.10730100     4.82793900
H       -3.36816400    -2.57958300     4.82817300
H       -1.41248600    -4.00023700     4.82769900
H        1.41248600    -4.00023700     4.82769900
H        3.36816400    -2.57958300     4.82817300
H        4.24094400     0.10730100     4.82793900
H        3.49401500     2.40602700     4.82800700
units angstrom
}
""")

CFLOW_BkybowlBkybowl_3p74 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97212500     2.41862700
C       -0.70622800     0.97212500     2.41862700
C       -1.14280400    -0.37137200     2.41849100
C        0.00000000    -1.20165400     2.41827400
C        1.14280400    -0.37137200     2.41849100
C        1.45779000     2.00650700     1.89581300
C       -1.45779000     2.00650700     1.89581300
C       -2.35873800    -0.76639200     1.89565100
C        0.00000000    -2.48003500     1.89534800
C        2.35873800    -0.76639200     1.89565100
C        0.69261800     3.17924500     1.54846400
C       -0.69261800     3.17924500     1.54846400
C       -2.80958100     1.64120300     1.54875100
C       -3.23765700     0.32374300     1.54864100
C       -2.42918200    -2.16498400     1.54865300
C       -1.30841500    -2.97916300     1.54840200
C        1.30841500    -2.97916300     1.54840200
C        2.42918200    -2.16498400     1.54865300
C        3.23765700     0.32374300     1.54864100
C        2.80958100     1.64120300     1.54875100
H        1.20851300     4.06642600     1.18749100
H       -1.20851300     4.06642600     1.18749100
H       -3.49401500     2.40602700     1.18800700
H       -4.24094400     0.10730100     1.18793900
H       -3.36816400    -2.57958300     1.18817300
H       -1.41248600    -4.00023700     1.18769900
H        1.41248600    -4.00023700     1.18769900
H        3.36816400    -2.57958300     1.18817300
H        4.24094400     0.10730100     1.18793900
H        3.49401500     2.40602700     1.18800700
--
0 1
C        0.70622800     0.97212500     6.15862700
C       -0.70622800     0.97212500     6.15862700
C       -1.14280400    -0.37137200     6.15849100
C        0.00000000    -1.20165400     6.15827400
C        1.14280400    -0.37137200     6.15849100
C        1.45779000     2.00650700     5.63581300
C       -1.45779000     2.00650700     5.63581300
C       -2.35873800    -0.76639200     5.63565100
C        0.00000000    -2.48003500     5.63534800
C        2.35873800    -0.76639200     5.63565100
C        0.69261800     3.17924500     5.28846400
C       -0.69261800     3.17924500     5.28846400
C       -2.80958100     1.64120300     5.28875100
C       -3.23765700     0.32374300     5.28864100
C       -2.42918200    -2.16498400     5.28865300
C       -1.30841500    -2.97916300     5.28840200
C        1.30841500    -2.97916300     5.28840200
C        2.42918200    -2.16498400     5.28865300
C        3.23765700     0.32374300     5.28864100
C        2.80958100     1.64120300     5.28875100
H        1.20851300     4.06642600     4.92749100
H       -1.20851300     4.06642600     4.92749100
H       -3.49401500     2.40602700     4.92800700
H       -4.24094400     0.10730100     4.92793900
H       -3.36816400    -2.57958300     4.92817300
H       -1.41248600    -4.00023700     4.92769900
H        1.41248600    -4.00023700     4.92769900
H        3.36816400    -2.57958300     4.92817300
H        4.24094400     0.10730100     4.92793900
H        3.49401500     2.40602700     4.92800700
units angstrom
}
""")

CFLOW_BkybowlBkybowl_3p84 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97212500     2.41862700
C       -0.70622800     0.97212500     2.41862700
C       -1.14280400    -0.37137200     2.41849100
C        0.00000000    -1.20165400     2.41827400
C        1.14280400    -0.37137200     2.41849100
C        1.45779000     2.00650700     1.89581300
C       -1.45779000     2.00650700     1.89581300
C       -2.35873800    -0.76639200     1.89565100
C        0.00000000    -2.48003500     1.89534800
C        2.35873800    -0.76639200     1.89565100
C        0.69261800     3.17924500     1.54846400
C       -0.69261800     3.17924500     1.54846400
C       -2.80958100     1.64120300     1.54875100
C       -3.23765700     0.32374300     1.54864100
C       -2.42918200    -2.16498400     1.54865300
C       -1.30841500    -2.97916300     1.54840200
C        1.30841500    -2.97916300     1.54840200
C        2.42918200    -2.16498400     1.54865300
C        3.23765700     0.32374300     1.54864100
C        2.80958100     1.64120300     1.54875100
H        1.20851300     4.06642600     1.18749100
H       -1.20851300     4.06642600     1.18749100
H       -3.49401500     2.40602700     1.18800700
H       -4.24094400     0.10730100     1.18793900
H       -3.36816400    -2.57958300     1.18817300
H       -1.41248600    -4.00023700     1.18769900
H        1.41248600    -4.00023700     1.18769900
H        3.36816400    -2.57958300     1.18817300
H        4.24094400     0.10730100     1.18793900
H        3.49401500     2.40602700     1.18800700
--
0 1
C        0.70622800     0.97212500     6.25862700
C       -0.70622800     0.97212500     6.25862700
C       -1.14280400    -0.37137200     6.25849100
C        0.00000000    -1.20165400     6.25827400
C        1.14280400    -0.37137200     6.25849100
C        1.45779000     2.00650700     5.73581300
C       -1.45779000     2.00650700     5.73581300
C       -2.35873800    -0.76639200     5.73565100
C        0.00000000    -2.48003500     5.73534800
C        2.35873800    -0.76639200     5.73565100
C        0.69261800     3.17924500     5.38846400
C       -0.69261800     3.17924500     5.38846400
C       -2.80958100     1.64120300     5.38875100
C       -3.23765700     0.32374300     5.38864100
C       -2.42918200    -2.16498400     5.38865300
C       -1.30841500    -2.97916300     5.38840200
C        1.30841500    -2.97916300     5.38840200
C        2.42918200    -2.16498400     5.38865300
C        3.23765700     0.32374300     5.38864100
C        2.80958100     1.64120300     5.38875100
H        1.20851300     4.06642600     5.02749100
H       -1.20851300     4.06642600     5.02749100
H       -3.49401500     2.40602700     5.02800700
H       -4.24094400     0.10730100     5.02793900
H       -3.36816400    -2.57958300     5.02817300
H       -1.41248600    -4.00023700     5.02769900
H        1.41248600    -4.00023700     5.02769900
H        3.36816400    -2.57958300     5.02817300
H        4.24094400     0.10730100     5.02793900
H        3.49401500     2.40602700     5.02800700
units angstrom
}
""")

CFLOW_BkybowlBkybowl_3p94 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97212500     2.41862700
C       -0.70622800     0.97212500     2.41862700
C       -1.14280400    -0.37137200     2.41849100
C        0.00000000    -1.20165400     2.41827400
C        1.14280400    -0.37137200     2.41849100
C        1.45779000     2.00650700     1.89581300
C       -1.45779000     2.00650700     1.89581300
C       -2.35873800    -0.76639200     1.89565100
C        0.00000000    -2.48003500     1.89534800
C        2.35873800    -0.76639200     1.89565100
C        0.69261800     3.17924500     1.54846400
C       -0.69261800     3.17924500     1.54846400
C       -2.80958100     1.64120300     1.54875100
C       -3.23765700     0.32374300     1.54864100
C       -2.42918200    -2.16498400     1.54865300
C       -1.30841500    -2.97916300     1.54840200
C        1.30841500    -2.97916300     1.54840200
C        2.42918200    -2.16498400     1.54865300
C        3.23765700     0.32374300     1.54864100
C        2.80958100     1.64120300     1.54875100
H        1.20851300     4.06642600     1.18749100
H       -1.20851300     4.06642600     1.18749100
H       -3.49401500     2.40602700     1.18800700
H       -4.24094400     0.10730100     1.18793900
H       -3.36816400    -2.57958300     1.18817300
H       -1.41248600    -4.00023700     1.18769900
H        1.41248600    -4.00023700     1.18769900
H        3.36816400    -2.57958300     1.18817300
H        4.24094400     0.10730100     1.18793900
H        3.49401500     2.40602700     1.18800700
--
0 1
C        0.70622800     0.97212500     6.35862700
C       -0.70622800     0.97212500     6.35862700
C       -1.14280400    -0.37137200     6.35849100
C        0.00000000    -1.20165400     6.35827400
C        1.14280400    -0.37137200     6.35849100
C        1.45779000     2.00650700     5.83581300
C       -1.45779000     2.00650700     5.83581300
C       -2.35873800    -0.76639200     5.83565100
C        0.00000000    -2.48003500     5.83534800
C        2.35873800    -0.76639200     5.83565100
C        0.69261800     3.17924500     5.48846400
C       -0.69261800     3.17924500     5.48846400
C       -2.80958100     1.64120300     5.48875100
C       -3.23765700     0.32374300     5.48864100
C       -2.42918200    -2.16498400     5.48865300
C       -1.30841500    -2.97916300     5.48840200
C        1.30841500    -2.97916300     5.48840200
C        2.42918200    -2.16498400     5.48865300
C        3.23765700     0.32374300     5.48864100
C        2.80958100     1.64120300     5.48875100
H        1.20851300     4.06642600     5.12749100
H       -1.20851300     4.06642600     5.12749100
H       -3.49401500     2.40602700     5.12800700
H       -4.24094400     0.10730100     5.12793900
H       -3.36816400    -2.57958300     5.12817300
H       -1.41248600    -4.00023700     5.12769900
H        1.41248600    -4.00023700     5.12769900
H        3.36816400    -2.57958300     5.12817300
H        4.24094400     0.10730100     5.12793900
H        3.49401500     2.40602700     5.12800700
units angstrom
}
""")

CFLOW_BzBz_S_3p9 = input.process_input("""
molecule dimer {
0 1
C        1.00000000     0.00000000    -0.39150000
C        1.00000000     1.20507400     0.30425000
C        1.00000000     1.20507400     1.69575000
C        1.00000000     0.00000000     2.39150000
C        1.00000000    -1.20507400     1.69575000
C        1.00000000    -1.20507400     0.30425000
H        1.00000000     0.00000000    -1.47150000
H        1.00000000     2.14038200    -0.23575000
H        1.00000000     2.14038200     2.23575000
H        1.00000000     0.00000000     3.47150000
H        1.00000000    -2.14038200     2.23575000
H        1.00000000    -2.14038200    -0.23575000
--
0 1
C        4.90000000     0.00000000    -0.39150000
C        4.90000000    -1.20507400     0.30425000
C        4.90000000    -1.20507400     1.69575000
C        4.90000000     0.00000000     2.39150000
C        4.90000000     1.20507400     1.69575000
C        4.90000000     1.20507400     0.30425000
H        4.90000000     0.00000000    -1.47150000
H        4.90000000    -2.14038200    -0.23575000
H        4.90000000    -2.14038200     2.23575000
H        4.90000000     0.00000000     3.47150000
H        4.90000000     2.14038200     2.23575000
H        4.90000000     2.14038200    -0.23575000
units angstrom
}
""")

CFLOW_C60Bkybowl = input.process_input("""
molecule dimer {
0 1
C       -2.44476100    -3.76629900     1.42580800
C       -1.71041300    -3.33049500     2.60053500
C       -1.61107700    -4.70789600     0.69932800
C       -3.24778300    -2.85870500     0.72644000
C       -1.80581100    -2.00338100     3.03325100
C       -0.42314300    -4.00342200     2.60056800
C       -0.36181500    -4.85453300     1.42540100
C       -1.61107700    -4.70789600    -0.69932800
C       -3.24778300    -2.85870500    -0.72644000
C       -3.34824600    -1.48121500     1.17593300
C       -0.61713400    -1.29882700     3.48146300
C       -2.63977500    -1.06122400     2.30729100
C        0.72106400    -3.32482700     3.03281100
C        0.84121400    -4.99627900     0.72619800
C       -2.44476100    -3.76629900    -1.42580800
C       -0.36181500    -4.85453300    -1.42540100
C       -3.34824600    -1.48121500    -1.17593300
C       -3.40856700    -0.62974400     0.00000000
C       -0.71652100     0.07949700     3.03452100
C        0.62236800    -1.94724500     3.48199500
C       -1.96672300     0.22646500     2.30843400
C        1.97091600    -3.47212300     2.30712300
C        2.02938800    -4.29196400     1.17524900
C        0.84121400    -4.99627900    -0.72619800
C       -1.71041300    -3.33049500    -2.60053500
C       -0.42314300    -4.00342200    -2.60056800
C       -2.63977500    -1.06122400    -2.30729100
C       -2.75699600     0.60742000     0.00000000
C        0.42766000     0.75551000     2.59901600
C        1.81124200    -1.24340900     3.03441800
C       -2.02080700     1.04039200     1.17358000
C        2.64493700    -2.18551800     2.30759300
C        2.76375600    -3.85694600     0.00000000
C        2.02938800    -4.29196400    -1.17524900
C       -1.80581100    -2.00338100    -3.03325100
C        0.72106400    -3.32482700    -3.03281100
C       -1.96672300     0.22646500    -2.30843400
C       -2.02080700     1.04039200    -1.17358000
C        1.71530400     0.08347500     2.60120600
C        0.36629300     1.60016300     1.42136800
C       -0.83328100     1.74007300     0.72347400
C        3.35260400    -1.76620300     1.17561900
C        3.41301600    -2.61807500     0.00000000
C        1.97091600    -3.47212300    -2.30712300
C       -0.61713400    -1.29882700    -3.48146300
C        0.62236800    -1.94724500    -3.48199500
C       -0.71652100     0.07949700    -3.03452100
C       -0.83328100     1.74007300    -0.72347400
C        2.44815900     0.51675300     1.42446800
C        1.61170500     1.45224300     0.69742700
C        3.25434100    -0.38835500     0.72661400
C        3.35260400    -1.76620300    -1.17561900
C        2.64493700    -2.18551800    -2.30759300
C        1.81124200    -1.24340900    -3.03441800
C        0.42766000     0.75551000    -2.59901600
C        0.36629300     1.60016300    -1.42136800
C        1.61170500     1.45224300    -0.69742700
C        3.25434100    -0.38835500    -0.72661400
C        1.71530400     0.08347500    -2.60120600
C        2.44815900     0.51675300    -1.42446800
--
0 1
C       -3.19032700     4.20143700    -0.69626700
C       -3.19032700     4.20143700     0.69626700
C       -2.02031900     4.57636100    -1.46287500
C       -1.65204100     4.22202600    -2.81770800
C       -0.32763200     4.23532900    -3.24742400
C        0.76021100     4.60662800    -2.36621300
C        2.16499200     4.26248900    -2.43659000
C        2.98394600     4.27348900    -1.30950700
C        2.48043100     4.63010100     0.00000000
C        2.98394600     4.27348900     1.30950700
C        2.16499200     4.26248900     2.43659000
C        0.76021100     4.60662800     2.36621300
C       -0.32763200     4.23532900     3.24742400
C       -1.65204100     4.22202600     2.81770800
C       -2.02031900     4.57636100     1.46287500
C       -0.98839100     5.12111700    -0.70843400
C        0.35840400     5.13651000    -1.14620400
C        1.19094000     5.14556000     0.00000000
C        0.35840400     5.13651000     1.14620400
C       -0.98839100     5.12111700     0.70843400
H       -4.06441800     3.80879000    -1.21377600
H       -4.06441800     3.80879000     1.21377600
H       -2.41195000     3.83189900    -3.49341900
H       -0.10220000     3.85419200    -4.24246700
H        2.57558400     3.88122400    -3.37072200
H        4.00276200     3.90074000    -1.40663400
H        4.00276200     3.90074000     1.40663400
H        2.57558400     3.88122400     3.37072200
H       -0.10220000     3.85419200     4.24246700
H       -2.41195000     3.83189900     3.49341900
units angstrom
}
""")

CFLOW_C60Bkycatch = input.process_input("""
molecule dimer {
0 1
C        0.66331015     2.47224543    -3.46907102
C        1.42015617     1.32292559    -3.00396243
C        2.57041166     1.81371387    -2.26848602
C        2.52992682     3.26434443    -2.28122351
C        3.00401714     1.13997115    -1.12532844
C        0.75092534     0.17313602    -2.56896973
C        1.34970151     3.67265183    -3.02198754
C       -0.73487266     2.43006841    -3.48162145
C        0.61293558     4.78597150    -2.60443253
C        2.93104019     3.98505989    -1.15109519
C        3.41951893     1.88565973     0.04424416
C        2.31290756    -0.04997921    -0.67658403
C        3.38580796     3.28115949     0.03353679
C        2.16608088     5.14067139    -0.71795451
C        1.02873635     5.53386529    -1.43083006
C       -0.83895727     4.74223871    -2.61750965
C       -1.50000158     3.58706631    -3.04798368
C       -1.42958972     1.23720713    -3.02984494
C       -0.70079471     0.12932438    -2.58209598
C        1.20740890    -0.52784143    -1.38299523
C        2.98266916     1.15706337     1.21677214
C        2.29952642    -0.03958154     0.77227513
C       -1.13558984    -0.59857808    -1.40416302
C        0.04303279    -1.01237272    -0.66599530
C       -2.66749383     3.10869624    -2.32900684
C       -2.62049531     1.65808299    -2.31588231
C       -1.31985524     5.46298442    -1.45206864
C       -0.16581543     5.95236463    -0.71840803
C        2.15286041     5.15191243     0.73517392
C        2.90920325     4.00251188     1.19917880
C        2.48763349     3.29814622     2.33234032
C        2.52780567     1.84727561     2.34162038
C        1.18111880    -0.50768798     1.46504110
C        0.02993604    -1.00226515     0.73319752
C       -3.03322065     0.95892578    -1.18037564
C       -2.28010505    -0.18753743    -0.71828150
C       -2.44484538     5.00276677    -0.75992948
C       -3.13123381     3.80331865    -1.20620187
C       -0.17891833     5.96313993     0.68028494
C        1.00200512     5.55454101     1.42075848
C        1.36386624     1.36716118     3.06263657
C        0.70261283     0.21108629     2.63199472
C       -1.16182217    -0.57752890     1.44348217
C       -2.29366335    -0.17710274     0.73063368
C       -3.51391707     1.67755719    -0.01903378
C       -3.56435945     3.07266849    -0.02959753
C       -2.45823180     5.01338579     0.69323184
C       -1.34661889     5.48448680     1.39963613
C        0.56463339     4.82416266     2.59737931
C        1.29339986     3.71692550     3.04448233
C       -0.88724635     4.78056517     2.58420771
C        0.59866875     2.52329539     3.49615878
C       -0.74897097     0.16745545     2.61860472
C       -3.05485921     0.97622590     1.16151858
C       -3.15281559     3.82063126     1.14401920
C       -2.70996852     3.14249831     2.28491906
C       -1.55640556     3.63153104     3.01887107
C       -0.79950548     2.48140954     3.48332736
C       -1.48578247     1.28191228     3.03625100
C       -2.66300775     1.69198328     2.29408228
--
0 1
C       -1.19705112    -6.78248429     2.35733806
C       -0.54183599    -5.86353674     1.51759331
C        0.86877714    -5.82567057     1.51551975
C        1.57530835    -6.70795878     2.35282537
C        0.91318038    -7.62330199     3.17160076
C       -0.48391191    -7.66068296     3.17395987
C        1.65636644    -4.85143802     0.71258456
C        2.53466306    -3.99607419     1.38230751
C        3.44173404    -3.15497768     0.71952324
C        3.43347840    -3.15627326    -0.73434261
C        2.52308649    -4.00274019    -1.38517267
C        1.65254533    -4.85662969    -0.70389877
C        4.33866685    -2.25977428     1.46268916
C        5.29893981    -1.63147919     0.69353296
C        5.28849745    -1.63022409    -0.73419220
C        4.31763559    -2.25848574    -1.48948521
C        5.96639205    -0.46906645    -1.17165663
C        6.39270574     0.24750766    -0.02707703
C        5.98483765    -0.47201511     1.12242424
C        6.58728953     1.62256504    -0.02651167
C        6.52216991     2.23829035    -1.33564482
C        6.10796161     1.53401343    -2.46244049
C        5.71333201     0.14267748    -2.39453033
C        4.84959509    -0.61209596    -3.27015462
C        4.18068986    -1.75590298    -2.83519980
C        6.54676159     2.23413884     1.28569334
C        6.15106671     1.52676026     2.41717370
C        5.75355853     0.13624687     2.35142539
C        4.90510306    -0.62072419     3.23997133
C        4.22670402    -1.76165226     2.81251158
C        0.86571310    -5.84104641    -1.49472383
C       -0.54438716    -5.88234817    -1.49266002
C       -1.20011445    -6.81314515    -2.31874879
C       -0.48685981    -7.70048225    -3.12531442
C        0.91022067    -7.65955460    -3.12749817
C        1.57257682    -6.73206950    -2.32259123
C       -1.38387740    -4.93529874     0.71538831
C       -2.30574161    -4.12797099     1.38643934
C       -3.26027551    -3.34024969     0.72478778
C       -3.25759718    -3.34562083    -0.72918632
C       -2.30370988    -4.14145572    -1.38152079
C       -1.38460169    -4.94383955    -0.70103552
C       -4.19374423    -2.50128573    -1.48319443
C       -5.19496250    -1.92539082    -0.72542320
C       -5.19945096    -1.92195629     0.70223145
C       -4.20205710    -2.49296788     1.46860355
C       -5.94731322    -0.80127853     1.13065456
C       -6.40010694    -0.11009739    -0.01972456
C       -5.93965120    -0.80649769    -1.16394564
C       -4.09234542    -1.99803930    -2.83168952
C       -4.82592114    -0.89494197    -3.26741384
C       -5.72605993    -0.18609697    -2.39018661
C       -6.19636273     1.18118633    -2.46124119
C       -6.64449582     1.86563842    -1.33531292
C       -6.66913664     1.25233479    -0.02360617
C       -4.11043606    -1.98413115     2.81583531
C       -4.84802453    -0.88008772     3.24208403
C       -5.74247476    -0.17547413     2.35567384
C       -6.21393303     1.19185429     2.41759528
C       -6.65302998     1.87152813     1.28535514
H        4.62000664    -0.22505906    -4.26169882
H        3.45469350    -2.20961902    -3.50677668
H        3.51220909    -2.21687142     3.49537998
H        4.69432556    -0.23734674     4.23720302
H        6.04077435     2.06950397     3.35486993
H        6.73054921     3.30332902     1.38169553
H        6.70293892     3.30801682    -1.43139928
H        5.98002372     2.07988699    -3.39606436
H       -6.10314216     1.72961497    -3.39750013
H       -6.88409228     2.92342563    -1.43388599
H       -6.89314925     2.92981266     1.37736599
H       -6.12671465     1.74432777     3.35205891
H       -4.65256863    -0.48054600     4.23598684
H       -3.36791024    -2.39592116     3.49632750
H       -3.34579680    -2.41359900    -3.50550755
H       -4.62325143    -0.49964918    -4.26160501
H        2.53813733    -4.01188964     2.46976786
H        2.51823415    -4.02462365    -2.47254584
H       -2.30437970    -4.14092364     2.47394371
H       -2.30081967    -4.16511905    -2.46886053
H        2.65994015    -6.69147901    -2.32006670
H        1.48177936    -8.34528562    -3.74980546
H       -1.01911751    -8.41865672    -3.74589702
H       -2.28795954    -6.83617949    -2.31329457
H        2.66270620    -6.67021288     2.34680032
H        1.48490256    -8.30236320     3.80101454
H       -1.01634978    -8.36944419     3.80511677
H       -2.28490570    -6.80366055     2.35468970
units angstrom
}
""")

# <<< Geometry Specification Strings >>>
GEOS = {}
for rxn in HRXN:

   if rxn is 'C60Bkybowl' or rxn is 'C60Bkycatch':
      GEOS['%s-%s-dimer'    % (dbse, rxn)] = eval('%s_%s' % (dbse, rxn ))
      GEOS['%s-%s-monoA-CP' % (dbse, rxn)] = eval('%s_%s' % (dbse, rxn )) + monoA_CP
      GEOS['%s-%s-monoB-CP' % (dbse, rxn)] = eval('%s_%s' % (dbse, rxn )) + monoB_CP
      GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = eval('%s_%s' % (dbse, rxn )) + monoA_unCP
      GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = eval('%s_%s' % (dbse, rxn )) + monoB_unCP

   else:
      distance = rxnpattern.match(rxn)
      GEOS['%s-%s-dimer'    % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))
      GEOS['%s-%s-monoA-CP' % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) ))) + monoA_CP
      GEOS['%s-%s-monoB-CP' % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) ))) + monoB_CP
      GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) ))) + monoA_unCP
      GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) ))) + monoB_unCP

