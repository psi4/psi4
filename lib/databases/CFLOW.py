"""
| Database of extended conjugated bimolecular systems.
| Geometries and Reference interaction energies from the following articles:
|   Polyene geometries from Marshall et al. JCTC 6 3681 (2010).
|   Polyene reference interaction energies from Sherrill group by ccsd(t**)-f12b/heavy-aug-cc-pvdz.
|   Acene geometries (except benzene) from Sherrill group by df-mp2/cc-pvtz c.2011.
|   Benzene geometry from NBC10 database and citations therein.
|   Acene reference interaction energies (incl. benzene dimer) from Sherrill group by ccsd(t**)-f12b/aug-cc-pvdz.
|   Buckybowl (Pulay-labeled) geometries from Sherrill group by PBE1PBE/6-31G*, following Pulay's instructions in Janowski et al. CPL 512 155 (2011).
|   Buckybowl (Pulay-labeled) reference interaction energies from Janowski et al. CPL 512 155 (2011).
|   Buckyware (Grimme-labeled) geometries from Grimme PCCP 12 7091 (2010).
|   Buckyware (Grimme-labeled) reference interaction energies from Grimme PCCP 12 7091 (2010) by B97-D2/TZVP.
|   Collection into CFLOW by Parrish et al. XXX XXX XXXXXX (2012).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'equilibrium'``
  - ``'Polyenes'`` equilibrium for linear polyene dimers for 2 through 16 monomer carbons
  - ``'cBzBz'`` 5-point dissociation curve for benzene dimer
  - ``'c2BzBz'`` 5-point dissociation curve for napthalene-benzene complex
  - ``'c2Bz2Bz'`` 5-point dissociation curve for napthalene dimer
  - ``'c3Bz2Bz'`` 5-point dissociation curve for anthracene-napthalene complex
  - ``'c3Bz3Bz'`` 5-point dissociation curve for anthracene dimer
  - ``'c4Bz3Bz'`` 5-point dissociation curve for tetracene-anthracene complex
  - ``'Arenes'`` equilibrium for benzene dimer through tetracene-anthracene complex linear arenes
  - ``'cArenes'`` 5-point curves around benzene dimer through tetracene-anthracene complex linear arenes
  - ``'cPulay'`` 4-point dissociation curve for bowl-in-bowl corannulene dimer
  - ``'Pulay'`` Pulay bowl-in-bowl corannulene dimer dissociation curve and extra point
  - ``'Grimme60'`` Grimme corannulene dimer, C60 @ buckybowl, and C60 @ buckycatcher
  - ``'Grimme70'`` Grimme C70 @ buckycatcher at three orientations
  - ``'Paper'`` linear polyene dimers, equilibrium arene complexes, Pulay corannulene dimer curve, and Grimme corannulene dimer and C60 complexes
  - ``'cPaper'`` linear polyene dimers, arene complex curves, Pulay corannulene dimer curve, and Grimme corannulene dimer and C60 complexes

"""
import re
import input

# <<< CFLOW Database Module >>>
dbse = 'CFLOW'

# <<< Database Members >>>
HRXN_SM = ['2Ae2Ae-3.8', 'BzBz_S-3.9']
HRXN_LG = ['C70Bkycatch']

# Polyenes
Polyenes = ['2Ae2Ae-3.8', '4Ae4Ae-3.8', '6Ae6Ae-3.8', '8Ae8Ae-3.8', '10Ae10Ae-3.8', '12Ae12Ae-3.8', '14Ae14Ae-3.8', '16Ae16Ae-3.8']

# Arenes
# geometries are flexible; dist arrays may be filled with any positive intermonomer distance
cBzBz = []
dist = [3.7, 3.8, 3.9, 4.0, 4.1]
for d in dist:
    cBzBz.append('BzBz_S-' + str(d))
c2BzBz = []
dist = [3.3, 3.4, 3.5, 3.6, 3.7]
for d in dist:
    c2BzBz.append('2BzBz_S-' + str(d))
c2Bz2Bz = []
dist = [3.6, 3.7, 3.8, 3.9, 4.0]
for d in dist:
    c2Bz2Bz.append('2Bz2Bz_S-' + str(d))
c3Bz2Bz = []
dist = [3.3, 3.4, 3.5, 3.6, 3.7]
for d in dist:
    c3Bz2Bz.append('3Bz2Bz_S-' + str(d))
c3Bz3Bz = []
dist = [3.5, 3.6, 3.7, 3.8, 3.9]
for d in dist:
    c3Bz3Bz.append('3Bz3Bz_S-' + str(d))
c4Bz3Bz = []
dist = [3.3, 3.4, 3.5, 3.6, 3.7]
for d in dist:
    c4Bz3Bz.append('4Bz3Bz_S-' + str(d))
Arenes = ['BzBz_S-3.9', '2BzBz_S-3.5', '2Bz2Bz_S-3.8', '3Bz2Bz_S-3.5', '3Bz3Bz_S-3.7', '4Bz3Bz_S-3.5']
temp = [cBzBz, c2BzBz, c2Bz2Bz, c3Bz2Bz, c3Bz3Bz, c4Bz3Bz]
cArenes = sum(temp, [])

# Pulay Buckyware
cPulay = ['BkybowlBkybowl-3.54', 'BkybowlBkybowl-3.64', 'BkybowlBkybowl-3.74', 'BkybowlBkybowl-3.84']
temp = [cPulay, ['BkybowlBkybowl-3.73']]
Pulay = sum(temp, [])

# Grimme Buckyware
Grimme60 = ['BkybowlBkybowl-3.63', 'C60Bkybowl', 'C60Bkycatch', ]
Grimme70 = ['C70Bkycatch', 'C70Bkycatch_W', 'C70Bkycatch_T', ]

# Aggregates
temp = [Polyenes, Arenes, ['BkybowlBkybowl-3.73'], Grimme60, ['C70Bkycatch']]
HRXN_EQ = sum(temp, [])
temp = [Polyenes, Arenes, Pulay, Grimme60]
Paper = sum(temp, [])
temp = [Polyenes, cArenes, Pulay, Grimme60]
cPaper = sum(temp, [])
temp = [Polyenes, cArenes, Pulay, Grimme60, Grimme70]
HRXN = sum(temp, [])

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supermolecular calculations
for rxn in HRXN:

    if rxn in Polyenes:  # homomonomer, CP-symmetric, w/o db monomer reuse
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'            % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'         % (dbse, rxn) : -2,
                                          '%s-%s-monoA-unCP'       % (dbse, rxn) : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-CP'         % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-unCP'       % (dbse, rxn) ]

    elif (rxn in cBzBz) or (rxn in c2Bz2Bz) or (rxn in c3Bz3Bz):  # homomonomer, CP-symmetric, w/ db monomer reuse
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'            % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'         % (dbse, rxn) : -2,
                                          '%s-Bz-mono-unCP'        % (dbse)      : -2,
                                          '%s-2Bz-mono-unCP'       % (dbse)      : -2,
                                          '%s-3Bz-mono-unCP'       % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-CP'         % (dbse, rxn) ]

        if rxn in cBzBz:
            ACTV['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'           % (dbse, rxn),
                                           '%s-Bz-mono-unCP'       % (dbse) ]
        elif rxn in c2Bz2Bz:
            ACTV['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'           % (dbse, rxn),
                                           '%s-2Bz-mono-unCP'      % (dbse) ]
        elif rxn in c3Bz3Bz:
            ACTV['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'           % (dbse, rxn),
                                           '%s-3Bz-mono-unCP'      % (dbse) ]

    elif (rxn in c2BzBz) or (rxn in c3Bz2Bz) or (rxn in c4Bz3Bz):  # heteromonomer, w/ db monomer reuse
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'            % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'         % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'         % (dbse, rxn) : -1,
                                          '%s-Bz-mono-unCP'        % (dbse)      : -1,
                                          '%s-2Bz-mono-unCP'       % (dbse)      : -1,
                                          '%s-3Bz-mono-unCP'       % (dbse)      : -1,
                                          '%s-4Bz-mono-unCP'       % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-CP'         % (dbse, rxn),
                                          '%s-%s-monoB-CP'         % (dbse, rxn) ]

        if rxn in c2BzBz:
            ACTV['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'           % (dbse, rxn),
                                           '%s-2Bz-mono-unCP'      % (dbse),
                                           '%s-Bz-mono-unCP'       % (dbse) ]
        elif rxn in c3Bz2Bz:
            ACTV['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'           % (dbse, rxn),
                                           '%s-3Bz-mono-unCP'      % (dbse),
                                           '%s-2Bz-mono-unCP'      % (dbse) ]
        elif rxn in c4Bz3Bz:
            ACTV['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'           % (dbse, rxn),
                                           '%s-4Bz-mono-unCP'      % (dbse),
                                           '%s-3Bz-mono-unCP'      % (dbse) ]

    elif rxn in Pulay:  # homomonomer, not-CP-symmetric, w/ db monomer reuse
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'            % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'         % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'         % (dbse, rxn) : -1,
                                          '%s-Bkybowl_P-mono-unCP' % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-CP'         % (dbse, rxn),
                                          '%s-%s-monoB-CP'         % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-Bkybowl_P-mono-unCP' % (dbse) ]

    elif rxn is 'BkybowlBkybowl-3.63':  # homomonomer, not-CP-symmetric, w/o db monomer reuse
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'            % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'         % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'         % (dbse, rxn) : -1,
                                          '%s-%s-monoA-unCP'       % (dbse, rxn) : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-CP'         % (dbse, rxn),
                                          '%s-%s-monoB-CP'         % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-unCP'       % (dbse, rxn) ]

    elif (rxn in Grimme60) or (rxn in Grimme70):  # heteromonomer, w/o db monomer reuse
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'            % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'         % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'         % (dbse, rxn) : -1,
                                          '%s-%s-monoA-unCP'       % (dbse, rxn) : -1,
                                          '%s-%s-monoB-unCP'       % (dbse, rxn) : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-CP'         % (dbse, rxn),
                                          '%s-%s-monoB-CP'         % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'            % (dbse, rxn),
                                          '%s-%s-monoA-unCP'       % (dbse, rxn),
                                          '%s-%s-monoB-unCP'       % (dbse, rxn) ]

# <<< Reference Values [kcal/mol] >>>
BIND = {}
nan = float('NaN')
# Polyenes
BIND['%s-%s'            % (dbse, '2Ae2Ae-3.8'            )] =    0.266
BIND['%s-%s'            % (dbse, '4Ae4Ae-3.8'            )] =   -0.373
BIND['%s-%s'            % (dbse, '6Ae6Ae-3.8'            )] =   -1.078
BIND['%s-%s'            % (dbse, '8Ae8Ae-3.8'            )] =   -1.796
BIND['%s-%s'            % (dbse, '10Ae10Ae-3.8'          )] =   -2.519
BIND['%s-%s'            % (dbse, '12Ae12Ae-3.8'          )] =    nan    # may exist, check with MSM
BIND['%s-%s'            % (dbse, '14Ae14Ae-3.8'          )] =    nan
BIND['%s-%s'            % (dbse, '16Ae16Ae-3.8'          )] =    nan
# Arenes
for item in cArenes:  # backstop to allow arbitrary intermonomer geometries
    BIND['%s-%s' % (dbse, item)] = nan
BIND['%s-%s'            % (dbse, 'BzBz_S-3.7'            )] =   -1.623
BIND['%s-%s'            % (dbse, 'BzBz_S-3.8'            )] =   -1.765
BIND['%s-%s'            % (dbse, 'BzBz_S-3.9'            )] =   -1.793  # BzBz minimum
BIND['%s-%s'            % (dbse, 'BzBz_S-4.0'            )] =   -1.746
BIND['%s-%s'            % (dbse, 'BzBz_S-4.1'            )] =   -1.651
BIND['%s-%s'            % (dbse, '2BzBz_S-3.3'           )] =    nan
BIND['%s-%s'            % (dbse, '2BzBz_S-3.4'           )] =    nan
BIND['%s-%s'            % (dbse, '2BzBz_S-3.5'           )] =    nan    # 2BzBz unconfirmed minimum
BIND['%s-%s'            % (dbse, '2BzBz_S-3.6'           )] =    nan
BIND['%s-%s'            % (dbse, '2BzBz_S-3.7'           )] =    nan
BIND['%s-%s'            % (dbse, '2Bz2Bz_S-3.6'          )] =    nan
BIND['%s-%s'            % (dbse, '2Bz2Bz_S-3.7'          )] =   -4.314
BIND['%s-%s'            % (dbse, '2Bz2Bz_S-3.8'          )] =   -4.384  # 2Bz2Bz minimum
BIND['%s-%s'            % (dbse, '2Bz2Bz_S-3.9'          )] =   -4.283
BIND['%s-%s'            % (dbse, '2Bz2Bz_S-4.0'          )] =    nan
BIND['%s-%s'            % (dbse, '3Bz2Bz_S-3.3'          )] =    nan
BIND['%s-%s'            % (dbse, '3Bz2Bz_S-3.4'          )] =   -7.607
BIND['%s-%s'            % (dbse, '3Bz2Bz_S-3.5'          )] =   -7.802  # 3Bz2Bz minimum
BIND['%s-%s'            % (dbse, '3Bz2Bz_S-3.6'          )] =   -7.680
BIND['%s-%s'            % (dbse, '3Bz2Bz_S-3.7'          )] =    nan
BIND['%s-%s'            % (dbse, '3Bz3Bz_S-3.5'          )] =    nan
BIND['%s-%s'            % (dbse, '3Bz3Bz_S-3.6'          )] =    nan
BIND['%s-%s'            % (dbse, '3Bz3Bz_S-3.7'          )] =    nan    # 3Bz3Bz unconfirmed minimum
BIND['%s-%s'            % (dbse, '3Bz3Bz_S-3.8'          )] =    nan
BIND['%s-%s'            % (dbse, '3Bz3Bz_S-3.9'          )] =    nan
BIND['%s-%s'            % (dbse, '4Bz3Bz_S-3.3'          )] =    nan
BIND['%s-%s'            % (dbse, '4Bz3Bz_S-3.4'          )] =    nan
BIND['%s-%s'            % (dbse, '4Bz3Bz_S-3.5'          )] =    nan    # 4Bz3Bz unconfirmed minimum
BIND['%s-%s'            % (dbse, '4Bz3Bz_S-3.6'          )] =    nan
BIND['%s-%s'            % (dbse, '4Bz3Bz_S-3.7'          )] =    nan
# Pulay Buckyware
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.54'   )] =  -14.8
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.64'   )] =  -15.4
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.74'   )] =  -15.4
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.84'   )] =  -15.0
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.73'   )] =  -15.6  # bootstrapped, Pulay does not report
# Grimme Buckyware
BIND['%s-%s'            % (dbse, 'BkybowlBkybowl-3.63'   )] =  -17.0
BIND['%s-%s'            % (dbse, 'C60Bkybowl'            )] =  -19.5
BIND['%s-%s'            % (dbse, 'C60Bkycatch'           )] =  -43.1
BIND['%s-%s'            % (dbse, 'C70Bkycatch'           )] =  -45.1
BIND['%s-%s'            % (dbse, 'C70Bkycatch_W'         )] =  -44.0
BIND['%s-%s'            % (dbse, 'C70Bkycatch_T'         )] =  -40.3

# <<< Comment Lines >>>
TAGL = {}
rxnpattern = re.compile(r'^(.+)-(.+)$')

# Polyenes
TAGL['%s-%s'            % (dbse,  '2Ae2Ae-3.8'  )] = 'Ethene (C2H4) Dimer, stacked at 3.8 A'
TAGL['%s-%s-dimer'      % (dbse,  '2Ae2Ae-3.8'  )] = 'Ethene (C2H4) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-CP'   % (dbse,  '2Ae2Ae-3.8'  )] = 'Ethene from Ethene (C2H4) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-unCP' % (dbse,  '2Ae2Ae-3.8'  )] = 'Ethene from Ethene (C2H4) Dimer, stacked at 3.8 A'
TAGL['%s-%s'            % (dbse,  '4Ae4Ae-3.8'  )] = 'Butadiene (C4H6) Dimer, stacked at 3.8 A'
TAGL['%s-%s-dimer'      % (dbse,  '4Ae4Ae-3.8'  )] = 'Butadiene (C4H6) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-CP'   % (dbse,  '4Ae4Ae-3.8'  )] = 'Butadiene from Butadiene (C4H6) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-unCP' % (dbse,  '4Ae4Ae-3.8'  )] = 'Butadiene from Butadiene (C4H6) Dimer, stacked at 3.8 A'
TAGL['%s-%s'            % (dbse,  '6Ae6Ae-3.8'  )] = 'Hexatriene (C6H8) Dimer, stacked at 3.8 A'
TAGL['%s-%s-dimer'      % (dbse,  '6Ae6Ae-3.8'  )] = 'Hexatriene (C6H8) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-CP'   % (dbse,  '6Ae6Ae-3.8'  )] = 'Hexatriene from Hexatriene (C6H8) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-unCP' % (dbse,  '6Ae6Ae-3.8'  )] = 'Hexatriene from Hexatriene (C6H8) Dimer, stacked at 3.8 A'
TAGL['%s-%s'            % (dbse,  '8Ae8Ae-3.8'  )] = 'Octatetraene (C8H10) Dimer, stacked at 3.8 A'
TAGL['%s-%s-dimer'      % (dbse,  '8Ae8Ae-3.8'  )] = 'Octatetraene (C8H10) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-CP'   % (dbse,  '8Ae8Ae-3.8'  )] = 'Octatetraene from Octatetraene (C8H10) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-unCP' % (dbse,  '8Ae8Ae-3.8'  )] = 'Octatetraene from Octatetraene (C8H10) Dimer, stacked at 3.8 A'
TAGL['%s-%s'            % (dbse,  '10Ae10Ae-3.8')] = 'Decapentaene (C10H12) Dimer, stacked at 3.8 A'
TAGL['%s-%s-dimer'      % (dbse,  '10Ae10Ae-3.8')] = 'Decapentaene (C10H12) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-CP'   % (dbse,  '10Ae10Ae-3.8')] = 'Decapentaene from Decapentaene (C10H12) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-unCP' % (dbse,  '10Ae10Ae-3.8')] = 'Decapentaene from Decapentaene (C10H12) Dimer, stacked at 3.8 A'
TAGL['%s-%s'            % (dbse,  '12Ae12Ae-3.8')] = 'Dodecahexaene (C12H14) Dimer, stacked at 3.8 A'
TAGL['%s-%s-dimer'      % (dbse,  '12Ae12Ae-3.8')] = 'Dodecahexaene (C12H14) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-CP'   % (dbse,  '12Ae12Ae-3.8')] = 'Dodecahexaene from Dodecahexaene (C12H14) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-unCP' % (dbse,  '12Ae12Ae-3.8')] = 'Dodecahexaene from Dodecahexaene (C12H14) Dimer, stacked at 3.8 A'
TAGL['%s-%s'            % (dbse,  '14Ae14Ae-3.8')] = 'Tetradecaheptaene (C14H16) Dimer, stacked at 3.8 A'
TAGL['%s-%s-dimer'      % (dbse,  '14Ae14Ae-3.8')] = 'Tetradecaheptaene (C14H16) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-CP'   % (dbse,  '14Ae14Ae-3.8')] = 'Tetradecaheptaene from Tetradecaheptaene (C14H16) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-unCP' % (dbse,  '14Ae14Ae-3.8')] = 'Tetradecaheptaene from Tetradecaheptaene (C14H16) Dimer, stacked at 3.8 A'
TAGL['%s-%s'            % (dbse,  '16Ae16Ae-3.8')] = 'Hexadecaoctaene (C16H18) Dimer, stacked at 3.8 A'
TAGL['%s-%s-dimer'      % (dbse,  '16Ae16Ae-3.8')] = 'Hexadecaoctaene (C16H18) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-CP'   % (dbse,  '16Ae16Ae-3.8')] = 'Hexadecaoctaene from Hexadecaoctaene (C16H18) Dimer, stacked at 3.8 A'
TAGL['%s-%s-monoA-unCP' % (dbse,  '16Ae16Ae-3.8')] = 'Hexadecaoctaene from Hexadecaoctaene (C16H18) Dimer, stacked at 3.8 A'

# Arenes
for item in cBzBz:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Benzene Dimer, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Benzene Dimer, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Benzene Dimer, stacked at %s A' % (distance.group(2))

for item in c2BzBz:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Napthalene-Benzene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Napthalene-Benzene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Napthalene from Napthalene-Benzene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Benzene from Napthalene-Benzene Complex, stacked at %s A' % (distance.group(2))

for item in c2Bz2Bz:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Napthalene Dimer, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Napthalene Dimer, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Napthalene from Napthalene Dimer, stacked at %s A' % (distance.group(2))

for item in c3Bz2Bz:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Anthracene-Napthalene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Anthracene-Napthalene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Anthracene from Anthracene-Napthalene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Napthalene from Anthracene-Napthalene Complex, stacked at %s A' % (distance.group(2))

for item in c3Bz3Bz:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Anthracene Dimer, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Anthracene Dimer, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Anthracene from Anthracene Dimer, stacked at %s A' % (distance.group(2))

for item in c4Bz3Bz:  # DON'T YET HAVE GEOMETRY FOR 4Bz
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Tetracene-Anthracene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Tetracene-Anthracene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Tetracene from Tetracene-Anthracene Complex, stacked at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Anthracene from Tetracene-Anthracene Complex, stacked at %s A' % (distance.group(2))

# Pulay Buckyware
for item in Pulay:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Corannulene Dimer, stacked at %s A, Pulay geometry' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Corannulene Dimer, stacked at %s A, Pulay geometry' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Corannulene from Corannulene Dimer, stacked at %s A, Pulay geometry' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Corannulene from Corannulene Dimer, stacked at %s A, Pulay geometry' % (distance.group(2))

# Grimme Buckyware
TAGL['%s-%s'            % (dbse, 'BkybowlBkybowl-3.63'   )] = """Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-dimer'      % (dbse, 'BkybowlBkybowl-3.63'   )] = """Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-monoA-CP'   % (dbse, 'BkybowlBkybowl-3.63'   )] = """Corannulene from Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-monoB-CP'   % (dbse, 'BkybowlBkybowl-3.63'   )] = """Corannulene from Corannulene Dimer, stacked at ~3.63 A, Grimme geometry """
TAGL['%s-%s-monoA-unCP' % (dbse, 'BkybowlBkybowl-3.63'   )] = """Corannulene, Grimme geometry """

TAGL['%s-%s'            % (dbse, 'C60Bkybowl'            )] = """C60 @ Corannulene Buckybowl """
TAGL['%s-%s-dimer'      % (dbse, 'C60Bkybowl'            )] = """C60 @ Corannulene Buckybowl """
TAGL['%s-%s-monoA-CP'   % (dbse, 'C60Bkybowl'            )] = """Buckyball from C60 @ Corannulene Buckybowl """
TAGL['%s-%s-monoB-CP'   % (dbse, 'C60Bkybowl'            )] = """Corannulene from C60 @ Corannulene Buckybowl """
TAGL['%s-%s-monoA-unCP' % (dbse, 'C60Bkybowl'            )] = """Buckyball from C60 @ Corannulene Buckybowl """
TAGL['%s-%s-monoB-unCP' % (dbse, 'C60Bkybowl'            )] = """Corannulene from C60 @ Corannulene Buckybowl """

TAGL['%s-%s'            % (dbse, 'C60Bkycatch'           )] = """C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-dimer'      % (dbse, 'C60Bkycatch'           )] = """C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-monoA-CP'   % (dbse, 'C60Bkycatch'           )] = """Buckyball from C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-monoB-CP'   % (dbse, 'C60Bkycatch'           )] = """Buckycatcher from C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-monoA-unCP' % (dbse, 'C60Bkycatch'           )] = """Buckyball from C60 @ C60H28 Buckycatcher """
TAGL['%s-%s-monoB-unCP' % (dbse, 'C60Bkycatch'           )] = """Buckycatcher from C60 @ C60H28 Buckycatcher """

TAGL['%s-%s'            % (dbse, 'C70Bkycatch'           )] = """C70 @ C60H28 Buckycatcher, ball tilted out of maw (min) """
TAGL['%s-%s-dimer'      % (dbse, 'C70Bkycatch'           )] = """C70 @ C60H28 Buckycatcher, ball tilted out of maw (min) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'C70Bkycatch'           )] = """Buckyball from C70 @ C60H28 Buckycatcher, ball tilted out of maw (min) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'C70Bkycatch'           )] = """Buckycatcher from C70 @ C60H28 Buckycatcher, ball tilted out of maw (min) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'C70Bkycatch'           )] = """Buckyball from C70 @ C60H28 Buckycatcher, ball tilted out of maw (min) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'C70Bkycatch'           )] = """Buckycatcher from C70 @ C60H28 Buckycatcher, ball tilted out of maw (min) """

TAGL['%s-%s'            % (dbse, 'C70Bkycatch_W'         )] = """C70 @ C60H28 Buckycatcher, ball wide in maw """
TAGL['%s-%s-dimer'      % (dbse, 'C70Bkycatch_W'         )] = """C70 @ C60H28 Buckycatcher, ball wide in maw """
TAGL['%s-%s-monoA-CP'   % (dbse, 'C70Bkycatch_W'         )] = """Buckyball from C70 @ C60H28 Buckycatcher, ball wide in maw """
TAGL['%s-%s-monoB-CP'   % (dbse, 'C70Bkycatch_W'         )] = """Buckycatcher from C70 @ C60H28 Buckycatcher, ball wide in maw """
TAGL['%s-%s-monoA-unCP' % (dbse, 'C70Bkycatch_W'         )] = """Buckyball from C70 @ C60H28 Buckycatcher, ball wide in maw """
TAGL['%s-%s-monoB-unCP' % (dbse, 'C70Bkycatch_W'         )] = """Buckycatcher from C70 @ C60H28 Buckycatcher, ball wide in maw """

TAGL['%s-%s'            % (dbse, 'C70Bkycatch_T'         )] = """C70 @ C60H28 Buckycatcher, ball tall along maw """
TAGL['%s-%s-dimer'      % (dbse, 'C70Bkycatch_T'         )] = """C70 @ C60H28 Buckycatcher, ball tall along maw """
TAGL['%s-%s-monoA-CP'   % (dbse, 'C70Bkycatch_T'         )] = """Buckyball from C70 @ C60H28 Buckycatcher, ball tall along maw """
TAGL['%s-%s-monoB-CP'   % (dbse, 'C70Bkycatch_T'         )] = """Buckycatcher from C70 @ C60H28 Buckycatcher, ball tall along maw """
TAGL['%s-%s-monoA-unCP' % (dbse, 'C70Bkycatch_T'         )] = """Buckyball from C70 @ C60H28 Buckycatcher, ball tall along maw """
TAGL['%s-%s-monoB-unCP' % (dbse, 'C70Bkycatch_T'         )] = """Buckycatcher from C70 @ C60H28 Buckycatcher, ball tall along maw """

# Shared monomers
TAGL['%s-Bz-mono-unCP'  % (dbse)] = 'Benzene'
TAGL['%s-2Bz-mono-unCP' % (dbse)] = 'Napthalene'
TAGL['%s-3Bz-mono-unCP' % (dbse)] = 'Anthracene'
TAGL['%s-4Bz-mono-unCP'  % (dbse)] = 'Tetracene'
TAGL['%s-Bkybowl_P-mono-unCP'  % (dbse)] = 'Corannulene, Pulay geometry'

# <<< Molecule Specifications >>>
monoA_unCP = 'monoA = dimer.extract_subsets(1)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_unCP = 'monoB = dimer.extract_subsets(2)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'
monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\nmonoA.set_name("monoA")\nPsiMod.set_active_molecule(monoA)\nPsiMod.IO.set_default_namespace("monoA")\n'
monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\nmonoB.set_name("monoB")\nPsiMod.set_active_molecule(monoB)\nPsiMod.IO.set_default_namespace("monoB")\n'

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
""", 0)

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
""", 0)

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
""", 0)

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
""", 0)

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
""", 0)

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
""", 0)

CFLOW_14Ae14Ae_3p8 = input.process_input("""
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
H        0.00000000     3.02849000     7.92041100
C        0.00000000     5.30353900     8.13791100
H        0.00000000     6.25270300     7.58991100
C        0.00000000     5.30353900     9.46988100
H        0.00000000     6.25270300    10.01788100
H        0.00000000     4.35437500    10.01788100
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
H        3.80000000     3.02849000     7.92041100
C        3.80000000     5.30353900     8.13791100
H        3.80000000     6.25270300     7.58991100
C        3.80000000     5.30353900     9.46988100
H        3.80000000     6.25270300    10.01788100
H        3.80000000     4.35437500    10.01788100
units angstrom
}
""", 0)

CFLOW_16Ae16Ae_3p8 = input.process_input("""
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
H        0.00000000    -1.70260600    -4.99491100
C        0.00000000     3.97765400     6.04044100
H        0.00000000     4.92681800     5.49244100
C        0.00000000     3.97765400     7.37241100
H        0.00000000     3.02849000     7.92041100
C        0.00000000     5.30353900     8.13791100
H        0.00000000     6.25270300     7.58991100
C        0.00000000     5.30353900     9.46988100
H        0.00000000     6.25270300    10.01788100
H        0.00000000     4.35437500    10.01788100
C        0.00000000    -3.97765400    -5.21241100
H        0.00000000    -4.92681800    -4.66441100
C        0.00000000    -3.97765400    -6.54438100
H        0.00000000    -4.92681800    -7.09238100
H        0.00000000    -3.02849000    -7.09238100
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
H        3.80000000    -1.70260600    -4.99491100
C        3.80000000     3.97765400     6.04044100
H        3.80000000     4.92681800     5.49244100
C        3.80000000     3.97765400     7.37241100
H        3.80000000     3.02849000     7.92041100
C        3.80000000     5.30353900     8.13791100
H        3.80000000     6.25270300     7.58991100
C        3.80000000     5.30353900     9.46988100
H        3.80000000     6.25270300    10.01788100
H        3.80000000     4.35437500    10.01788100
C        3.80000000    -3.97765400    -5.21241100
H        3.80000000    -4.92681800    -4.66441100
C        3.80000000    -3.97765400    -6.54438100
H        3.80000000    -4.92681800    -7.09238100
H        3.80000000    -3.02849000    -7.09238100
units angstrom
}
""", 0)

for item in cBzBz:
    distance = rxnpattern.match(item)
    ffdistance = '%14.8f' % (float(distance.group(2)))
    itemclean = dbse + '_' + re.sub('-', '_', re.sub(r'\.', 'p', item ))
    vars()[itemclean] = input.process_input("""
molecule dimer {
0 1
C        0.00000000    -0.39150000     0.00000000
C        1.20507400     0.30425000     0.00000000
C        1.20507400     1.69575000     0.00000000
C        0.00000000     2.39150000     0.00000000
C       -1.20507400     1.69575000     0.00000000
C       -1.20507400     0.30425000     0.00000000
H        0.00000000    -1.47150000     0.00000000
H        2.14038200    -0.23575000     0.00000000
H        2.14038200     2.23575000     0.00000000
H        0.00000000     3.47150000     0.00000000
H       -2.14038200     2.23575000     0.00000000
H       -2.14038200    -0.23575000     0.00000000
--
0 1
C        0.00000000    -0.39150000 %(ffdistance)s
C       -1.20507400     0.30425000 %(ffdistance)s
C       -1.20507400     1.69575000 %(ffdistance)s
C        0.00000000     2.39150000 %(ffdistance)s
C        1.20507400     1.69575000 %(ffdistance)s
C        1.20507400     0.30425000 %(ffdistance)s
H        0.00000000    -1.47150000 %(ffdistance)s
H       -2.14038200    -0.23575000 %(ffdistance)s
H       -2.14038200     2.23575000 %(ffdistance)s
H        0.00000000     3.47150000 %(ffdistance)s
H        2.14038200     2.23575000 %(ffdistance)s
H        2.14038200    -0.23575000 %(ffdistance)s
units angstrom
}
""" % vars() )

for item in c2BzBz:
    distance = rxnpattern.match(item)
    ffdistance = '%14.8f' % (float(distance.group(2)))
    itemclean = dbse + '_' + re.sub('-', '_', re.sub(r'\.', 'p', item ))
    vars()[itemclean] = input.process_input("""
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
C        0.00000000    -1.39150000 %(ffdistance)s
C        1.20507400    -0.69575000 %(ffdistance)s
C        1.20507400     0.69575000 %(ffdistance)s
C        0.00000000     1.39150000 %(ffdistance)s
C       -1.20507400     0.69575000 %(ffdistance)s
C       -1.20507400    -0.69575000 %(ffdistance)s
H        0.00000000    -2.47150000 %(ffdistance)s
H        2.14038200    -1.23575000 %(ffdistance)s
H        2.14038200     1.23575000 %(ffdistance)s
H        0.00000000     2.47150000 %(ffdistance)s
H       -2.14038200     1.23575000 %(ffdistance)s
H       -2.14038200    -1.23575000 %(ffdistance)s
units angstrom
}
""" % vars() )

for item in c2Bz2Bz:
    distance = rxnpattern.match(item)
    ffdistance = '%14.8f' % (float(distance.group(2)))
    itemclean = dbse + '_' + re.sub('-', '_', re.sub(r'\.', 'p', item ))
    vars()[itemclean] = input.process_input("""
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
C       -2.41591300    -0.70276500 %(ffdistance)s
C       -1.23100600    -1.39300200 %(ffdistance)s
C        0.00000000    -0.71089400 %(ffdistance)s
C        0.00000000     0.71089400 %(ffdistance)s
C       -1.23100600     1.39300200 %(ffdistance)s
C       -2.41591300     0.70276500 %(ffdistance)s
C        1.23100600    -1.39300200 %(ffdistance)s
C        1.23100600     1.39300200 %(ffdistance)s
C        2.41591300     0.70276500 %(ffdistance)s
C        2.41591300    -0.70276500 %(ffdistance)s
H        1.22889800    -2.47183700 %(ffdistance)s
H       -3.35016300    -1.23782500 %(ffdistance)s
H       -1.22889800    -2.47183700 %(ffdistance)s
H       -1.22889800     2.47183700 %(ffdistance)s
H       -3.35016300     1.23782500 %(ffdistance)s
H        1.22889800     2.47183700 %(ffdistance)s
H        3.35016300     1.23782500 %(ffdistance)s
H        3.35016300    -1.23782500 %(ffdistance)s
units angstrom
}
""" % vars() )

for item in c3Bz2Bz:
    distance = rxnpattern.match(item)
    ffdistance = '%14.8f' % (float(distance.group(2)))
    itemclean = dbse + '_' + re.sub('-', '_', re.sub(r'\.', 'p', item ))
    vars()[itemclean] = input.process_input("""
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
C       -2.41591300    -0.70276500 %(ffdistance)s
C       -1.23100600    -1.39300200 %(ffdistance)s
C        0.00000000    -0.71089400 %(ffdistance)s
C        0.00000000     0.71089400 %(ffdistance)s
C       -1.23100600     1.39300200 %(ffdistance)s
C       -2.41591300     0.70276500 %(ffdistance)s
C        1.23100600    -1.39300200 %(ffdistance)s
C        1.23100600     1.39300200 %(ffdistance)s
C        2.41591300     0.70276500 %(ffdistance)s
C        2.41591300    -0.70276500 %(ffdistance)s
H        1.22889800    -2.47183700 %(ffdistance)s
H       -3.35016300    -1.23782500 %(ffdistance)s
H       -1.22889800    -2.47183700 %(ffdistance)s
H       -1.22889800     2.47183700 %(ffdistance)s
H       -3.35016300     1.23782500 %(ffdistance)s
H        1.22889800     2.47183700 %(ffdistance)s
H        3.35016300     1.23782500 %(ffdistance)s
H        3.35016300    -1.23782500 %(ffdistance)s
units angstrom
}
""" % vars() )

for item in c3Bz3Bz:
    distance = rxnpattern.match(item)
    ffdistance = '%14.8f' % (float(distance.group(2)))
    itemclean = dbse + '_' + re.sub('-', '_', re.sub(r'\.', 'p', item ))
    vars()[itemclean] = input.process_input("""
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
C       -3.63206100    -0.70629600 %(ffdistance)s
C       -2.45292600    -1.39698600 %(ffdistance)s
C       -1.21398000    -0.71604800 %(ffdistance)s
C       -1.21398000     0.71604800 %(ffdistance)s
C       -2.45292600     1.39698600 %(ffdistance)s
C       -3.63206100     0.70629600 %(ffdistance)s
C        0.00000000    -1.39478500 %(ffdistance)s
C        0.00000000     1.39478500 %(ffdistance)s
C        1.21398000     0.71604800 %(ffdistance)s
C        1.21398000    -0.71604800 %(ffdistance)s
C        2.45292600    -1.39698600 %(ffdistance)s
H        2.45137100    -2.47594200 %(ffdistance)s
C        3.63206100    -0.70629600 %(ffdistance)s
C        3.63206100     0.70629600 %(ffdistance)s
C        2.45292600     1.39698600 %(ffdistance)s
H        0.00000000    -2.47562800 %(ffdistance)s
H       -4.56737900    -1.23960400 %(ffdistance)s
H       -2.45137100    -2.47594200 %(ffdistance)s
H       -2.45137100     2.47594200 %(ffdistance)s
H       -4.56737900     1.23960400 %(ffdistance)s
H        0.00000000     2.47562800 %(ffdistance)s
H        4.56737900    -1.23960400 %(ffdistance)s
H        4.56737900     1.23960400 %(ffdistance)s
H        2.45137100     2.47594200 %(ffdistance)s
units angstrom
}
""" % vars() )

for item in c4Bz3Bz:
    distance = rxnpattern.match(item)
    ffdistance = '%14.8f' % (float(distance.group(2)))
    itemclean = dbse + '_' + re.sub('-', '_', re.sub(r'\.', 'p', item ))
    vars()[itemclean] = input.process_input("""
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
C       -3.63206100    -0.70629600 %(ffdistance)s
C       -2.45292600    -1.39698600 %(ffdistance)s
C       -1.21398000    -0.71604800 %(ffdistance)s
C       -1.21398000     0.71604800 %(ffdistance)s
C       -2.45292600     1.39698600 %(ffdistance)s
C       -3.63206100     0.70629600 %(ffdistance)s
C        0.00000000    -1.39478500 %(ffdistance)s
C        0.00000000     1.39478500 %(ffdistance)s
C        1.21398000     0.71604800 %(ffdistance)s
C        1.21398000    -0.71604800 %(ffdistance)s
C        2.45292600    -1.39698600 %(ffdistance)s
H        2.45137100    -2.47594200 %(ffdistance)s
C        3.63206100    -0.70629600 %(ffdistance)s
C        3.63206100     0.70629600 %(ffdistance)s
C        2.45292600     1.39698600 %(ffdistance)s
H        0.00000000    -2.47562800 %(ffdistance)s
H       -4.56737900    -1.23960400 %(ffdistance)s
H       -2.45137100    -2.47594200 %(ffdistance)s
H       -2.45137100     2.47594200 %(ffdistance)s
H       -4.56737900     1.23960400 %(ffdistance)s
H        0.00000000     2.47562800 %(ffdistance)s
H        4.56737900    -1.23960400 %(ffdistance)s
H        4.56737900     1.23960400 %(ffdistance)s
H        2.45137100     2.47594200 %(ffdistance)s
units angstrom
}
""" % vars() )

CFLOW_BkybowlBkybowl_3p54 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97211978     0.61694803
C       -0.70622800     0.97211978     0.61694803
C       -1.14280400    -0.37137722     0.61681203
C        0.00000000    -1.20165922     0.61659503
C        1.14280400    -0.37137722     0.61681203
C        1.45779000     2.00650178     0.09413403
C       -1.45779000     2.00650178     0.09413403
C       -2.35873800    -0.76639722     0.09397203
C        0.00000000    -2.48004022     0.09366903
C        2.35873800    -0.76639722     0.09397203
C        0.69261800     3.17923978    -0.25321497
C       -0.69261800     3.17923978    -0.25321497
C       -2.80958100     1.64119778    -0.25292797
C       -3.23765700     0.32373778    -0.25303797
C       -2.42918200    -2.16498922    -0.25302597
C       -1.30841500    -2.97916822    -0.25327697
C        1.30841500    -2.97916822    -0.25327697
C        2.42918200    -2.16498922    -0.25302597
C        3.23765700     0.32373778    -0.25303797
C        2.80958100     1.64119778    -0.25292797
H        1.20851300     4.06642078    -0.61418797
H       -1.20851300     4.06642078    -0.61418797
H       -3.49401500     2.40602178    -0.61367197
H       -4.24094400     0.10729578    -0.61373997
H       -3.36816400    -2.57958822    -0.61350597
H       -1.41248600    -4.00024222    -0.61397997
H        1.41248600    -4.00024222    -0.61397997
H        3.36816400    -2.57958822    -0.61350597
H        4.24094400     0.10729578    -0.61373997
H        3.49401500     2.40602178    -0.61367197
--
0 1
C        0.70622800     0.97211978     4.15694803
C       -0.70622800     0.97211978     4.15694803
C       -1.14280400    -0.37137722     4.15681203
C        0.00000000    -1.20165922     4.15659503
C        1.14280400    -0.37137722     4.15681203
C        1.45779000     2.00650178     3.63413403
C       -1.45779000     2.00650178     3.63413403
C       -2.35873800    -0.76639722     3.63397203
C        0.00000000    -2.48004022     3.63366903
C        2.35873800    -0.76639722     3.63397203
C        0.69261800     3.17923978     3.28678503
C       -0.69261800     3.17923978     3.28678503
C       -2.80958100     1.64119778     3.28707203
C       -3.23765700     0.32373778     3.28696203
C       -2.42918200    -2.16498922     3.28697403
C       -1.30841500    -2.97916822     3.28672303
C        1.30841500    -2.97916822     3.28672303
C        2.42918200    -2.16498922     3.28697403
C        3.23765700     0.32373778     3.28696203
C        2.80958100     1.64119778     3.28707203
H        1.20851300     4.06642078     2.92581203
H       -1.20851300     4.06642078     2.92581203
H       -3.49401500     2.40602178     2.92632803
H       -4.24094400     0.10729578     2.92626003
H       -3.36816400    -2.57958822     2.92649403
H       -1.41248600    -4.00024222     2.92602003
H        1.41248600    -4.00024222     2.92602003
H        3.36816400    -2.57958822     2.92649403
H        4.24094400     0.10729578     2.92626003
H        3.49401500     2.40602178     2.92632803
units angstrom
}
""", 0)

CFLOW_BkybowlBkybowl_3p64 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97211978     0.61694803
C       -0.70622800     0.97211978     0.61694803
C       -1.14280400    -0.37137722     0.61681203
C        0.00000000    -1.20165922     0.61659503
C        1.14280400    -0.37137722     0.61681203
C        1.45779000     2.00650178     0.09413403
C       -1.45779000     2.00650178     0.09413403
C       -2.35873800    -0.76639722     0.09397203
C        0.00000000    -2.48004022     0.09366903
C        2.35873800    -0.76639722     0.09397203
C        0.69261800     3.17923978    -0.25321497
C       -0.69261800     3.17923978    -0.25321497
C       -2.80958100     1.64119778    -0.25292797
C       -3.23765700     0.32373778    -0.25303797
C       -2.42918200    -2.16498922    -0.25302597
C       -1.30841500    -2.97916822    -0.25327697
C        1.30841500    -2.97916822    -0.25327697
C        2.42918200    -2.16498922    -0.25302597
C        3.23765700     0.32373778    -0.25303797
C        2.80958100     1.64119778    -0.25292797
H        1.20851300     4.06642078    -0.61418797
H       -1.20851300     4.06642078    -0.61418797
H       -3.49401500     2.40602178    -0.61367197
H       -4.24094400     0.10729578    -0.61373997
H       -3.36816400    -2.57958822    -0.61350597
H       -1.41248600    -4.00024222    -0.61397997
H        1.41248600    -4.00024222    -0.61397997
H        3.36816400    -2.57958822    -0.61350597
H        4.24094400     0.10729578    -0.61373997
H        3.49401500     2.40602178    -0.61367197
--
0 1
C        0.70622800     0.97211978     4.25694803
C       -0.70622800     0.97211978     4.25694803
C       -1.14280400    -0.37137722     4.25681203
C        0.00000000    -1.20165922     4.25659503
C        1.14280400    -0.37137722     4.25681203
C        1.45779000     2.00650178     3.73413403
C       -1.45779000     2.00650178     3.73413403
C       -2.35873800    -0.76639722     3.73397203
C        0.00000000    -2.48004022     3.73366903
C        2.35873800    -0.76639722     3.73397203
C        0.69261800     3.17923978     3.38678503
C       -0.69261800     3.17923978     3.38678503
C       -2.80958100     1.64119778     3.38707203
C       -3.23765700     0.32373778     3.38696203
C       -2.42918200    -2.16498922     3.38697403
C       -1.30841500    -2.97916822     3.38672303
C        1.30841500    -2.97916822     3.38672303
C        2.42918200    -2.16498922     3.38697403
C        3.23765700     0.32373778     3.38696203
C        2.80958100     1.64119778     3.38707203
H        1.20851300     4.06642078     3.02581203
H       -1.20851300     4.06642078     3.02581203
H       -3.49401500     2.40602178     3.02632803
H       -4.24094400     0.10729578     3.02626003
H       -3.36816400    -2.57958822     3.02649403
H       -1.41248600    -4.00024222     3.02602003
H        1.41248600    -4.00024222     3.02602003
H        3.36816400    -2.57958822     3.02649403
H        4.24094400     0.10729578     3.02626003
H        3.49401500     2.40602178     3.02632803
units angstrom
}
""", 0)

CFLOW_BkybowlBkybowl_3p73 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97211978     0.61694803
C       -0.70622800     0.97211978     0.61694803
C       -1.14280400    -0.37137722     0.61681203
C        0.00000000    -1.20165922     0.61659503
C        1.14280400    -0.37137722     0.61681203
C        1.45779000     2.00650178     0.09413403
C       -1.45779000     2.00650178     0.09413403
C       -2.35873800    -0.76639722     0.09397203
C        0.00000000    -2.48004022     0.09366903
C        2.35873800    -0.76639722     0.09397203
C        0.69261800     3.17923978    -0.25321497
C       -0.69261800     3.17923978    -0.25321497
C       -2.80958100     1.64119778    -0.25292797
C       -3.23765700     0.32373778    -0.25303797
C       -2.42918200    -2.16498922    -0.25302597
C       -1.30841500    -2.97916822    -0.25327697
C        1.30841500    -2.97916822    -0.25327697
C        2.42918200    -2.16498922    -0.25302597
C        3.23765700     0.32373778    -0.25303797
C        2.80958100     1.64119778    -0.25292797
H        1.20851300     4.06642078    -0.61418797
H       -1.20851300     4.06642078    -0.61418797
H       -3.49401500     2.40602178    -0.61367197
H       -4.24094400     0.10729578    -0.61373997
H       -3.36816400    -2.57958822    -0.61350597
H       -1.41248600    -4.00024222    -0.61397997
H        1.41248600    -4.00024222    -0.61397997
H        3.36816400    -2.57958822    -0.61350597
H        4.24094400     0.10729578    -0.61373997
H        3.49401500     2.40602178    -0.61367197
--
0 1
C        0.70622800     0.97211978     4.34694803
C       -0.70622800     0.97211978     4.34694803
C       -1.14280400    -0.37137722     4.34681203
C        0.00000000    -1.20165922     4.34659503
C        1.14280400    -0.37137722     4.34681203
C        1.45779000     2.00650178     3.82413403
C       -1.45779000     2.00650178     3.82413403
C       -2.35873800    -0.76639722     3.82397203
C        0.00000000    -2.48004022     3.82366903
C        2.35873800    -0.76639722     3.82397203
C        0.69261800     3.17923978     3.47678503
C       -0.69261800     3.17923978     3.47678503
C       -2.80958100     1.64119778     3.47707203
C       -3.23765700     0.32373778     3.47696203
C       -2.42918200    -2.16498922     3.47697403
C       -1.30841500    -2.97916822     3.47672303
C        1.30841500    -2.97916822     3.47672303
C        2.42918200    -2.16498922     3.47697403
C        3.23765700     0.32373778     3.47696203
C        2.80958100     1.64119778     3.47707203
H        1.20851300     4.06642078     3.11581203
H       -1.20851300     4.06642078     3.11581203
H       -3.49401500     2.40602178     3.11632803
H       -4.24094400     0.10729578     3.11626003
H       -3.36816400    -2.57958822     3.11649403
H       -1.41248600    -4.00024222     3.11602003
H        1.41248600    -4.00024222     3.11602003
H        3.36816400    -2.57958822     3.11649403
H        4.24094400     0.10729578     3.11626003
H        3.49401500     2.40602178     3.11632803
units angstrom
}
""", 0)

CFLOW_BkybowlBkybowl_3p74 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97211978     0.61694803
C       -0.70622800     0.97211978     0.61694803
C       -1.14280400    -0.37137722     0.61681203
C        0.00000000    -1.20165922     0.61659503
C        1.14280400    -0.37137722     0.61681203
C        1.45779000     2.00650178     0.09413403
C       -1.45779000     2.00650178     0.09413403
C       -2.35873800    -0.76639722     0.09397203
C        0.00000000    -2.48004022     0.09366903
C        2.35873800    -0.76639722     0.09397203
C        0.69261800     3.17923978    -0.25321497
C       -0.69261800     3.17923978    -0.25321497
C       -2.80958100     1.64119778    -0.25292797
C       -3.23765700     0.32373778    -0.25303797
C       -2.42918200    -2.16498922    -0.25302597
C       -1.30841500    -2.97916822    -0.25327697
C        1.30841500    -2.97916822    -0.25327697
C        2.42918200    -2.16498922    -0.25302597
C        3.23765700     0.32373778    -0.25303797
C        2.80958100     1.64119778    -0.25292797
H        1.20851300     4.06642078    -0.61418797
H       -1.20851300     4.06642078    -0.61418797
H       -3.49401500     2.40602178    -0.61367197
H       -4.24094400     0.10729578    -0.61373997
H       -3.36816400    -2.57958822    -0.61350597
H       -1.41248600    -4.00024222    -0.61397997
H        1.41248600    -4.00024222    -0.61397997
H        3.36816400    -2.57958822    -0.61350597
H        4.24094400     0.10729578    -0.61373997
H        3.49401500     2.40602178    -0.61367197
--
0 1
C        0.70622800     0.97211978     4.35694803
C       -0.70622800     0.97211978     4.35694803
C       -1.14280400    -0.37137722     4.35681203
C        0.00000000    -1.20165922     4.35659503
C        1.14280400    -0.37137722     4.35681203
C        1.45779000     2.00650178     3.83413403
C       -1.45779000     2.00650178     3.83413403
C       -2.35873800    -0.76639722     3.83397203
C        0.00000000    -2.48004022     3.83366903
C        2.35873800    -0.76639722     3.83397203
C        0.69261800     3.17923978     3.48678503
C       -0.69261800     3.17923978     3.48678503
C       -2.80958100     1.64119778     3.48707203
C       -3.23765700     0.32373778     3.48696203
C       -2.42918200    -2.16498922     3.48697403
C       -1.30841500    -2.97916822     3.48672303
C        1.30841500    -2.97916822     3.48672303
C        2.42918200    -2.16498922     3.48697403
C        3.23765700     0.32373778     3.48696203
C        2.80958100     1.64119778     3.48707203
H        1.20851300     4.06642078     3.12581203
H       -1.20851300     4.06642078     3.12581203
H       -3.49401500     2.40602178     3.12632803
H       -4.24094400     0.10729578     3.12626003
H       -3.36816400    -2.57958822     3.12649403
H       -1.41248600    -4.00024222     3.12602003
H        1.41248600    -4.00024222     3.12602003
H        3.36816400    -2.57958822     3.12649403
H        4.24094400     0.10729578     3.12626003
H        3.49401500     2.40602178     3.12632803
units angstrom
}
""", 0)

CFLOW_BkybowlBkybowl_3p84 = input.process_input("""
molecule dimer {
0 1
C        0.70622800     0.97211978     0.61694803
C       -0.70622800     0.97211978     0.61694803
C       -1.14280400    -0.37137722     0.61681203
C        0.00000000    -1.20165922     0.61659503
C        1.14280400    -0.37137722     0.61681203
C        1.45779000     2.00650178     0.09413403
C       -1.45779000     2.00650178     0.09413403
C       -2.35873800    -0.76639722     0.09397203
C        0.00000000    -2.48004022     0.09366903
C        2.35873800    -0.76639722     0.09397203
C        0.69261800     3.17923978    -0.25321497
C       -0.69261800     3.17923978    -0.25321497
C       -2.80958100     1.64119778    -0.25292797
C       -3.23765700     0.32373778    -0.25303797
C       -2.42918200    -2.16498922    -0.25302597
C       -1.30841500    -2.97916822    -0.25327697
C        1.30841500    -2.97916822    -0.25327697
C        2.42918200    -2.16498922    -0.25302597
C        3.23765700     0.32373778    -0.25303797
C        2.80958100     1.64119778    -0.25292797
H        1.20851300     4.06642078    -0.61418797
H       -1.20851300     4.06642078    -0.61418797
H       -3.49401500     2.40602178    -0.61367197
H       -4.24094400     0.10729578    -0.61373997
H       -3.36816400    -2.57958822    -0.61350597
H       -1.41248600    -4.00024222    -0.61397997
H        1.41248600    -4.00024222    -0.61397997
H        3.36816400    -2.57958822    -0.61350597
H        4.24094400     0.10729578    -0.61373997
H        3.49401500     2.40602178    -0.61367197
--
0 1
C        0.70622800     0.97211978     4.45694803
C       -0.70622800     0.97211978     4.45694803
C       -1.14280400    -0.37137722     4.45681203
C        0.00000000    -1.20165922     4.45659503
C        1.14280400    -0.37137722     4.45681203
C        1.45779000     2.00650178     3.93413403
C       -1.45779000     2.00650178     3.93413403
C       -2.35873800    -0.76639722     3.93397203
C        0.00000000    -2.48004022     3.93366903
C        2.35873800    -0.76639722     3.93397203
C        0.69261800     3.17923978     3.58678503
C       -0.69261800     3.17923978     3.58678503
C       -2.80958100     1.64119778     3.58707203
C       -3.23765700     0.32373778     3.58696203
C       -2.42918200    -2.16498922     3.58697403
C       -1.30841500    -2.97916822     3.58672303
C        1.30841500    -2.97916822     3.58672303
C        2.42918200    -2.16498922     3.58697403
C        3.23765700     0.32373778     3.58696203
C        2.80958100     1.64119778     3.58707203
H        1.20851300     4.06642078     3.22581203
H       -1.20851300     4.06642078     3.22581203
H       -3.49401500     2.40602178     3.22632803
H       -4.24094400     0.10729578     3.22626003
H       -3.36816400    -2.57958822     3.22649403
H       -1.41248600    -4.00024222     3.22602003
H        1.41248600    -4.00024222     3.22602003
H        3.36816400    -2.57958822     3.22649403
H        4.24094400     0.10729578     3.22626003
H        3.49401500     2.40602178     3.22632803
units angstrom
}
""", 0)

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
""", 0)

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
""", 0)

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
""", 0)

CFLOW_C70Bkycatch = input.process_input("""
molecule dimer {
0 1
C        5.06797524     2.29936385    -2.08110191
C        5.89897484     1.17307146    -2.09318976
C        6.52371734     0.72403427    -0.86283944
C        6.52298902    -0.72680450    -0.86336685
C        4.82422217     3.00948035    -0.84308331
C        6.29135276     1.41952156     0.32994032
C        3.82179516     2.29178550    -2.81832071
C        5.51176916    -0.00025391    -2.85399954
C        5.89789027    -1.17446402    -2.09399943
C        4.30966839     0.00055758    -3.57145104
C        6.29011693    -1.42294291     0.32897281
C        5.42470772     2.57891710     0.33895272
C        6.04330495     0.69484063     1.55837520
C        5.06591606    -2.30004878    -2.08272685
C        6.04268024    -0.69882915     1.55788732
C        5.42243943    -2.58156480     0.33712048
C        3.44832217     3.45034855    -0.82960580
C        3.44988130     1.16508719    -3.55024381
C        2.82895255     3.00645603    -2.04921429
C        4.66725174     2.57712839     1.57044432
C        5.04939546     1.41310607     2.32396359
C        5.04810979    -1.41672651     2.32293081
C        4.66493768    -2.57988475     1.56857265
C        4.82158603    -3.01082984    -0.84519530
C        2.70895392     3.53703254     0.38071126
C        3.34127529     3.08846007     1.62819035
C        4.11489123     0.73371052     3.15446026
C        4.11426399    -0.73708517     3.15391512
C        3.33853737    -3.09005761     1.62599817
C        3.44526480    -3.45051077    -0.83202815
C        1.29743241     3.43959955     0.26310748
C        1.45628772     2.63991265    -2.08630639
C        2.51079522     2.57621272     2.66402205
C        2.89265937     1.41354883     3.41790352
C        2.89143447    -1.41605179     3.41687270
C        2.50851935    -2.57780057     2.66217162
C        2.70579418    -3.53723914     0.37821857
C        0.47621507     2.99837080     1.36291378
C        0.68026831     2.99610084    -0.95227396
C        1.07080918     2.57550069     2.54650265
C        1.68679031     0.69590182     3.76914407
C        1.68618267    -0.69763588     3.76869291
C        1.06855188    -2.57584146     2.54471320
C        1.29437200    -3.43841716     0.26065200
C        3.44875717    -1.16320657    -3.55093615
C        3.81970105    -2.29077144    -2.81984491
C        2.07403888     0.72077706    -3.53734415
C        2.82623686    -3.00516973    -2.05131774
C        1.06126875     1.45650863    -2.86129493
C       -0.66398842     2.29137013     0.83100892
C       -0.52458179     2.28389863    -0.60577871
C        0.54934101     1.42069400     3.24361876
C        0.54804539    -1.42107867     3.24264168
C        0.47356131    -2.99721728     1.36075992
C        0.67755844    -2.99355415    -0.95440564
C        2.07336548    -0.71768277    -3.53775270
C        1.45390244    -2.63728192    -2.08814605
C       -0.08290894     0.72105199    -2.44133621
C       -1.17309582     1.17226915     1.49464872
C       -0.90371109     1.16363667    -1.33779056
C       -0.55193239     0.72570549     2.72677805
C       -0.55261549    -0.72475925     2.72626074
C       -0.66605016    -2.28888993     0.82939921
C       -0.52671052    -2.28057480    -0.60744776
C        1.05992553    -1.45297180    -2.86232717
C       -0.08361305    -0.71677179    -2.44188719
C       -1.56544918     0.00170157     0.73502846
C       -1.43501618     0.00209330    -0.65904002
C       -1.17423771    -1.16979831     1.49380312
C       -0.90486484    -1.15945883    -1.33867218
--
0 1
C       -7.40282719     1.38705903     1.93015402
C       -6.43174906     0.70562386     1.17413928
C       -6.43169916    -0.70519539     1.17414003
C       -7.40275226    -1.38668799     1.93010957
C       -8.36853495    -0.69871005     2.66530912
C       -8.36854539     0.69898155     2.66534829
C       -5.41310213    -1.51929680     0.45760552
C       -4.63963081    -2.41637306     1.19827021
C       -3.76921098    -3.34625432     0.60970038
C       -3.65259791    -3.34732804    -0.83951400
C       -4.41678917    -2.41707422    -1.56061829
C       -5.29898154    -1.51959482    -0.95426824
C       -2.95693866    -4.25760533     1.42700325
C       -2.29489905    -5.23998494     0.71681149
C       -2.18285711    -5.24315227    -0.70653413
C       -2.72335425    -4.26251856    -1.51561186
C       -1.01136268    -5.95754614    -1.04754748
C       -0.39541059    -6.38900992     0.15230329
C       -1.18900026    -5.94920422     1.23961635
C        0.97082709    -6.61482652     0.25949369
C        1.68544147    -6.58301855    -0.99923794
C        1.08011363    -6.16437233    -2.18021740
C       -0.30058099    -5.73208286    -2.22175342
C       -0.96196950    -4.85841411    -3.16063317
C       -2.11638769    -4.15390053    -2.82007590
C        1.48320537    -6.57021747     1.61327170
C        0.70222080    -6.14405426     2.68364968
C       -0.66956541    -5.71632990     2.50843723
C       -1.47037063    -4.83949666     3.32828362
C       -2.56009673    -4.14107688     2.80945524
C       -6.18868895    -0.70522471    -1.82559182
C       -6.18871027     0.70572034    -1.82568751
C       -7.02482020     1.38698780    -2.72890173
C       -7.85837343     0.69895754    -3.61121594
C       -7.85834760    -0.69864423    -3.61118479
C       -7.02479039    -1.38654286    -2.72878356
C       -5.41326416     1.51980330     0.45755526
C       -4.63983861     2.41692463     1.19824812
C       -3.76918675     3.34660343     0.60971790
C       -3.65271120     3.34788783    -0.83950039
C       -4.41687746     2.41759136    -1.56060983
C       -5.29914943     1.52015400    -0.95430570
C       -2.72349367     4.26305124    -1.51555493
C       -2.18267268     5.24349468    -0.70643704
C       -2.29462293     5.24020773     0.71693380
C       -2.95682526     4.25784536     1.42703012
C       -1.18831948     5.94886912     1.23975344
C       -0.39462729     6.38839689     0.15241625
C       -1.01086572     5.95737367    -1.04748000
C       -2.11652987     4.15418557    -2.82002585
C       -0.96188010     4.85826802    -3.16061178
C       -0.30024170     5.73169669    -2.22174698
C        1.08064619     6.16338327    -2.18029293
C        1.68631879     6.58125160    -0.99924975
C        0.97177546     6.61324897     0.25952133
C       -2.56006900     4.14111534     2.80948951
C       -1.46994597     4.83887788     3.32828060
C       -0.66884844     5.71552491     2.50850096
C        0.70325213     6.14213626     2.68356052
C        1.48445470     6.56764025     1.61310560
H       -0.49452211    -4.65119484    -4.12204019
H       -2.49757409    -3.42346957    -3.53080365
H       -3.04852812    -3.40829515     3.44863134
H       -1.15740450    -4.62511667     4.34901029
H        1.17531918    -6.03336827     3.65833201
H        2.53722607    -6.77805717     1.79170491
H        2.75414559    -6.79110843    -1.01131275
H        1.69862775    -6.06129033    -3.07045038
H        1.69910712     6.06009438    -3.07053950
H        2.75515604     6.78864211    -1.01132691
H        2.53880160     6.77412731     1.79128346
H        1.17630534     6.03075868     3.65815555
H       -1.15699754     4.62430526     4.34897406
H       -3.04877949     3.40852960     3.44852522
H       -2.49798014     3.42386231    -3.53074179
H       -0.49467880     4.65101420    -4.12213549
H       -4.74557626    -2.41475495     2.28076588
H       -4.34825009    -2.41602943    -2.64608584
H       -4.74578607     2.41525372     2.28075286
H       -4.34814581     2.41640118    -2.64605832
H       -7.01680686    -2.47454834    -2.72375389
H       -8.50349987    -1.25092834    -4.29157316
H       -8.50354436     1.25115748    -4.29165732
H       -7.01688631     2.47502594    -2.72393185
H       -7.39376967    -2.47475553     1.92672947
H       -9.11566516    -1.25081347     3.23202192
H       -9.11570303     1.25103930     3.23207489
H       -7.39394496     2.47512523     1.92680562
units angstrom
}
""", 0)

CFLOW_C70Bkycatch_W = input.process_input("""
molecule dimer {
0 1
C       -0.10781880     0.97213108    -3.21727958
C       -1.21658377     0.56016424    -3.96683662
C       -1.60606545    -0.83787885    -3.96741791
C       -3.05611231    -0.89980546    -3.96799581
C        0.65010364     0.01133024    -2.44125693
C       -0.87042029    -1.76232215    -3.21687823
C       -0.17483466     2.18523992    -2.43423388
C       -2.42611638     1.36142648    -3.96715163
C       -3.56286502     0.45940653    -3.96775979
C       -2.47717727     2.54247415    -3.21712483
C       -3.71077366    -1.88615205    -3.21929949
C        0.27185059    -1.32865763    -2.43948761
C       -1.55443287    -2.76896301    -2.43596133
C       -4.70364950     0.77724639    -3.21974474
C       -2.94476283    -2.83530034    -2.44063531
C       -4.88527646    -1.55093347    -2.44217657
C        1.06632329     0.64189367    -1.20764345
C       -1.33314390     2.95519242    -2.43574521
C        0.55300326     1.98241819    -1.20441129
C        0.29411633    -2.07660796    -1.20435385
C       -0.83059707    -2.96525534    -1.20379370
C       -3.65303149    -3.09845992    -1.20815628
C       -4.85208388    -2.30494433    -1.20895454
C       -5.37162429    -0.24481259    -2.44208222
C        1.19298611    -0.09810670     0.00000000
C        0.78781226    -1.50971043     0.00000000
C       -1.49412360    -3.31477774     0.00000000
C       -2.96060313    -3.38359060     0.00000000
C       -5.38776660    -1.77934819     0.00000000
C       -5.83952228     0.34709243    -1.20894189
C        1.06632329     0.64189367     1.20764345
C        0.16294800     2.62299453     0.00000000
C        0.29411633    -2.07660796     1.20435385
C       -0.83059707    -2.96525534     1.20379370
C       -3.65303149    -3.09845992     1.20815628
C       -4.85208388    -2.30494433     1.20895454
C       -5.89999088    -0.40081285     0.00000000
C        0.65010364     0.01133024     2.44125693
C        0.55300326     1.98241819     1.20441129
C        0.27185059    -1.32865763     2.43948761
C       -1.55443287    -2.76896301     2.43596133
C       -2.94476283    -2.83530034     2.44063531
C       -4.88527646    -1.55093347     2.44217657
C       -5.83952228     0.34709243     1.20894189
C       -3.65660274     2.86379503    -2.44085073
C       -4.74846628     1.99799634    -2.44246228
C       -1.80286672     3.54408474    -1.20473355
C       -5.45436091     1.73254453    -1.20900272
C       -1.06021136     3.43587667     0.00000000
C       -0.10781880     0.97213108     3.21727958
C       -0.17483466     2.18523992     2.43423388
C       -0.87042029    -1.76232215     3.21687823
C       -3.71077366    -1.88615205     3.21929949
C       -5.37162429    -0.24481259     2.44208222
C       -5.45436091     1.73254453     1.20900272
C       -3.23828899     3.49132855    -1.20795017
C       -5.12070802     2.40486528     0.00000000
C       -1.80286672     3.54408474     1.20473355
C       -1.21658377     0.56016424     3.96683662
C       -1.33314390     2.95519242     2.43574521
C       -1.60606545    -0.83787885     3.96741791
C       -3.05611231    -0.89980546     3.96799581
C       -4.70364950     0.77724639     3.21974474
C       -4.74846628     1.99799634     2.44246228
C       -3.96933186     3.31952121     0.00000000
C       -3.23828899     3.49132855     1.20795017
C       -2.42611638     1.36142648     3.96715163
C       -2.47717727     2.54247415     3.21712483
C       -3.56286502     0.45940653     3.96775979
C       -3.65660274     2.86379503     2.44085073
--
0 1
C        6.90570687    -1.79603453     2.33940689
C        6.04083719    -1.06338601     1.50569884
C        6.12246354     0.34506609     1.50601979
C        7.06625335     0.97211755     2.34017719
C        7.92605173     0.23321429     3.15311790
C        7.84520480    -1.16204218     3.15271259
C        5.21867846     1.21718207     0.70841335
C        4.45562120     2.17276874     1.38379919
C        3.71001944     3.16414605     0.72780887
C        3.71001944     3.16414605    -0.72780887
C        4.45562120     2.17276874    -1.38379919
C        5.21867846     1.21718207    -0.70841335
C        2.92506029     4.15311218     1.47954890
C        2.39597164     5.17208365     0.71316580
C        2.39597164     5.17208365    -0.71316580
C        2.92506029     4.15311218    -1.47954890
C        1.31510289     5.97185277    -1.14615986
C        0.64187080     6.45566387     0.00000000
C        1.31510289     5.97185277     1.14615986
C       -0.70357855     6.80011635     0.00000000
C       -1.31136033     6.82425074    -1.31439107
C       -0.64503868     6.35873108    -2.44656282
C        0.69450315     5.81377382    -2.37982472
C        1.37566404     4.90527922    -3.26985567
C        2.43718833     4.10925224    -2.83661911
C       -1.31136033     6.82425074     1.31439107
C       -0.64503868     6.35873108     2.44656282
C        0.69450315     5.81377382     2.37982472
C        1.37566404     4.90527922     3.26985567
C        2.43718833     4.10925224     2.83661911
C        6.12246354     0.34506609    -1.50601979
C        6.04083719    -1.06338601    -1.50569884
C        6.90570687    -1.79603453    -2.33940689
C        7.84520480    -1.16204218    -3.15271259
C        7.92605173     0.23321429    -3.15311790
C        7.06625335     0.97211755    -2.34017719
C        5.04147517    -1.82425591     0.70843700
C        4.17151913    -2.68356989     1.38393595
C        3.31296095    -3.57829123     0.72728601
C        3.31296095    -3.57829123    -0.72728601
C        4.17151913    -2.68356989    -1.38393595
C        5.04147517    -1.82425591    -0.70843700
C        2.40407610    -4.45524164    -1.47787004
C        1.75156130    -5.40255934    -0.71382632
C        1.75156130    -5.40255934     0.71382632
C        2.40407610    -4.45524164     1.47787004
C        0.57815418    -6.05993979     1.14701113
C       -0.15064668    -6.45653706     0.00000000
C        0.57815418    -6.05993979    -1.14701113
C        1.91731357    -4.33365482    -2.83061445
C        0.76170346    -4.98432977    -3.26336797
C       -0.02237763    -5.80907271    -2.37620974
C       -1.42354491    -6.16582790    -2.44250520
C       -2.14239143    -6.54130620    -1.31101719
C       -1.53106269    -6.60924580     0.00000000
C        1.91731357    -4.33365482     2.83061445
C        0.76170346    -4.98432977     3.26336797
C       -0.02237763    -5.80907271     2.37620974
C       -1.42354491    -6.16582790     2.44250520
C       -2.14239143    -6.54130620     1.31101719
H        0.98822457     4.74441198    -4.27449004
H        2.82463078     3.35979862    -3.52388738
H        2.82463078     3.35979862     3.52388738
H        0.98822457     4.74441198     4.27449004
H       -1.19454305     6.31196301     3.38564845
H       -2.35403703     7.12194754     1.41556147
H       -2.35403703     7.12194754    -1.41556147
H       -1.19454305     6.31196301    -3.38564845
H       -1.96346759    -6.03668769    -3.37907240
H       -3.21666968    -6.69178093    -1.40715986
H       -3.21666968    -6.69178093     1.40715986
H       -1.96346759    -6.03668769     3.37907240
H        0.39006180    -4.76624715     4.26315989
H        2.39353638    -3.63284792     3.51302759
H        2.39353638    -3.63284792    -3.51302759
H        0.39006180    -4.76624715    -4.26315989
H        4.47498141     2.17041498     2.47122669
H        4.47498141     2.17041498    -2.47122669
H        4.19101455    -2.68369728     2.47137160
H        4.19101455    -2.68369728    -2.47137160
H        7.12160205     2.05874398    -2.33571188
H        8.65439980     0.74410823    -3.77974819
H        8.50950791    -1.75406635    -3.77909453
H        6.83476091    -2.88178047    -2.33438053
H        7.12160205     2.05874398     2.33571188
H        8.65439980     0.74410823     3.77974819
H        8.50950791    -1.75406635     3.77909453
H        6.83476091    -2.88178047     2.33438053
units angstrom
}
""", 0)

CFLOW_C70Bkycatch_T = input.process_input("""
molecule dimer {
0 1
C        0.57719530     1.41260349     3.20966236
C        1.52957617     0.71536970     3.95340794
C        1.52341728    -0.72993399     3.95278030
C        2.89812010    -1.18585626     3.96075656
C       -0.40966062     0.69680662     2.43743626
C        0.56534991    -1.41818467     3.20852920
C        0.96313082     2.57233822     2.43862832
C        2.90834840     1.15899686     3.96106572
C        3.75619807    -0.01716596     3.96474958
C        3.27934042     2.28567333     3.21938256
C        3.25916887    -2.31580260     3.21924696
C       -0.41541659    -0.69367716     2.43698313
C        0.94062341    -2.58149910     2.43801192
C        4.94177554    -0.02241586     3.21939933
C        2.26197530    -3.02218407     2.44150203
C        4.47981288    -2.31203140     2.44172906
C       -0.64549333     1.41693343     1.20789974
C        2.28846728     3.00091312     2.44177690
C        0.20426813     2.57766514     1.20847235
C       -0.65792380    -1.41181580     1.20773964
C        0.18139487    -2.58022583     1.20826884
C        2.87125346    -3.46599697     1.20827994
C        4.24124564    -3.02792521     1.20854561
C        5.30512139    -1.18828847     2.44235511
C       -0.97096884     0.73981891     0.00000000
C       -0.97752486    -0.73186485     0.00000000
C        0.72707542    -3.09635659     0.00000000
C        2.12500277    -3.55575480     0.00000000
C        4.89941243    -2.66744630     0.00000000
C        5.91523627    -0.74577536     1.20873695
C       -0.64549333     1.41693343    -1.20789974
C        0.75449827     3.08875686     0.00000000
C       -0.65792380    -1.41181580    -1.20773964
C        0.18139487    -2.58022583    -1.20826884
C        2.87125346    -3.46599697    -1.20827994
C        4.24124564    -3.02792521    -1.20854561
C        5.76932895    -1.48168053     0.00000000
C       -0.40966062     0.69680662    -2.43743626
C        0.20426813     2.57766514    -1.20847235
C       -0.41541659    -0.69367716    -2.43698313
C        0.94062341    -2.58149910    -2.43801192
C        2.26197530    -3.02218407    -2.44150203
C        4.47981288    -2.31203140    -2.44172906
C        5.91523627    -0.74577536    -1.20873695
C        4.49983898     2.27104449     2.44175703
C        5.31528705     1.14018198     2.44231634
C        2.90165291     3.43917289     1.20837050
C        5.92147887     0.69230566     1.20870813
C        2.15639371     3.53564498     0.00000000
C        0.57719530     1.41260349    -3.20966236
C        0.96313082     2.57233822    -2.43862832
C        0.56534991    -1.41818467    -3.20852920
C        3.25916887    -2.31580260    -3.21924696
C        5.30512139    -1.18828847    -2.44235511
C        5.92147887     0.69230566    -1.20870813
C        4.26761240     2.98899299     1.20857391
C        5.78202789     1.42944850     0.00000000
C        2.90165291     3.43917289    -1.20837050
C        1.52957617     0.71536970    -3.95340794
C        2.28846728     3.00091312    -2.44177690
C        1.52341728    -0.72993399    -3.95278030
C        2.89812010    -1.18585626    -3.96075656
C        4.94177554    -0.02241586    -3.21939933
C        5.31528705     1.14018198    -2.44231634
C        4.92256137     2.62278876     0.00000000
C        4.26761240     2.98899299    -1.20857391
C        2.90834840     1.15899686    -3.96106572
C        3.27934042     2.28567333    -3.21938256
C        3.75619807    -0.01716596    -3.96474958
C        4.49983898     2.27104449    -2.44175703
--
0 1
C       -7.01793526    -2.37728744     1.38614650
C       -6.12012591    -1.53409251     0.70549270
C       -6.12012591    -1.53409251    -0.70549270
C       -7.01793526    -2.37728744    -1.38614650
C       -7.91089615    -3.19956318    -0.69879261
C       -7.91089615    -3.19956318     0.69879261
C       -5.17785030    -0.72507011    -1.52494216
C       -4.36310761    -1.38976568    -2.44588892
C       -3.58623860    -0.72151453    -3.40367410
C       -3.59438859     0.73310546    -3.39294709
C       -4.38580987     1.37882607    -2.43207640
C       -5.19248222     0.69146053    -1.52133143
C       -2.77995584    -1.45863485    -4.38745582
C       -2.27829826    -0.68396421    -5.41450897
C       -2.28301127     0.74246168    -5.40109696
C       -2.79157438     1.49415240    -4.36080444
C       -1.22676695     1.18724256    -6.22627270
C       -0.56838807     0.04798403    -6.74806679
C       -1.21969845    -1.10560660    -6.24901765
C        0.76549781     0.05610751    -7.13269927
C        1.37244130     1.37142255    -7.15434080
C        0.72312375     2.49439362    -6.64612635
C       -0.59779416     2.41625096    -6.05930820
C       -1.23911789     3.28483581    -5.10153284
C       -2.28619584     2.84276337    -4.29021637
C        1.38053185    -1.25455913    -7.18122242
C        0.73814235    -2.39159678    -6.69582606
C       -0.58313594    -2.33375989    -6.10739626
C       -1.21814452    -3.22513820    -5.16663408
C       -2.26621218    -2.80533066    -4.34473493
C       -6.15717387     1.47619172    -0.70500458
C       -6.15717387     1.47619172     0.70500458
C       -7.07613186     2.29522026     1.38658909
C       -7.99050720     3.09342876     0.69883936
C       -7.99050720     3.09342876    -0.69883936
C       -7.07613186     2.29522026    -1.38658909
C       -5.17785030    -0.72507011     1.52494216
C       -4.36310761    -1.38976568     2.44588892
C       -3.58623860    -0.72151453     3.40367410
C       -3.59438859     0.73310546     3.39294709
C       -4.38580987     1.37882607     2.43207640
C       -5.19248222     0.69146053     1.52133143
C       -2.79157438     1.49415240     4.36080444
C       -2.28301127     0.74246168     5.40109696
C       -2.27829826    -0.68396421     5.41450897
C       -2.77995584    -1.45863485     4.38745582
C       -1.21969845    -1.10560660     6.24901765
C       -0.56838807     0.04798403     6.74806679
C       -1.22676695     1.18724256     6.22627270
C       -2.28619584     2.84276337     4.29021637
C       -1.23911789     3.28483581     5.10153284
C       -0.59779416     2.41625096     6.05930820
C        0.72312375     2.49439362     6.64612635
C        1.37244130     1.37142255     7.15434080
C        0.76549781     0.05610751     7.13269927
C       -2.26621218    -2.80533066     4.34473493
C       -1.21814452    -3.22513820     5.16663408
C       -0.58313594    -2.33375989     6.10739626
C        0.73814235    -2.39159678     6.69582606
C        1.38053185    -1.25455913     7.18122242
H       -0.83403138     4.27994472    -4.92401873
H       -2.64820499     3.51289365    -3.51289465
H       -2.62254343    -3.49306260    -3.58034486
H       -0.80715384    -4.22124181    -5.00955636
H        1.29624156    -3.32649408    -6.66589095
H        2.41398078    -1.34747453    -7.51215834
H        2.40524642     1.47731570    -7.48334724
H        1.27533175     3.43197205    -6.59722947
H        1.27533175     3.43197205     6.59722947
H        2.40524642     1.47731570     7.48334724
H        2.41398078    -1.34747453     7.51215834
H        1.29624156    -3.32649408     6.66589095
H       -0.80715384    -4.22124181     5.00955636
H       -2.62254343    -3.49306260     3.58034486
H       -2.64820499     3.51289365     3.51289465
H       -0.83403138     4.27994472     4.92401873
H       -4.37587012    -2.47734670    -2.45283752
H       -4.41722929     2.46606700    -2.42897203
H       -4.37587012    -2.47734670     2.45283752
H       -4.41722929     2.46606700     2.42897203
H       -7.06812939     2.29040531    -2.47459918
H       -8.69824914     3.70841580    -1.25116110
H       -8.69824914     3.70841580     1.25116110
H       -7.06812939     2.29040531     2.47459918
H       -7.01096558    -2.37180000    -2.47413143
H       -8.60234201    -3.83259883    -1.25136916
H       -8.60234201    -3.83259883     1.25136916
H       -7.01096558    -2.37180000     2.47413143
units angstrom
}
""", 0)

CFLOW_Benzene_monomer = input.process_input("""
molecule monomer {
0 1
C        0.00000000    -1.39150000     0.00000000
C        1.20507400    -0.69575000     0.00000000
C        1.20507400     0.69575000     0.00000000
C        0.00000000     1.39150000     0.00000000
C       -1.20507400     0.69575000     0.00000000
C       -1.20507400    -0.69575000     0.00000000
H        0.00000000    -2.47150000     0.00000000
H        2.14038200    -1.23575000     0.00000000
H        2.14038200     1.23575000     0.00000000
H        0.00000000     2.47150000     0.00000000
H       -2.14038200     1.23575000     0.00000000
H       -2.14038200    -1.23575000     0.00000000
units angstrom
}
""", 0)

CFLOW_Napthalene_monomer = input.process_input("""
molecule monomer {
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
units angstrom
}
""", 0)

CFLOW_Anthracene_monomer = input.process_input("""
molecule monomer {
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
units angstrom
}
""", 0)

# 4Bz geometry NOT finalized
CFLOW_Tetracene_monomer = input.process_input("""
molecule monomer {
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
units angstrom
}
""", 0)

CFLOW_Bkybowl_Pulay_monomer = input.process_input("""
molecule monomer {
0 1
C        0.70622800     0.97211978     0.61694803
C       -0.70622800     0.97211978     0.61694803
C       -1.14280400    -0.37137722     0.61681203
C        0.00000000    -1.20165922     0.61659503
C        1.14280400    -0.37137722     0.61681203
C        1.45779000     2.00650178     0.09413403
C       -1.45779000     2.00650178     0.09413403
C       -2.35873800    -0.76639722     0.09397203
C        0.00000000    -2.48004022     0.09366903
C        2.35873800    -0.76639722     0.09397203
C        0.69261800     3.17923978    -0.25321497
C       -0.69261800     3.17923978    -0.25321497
C       -2.80958100     1.64119778    -0.25292797
C       -3.23765700     0.32373778    -0.25303797
C       -2.42918200    -2.16498922    -0.25302597
C       -1.30841500    -2.97916822    -0.25327697
C        1.30841500    -2.97916822    -0.25327697
C        2.42918200    -2.16498922    -0.25302597
C        3.23765700     0.32373778    -0.25303797
C        2.80958100     1.64119778    -0.25292797
H        1.20851300     4.06642078    -0.61418797
H       -1.20851300     4.06642078    -0.61418797
H       -3.49401500     2.40602178    -0.61367197
H       -4.24094400     0.10729578    -0.61373997
H       -3.36816400    -2.57958822    -0.61350597
H       -1.41248600    -4.00024222    -0.61397997
H        1.41248600    -4.00024222    -0.61397997
H        3.36816400    -2.57958822    -0.61350597
H        4.24094400     0.10729578    -0.61373997
H        3.49401500     2.40602178    -0.61367197
units angstrom
}
""", 0)

# <<< Geometry Specification Strings >>>
GEOS = {}
for rxn in HRXN:

    if (rxn is 'C60Bkybowl') or (rxn is 'C60Bkycatch') or (rxn in Grimme70):
        GEOS['%s-%s-dimer'    % (dbse, rxn)] = eval('%s_%s' % (dbse, rxn ))
        GEOS['%s-%s-monoA-CP' % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn ))) + monoA_CP
        GEOS['%s-%s-monoB-CP' % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn ))) + monoB_CP
        GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn ))) + monoA_unCP
        GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn ))) + monoB_unCP

    else:
        distance = rxnpattern.match(rxn)
        GEOS['%s-%s-dimer'    % (dbse, rxn)] = eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))
        GEOS['%s-%s-monoA-CP' % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))) + monoA_CP
        GEOS['%s-%s-monoB-CP' % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))) + monoB_CP
        GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))) + monoA_unCP
        GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = str(eval('%s_%s_%s' % (dbse, distance.group(1), re.sub(r'\.', 'p', distance.group(2) )))) + monoB_unCP


GEOS['%s-Bz-mono-unCP'         % (dbse)] = eval('%s_Benzene_monomer'        % (dbse))
GEOS['%s-2Bz-mono-unCP'        % (dbse)] = eval('%s_Napthalene_monomer'     % (dbse))
GEOS['%s-3Bz-mono-unCP'        % (dbse)] = eval('%s_Anthracene_monomer'     % (dbse))
GEOS['%s-4Bz-mono-unCP'        % (dbse)] = eval('%s_Tetracene_monomer'      % (dbse))
GEOS['%s-Bkybowl_P-mono-unCP'  % (dbse)] = eval('%s_Bkybowl_Pulay_monomer'  % (dbse))
