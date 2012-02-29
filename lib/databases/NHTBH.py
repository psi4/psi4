"""
**NHTBH**

| Database (Truhlar) of non-hydrogen-transfer barrier height reactions.
| Geometries and Reaction energies from Truhlar and coworkers at site http://t1.chem.umn.edu/misc/database_group/database_therm_bh/non_H.htm.

- **cp**  ``'off'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``

----

"""
import re
import input

# <<< NHTBH Database Module >>>
dbse = 'NHTBH'
isOS = 'true'

# <<< Database Members >>>
HRXN = range(1, 39)
HRXN_SM = [3, 4, 31, 32]
HRXN_LG = [36]

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s' % (dbse,  1)] = ['%s-%s-reagent' % (dbse, 'H'     ),
                              '%s-%s-reagent' % (dbse, 'N2O'   ),
                              '%s-%s-reagent' % (dbse, 'N2OHts') ]
RXNM['%s-%s' % (dbse,  1)] = dict(zip(ACTV['%s-%s' % (dbse,  1)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse,  2)] = ['%s-%s-reagent' % (dbse, 'OH'    ),
                              '%s-%s-reagent' % (dbse, 'N2'    ),
                              '%s-%s-reagent' % (dbse, 'N2OHts') ]
RXNM['%s-%s' % (dbse,  2)] = dict(zip(ACTV['%s-%s' % (dbse,  2)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse,  3)] = ['%s-%s-reagent' % (dbse, 'H'    ),
                              '%s-%s-reagent' % (dbse, 'HF'   ),
                              '%s-%s-reagent' % (dbse, 'HFHts') ]
RXNM['%s-%s' % (dbse,  3)] = dict(zip(ACTV['%s-%s' % (dbse,  3)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse,  4)] = ['%s-%s-reagent' % (dbse, 'H'    ),
                              '%s-%s-reagent' % (dbse, 'HF'   ),
                              '%s-%s-reagent' % (dbse, 'HFHts') ]
RXNM['%s-%s' % (dbse,  4)] = dict(zip(ACTV['%s-%s' % (dbse,  4)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse,  5)] = ['%s-%s-reagent' % (dbse, 'H'     ),
                              '%s-%s-reagent' % (dbse, 'HCl'   ),
                              '%s-%s-reagent' % (dbse, 'HClHts') ]
RXNM['%s-%s' % (dbse,  5)] = dict(zip(ACTV['%s-%s' % (dbse,  5)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse,  6)] = ['%s-%s-reagent' % (dbse, 'H'     ),
                              '%s-%s-reagent' % (dbse, 'HCl'   ),
                              '%s-%s-reagent' % (dbse, 'HClHts') ]
RXNM['%s-%s' % (dbse,  6)] = dict(zip(ACTV['%s-%s' % (dbse,  6)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse,  7)] = ['%s-%s-reagent' % (dbse, 'H'      ),
                              '%s-%s-reagent' % (dbse, 'CH3F'   ),
                              '%s-%s-reagent' % (dbse, 'HFCH3ts') ]
RXNM['%s-%s' % (dbse,  7)] = dict(zip(ACTV['%s-%s' % (dbse,  7)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse,  8)] = ['%s-%s-reagent' % (dbse, 'HF'     ),
                              '%s-%s-reagent' % (dbse, 'CH3'    ),
                              '%s-%s-reagent' % (dbse, 'HFCH3ts') ]
RXNM['%s-%s' % (dbse,  8)] = dict(zip(ACTV['%s-%s' % (dbse,  8)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse,  9)] = ['%s-%s-reagent' % (dbse, 'H'    ),
                              '%s-%s-reagent' % (dbse, 'F2'   ),
                              '%s-%s-reagent' % (dbse, 'HF2ts') ]
RXNM['%s-%s' % (dbse,  9)] = dict(zip(ACTV['%s-%s' % (dbse,  9)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 10)] = ['%s-%s-reagent' % (dbse, 'HF'   ),
                              '%s-%s-reagent' % (dbse, 'F'    ),
                              '%s-%s-reagent' % (dbse, 'HF2ts') ]
RXNM['%s-%s' % (dbse, 10)] = dict(zip(ACTV['%s-%s' % (dbse, 10)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 11)] = ['%s-%s-reagent' % (dbse, 'CH3'     ),
                              '%s-%s-reagent' % (dbse, 'ClF'     ),
                              '%s-%s-reagent' % (dbse, 'CH3FClts') ]
RXNM['%s-%s' % (dbse, 11)] = dict(zip(ACTV['%s-%s' % (dbse, 11)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 12)] = ['%s-%s-reagent' % (dbse, 'CH3F'    ),
                              '%s-%s-reagent' % (dbse, 'Cl'      ),
                              '%s-%s-reagent' % (dbse, 'CH3FClts') ]
RXNM['%s-%s' % (dbse, 12)] = dict(zip(ACTV['%s-%s' % (dbse, 12)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 13)] = ['%s-%s-reagent' % (dbse, 'F_anion'),
                              '%s-%s-reagent' % (dbse, 'CH3F'   ),
                              '%s-%s-reagent' % (dbse, 'FCH3Fts') ]
RXNM['%s-%s' % (dbse, 13)] = dict(zip(ACTV['%s-%s' % (dbse, 13)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 14)] = ['%s-%s-reagent' % (dbse, 'F_anion'),
                              '%s-%s-reagent' % (dbse, 'CH3F'   ),
                              '%s-%s-reagent' % (dbse, 'FCH3Fts') ]
RXNM['%s-%s' % (dbse, 14)] = dict(zip(ACTV['%s-%s' % (dbse, 14)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 15)] = ['%s-%s-reagent' % (dbse, 'FCH3Fcomp'),
                              '%s-%s-reagent' % (dbse, 'FCH3Fts'  ) ]
RXNM['%s-%s' % (dbse, 15)] = dict(zip(ACTV['%s-%s' % (dbse, 15)], [-1, +1]))

ACTV['%s-%s' % (dbse, 16)] = ['%s-%s-reagent' % (dbse, 'FCH3Fcomp'),
                              '%s-%s-reagent' % (dbse, 'FCH3Fts'  ) ]
RXNM['%s-%s' % (dbse, 16)] = dict(zip(ACTV['%s-%s' % (dbse, 16)], [-1, +1]))

ACTV['%s-%s' % (dbse, 17)] = ['%s-%s-reagent' % (dbse, 'Cl_anion' ),
                              '%s-%s-reagent' % (dbse, 'CH3Cl'    ),
                              '%s-%s-reagent' % (dbse, 'ClCH3Clts') ]
RXNM['%s-%s' % (dbse, 17)] = dict(zip(ACTV['%s-%s' % (dbse, 17)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 18)] = ['%s-%s-reagent' % (dbse, 'Cl_anion' ),
                              '%s-%s-reagent' % (dbse, 'CH3Cl'    ),
                              '%s-%s-reagent' % (dbse, 'ClCH3Clts') ]
RXNM['%s-%s' % (dbse, 18)] = dict(zip(ACTV['%s-%s' % (dbse, 18)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 19)] = ['%s-%s-reagent' % (dbse, 'ClCH3Clcomp'),
                              '%s-%s-reagent' % (dbse, 'ClCH3Clts'  ) ]
RXNM['%s-%s' % (dbse, 19)] = dict(zip(ACTV['%s-%s' % (dbse, 19)], [-1, +1]))

ACTV['%s-%s' % (dbse, 20)] = ['%s-%s-reagent' % (dbse, 'ClCH3Clcomp'),
                              '%s-%s-reagent' % (dbse, 'ClCH3Clts'  ) ]
RXNM['%s-%s' % (dbse, 20)] = dict(zip(ACTV['%s-%s' % (dbse, 20)], [-1, +1]))

ACTV['%s-%s' % (dbse, 21)] = ['%s-%s-reagent' % (dbse, 'F_anion' ),
                              '%s-%s-reagent' % (dbse, 'CH3Cl'   ),
                              '%s-%s-reagent' % (dbse, 'FCH3Clts') ]
RXNM['%s-%s' % (dbse, 21)] = dict(zip(ACTV['%s-%s' % (dbse, 21)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 22)] = ['%s-%s-reagent' % (dbse, 'CH3F'),
                              '%s-%s-reagent' % (dbse, 'Cl_anion'),
                              '%s-%s-reagent' % (dbse, 'FCH3Clts') ]
RXNM['%s-%s' % (dbse, 22)] = dict(zip(ACTV['%s-%s' % (dbse, 22)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 23)] = ['%s-%s-reagent' % (dbse, 'FCH3Clcomp1'),
                              '%s-%s-reagent' % (dbse, 'FCH3Clts'   ) ]
RXNM['%s-%s' % (dbse, 23)] = dict(zip(ACTV['%s-%s' % (dbse, 23)], [-1, +1]))

ACTV['%s-%s' % (dbse, 24)] = ['%s-%s-reagent' % (dbse, 'FCH3Clcomp2'),
                              '%s-%s-reagent' % (dbse, 'FCH3Clts'   ) ]
RXNM['%s-%s' % (dbse, 24)] = dict(zip(ACTV['%s-%s' % (dbse, 24)], [-1, +1]))

ACTV['%s-%s' % (dbse, 25)] = ['%s-%s-reagent' % (dbse, 'OH_anion'),
                              '%s-%s-reagent' % (dbse, 'CH3F'    ),
                              '%s-%s-reagent' % (dbse, 'HOCH3Fts') ]
RXNM['%s-%s' % (dbse, 25)] = dict(zip(ACTV['%s-%s' % (dbse, 25)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 26)] = ['%s-%s-reagent' % (dbse, 'CH3OH'   ),
                              '%s-%s-reagent' % (dbse, 'F_anion' ),
                              '%s-%s-reagent' % (dbse, 'HOCH3Fts') ]
RXNM['%s-%s' % (dbse, 26)] = dict(zip(ACTV['%s-%s' % (dbse, 26)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 27)] = ['%s-%s-reagent' % (dbse, 'HOCH3Fcomp2'),
                              '%s-%s-reagent' % (dbse, 'HOCH3Fts'   ) ]
RXNM['%s-%s' % (dbse, 27)] = dict(zip(ACTV['%s-%s' % (dbse, 27)], [-1, +1]))

ACTV['%s-%s' % (dbse, 28)] = ['%s-%s-reagent' % (dbse, 'HOCH3Fcomp1'),
                              '%s-%s-reagent' % (dbse, 'HOCH3Fts'   ) ]
RXNM['%s-%s' % (dbse, 28)] = dict(zip(ACTV['%s-%s' % (dbse, 28)], [-1, +1]))

ACTV['%s-%s' % (dbse, 29)] = ['%s-%s-reagent' % (dbse, 'H'    ),
                              '%s-%s-reagent' % (dbse, 'N2'   ),
                              '%s-%s-reagent' % (dbse, 'HN2ts') ]
RXNM['%s-%s' % (dbse, 29)] = dict(zip(ACTV['%s-%s' % (dbse, 29)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 30)] = ['%s-%s-reagent' % (dbse, 'HN2'  ),
                              '%s-%s-reagent' % (dbse, 'HN2ts') ]
RXNM['%s-%s' % (dbse, 30)] = dict(zip(ACTV['%s-%s' % (dbse, 30)], [-1, +1]))

ACTV['%s-%s' % (dbse, 31)] = ['%s-%s-reagent' % (dbse, 'H'    ),
                              '%s-%s-reagent' % (dbse, 'CO'   ),
                              '%s-%s-reagent' % (dbse, 'HCOts') ]
RXNM['%s-%s' % (dbse, 31)] = dict(zip(ACTV['%s-%s' % (dbse, 31)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 32)] = ['%s-%s-reagent' % (dbse, 'HCO'  ),
                              '%s-%s-reagent' % (dbse, 'HCOts') ]
RXNM['%s-%s' % (dbse, 32)] = dict(zip(ACTV['%s-%s' % (dbse, 32)], [-1, +1]))

ACTV['%s-%s' % (dbse, 33)] = ['%s-%s-reagent' % (dbse, 'H'     ),
                              '%s-%s-reagent' % (dbse, 'C2H4'  ),
                              '%s-%s-reagent' % (dbse, 'C2H5ts') ]
RXNM['%s-%s' % (dbse, 33)] = dict(zip(ACTV['%s-%s' % (dbse, 33)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 34)] = ['%s-%s-reagent' % (dbse, 'C2H5'  ),
                              '%s-%s-reagent' % (dbse, 'C2H5ts') ]
RXNM['%s-%s' % (dbse, 34)] = dict(zip(ACTV['%s-%s' % (dbse, 34)], [-1, +1]))

ACTV['%s-%s' % (dbse, 35)] = ['%s-%s-reagent' % (dbse, 'CH3'   ),
                              '%s-%s-reagent' % (dbse, 'C2H4'  ),
                              '%s-%s-reagent' % (dbse, 'C3H7ts') ]
RXNM['%s-%s' % (dbse, 35)] = dict(zip(ACTV['%s-%s' % (dbse, 35)], [-1, -1, +1]))

ACTV['%s-%s' % (dbse, 36)] = ['%s-%s-reagent' % (dbse, 'C3H7'  ),
                              '%s-%s-reagent' % (dbse, 'C3H7ts') ]
RXNM['%s-%s' % (dbse, 36)] = dict(zip(ACTV['%s-%s' % (dbse, 36)], [-1, +1]))

ACTV['%s-%s' % (dbse, 37)] = ['%s-%s-reagent' % (dbse, 'HCN'  ),
                              '%s-%s-reagent' % (dbse, 'HCNts') ]
RXNM['%s-%s' % (dbse, 37)] = dict(zip(ACTV['%s-%s' % (dbse, 37)], [-1, +1]))

ACTV['%s-%s' % (dbse, 38)] = ['%s-%s-reagent' % (dbse, 'HNC'  ),
                              '%s-%s-reagent' % (dbse, 'HCNts') ]
RXNM['%s-%s' % (dbse, 38)] = dict(zip(ACTV['%s-%s' % (dbse, 38)], [-1, +1]))

# <<< Reference Values >>>
BIND = {}
BIND['%s-%s' % (dbse,  1)] =  18.14
BIND['%s-%s' % (dbse,  2)] =  83.22
BIND['%s-%s' % (dbse,  3)] =  42.18
BIND['%s-%s' % (dbse,  4)] =  42.18
BIND['%s-%s' % (dbse,  5)] =  18.00
BIND['%s-%s' % (dbse,  6)] =  18.00
BIND['%s-%s' % (dbse,  7)] =  30.38
BIND['%s-%s' % (dbse,  8)] =  57.02
BIND['%s-%s' % (dbse,  9)] =   2.27
BIND['%s-%s' % (dbse, 10)] = 106.18
BIND['%s-%s' % (dbse, 11)] =   7.43
BIND['%s-%s' % (dbse, 12)] =  60.17
BIND['%s-%s' % (dbse, 13)] =  -0.34
BIND['%s-%s' % (dbse, 14)] =  -0.34
BIND['%s-%s' % (dbse, 15)] =  13.38
BIND['%s-%s' % (dbse, 16)] =  13.38
BIND['%s-%s' % (dbse, 17)] =   3.10
BIND['%s-%s' % (dbse, 18)] =   3.10
BIND['%s-%s' % (dbse, 19)] =  13.61
BIND['%s-%s' % (dbse, 20)] =  13.61
BIND['%s-%s' % (dbse, 21)] = -12.54
BIND['%s-%s' % (dbse, 22)] =  20.11
BIND['%s-%s' % (dbse, 23)] =   2.89
BIND['%s-%s' % (dbse, 24)] =  29.62
BIND['%s-%s' % (dbse, 25)] =  -2.78
BIND['%s-%s' % (dbse, 26)] =  17.33
BIND['%s-%s' % (dbse, 27)] =  10.96
BIND['%s-%s' % (dbse, 28)] =  47.20
BIND['%s-%s' % (dbse, 29)] =  14.69
BIND['%s-%s' % (dbse, 30)] =  10.72
BIND['%s-%s' % (dbse, 31)] =   3.17
BIND['%s-%s' % (dbse, 32)] =  22.68
BIND['%s-%s' % (dbse, 33)] =   1.72
BIND['%s-%s' % (dbse, 34)] =  41.75
BIND['%s-%s' % (dbse, 35)] =   6.85
BIND['%s-%s' % (dbse, 36)] =  32.97
BIND['%s-%s' % (dbse, 37)] =  48.16
BIND['%s-%s' % (dbse, 38)] =  33.11

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s' % (dbse,  1)] = '{ H + N2O <-- [HN2O] } --> OH + N2'
TAGL['%s-%s' % (dbse,  2)] = 'H + N2O <-- { [HN2O] --> OH + N2 }'
TAGL['%s-%s' % (dbse,  3)] = '{ H + FH <-- [HFH] } --> HF + H'
TAGL['%s-%s' % (dbse,  4)] = 'H + FH <-- { [HFH] --> HF + H }'
TAGL['%s-%s' % (dbse,  5)] = '{ H + ClH <-- [HClH] } --> HCl + H'
TAGL['%s-%s' % (dbse,  6)] = 'H + ClH <-- { [HClH] --> HCl + H }'
TAGL['%s-%s' % (dbse,  7)] = '{ H + FCH3 <-- [HFCH3] } --> HF + CH3'
TAGL['%s-%s' % (dbse,  8)] = 'H + FCH3 <-- { [HFCH3] --> HF + CH3 }'
TAGL['%s-%s' % (dbse,  9)] = '{ H + F2 <-- [HF2] } --> HF + F'
TAGL['%s-%s' % (dbse, 10)] = 'H + F2 <-- { [HF2] --> HF + F }'
TAGL['%s-%s' % (dbse, 11)] = '{ CH3 + FCl <-- [CH3FCl] } --> CH3F + Cl'
TAGL['%s-%s' % (dbse, 12)] = 'CH3 + FCl <-- { [CH3FCl] --> CH3F + Cl }'
TAGL['%s-%s' % (dbse, 13)] = '{ F- + CH3F <-- [FCH3F-] } --> FCH3 + F-'
TAGL['%s-%s' % (dbse, 14)] = 'F- + CH3F <-- { [FCH3F-] --> FCH3 + F- }'
TAGL['%s-%s' % (dbse, 15)] = '{ F- ... CH3F <-- [FCH3F-] } --> FCH3 ... F-'
TAGL['%s-%s' % (dbse, 16)] = 'F- ... CH3F <-- { [FCH3F-] --> FCH3 ... F- }'
TAGL['%s-%s' % (dbse, 17)] = '{ Cl- + CH3Cl <-- [ClCH3Cl-] } --> ClCH3 + Cl-'
TAGL['%s-%s' % (dbse, 18)] = 'Cl- + CH3Cl <-- { [ClCH3Cl-] --> ClCH3 + Cl- }'
TAGL['%s-%s' % (dbse, 19)] = '{ Cl- ... CH3Cl <-- [ClCH3Cl-] } --> ClCH3 ... Cl-'
TAGL['%s-%s' % (dbse, 20)] = 'Cl- ... CH3Cl <-- { [ClCH3Cl-] --> ClCH3 ... Cl- }'
TAGL['%s-%s' % (dbse, 21)] = '{ F- + CH3Cl <-- [FCH3Cl-] } --> FCH3 + Cl-'
TAGL['%s-%s' % (dbse, 22)] = 'F- + CH3Cl <-- { [FCH3Cl-] --> FCH3 + Cl- }'
TAGL['%s-%s' % (dbse, 23)] = '{ F- ... CH3Cl <-- [FCH3Cl-] } --> FCH3 ... Cl-'
TAGL['%s-%s' % (dbse, 24)] = 'F- ... CH3Cl <-- { [FCH3Cl-] --> FCH3 ... Cl- }'
TAGL['%s-%s' % (dbse, 25)] = '{ OH- + CH3F <-- [OHCH3F-] } --> HOCH3 + F-'
TAGL['%s-%s' % (dbse, 26)] = 'OH- + CH3F <-- { [OHCH3F-] --> HOCH3 + F- }'
TAGL['%s-%s' % (dbse, 27)] = '{ OH- ... CH3F <-- [OHCH3F-] } --> HOCH3 ... F-'
TAGL['%s-%s' % (dbse, 28)] = 'OH- ... CH3F <-- { [OHCH3F-] --> HOCH3 ... F- }'
TAGL['%s-%s' % (dbse, 29)] = '{ H + N2 <-- [HN2] } --> HN2'
TAGL['%s-%s' % (dbse, 30)] = 'H + N2 <-- { [HN2] --> HN2 }'
TAGL['%s-%s' % (dbse, 31)] = '{ H + CO <-- [HCO] } --> HCO'
TAGL['%s-%s' % (dbse, 32)] = 'H + CO <-- { [HCO] --> HCO }'
TAGL['%s-%s' % (dbse, 33)] = '{ H + C2H4 <-- [HC2H4] } --> CH3CH2'
TAGL['%s-%s' % (dbse, 34)] = 'H + C2H4 <-- { [HC2H4] --> CH3CH2 }'
TAGL['%s-%s' % (dbse, 35)] = '{ CH3 + C2H4 <-- [CH3C2H4] } --> CH3CH2CH2'
TAGL['%s-%s' % (dbse, 36)] = 'CH3 + C2H4 <-- { [CH3C2H4] --> CH3CH2CH2 }'
TAGL['%s-%s' % (dbse, 37)] = '{ HCN <-- [HCN] } --> HNC'
TAGL['%s-%s' % (dbse, 38)] = 'HCN <-- { [HCN] --> HNC }'
TAGL['%s-%s-reagent' % (dbse, 'C2H4'       )] = 'Ethene'
TAGL['%s-%s-reagent' % (dbse, 'C2H5ts'     )] = 'Transition State of H + C2H4 <--> CH3CH2'
TAGL['%s-%s-reagent' % (dbse, 'C2H5'       )] = 'C2H5'
TAGL['%s-%s-reagent' % (dbse, 'C3H7ts'     )] = 'Transition State of CH3 + C2H4 <--> CH3CH2CH2'
TAGL['%s-%s-reagent' % (dbse, 'C3H7'       )] = 'C3H7'
TAGL['%s-%s-reagent' % (dbse, 'CH3Cl'      )] = 'CH3Cl'
TAGL['%s-%s-reagent' % (dbse, 'CH3FClts'   )] = 'Transition State of CH3 + FCL <--> CH3F + Cl'
TAGL['%s-%s-reagent' % (dbse, 'CH3F'       )] = 'CH3F'
TAGL['%s-%s-reagent' % (dbse, 'CH3OH'      )] = 'Methanol'
TAGL['%s-%s-reagent' % (dbse, 'CH3'        )] = 'CH3'
TAGL['%s-%s-reagent' % (dbse, 'ClCH3Clcomp')] = 'Complex of Cl- + CH3Cl'
TAGL['%s-%s-reagent' % (dbse, 'ClCH3Clts'  )] = 'Transition State of Cl- + CH3Cl <--> ClCH3 + Cl-'
TAGL['%s-%s-reagent' % (dbse, 'ClF'        )] = 'ClF'
TAGL['%s-%s-reagent' % (dbse, 'Cl_anion'   )] = 'Chloride Anion'
TAGL['%s-%s-reagent' % (dbse, 'Cl'         )] = 'Chlorine Atom'
TAGL['%s-%s-reagent' % (dbse, 'CO'         )] = 'Carbon Monoxide'
TAGL['%s-%s-reagent' % (dbse, 'F2'         )] = 'Fluorine Molecule'
TAGL['%s-%s-reagent' % (dbse, 'FCH3Clcomp1')] = 'Complex of F- + CH3Cl'
TAGL['%s-%s-reagent' % (dbse, 'FCH3Clcomp2')] = 'Complex of FCH3 + Cl-'
TAGL['%s-%s-reagent' % (dbse, 'FCH3Clts'   )] = 'Transition State of F- + CH3Cl <--> FCH3 + Cl-'
TAGL['%s-%s-reagent' % (dbse, 'FCH3Fcomp'  )] = 'Complex of F- + CH3F'
TAGL['%s-%s-reagent' % (dbse, 'FCH3Fts'    )] = 'Transition State of F- CH3F <--> FCH3 + F-'
TAGL['%s-%s-reagent' % (dbse, 'F_anion'    )] = 'Fluoride Anion'
TAGL['%s-%s-reagent' % (dbse, 'F'          )] = 'Fluorine Atom'
TAGL['%s-%s-reagent' % (dbse, 'HClHts'     )] = 'Transition State of H + ClH <--> HCl + H'
TAGL['%s-%s-reagent' % (dbse, 'HCl'        )] = 'Hydrogen Chloride'
TAGL['%s-%s-reagent' % (dbse, 'HCNts'      )] = 'Transition State of HCN <--> HNC'
TAGL['%s-%s-reagent' % (dbse, 'HCN'        )] = 'Hydrogen Cyanide'
TAGL['%s-%s-reagent' % (dbse, 'HCOts'      )] = 'Transition State of H + CO <--> HCO'
TAGL['%s-%s-reagent' % (dbse, 'HCO'        )] = 'HCO'
TAGL['%s-%s-reagent' % (dbse, 'HF2ts'      )] = 'Transition State of H + F2 <--> HF + F'
TAGL['%s-%s-reagent' % (dbse, 'HFCH3ts'    )] = 'Transition State of H + FCH3 <--> HF + CH3'
TAGL['%s-%s-reagent' % (dbse, 'HFHts'      )] = 'Transition State of H + FH <--> HF + H'
TAGL['%s-%s-reagent' % (dbse, 'HF'         )] = 'Hydrogen Fluoride'
TAGL['%s-%s-reagent' % (dbse, 'HN2ts'      )] = 'Transition State of H + N2 <--> HN2'
TAGL['%s-%s-reagent' % (dbse, 'HN2'        )] = 'HN2'
TAGL['%s-%s-reagent' % (dbse, 'HNC'        )] = 'HNC'
TAGL['%s-%s-reagent' % (dbse, 'HOCH3Fcomp1')] = 'Complex of HOCH3 + F-'
TAGL['%s-%s-reagent' % (dbse, 'HOCH3Fcomp2')] = 'Complex of OH- + CH3F'
TAGL['%s-%s-reagent' % (dbse, 'HOCH3Fts'   )] = 'Transition State of OH- + CH3F <--> HOCH3 + F-'
TAGL['%s-%s-reagent' % (dbse, 'H'          )] = 'Hydrogen Atom'
TAGL['%s-%s-reagent' % (dbse, 'N2OHts'     )] = 'Transition State of H + N2O <--> OH + N2'
TAGL['%s-%s-reagent' % (dbse, 'N2O'        )] = 'N2O'
TAGL['%s-%s-reagent' % (dbse, 'N2'         )] = 'Nitrogen Molecule'
TAGL['%s-%s-reagent' % (dbse, 'OH_anion'   )] = 'Hydroxide Anion'
TAGL['%s-%s-reagent' % (dbse, 'OH'         )] = 'OH'

# <<< Molecule Specifications >>>
NHTBH_C2H4 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.66559300
C        0.00000000    -0.00000000    -0.66559300
H        0.00000000     0.92149500     1.23166800
H        0.00000000    -0.92149500     1.23166800
H        0.00000000     0.92149500    -1.23166800
H        0.00000000    -0.92149500    -1.23166800
units angstrom
}
""", 0)

NHTBH_C2H5ts = input.process_input("""
molecule dimer {
0 2
C       -0.56787700     0.00005100    -0.21895800
C        0.75113900    -0.00003600     0.04193200
H       -1.49388400    -0.00048800     1.53176500
H       -1.10169100     0.92065100    -0.40862600
H       -1.10202200    -0.92023400    -0.40911000
H        1.29912800    -0.92234400     0.17376300
H        1.29889900     0.92232500     0.17436300
units angstrom
}
""", 0)

NHTBH_C2H5 = input.process_input("""
molecule dimer {
0 2
C       -0.25871900    -0.81682900     0.00000000
C       -0.25098700     0.67419100     0.00000000
H        0.75883000    -1.22593900     0.00000000
H       -0.75883000    -1.21386600     0.88341900
H       -0.75883000    -1.21386600    -0.88341900
H       -0.17002100     1.22593900    -0.92432000
H       -0.17002100     1.22593900     0.92432000
units angstrom
}
""", 0)

NHTBH_C3H7ts = input.process_input("""
molecule dimer {
0 2
C       -0.47213200     0.64593300    -0.00004300
C       -1.38261700    -0.36388500    -0.00000200
H       -0.23204400     1.16457500    -0.91726400
H       -0.23234200     1.16475900     0.91716900
H       -1.72712800    -0.80981000     0.92251900
H       -1.72693600    -0.81013100    -0.92243500
C        1.61201500    -0.24218900     0.00003500
H        2.19518200     0.66867100    -0.00126900
H        1.58942300    -0.80961900    -0.91863200
H        1.59024500    -0.80759800     0.91996900
units angstrom
}
""", 0)

NHTBH_C3H7 = input.process_input("""
molecule dimer {
0 2
C        1.20844000    -0.28718900     0.00005700
C       -0.06535900     0.57613200    -0.00005700
C       -1.31478700    -0.23951800    -0.00001100
H        1.24136900    -0.92839500     0.88123400
H        1.24139400    -0.92858600    -0.88098000
H        2.10187100     0.33872700     0.00000000
H       -0.04821800     1.22685100    -0.87708900
H       -0.04827200     1.22703700     0.87683400
H       -1.72914600    -0.61577100     0.92443500
H       -1.72876300    -0.61641500    -0.92436900
units angstrom
}
""", 0)

NHTBH_CH3Cl = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -1.12588600
Cl       0.00000000     0.00000000     0.65683000
H        0.00000000     1.02799300    -1.47026400
H        0.89026800    -0.51399700    -1.47026400
H       -0.89026800    -0.51399700    -1.47026400
units angstrom
}
""", 0)

NHTBH_CH3FClts = input.process_input("""
molecule dimer {
0 2
Cl       1.45474900    -0.00123700    -0.00004000
F       -0.32358700     0.00463100     0.00012400
C       -2.38741800    -0.00214700    -0.00007300
H       -2.49508600    -0.85536100    -0.64940400
H       -2.49731300    -0.13867300     1.06313900
H       -2.50153700     0.98626900    -0.41373400
units angstrom
}
""", 0)

NHTBH_CH3F = input.process_input("""
molecule dimer {
0 1
C       -0.63207400     0.00000100    -0.00000000
F        0.74911700     0.00000200    -0.00000200
H       -0.98318200    -0.33848900     0.97262500
H       -0.98322200     1.01155300    -0.19317200
H       -0.98320300    -0.67308400    -0.77943700
units angstrom
}
""", 0)

NHTBH_CH3OH = input.process_input("""
molecule dimer {
0 1
C       -0.04642300     0.66306900     0.00000000
O       -0.04642300    -0.75506300     0.00000000
H       -1.08695600     0.97593800     0.00000000
H        0.86059200    -1.05703900     0.00000000
H        0.43814500     1.07159400     0.88953900
H        0.43814500     1.07159400    -0.88953900
units angstrom
}
""", 0)

NHTBH_CH3 = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000     0.00000000
H        1.07731727     0.00000000     0.00000000
H       -0.53865863     0.93298412     0.00000000
H       -0.53865863    -0.93298412    -0.00000000
units angstrom
}
""", 0)

NHTBH_ClCH3Clcomp = input.process_input("""
molecule dimer {
-1 1
Cl       0.00000000     0.00000000    -2.38473500
C        0.00000000     0.00000000    -0.56633100
H        0.00000000     1.02506600    -0.22437900
H       -0.88773400    -0.51253300    -0.22437900
H        0.88773400    -0.51253300    -0.22437900
Cl       0.00000000     0.00000000     2.62421300
units angstrom
}
""", 0)

NHTBH_ClCH3Clts = input.process_input("""
molecule dimer {
-1 1
Cl       2.32258100    -0.00013200     0.00014000
C       -0.00008500     0.00049100    -0.00050900
H        0.00007700    -0.74429000    -0.76760500
H       -0.00032000    -0.29144300     1.02802100
H        0.00008100     1.03721800    -0.26195900
Cl      -2.32254200    -0.00012900     0.00013000
units angstrom
}
""", 0)

NHTBH_ClF = input.process_input("""
molecule dimer {
0 1
F        0.00000000     0.00000000     0.00000000
Cl       1.63033021     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_Cl_anion = input.process_input("""
molecule dimer {
-1 1
Cl       0.00000000     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_Cl = input.process_input("""
molecule dimer {
0 2
Cl       0.00000000     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_CO = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000     0.00000000
C        1.12960815     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_F2 = input.process_input("""
molecule dimer {
0 1
F        0.00000000     0.00000000     0.00000000
F        1.39520410     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_FCH3Clcomp1 = input.process_input("""
molecule dimer {
-1 1
Cl       0.00000000     0.00000000     1.62313800
C        0.00000000     0.00000000    -0.22735800
H        0.00000000     1.02632100    -0.55514100
H        0.88882000    -0.51316000    -0.55514100
H       -0.88882000    -0.51316000    -0.55514100
F        0.00000000     0.00000000    -2.72930800
units angstrom
}
""", 0)

NHTBH_FCH3Clcomp2 = input.process_input("""
molecule dimer {
-1 1
F        0.00000000     0.00000000    -2.64853900
C        0.00000000     0.00000000    -1.24017000
H        0.00000000     1.02471900    -0.88640600
H       -0.88743200    -0.51235900    -0.88640600
H        0.88743200    -0.51235900    -0.88640600
Cl       0.00000000     0.00000000     1.99629900
units angstrom
}
""", 0)

NHTBH_FCH3Clts = input.process_input("""
molecule dimer {
-1 1
F        0.00000000     0.00000000    -2.53792900
C        0.00000000     0.00000000    -0.48837200
H        0.00000000     1.06208700    -0.61497200
H       -0.91979500    -0.53104400    -0.61497200
H        0.91979500    -0.53104400    -0.61497200
Cl       0.00000000     0.00000000     1.62450100
units angstrom
}
""", 0)

NHTBH_FCH3Fcomp = input.process_input("""
molecule dimer {
-1 1
F        0.00000000     0.00000000    -1.84762600
C        0.00000000     0.00000000    -0.42187300
H        0.00000000     1.02358100    -0.07384300
H       -0.88644700    -0.51179100    -0.07384300
H        0.88644700    -0.51179100    -0.07384300
F        0.00000000     0.00000000     2.15348900
units angstrom
}
""", 0)

NHTBH_FCH3Fts = input.process_input("""
molecule dimer {
-1 1
F        0.00309800    -0.01889200    -0.01545600
C       -0.00014900    -0.00014000     1.80785700
H        1.06944900     0.00170800     1.80976100
H       -0.53660700     0.92513300     1.79693500
H       -0.53260100    -0.92778300     1.81705800
F       -0.00319100     0.01997400     3.63184500
units angstrom
}
""", 0)

NHTBH_F_anion = input.process_input("""
molecule dimer {
-1 1
F        0.00000000     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_F = input.process_input("""
molecule dimer {
0 2
F        0.00000000     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_HClHts = input.process_input("""
molecule dimer {
0 2
H        0.00000000     0.00000000     1.48580000
Cl       0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000    -1.48580000
units angstrom
}
""", 0)

NHTBH_HCl = input.process_input("""
molecule dimer {
0 1
Cl       0.00000000     0.00000000     0.00000000
H        1.27444789     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_HCNts = input.process_input("""
molecule dimer {
0 1
C        0.08031900     0.62025800     0.00000000
N        0.08031900    -0.56809500     0.00000000
H       -1.04414800     0.25512100     0.00000000
units angstrom
}
""", 0)

NHTBH_HCN = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -0.50036500
N        0.00000000     0.00000000     0.65264000
H        0.00000000     0.00000000    -1.56629100
units angstrom
}
""", 0)

NHTBH_HCOts = input.process_input("""
molecule dimer {
0 2
H       -1.52086400     1.38882900     0.00000000
C        0.10863300     0.54932900     0.00000000
O        0.10863300    -0.58560100     0.00000000
units angstrom
}
""", 0)

NHTBH_HCO = input.process_input("""
molecule dimer {
0 2
H       -0.00905700     0.00000000    -0.00708600
C       -0.00703500     0.00000000     1.10967800
O        0.95604000     0.00000000     1.78565600
units angstrom
}
""", 0)

NHTBH_HF2ts = input.process_input("""
molecule dimer {
0 2
H        0.00000000     0.00000000    -2.23127300
F        0.00000000     0.00000000    -0.61621800
F        0.00000000     0.00000000     0.86413800
units angstrom
}
""", 0)

NHTBH_HFCH3ts = input.process_input("""
molecule dimer {
0 2
H       -0.03976400     0.00000000     0.04410600
F       -0.04932100     0.00000000     1.28255400
C       -0.06154400     0.00000000     2.95115700
H        0.99049700     0.00000000     3.19427500
H       -0.59007000     0.91235500     3.18348100
H       -0.59007000    -0.91235500     3.18348100
units angstrom
}
""", 0)

NHTBH_HFHts = input.process_input("""
molecule dimer {
0 2
H        0.00000000     0.00000000     1.13721700
F        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000    -1.13721700
units angstrom
}
""", 0)

NHTBH_HF = input.process_input("""
molecule dimer {
0 1
F        0.00000000     0.00000000     0.00000000
H        0.91538107     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_HN2ts = input.process_input("""
molecule dimer {
0 2
N        0.00000000     0.00000000     0.00000000
N        1.12281100     0.00000000     0.00000000
H        1.78433286     1.26844651     0.00000000
units angstrom
}
""", 0)

NHTBH_HN2 = input.process_input("""
molecule dimer {
0 2
N        0.00000000     0.00000000     0.00000000
N        1.17820000     0.00000000     0.00000000
H        1.64496947     0.93663681     0.00000000
units angstrom
}
""", 0)

NHTBH_HNC = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000    -0.73724800
N        0.00000000     0.00000000     0.43208900
H        0.00000000     0.00000000     1.42696000
units angstrom
}
""", 0)

NHTBH_HOCH3Fcomp1 = input.process_input("""
molecule dimer {
-1 1
C       -1.29799700    -0.38951800    -0.00003400
O       -0.47722300     0.72802100     0.00005400
H       -2.35192200    -0.08023200    -0.00863900
H       -1.14085300    -1.03582100    -0.87810100
H       -1.15317800    -1.02751300     0.88635900
H        0.51058000     0.37116000     0.00024300
F        1.74901600    -0.19051700    -0.00001000
units angstrom
}
""", 0)

NHTBH_HOCH3Fcomp2 = input.process_input("""
molecule dimer {
-1 1
F        0.00037100    -2.46834000     0.02139000
C       -0.27664200    -1.07441800    -0.00269000
H        0.64929000    -0.51650000    -0.00901600
H       -0.84198900    -0.84711900    -0.89707500
H       -0.85102800    -0.82658900     0.88141700
O       -0.30171300     1.58252400    -0.20654400
H       -0.60511200     2.49243400    -0.16430500
units angstrom
}
""", 0)

NHTBH_HOCH3Fts = input.process_input("""
molecule dimer {
-1 1
F        0.02253600    -0.00745300     0.00552900
C       -0.01842000     0.00503700     1.76492500
H        1.04805000     0.00524000     1.85414600
H       -0.54781900     0.93470700     1.79222400
H       -0.54895500    -0.92343300     1.80576200
O        0.00126500     0.01920000     3.75059900
H       -0.92676300     0.03161500     3.99758100
units angstrom
}
""", 0)

NHTBH_H = input.process_input("""
molecule dimer {
0 2
H        0.00000000     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_N2OHts = input.process_input("""
molecule dimer {
0 2
H       -0.30328600    -1.93071200     0.00000000
O       -0.86100600    -0.62152600     0.00000000
N        0.00000000     0.25702700     0.00000000
N        1.02733300     0.72910400     0.00000000
units angstrom
}
""", 0)

NHTBH_N2O = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     0.00000000
N        1.12056262     0.00000000     0.00000000
O        2.30761092     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_N2 = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     0.00000000
N        1.09710935     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_OH_anion = input.process_input("""
molecule dimer {
-1 1
O        0.00000000     0.00000000     0.00000000
H        0.96204317     0.00000000     0.00000000
units angstrom
}
""", 0)

NHTBH_OH = input.process_input("""
molecule dimer {
0 2
O        0.00000000     0.00000000     0.00000000
H        0.96889819     0.00000000     0.00000000
units angstrom
}
""", 0)

# <<< Geometry Specification Strings >>>
rxnpattern = re.compile(r'^(.+)-(.+)-(.+)$')
GEOS = {}
for rxn in HRXN:
    for rgt in ACTV['%s-%s' % (dbse, rxn)]:

        molname = rxnpattern.match(rgt)
        GEOS['%s' % (rgt)] = eval('%s_%s' % (dbse, molname.group(2)))
