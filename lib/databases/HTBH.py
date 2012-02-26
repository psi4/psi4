"""
**HTBH**

| Database (Truhlar) of hydrogen-transfer barrier height reactions.
| Geometries from Truhlar and coworkers at site http://t1.chem.umn.edu/misc/database_group/database_therm_bh/raw_geom.cgi .
| Reference energies from Zhao et al. JPCA, 109 2012-2018 (2005) doi: 10.1021/jp045141s [in supporting information].

- **cp**  ``'off'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``

----

"""
import re
import input

# <<< HTBH Database Module >>>
dbse = 'HTBH'
isOS = 'true'

# <<< Database Members >>>
HRXN = range (1, 39)
HRXN_SM = ['5','6','9','10','23','24']
HRXN_LG = ['13','14','33','34','37','38']

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse,  1)] = ['%s-%s-reagent'      % (dbse, 'H'     ),
                                         '%s-%s-reagent'      % (dbse, 'HCl'   ),
                                         '%s-%s-reagent'      % (dbse, 'HHClts') ]
RXNM['%s-%s'            % (dbse,  1)] = dict(zip(ACTV['%s-%s' % (dbse,  1)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse,  2)] = ['%s-%s-reagent'      % (dbse, 'H2'    ),
                                         '%s-%s-reagent'      % (dbse, 'Cl'    ),
                                         '%s-%s-reagent'      % (dbse, 'HHClts') ]
RXNM['%s-%s'            % (dbse,  2)] = dict(zip(ACTV['%s-%s' % (dbse,  2)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse,  3)] = ['%s-%s-reagent'      % (dbse, 'OH'    ),
                                         '%s-%s-reagent'      % (dbse, 'H2'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHH2ts') ]
RXNM['%s-%s'            % (dbse,  3)] = dict(zip(ACTV['%s-%s' % (dbse,  3)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse,  4)] = ['%s-%s-reagent'      % (dbse, 'H'     ),
                                         '%s-%s-reagent'      % (dbse, 'H2O'   ),
                                         '%s-%s-reagent'      % (dbse, 'OHH2ts') ]
RXNM['%s-%s'            % (dbse,  4)] = dict(zip(ACTV['%s-%s' % (dbse,  4)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse,  5)] = ['%s-%s-reagent'      % (dbse, 'CH3'    ),
                                         '%s-%s-reagent'      % (dbse, 'H2'     ),
                                         '%s-%s-reagent'      % (dbse, 'CH3H2ts') ]
RXNM['%s-%s'            % (dbse,  5)] = dict(zip(ACTV['%s-%s' % (dbse,  5)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse,  6)] = ['%s-%s-reagent'      % (dbse, 'H'      ),
                                         '%s-%s-reagent'      % (dbse, 'CH4'    ),
                                         '%s-%s-reagent'      % (dbse, 'CH3H2ts') ]
RXNM['%s-%s'            % (dbse,  6)] = dict(zip(ACTV['%s-%s' % (dbse,  6)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse,  7)] = ['%s-%s-reagent'      % (dbse, 'OH'     ),
                                         '%s-%s-reagent'      % (dbse, 'CH4'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHCH4ts') ]
RXNM['%s-%s'            % (dbse,  7)] = dict(zip(ACTV['%s-%s' % (dbse,  7)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse,  8)] = ['%s-%s-reagent'      % (dbse, 'CH3'    ),
                                         '%s-%s-reagent'      % (dbse, 'H2O'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHCH4ts') ]
RXNM['%s-%s'            % (dbse,  8)] = dict(zip(ACTV['%s-%s' % (dbse,  8)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse,  9)] = ['%s-%s-reagent'      % (dbse, 'H'    ),
                                         '%s-%s-reagent'      % (dbse, 'H2'   ),
                                         '%s-%s-reagent'      % (dbse, 'HH2ts') ]
RXNM['%s-%s'            % (dbse,  9)] = dict(zip(ACTV['%s-%s' % (dbse,  9)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 10)] = ['%s-%s-reagent'      % (dbse, 'H2'   ),
                                         '%s-%s-reagent'      % (dbse, 'H'    ),
                                         '%s-%s-reagent'      % (dbse, 'HH2ts') ]
RXNM['%s-%s'            % (dbse, 10)] = dict(zip(ACTV['%s-%s' % (dbse, 10)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 11)] = ['%s-%s-reagent'      % (dbse, 'OH'     ),
                                         '%s-%s-reagent'      % (dbse, 'NH3'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHNH3ts') ]
RXNM['%s-%s'            % (dbse, 11)] = dict(zip(ACTV['%s-%s' % (dbse, 11)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 12)] = ['%s-%s-reagent'      % (dbse, 'H2O'    ),
                                         '%s-%s-reagent'      % (dbse, 'NH2'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHNH3ts') ]
RXNM['%s-%s'            % (dbse, 12)] = dict(zip(ACTV['%s-%s' % (dbse, 12)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 13)] = ['%s-%s-reagent'      % (dbse, 'HCl'     ),
                                         '%s-%s-reagent'      % (dbse, 'CH3'     ),
                                         '%s-%s-reagent'      % (dbse, 'HClCH3ts') ]
RXNM['%s-%s'            % (dbse, 13)] = dict(zip(ACTV['%s-%s' % (dbse, 13)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 14)] = ['%s-%s-reagent'      % (dbse, 'Cl'      ),
                                         '%s-%s-reagent'      % (dbse, 'CH4'     ),
                                         '%s-%s-reagent'      % (dbse, 'HClCH3ts') ]
RXNM['%s-%s'            % (dbse, 14)] = dict(zip(ACTV['%s-%s' % (dbse, 14)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 15)] = ['%s-%s-reagent'      % (dbse, 'OH'      ),
                                         '%s-%s-reagent'      % (dbse, 'C2H6'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHC2H6ts') ]
RXNM['%s-%s'            % (dbse, 15)] = dict(zip(ACTV['%s-%s' % (dbse, 15)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 16)] = ['%s-%s-reagent'      % (dbse, 'H2O'     ),
                                         '%s-%s-reagent'      % (dbse, 'C2H5'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHC2H6ts') ]
RXNM['%s-%s'            % (dbse, 16)] = dict(zip(ACTV['%s-%s' % (dbse, 16)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 17)] = ['%s-%s-reagent'      % (dbse, 'F'    ),
                                         '%s-%s-reagent'      % (dbse, 'H2'   ),
                                         '%s-%s-reagent'      % (dbse, 'FH2ts') ]
RXNM['%s-%s'            % (dbse, 17)] = dict(zip(ACTV['%s-%s' % (dbse, 17)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 18)] = ['%s-%s-reagent'      % (dbse, 'HF'   ),
                                         '%s-%s-reagent'      % (dbse, 'H'    ),
                                         '%s-%s-reagent'      % (dbse, 'FH2ts') ]
RXNM['%s-%s'            % (dbse, 18)] = dict(zip(ACTV['%s-%s' % (dbse, 18)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 19)] = ['%s-%s-reagent'      % (dbse, 'O'      ),
                                         '%s-%s-reagent'      % (dbse, 'CH4'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHCH3ts') ]
RXNM['%s-%s'            % (dbse, 19)] = dict(zip(ACTV['%s-%s' % (dbse, 19)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 20)] = ['%s-%s-reagent'      % (dbse, 'OH'     ),
                                         '%s-%s-reagent'      % (dbse, 'CH3'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHCH3ts') ]
RXNM['%s-%s'            % (dbse, 20)] = dict(zip(ACTV['%s-%s' % (dbse, 20)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 21)] = ['%s-%s-reagent'      % (dbse, 'H'     ),
                                         '%s-%s-reagent'      % (dbse, 'PH3'   ),
                                         '%s-%s-reagent'      % (dbse, 'HPH3ts') ]
RXNM['%s-%s'            % (dbse, 21)] = dict(zip(ACTV['%s-%s' % (dbse, 21)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 22)] = ['%s-%s-reagent'      % (dbse, 'PH2'   ),
                                         '%s-%s-reagent'      % (dbse, 'H2'    ),
                                         '%s-%s-reagent'      % (dbse, 'HPH3ts') ]
RXNM['%s-%s'            % (dbse, 22)] = dict(zip(ACTV['%s-%s' % (dbse, 22)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 23)] = ['%s-%s-reagent'      % (dbse, 'H'    ),
                                         '%s-%s-reagent'      % (dbse, 'OH'   ),
                                         '%s-%s-reagent'      % (dbse, 'OHHts') ]
RXNM['%s-%s'            % (dbse, 23)] = dict(zip(ACTV['%s-%s' % (dbse, 23)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 24)] = ['%s-%s-reagent'      % (dbse, 'H2'   ),
                                         '%s-%s-reagent'      % (dbse, 'O'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHHts') ]
RXNM['%s-%s'            % (dbse, 24)] = dict(zip(ACTV['%s-%s' % (dbse, 24)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 25)] = ['%s-%s-reagent'      % (dbse, 'H'     ),
                                         '%s-%s-reagent'      % (dbse, 'H2S'   ),
                                         '%s-%s-reagent'      % (dbse, 'HH2Sts') ]
RXNM['%s-%s'            % (dbse, 25)] = dict(zip(ACTV['%s-%s' % (dbse, 25)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 26)] = ['%s-%s-reagent'      % (dbse, 'H2'    ),
                                         '%s-%s-reagent'      % (dbse, 'HS'    ),
                                         '%s-%s-reagent'      % (dbse, 'HH2Sts') ]
RXNM['%s-%s'            % (dbse, 26)] = dict(zip(ACTV['%s-%s' % (dbse, 26)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 27)] = ['%s-%s-reagent'      % (dbse, 'O'     ),
                                         '%s-%s-reagent'      % (dbse, 'HCl'   ),
                                         '%s-%s-reagent'      % (dbse, 'OHClts') ]
RXNM['%s-%s'            % (dbse, 27)] = dict(zip(ACTV['%s-%s' % (dbse, 27)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 28)] = ['%s-%s-reagent'      % (dbse, 'OH'    ),
                                         '%s-%s-reagent'      % (dbse, 'Cl'    ),
                                         '%s-%s-reagent'      % (dbse, 'OHClts') ]
RXNM['%s-%s'            % (dbse, 28)] = dict(zip(ACTV['%s-%s' % (dbse, 28)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 29)] = ['%s-%s-reagent'      % (dbse, 'NH2'     ),
                                         '%s-%s-reagent'      % (dbse, 'CH3'     ),
                                         '%s-%s-reagent'      % (dbse, 'CH3NH2ts') ]
RXNM['%s-%s'            % (dbse, 29)] = dict(zip(ACTV['%s-%s' % (dbse, 29)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 30)] = ['%s-%s-reagent'      % (dbse, 'CH4'     ),
                                         '%s-%s-reagent'      % (dbse, 'NH'      ),
                                         '%s-%s-reagent'      % (dbse, 'CH3NH2ts') ]
RXNM['%s-%s'            % (dbse, 30)] = dict(zip(ACTV['%s-%s' % (dbse, 30)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 31)] = ['%s-%s-reagent'      % (dbse, 'NH2'      ),
                                         '%s-%s-reagent'      % (dbse, 'C2H5'     ),
                                         '%s-%s-reagent'      % (dbse, 'NH2C2H5ts') ]
RXNM['%s-%s'            % (dbse, 31)] = dict(zip(ACTV['%s-%s' % (dbse, 31)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 32)] = ['%s-%s-reagent'      % (dbse, 'C2H6'     ),
                                         '%s-%s-reagent'      % (dbse, 'NH'       ),
                                         '%s-%s-reagent'      % (dbse, 'NH2C2H5ts') ]
RXNM['%s-%s'            % (dbse, 32)] = dict(zip(ACTV['%s-%s' % (dbse, 32)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 33)] = ['%s-%s-reagent'      % (dbse, 'C2H6'     ),
                                         '%s-%s-reagent'      % (dbse, 'NH2'      ),
                                         '%s-%s-reagent'      % (dbse, 'C2H6NH2ts') ]
RXNM['%s-%s'            % (dbse, 33)] = dict(zip(ACTV['%s-%s' % (dbse, 33)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 34)] = ['%s-%s-reagent'      % (dbse, 'NH3'      ),
                                         '%s-%s-reagent'      % (dbse, 'C2H5'     ),
                                         '%s-%s-reagent'      % (dbse, 'C2H6NH2ts') ]
RXNM['%s-%s'            % (dbse, 34)] = dict(zip(ACTV['%s-%s' % (dbse, 34)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 35)] = ['%s-%s-reagent'      % (dbse, 'NH2'     ),
                                         '%s-%s-reagent'      % (dbse, 'CH4'     ),
                                         '%s-%s-reagent'      % (dbse, 'NH2CH4ts') ]
RXNM['%s-%s'            % (dbse, 35)] = dict(zip(ACTV['%s-%s' % (dbse, 35)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 36)] = ['%s-%s-reagent'      % (dbse, 'CH3'     ),
                                         '%s-%s-reagent'      % (dbse, 'NH3'     ),
                                         '%s-%s-reagent'      % (dbse, 'NH2CH4ts') ]
RXNM['%s-%s'            % (dbse, 36)] = dict(zip(ACTV['%s-%s' % (dbse, 36)], [-1, -1, +1]))

ACTV['%s-%s'            % (dbse, 37)] = ['%s-%s-reagent'      % (dbse, 'C5H8'  ),
                                         '%s-%s-reagent'      % (dbse, 'C5H8ts') ]
RXNM['%s-%s'            % (dbse, 37)] = dict(zip(ACTV['%s-%s' % (dbse, 37)], [-1, +1]))

ACTV['%s-%s'            % (dbse, 38)] = ['%s-%s-reagent'      % (dbse, 'C5H8'  ),
                                         '%s-%s-reagent'      % (dbse, 'C5H8ts') ]
RXNM['%s-%s'            % (dbse, 38)] = dict(zip(ACTV['%s-%s' % (dbse, 38)], [-1, +1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'% (dbse,  1)] =    5.7
BIND['%s-%s'% (dbse,  2)] =    8.7
BIND['%s-%s'% (dbse,  3)] =    5.1
BIND['%s-%s'% (dbse,  4)] =   21.2
BIND['%s-%s'% (dbse,  5)] =   12.1
BIND['%s-%s'% (dbse,  6)] =   15.3
BIND['%s-%s'% (dbse,  7)] =    6.7
BIND['%s-%s'% (dbse,  8)] =   19.6
BIND['%s-%s'% (dbse,  9)] =    9.6
BIND['%s-%s'% (dbse, 10)] =    9.6
BIND['%s-%s'% (dbse, 11)] =    3.2
BIND['%s-%s'% (dbse, 12)] =   12.7
BIND['%s-%s'% (dbse, 13)] =    1.7
BIND['%s-%s'% (dbse, 14)] =    7.9
BIND['%s-%s'% (dbse, 15)] =    3.4
BIND['%s-%s'% (dbse, 16)] =   19.9
BIND['%s-%s'% (dbse, 17)] =    1.8
BIND['%s-%s'% (dbse, 18)] =   33.4
BIND['%s-%s'% (dbse, 19)] =   13.7
BIND['%s-%s'% (dbse, 20)] =    8.1
BIND['%s-%s'% (dbse, 21)] =    3.1
BIND['%s-%s'% (dbse, 22)] =   23.2
BIND['%s-%s'% (dbse, 23)] =   10.7
BIND['%s-%s'% (dbse, 24)] =   13.1
BIND['%s-%s'% (dbse, 25)] =    3.5
BIND['%s-%s'% (dbse, 26)] =   17.3
BIND['%s-%s'% (dbse, 27)] =    9.8
BIND['%s-%s'% (dbse, 28)] =   10.4
BIND['%s-%s'% (dbse, 29)] =    8.0
BIND['%s-%s'% (dbse, 30)] =   22.4
BIND['%s-%s'% (dbse, 31)] =    7.5
BIND['%s-%s'% (dbse, 32)] =   18.3
BIND['%s-%s'% (dbse, 33)] =   10.4
BIND['%s-%s'% (dbse, 34)] =   17.4
BIND['%s-%s'% (dbse, 35)] =   14.5
BIND['%s-%s'% (dbse, 36)] =   17.8
BIND['%s-%s'% (dbse, 37)] =   38.4
BIND['%s-%s'% (dbse, 38)] =   38.4

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s' % (dbse,  1)] = '{ H + HCl <-- [HHCl] } --> H2 + Cl'
TAGL['%s-%s' % (dbse,  2)] = 'H + HCl <-- { [HHCl] --> H2 + Cl }' 
TAGL['%s-%s' % (dbse,  3)] = '{ OH + H2 <-- [OHH2] } --> H + H2O' 
TAGL['%s-%s' % (dbse,  4)] = 'OH + HCl <-- { [OHH2] --> H + H2O }'
TAGL['%s-%s' % (dbse,  5)] = '{ CH3 + H2 <-- [CH3H2] } --> H + CH4'
TAGL['%s-%s' % (dbse,  6)] = 'CH3 + H2 <-- { [CH3H2] --> H + CH4 }' 
TAGL['%s-%s' % (dbse,  7)] = '{ OH + CH4 <-- [OHCH4] } --> CH3 + H2O'
TAGL['%s-%s' % (dbse,  8)] = 'OH + CH4 <-- { [OHCH4] --> CH3 + H2O }' 
TAGL['%s-%s' % (dbse,  9)] = '{ H + H2 <-- [HH2] } --> H2 + H'
TAGL['%s-%s' % (dbse, 10)] =  'H + H2 <-- { [HH2] -- >H2 + H }'
TAGL['%s-%s' % (dbse, 11)] = '{ OH + NH3 <-- [OHNH3] } --> H2O + NH2'
TAGL['%s-%s' % (dbse, 12)] =  'OH + NH3 <-- { [OHNH3] --> H2O + NH2 }'
TAGL['%s-%s' % (dbse, 13)] = '{ HCl + CH3 <-- [HClCH3] } --> Cl + CH4'
TAGL['%s-%s' % (dbse, 14)] =  'HCl + CH3 <-- { [HClCH3] --> Cl + CH4 }'
TAGL['%s-%s' % (dbse, 15)] = '{ OH + C2H6 <-- [OHC2H6] } --> H2O + C2H5'
TAGL['%s-%s' % (dbse, 16)] =  'OH + C2H6 <-- { [OHC2H6] --> H2O + C2H5 }'
TAGL['%s-%s' % (dbse, 17)] = '{ F + H2 <-- [FH2] } --> HF + H'
TAGL['%s-%s' % (dbse, 18)] =  'F + H2 <-- { [FH2] --> HF + H}'
TAGL['%s-%s' % (dbse, 19)] = '{ O + CH4 <-- [OHCH3] } --> OH + CH3'
TAGL['%s-%s' % (dbse, 20)] =  'O + CH4 <-- { [OHCH3] --> OH + CH3 }'
TAGL['%s-%s' % (dbse, 21)] = '{ H + PH3 <-- [HPH3] } --> PH2 + H2'
TAGL['%s-%s' % (dbse, 22)] =  'H + PH3 <-- { [HPH3] --> PH2 + H2 }'
TAGL['%s-%s' % (dbse, 23)] = '{ H + OH <-- [OHH] } --> H2 + O'
TAGL['%s-%s' % (dbse, 24)] =  'H + OH <-- { [OHH] --> H2 + O }'
TAGL['%s-%s' % (dbse, 25)] = '{ H + H2S <-- [HH2S] } --> H2 + HS'
TAGL['%s-%s' % (dbse, 26)] =  'H + H2S <-- { [HH2S] --> H2 + HS}'
TAGL['%s-%s' % (dbse, 27)] = '{ O + HCl <-- [OHCl] } --> OH + Cl'
TAGL['%s-%s' % (dbse, 28)] =  'O + HCl <-- { [OHCl] --> OH + Cl}'
TAGL['%s-%s' % (dbse, 29)] = '{ NH2 + CH3 <-- [CH3NH2] } --> CH4 + NH'
TAGL['%s-%s' % (dbse, 30)] =  'NH2 + CH3 <-- { [CH3NH2] --> CH4 + NH }'
TAGL['%s-%s' % (dbse, 31)] = '{ NH2 + C2H5 <-- [NH2C2H5] } --> C2H6 + NH'
TAGL['%s-%s' % (dbse, 32)] =  'NH2 + C2H5 <-- { [NH2C2H5] --> C2H6 + NH }'
TAGL['%s-%s' % (dbse, 33)] = '{ C2H6 + NH2 <-- [C2H6NH2] } --> NH3 + C2H5'
TAGL['%s-%s' % (dbse, 34)] =  'C2H6 + NH2 <-- { [C2H6NH2] --> NH3 + C2H5 }'
TAGL['%s-%s' % (dbse, 35)] = '{ NH2 + CH4 <-- [NH2CH4] } --> CH3 + NH3'
TAGL['%s-%s' % (dbse, 36)] =  'NH2 + CH4 <-- { [NH2CH4] --> CH3 + NH3 }'
TAGL['%s-%s' % (dbse, 37)] = '{ C5H8 <-- [C5H8] } --> C5H8'
TAGL['%s-%s' % (dbse, 38)] =  'C5H8 <-- { [C5H8] --> C5H8 }'
TAGL['%s-%s-reagent'    % (dbse, 'C2H5'         )] = 'C2H5'
TAGL['%s-%s-reagent'    % (dbse, 'C2H6'         )] = 'Ethane'
TAGL['%s-%s-reagent'    % (dbse, 'C2H6NH2ts'    )] = 'Transition state of C2H6 + NH2 <--> NH3 + C2H5'
TAGL['%s-%s-reagent'    % (dbse, 'C5H8'         )] = 's-trans cis-C5H8'
TAGL['%s-%s-reagent'    % (dbse, 'C5H8ts'       )] = 'Transition state of s-trans cis-C5H8 <--> s-trans cis C5H8'
TAGL['%s-%s-reagent'    % (dbse, 'CH3'          )] = 'CH3'
TAGL['%s-%s-reagent'    % (dbse, 'CH3H2ts'      )] = 'Transition state of CH3 + H2 <--> H + CH4'
TAGL['%s-%s-reagent'    % (dbse, 'CH3NH2ts'     )] = 'Transition state of CH3 + NH2 <--> CH4 + NH'
TAGL['%s-%s-reagent'    % (dbse, 'CH4'          )] = 'Methane'
TAGL['%s-%s-reagent'    % (dbse, 'Cl'           )] = 'Chlorine atom'
TAGL['%s-%s-reagent'    % (dbse, 'F'            )] = 'Fluorine atom'
TAGL['%s-%s-reagent'    % (dbse, 'FH2ts'        )] = 'Transition state of F + H2 <--> HF + H'
TAGL['%s-%s-reagent'    % (dbse, 'H'            )] = 'Hydrogen atom'
TAGL['%s-%s-reagent'    % (dbse, 'H2'           )] = 'Hydrogen molecule'
TAGL['%s-%s-reagent'    % (dbse, 'H2O'          )] = 'Water'
TAGL['%s-%s-reagent'    % (dbse, 'H2S'          )] = 'Hydrogen Sulfide'
TAGL['%s-%s-reagent'    % (dbse, 'HCl'          )] = 'Hydrogen Chloride'
TAGL['%s-%s-reagent'    % (dbse, 'HClCH3ts'     )] = 'Transition state of HCl + CH3 <--> Cl + CH4'
TAGL['%s-%s-reagent'    % (dbse, 'HHClts'       )] = 'Transition state of H + HCl <--> H2 + Cl'
TAGL['%s-%s-reagent'    % (dbse, 'HF'           )] = 'Hydrogen Fluoride'
TAGL['%s-%s-reagent'    % (dbse, 'HH2Sts'       )] = 'Transition state of H + H2S <--> H2 + HS'
TAGL['%s-%s-reagent'    % (dbse, 'HH2ts'        )] = 'Transition state of H + H2 <--> H2 + H' 
TAGL['%s-%s-reagent'    % (dbse, 'NH'           )] = 'NH'
TAGL['%s-%s-reagent'    % (dbse, 'HPH3ts'       )] = 'Transition state of H + PH3 <--> PH2 + H2'
TAGL['%s-%s-reagent'    % (dbse, 'NH2'          )] = 'NH2'
TAGL['%s-%s-reagent'    % (dbse, 'NH2C2H5ts'    )] = 'Transition state of C2H5 + NH2 <--> NH + C2H6'
TAGL['%s-%s-reagent'    % (dbse, 'NH2CH4ts'     )] = 'Transition state of CH4 + NH2 <--> NH3 + CH3'
TAGL['%s-%s-reagent'    % (dbse, 'NH3'          )] = 'Ammonia'
TAGL['%s-%s-reagent'    % (dbse, 'O'            )] = 'Oxygen atom'
TAGL['%s-%s-reagent'    % (dbse, 'OH'           )] = 'OH'
TAGL['%s-%s-reagent'    % (dbse, 'OHC2H6ts'     )] = 'Transition state of C2H6 + OH <--> H2O + C2H5'
TAGL['%s-%s-reagent'    % (dbse, 'OHCH3ts'      )] = 'Transition state of O + CH4 <--> OH + CH3'
TAGL['%s-%s-reagent'    % (dbse, 'OHCH4ts'      )] = 'Transition state of OH + CH4 <--> CH3 + H2O'
TAGL['%s-%s-reagent'    % (dbse, 'OHClts'       )] = 'Transition state of O + HCl <--> OH + Cl'
TAGL['%s-%s-reagent'    % (dbse, 'OHH2ts'       )] = 'Transition state of OH + H2 <--> H + H2O'
TAGL['%s-%s-reagent'    % (dbse, 'OHHts'        )] = 'Transition state of OH + H <--> H2 + O'
TAGL['%s-%s-reagent'    % (dbse, 'OHNH3ts'      )] = 'Transition state of OH + NH3 <--> NH2 + H2O'
TAGL['%s-%s-reagent'    % (dbse, 'PH2'          )] = 'PH2'
TAGL['%s-%s-reagent'    % (dbse, 'PH3'          )] = 'Phosphine'
TAGL['%s-%s-reagent'    % (dbse, 'HS'           )] = 'HS'

# <<< Molecule Specifications >>>
HTBH_C2H5 = input.process_input("""
molecule dimer {
0 2
C        0.00550995    -0.00307714    -0.77443959
C        0.00550995    -0.00307714     0.71569982
H        0.00550995    -1.01684444     1.11670108
H        0.37964525     0.84547158    -1.32730429
H       -0.88217468     0.49798042     1.12141209
H        0.87299475     0.52193057     1.11660682
H       -0.50718726    -0.77526005    -1.32801142
units angstrom
}
""")

HTBH_C2H6 = input.process_input("""
molecule dimer {
0 1
C        0.00000020    -0.00000013    -0.76309187
C        0.00000020    -0.00000013     0.76309163
H        0.00000020    -1.01606691     1.15831231
H       -0.87903844    -0.50959541    -1.15830943
H       -0.87994508     0.50802887     1.15831013
H        0.87993813     0.50804049     1.15830883
H       -0.00180313     1.01606605    -1.15830975
H        0.88084363    -0.50646996    -1.15830912
units angstrom
}
""")

HTBH_C2H6NH2ts = input.process_input("""
molecule dimer {
0 2
C       -1.48570000    -0.44815600    -0.00001900
C       -0.50504200     0.70174000     0.00002900
N        1.86516100    -0.34016700    -0.00005700
H       -1.35419300    -1.07650500    -0.88050300
H       -1.35415900    -1.07661100     0.88038500
H       -2.51702500    -0.08617300     0.00002500
H       -0.52222400     1.31611800    -0.89721800
H       -0.52220500     1.31602900     0.89733800
H        0.66504700     0.14796100    -0.00003400
H        2.24664400     0.15971700    -0.80480600
H        2.24643900     0.15913300     0.80515100
units angstrom
}
""")

HTBH_C5H8 = input.process_input("""
molecule dimer {
0 1
C       -2.05563800    -0.61227200     0.00000700
C       -1.23109600     0.64044800     0.00004900
C        0.10563400     0.73427300     0.00002600
C        1.05755500    -0.37440700    -0.00004400
C        2.38358300    -0.19893600    -0.00003600
H       -2.70508500    -0.64159700     0.87713200
H       -2.70512900    -0.64150800    -0.87708900
H       -1.45133200    -1.51607900    -0.00005500
H       -1.79366500     1.56758600     0.00010300
H        0.54575600     1.72564300     0.00006400
H        0.66526200    -1.38324200    -0.00010500
H        3.06468900    -1.03771900    -0.00008800
H        2.81927500     0.79228500     0.00002300
units angstrom
}
""")

HTBH_C5H8ts = input.process_input("""
molecule dimer {
0 1
C       -1.29962300    -0.90485300    -0.02015500
C       -1.20594700     0.50581700    -0.01341400
C        0.00000000     1.18336100     0.15330100
C        1.20594800     0.50581400    -0.01342200
C        1.29962600    -0.90485100    -0.02014700
H        2.16879700    -1.32754900    -0.51569700
H        1.03204100    -1.45438500     0.87316600
H        2.03713000     1.08558300    -0.39850400
H        0.00000100     2.26291300     0.08590500
H       -2.03713300     1.08558700    -0.39848100
H       -2.16879600    -1.32754000    -0.51571600
H       -0.00001100    -1.18194200    -0.52080800
H       -1.03205900    -1.45439400     0.87315800
units angstrom
}
""")

HTBH_CH3 = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.00000000    -0.00000000
H        0.00000000     0.00000000     1.07731727
H       -0.00000000     0.93298412    -0.53865863
H        0.00000000    -0.93298412    -0.53865863
units angstrom
}
""")

HTBH_CH3H2ts = input.process_input("""
molecule dimer {
0 2
C        0.00000000     0.26481300     0.00000000
H        1.05342900     0.51666800     0.00000000
H       -0.52662700     0.51702500     0.91225000
H       -0.52662700     0.51702500    -0.91225000
H       -0.00026000    -1.11777100     0.00000000
H        0.00008400    -2.02182500     0.00000000
units angstrom
}
""")

HTBH_CH3NH2ts = input.process_input("""
molecule dimer {
0 3
C       -1.19957700    -0.01112600    -0.00003000
N        1.40071500     0.12986200     0.00001500
H       -1.42666000    -0.51293200     0.93305700
H       -1.41990700    -0.59138200    -0.88814300
H       -1.52023700     1.02280600    -0.04578300
H        0.18892600     0.12689600     0.00100100
H        1.57033800    -0.88766700    -0.00005300
units angstrom
}
""")

HTBH_CH4 = input.process_input("""
molecule dimer {
0 1
C        0.00000000     0.00000000     0.00000000
H        0.00000000     1.08744517     0.00000000
H       -0.51262657    -0.36248173     0.88789526
H       -0.51262657    -0.36248173    -0.88789526
H        1.02525314    -0.36248173     0.00000000
units angstrom
}
""")

HTBH_Cl = input.process_input("""
molecule dimer {
0 2
Cl       0.00000000     0.00000000     0.00000000
units angstrom
}
""")

HTBH_F = input.process_input("""
molecule dimer {
0 2
F        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

HTBH_FH2ts = input.process_input("""
molecule dimer {
0 2
H        0.14656800    -1.12839000     0.00000000
F        0.00000000     0.33042200     0.00000000
H       -0.14656800    -1.84541000     0.00000000
units angstrom
}
""")

HTBH_H = input.process_input("""
molecule dimer {
0 2
H        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

HTBH_H2 = input.process_input("""
molecule dimer {
0 1
H        0.00000000     0.00000000     0.00000000
H        0.74187646     0.00000000     0.00000000
units angstrom
}
""")

HTBH_H2O = input.process_input("""
molecule dimer {
0 1
O        0.00000000     0.00000000    -0.06555155
H        0.00000000    -0.75670946     0.52017534
H        0.00000000     0.75670946     0.52017534
units angstrom
}
""")

HTBH_H2S = input.process_input("""
molecule dimer {
0 1
S        0.00000000     0.00000000     0.10251900
H        0.00000000     0.96624900    -0.82015400
H        0.00000000    -0.96624900    -0.82015400
units angstrom
}
""")

HTBH_HCl = input.process_input("""
molecule dimer {
0 1
Cl       0.00000000     0.00000000     0.00000000
H        1.27444789     0.00000000     0.00000000
units angstrom
}
""")

HTBH_HClCH3ts = input.process_input("""
molecule dimer {
0 2
C        0.24411700     0.59991600     1.70242300
H       -0.67559700     0.27848200     2.17293900
H        0.35191000     1.66378600     1.53767200
H        1.14068600     0.06578700     1.98782200
H        0.05716300     0.13997300     0.39711200
Cl      -0.13758000    -0.33809000    -0.95941600
units angstrom
}
""")

HTBH_HHClts = input.process_input("""
molecule dimer {
0 2
H        0.00048000    -1.34062700     0.00000000
Cl       0.00000000     0.20325200     0.00000000
H       -0.00048000    -2.11465900     0.00000000
units angstrom
}
""")

HTBH_HF = input.process_input("""
molecule dimer {
0 1
F        0.00000000     0.00000000     0.00000000
H        0.91538107     0.00000000     0.00000000
units angstrom
}
""")

HTBH_HH2Sts = input.process_input("""
molecule dimer {
0 2
H        1.26209700    -0.22009700     0.00000000
S        0.00000000     0.22315300     0.00000000
H       -0.50057600    -1.11544500     0.00000000
H       -0.76152100    -2.23491300     0.00000000
units angstrom
}
""")

HTBH_HH2ts = input.process_input("""
molecule dimer {
0 2
H        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     0.92947400
H        0.00000000     0.00000000    -0.92947400
units angstrom
}
""")

HTBH_NH = input.process_input("""
molecule dimer {
0 3
N        0.00000000     0.00000000     0.00000000
H        1.03673136     0.00000000     0.00000000
units angstrom
}
""")

HTBH_HPH3ts = input.process_input("""
molecule dimer {
0 2
P        0.21742900     0.00008800    -0.11124900
H        0.24660900     1.03466800     0.85216400
H        0.26266100    -1.02505800     0.86162300
H       -1.26641800    -0.01095200    -0.15062600
H       -2.50429000     0.00002800     0.10557500
units angstrom
}
""")

HTBH_NH2 = input.process_input("""
molecule dimer {
0 2
N        0.00000000     0.00000000    -0.08007491
H        0.00000000    -0.80231373     0.55629442
H        0.00000000     0.80231373     0.55629442
units angstrom
}
""")

HTBH_NH2C2H5ts = input.process_input("""
molecule dimer {
0 3
C       -1.39498400    -0.44966100     0.00070300
C       -0.43574600     0.71406300     0.00202700
N        1.92757000    -0.37835200     0.00303600
H       -1.20008700    -1.12095100    -0.83568700
H       -1.32209500    -1.02788400     0.92177300
H       -2.42871300    -0.10535200    -0.08933400
H       -0.41768800     1.30848200    -0.90720100
H       -0.44112700     1.32909500     0.89746700
H        0.82850100     0.18059300    -0.02856100
H        2.47259200     0.49807300     0.00391000
units angstrom
}
""")

HTBH_NH2CH4ts = input.process_input("""
molecule dimer {
0 2
C       -1.26075000    -0.00000600     0.01229100
N        1.31325500    -0.00000500    -0.13678200
H       -1.58398700     0.90853800    -0.48474400
H       -1.46367200    -0.00457300     1.07730200
H       -1.58474800    -0.90388000    -0.49270000
H        0.04310800    -0.00006400    -0.15169200
H        1.48045900     0.80557700     0.46775100
H        1.48055700    -0.80552400     0.46780800
units angstrom
}
""")

HTBH_NH3 = input.process_input("""
molecule dimer {
0 1
N        0.00000000     0.00000000     0.11289000
H        0.00000000     0.93802400    -0.26340900
H        0.81235300    -0.46901200    -0.26340900
H       -0.81235300    -0.46901200    -0.26340900
units angstrom
}
""")

HTBH_O = input.process_input("""
molecule dimer {
0 3
O        0.00000000     0.00000000     0.00000000
units angstrom
}
""")

HTBH_OH = input.process_input("""
molecule dimer {
0 2
O        0.00000000     0.00000000     0.00000000
H        0.96889819     0.00000000     0.00000000
units angstrom
}
""")

HTBH_OHC2H6ts = input.process_input("""
molecule dimer {
0 2
C        1.45833400    -0.44636500     0.02547800
C        0.46942300     0.69742200    -0.02749300
O       -1.85303700    -0.31465900    -0.05305500
H        1.30176400    -1.06107900     0.91073700
H        1.36658500    -1.08618900    -0.85111800
H        2.48224500    -0.06687900     0.05715000
H        0.47106900     1.32544300     0.86103700
H        0.53352400     1.30349500    -0.92856000
H       -0.63023200     0.20781600    -0.07846500
H       -2.26720700     0.38832100     0.46575100
units angstrom
}
""")

HTBH_OHCH3ts = input.process_input("""
molecule dimer {
0 3
C        0.00029000    -1.14228900     0.00000000
H       -1.05595700    -1.38473500     0.00000000
H        0.52016700    -1.40738900     0.91244700
H        0.52016700    -1.40738900    -0.91244700
H        0.01156000     0.16009900     0.00000000
O        0.00029000     1.36164300     0.00000000
units angstrom
}
""")

HTBH_OHCH4ts = input.process_input("""
molecule dimer {
0 2
C       -1.21148700     0.00796800     0.00040700
O        1.29396500    -0.10869400     0.00013300
H        0.00947600    -0.11802000     0.00279900
H       -1.52552900    -0.23325000     1.01007000
H       -1.43066500     1.03323300    -0.27808200
H       -1.55271000    -0.71011400    -0.73770200
H        1.41663600     0.84989400    -0.00059100
units angstrom
}
""")

HTBH_OHClts = input.process_input("""
molecule dimer {
0 3
Cl       0.01882000    -0.81730100     0.00000000
H       -0.47048800     0.56948000     0.00000000
O        0.01882000     1.66557900     0.00000000
units angstrom
}
""")

HTBH_OHH2ts = input.process_input("""
molecule dimer {
0 2
O       -0.30106400    -0.10804900    -0.00000800
H       -0.42794500     0.85156900     0.00001600
H        1.01548600    -0.10036700     0.00011900
H        1.82096800     0.11318700    -0.00007300
units angstrom
}
""")

HTBH_OHHts = input.process_input("""
molecule dimer {
0 3
H        0.00000000     0.00000000    -0.86028700
O        0.00000000     0.00000000     0.32902400
H        0.00000000     0.00000000    -1.77190500
units angstrom
}
""")

HTBH_OHNH3ts = input.process_input("""
molecule dimer {
0 2
N       -1.15081600    -0.04393200    -0.10255900
O        1.17918600    -0.09269600    -0.01029000
H       -1.30318500    -0.54763800     0.76657100
H       -1.33891300     0.93580800     0.09185400
H       -0.03068700    -0.15383400    -0.35318400
H        1.29500900     0.81475300     0.29499100
units angstrom
}
""")

HTBH_PH2 = input.process_input("""
molecule dimer {
0 2
P        0.00000000     0.00000000    -0.11565700
H        1.02013000     0.00000000     0.86742700
H       -1.02013000     0.00000000     0.86742700
units angstrom
}
""")

HTBH_PH3 = input.process_input("""
molecule dimer {
0 1
P        0.00000000     0.00000000     0.12641100
H        1.19133900     0.00000000    -0.63205600
H       -0.59566900    -1.03173000    -0.63205600
H       -0.59566900     1.03173000    -0.63205600
units angstrom
}
""")

HTBH_HS = input.process_input("""
molecule dimer {
0 2
S        0.00000000     0.00000000     0.00000000
H        1.34020229     0.00000000     0.00000000
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

