"""
| Database of <description of members and reference energy type>.
| Geometries from <Reference>.
| Reference interaction energies from <Reference>.
| Taken from Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H J. Chem. Phys. 2010, 132, 154104.


- **benchmark**

  - ``'<benchmark_name>'`` <Reference>.
  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.

- **subset**

  - ``'small'`` <members_description>
  - ``'large'`` <members_description>
  - ``'<subset>'`` <members_description>

"""
import re
import qcdb

# <<< HEAVY28 Database Module >>>
dbse = 'HEAVY28'

# <<< Database Members >>>
#HRXN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', ]
HRXN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '22', '24', '25', '26', '27']
HRXN_SM = []
HRXN_LG = []

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, '1'                     )] = ['%s-%s-reagent'      % (dbse, 'BiH3_2'),
                                                               '%s-%s-reagent'      % (dbse, 'BiH3') ]
RXNM['%s-%s'            % (dbse, '1'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '1')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '2'                     )] = ['%s-%s-reagent'      % (dbse, 'BiH3_H2O'),
                                                               '%s-%s-reagent'      % (dbse, 'BiH3'),
                                                               '%s-%s-reagent'      % (dbse, 'H2O') ]
RXNM['%s-%s'            % (dbse, '2'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '2')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '3'                     )] = ['%s-%s-reagent'      % (dbse, 'BiH3_H2S'),
                                                               '%s-%s-reagent'      % (dbse, 'BiH3'),
                                                               '%s-%s-reagent'      % (dbse, 'H2S') ]
RXNM['%s-%s'            % (dbse, '3'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '3')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '4'                     )] = ['%s-%s-reagent'      % (dbse, 'BiH3_HCl'),
                                                               '%s-%s-reagent'      % (dbse, 'BiH3'),
                                                               '%s-%s-reagent'      % (dbse, 'HCl') ]
RXNM['%s-%s'            % (dbse, '4'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '4')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '5'                     )] = ['%s-%s-reagent'      % (dbse, 'BiH3_HBr'),
                                                               '%s-%s-reagent'      % (dbse, 'BiH3'),
                                                               '%s-%s-reagent'      % (dbse, 'HBr') ]
RXNM['%s-%s'            % (dbse, '5'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '5')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '6'                     )] = ['%s-%s-reagent'      % (dbse, 'BiH3_HI'),
                                                               '%s-%s-reagent'      % (dbse, 'BiH3'),
                                                               '%s-%s-reagent'      % (dbse, 'HI') ]
RXNM['%s-%s'            % (dbse, '6'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '6')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '7'                     )] = ['%s-%s-reagent'      % (dbse, 'BiH3_NH3'),
                                                               '%s-%s-reagent'      % (dbse, 'BiH3'),
                                                               '%s-%s-reagent'      % (dbse, 'NH3') ]
RXNM['%s-%s'            % (dbse, '7'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '7')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '8'                     )] = ['%s-%s-reagent'      % (dbse, 'PbH4_2'),
                                                               '%s-%s-reagent'      % (dbse, 'PbH4') ]
RXNM['%s-%s'            % (dbse, '8'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '8')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '9'                     )] = ['%s-%s-reagent'      % (dbse, 'PbH4_BiH3'),
                                                               '%s-%s-reagent'      % (dbse, 'PbH4'),
                                                               '%s-%s-reagent'      % (dbse, 'BiH3') ]
RXNM['%s-%s'            % (dbse, '9'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '9')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '10'                    )] = ['%s-%s-reagent'      % (dbse, 'PbH4_H2O'),
                                                               '%s-%s-reagent'      % (dbse, 'PbH4'),
                                                               '%s-%s-reagent'      % (dbse, 'H2O') ]
RXNM['%s-%s'            % (dbse, '10'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '10')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '11'                    )] = ['%s-%s-reagent'      % (dbse, 'PbH4_HCl'),
                                                               '%s-%s-reagent'      % (dbse, 'PbH4'),
                                                               '%s-%s-reagent'      % (dbse, 'HCl') ]
RXNM['%s-%s'            % (dbse, '11'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '11')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '12'                    )] = ['%s-%s-reagent'      % (dbse, 'PbH4_HBr'),
                                                               '%s-%s-reagent'      % (dbse, 'PbH4'),
                                                               '%s-%s-reagent'      % (dbse, 'HBr') ]
RXNM['%s-%s'            % (dbse, '12'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '12')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '13'                    )] = ['%s-%s-reagent'      % (dbse, 'PbH4_HI'),
                                                               '%s-%s-reagent'      % (dbse, 'PbH4'),
                                                               '%s-%s-reagent'      % (dbse, 'HI') ]
RXNM['%s-%s'            % (dbse, '13'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '13')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '14'                    )] = ['%s-%s-reagent'      % (dbse, 'PbH4_TeH2'),
                                                               '%s-%s-reagent'      % (dbse, 'PbH4'),
                                                               '%s-%s-reagent'      % (dbse, 'TeH2') ]
RXNM['%s-%s'            % (dbse, '14'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '14')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '15'                    )] = ['%s-%s-reagent'      % (dbse, 'SbH3_2'),
                                                               '%s-%s-reagent'      % (dbse, 'SbH3') ]
RXNM['%s-%s'            % (dbse, '15'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '15')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '16'                    )] = ['%s-%s-reagent'      % (dbse, 'SbH3_H2O'),
                                                               '%s-%s-reagent'      % (dbse, 'SbH3'),
                                                               '%s-%s-reagent'      % (dbse, 'H2O') ]
RXNM['%s-%s'            % (dbse, '16'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '16')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '17'                    )] = ['%s-%s-reagent'      % (dbse, 'SbH3_H2S'),
                                                               '%s-%s-reagent'      % (dbse, 'SbH3'),
                                                               '%s-%s-reagent'      % (dbse, 'H2S') ]
RXNM['%s-%s'            % (dbse, '17'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '17')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '18'                    )] = ['%s-%s-reagent'      % (dbse, 'SbH3_HCl'),
                                                               '%s-%s-reagent'      % (dbse, 'SbH3'),
                                                               '%s-%s-reagent'      % (dbse, 'HCl') ]
RXNM['%s-%s'            % (dbse, '18'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '18')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '19'                    )] = ['%s-%s-reagent'      % (dbse, 'SbH3_HBr'),
                                                               '%s-%s-reagent'      % (dbse, 'SbH3'),
                                                               '%s-%s-reagent'      % (dbse, 'HBr') ]
RXNM['%s-%s'            % (dbse, '19'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '19')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '20'                    )] = ['%s-%s-reagent'      % (dbse, 'SbH3_HI'),
                                                               '%s-%s-reagent'      % (dbse, 'SbH3'),
                                                               '%s-%s-reagent'      % (dbse, 'HI') ]
RXNM['%s-%s'            % (dbse, '20'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '20')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '21'                    )] = ['%s-%s-reagent'      % (dbse, 'SbH3_NH3'),  # proven BAD
                                                               '%s-%s-reagent'      % (dbse, 'SbH3'),
                                                               '%s-%s-reagent'      % (dbse, 'NH3') ]
RXNM['%s-%s'            % (dbse, '21'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '21')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '22'                    )] = ['%s-%s-reagent'      % (dbse, 'TeH2_2'),
                                                               '%s-%s-reagent'      % (dbse, 'TeH2') ]
RXNM['%s-%s'            % (dbse, '22'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '22')], [-1, +2]))

ACTV['%s-%s'            % (dbse, '23'                    )] = ['%s-%s-reagent'      % (dbse, 'TeH2_H2O'),  # proven BAD
                                                               '%s-%s-reagent'      % (dbse, 'TeH2'),
                                                               '%s-%s-reagent'      % (dbse, 'H2O') ]
RXNM['%s-%s'            % (dbse, '23'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '23')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '24'                    )] = ['%s-%s-reagent'      % (dbse, 'TeH2_H2S'),
                                                               '%s-%s-reagent'      % (dbse, 'TeH2'),
                                                               '%s-%s-reagent'      % (dbse, 'H2S') ]
RXNM['%s-%s'            % (dbse, '24'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '24')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '25'                    )] = ['%s-%s-reagent'      % (dbse, 'TeH2_HCl'),
                                                               '%s-%s-reagent'      % (dbse, 'TeH2'),
                                                               '%s-%s-reagent'      % (dbse, 'HCl') ]
RXNM['%s-%s'            % (dbse, '25'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '25')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '26'                    )] = ['%s-%s-reagent'      % (dbse, 'TeH2_HBr'),
                                                               '%s-%s-reagent'      % (dbse, 'TeH2'),
                                                               '%s-%s-reagent'      % (dbse, 'HBr') ]
RXNM['%s-%s'            % (dbse, '26'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '26')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '27'                    )] = ['%s-%s-reagent'      % (dbse, 'TeH2_HI'),
                                                               '%s-%s-reagent'      % (dbse, 'TeH2'),
                                                               '%s-%s-reagent'      % (dbse, 'HI') ]
RXNM['%s-%s'            % (dbse, '27'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '27')], [-1, +1, +1]))

ACTV['%s-%s'            % (dbse, '28'                    )] = ['%s-%s-reagent'      % (dbse, 'TeH2_NH3'),  # proven BAD
                                                               '%s-%s-reagent'      % (dbse, 'TeH2'),
                                                               '%s-%s-reagent'      % (dbse, 'NH3') ]
RXNM['%s-%s'            % (dbse, '28'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '28')], [-1, +1, +1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
# Original publication
# Current revision
BIND_HEAVY280 = {}
BIND_HEAVY280['%s-%s'            % (dbse, '1'                     )] =    1.29
BIND_HEAVY280['%s-%s'            % (dbse, '2'                     )] =    2.42
BIND_HEAVY280['%s-%s'            % (dbse, '3'                     )] =    1.40
BIND_HEAVY280['%s-%s'            % (dbse, '4'                     )] =    0.85
BIND_HEAVY280['%s-%s'            % (dbse, '5'                     )] =    1.16
BIND_HEAVY280['%s-%s'            % (dbse, '6'                     )] =    1.42
BIND_HEAVY280['%s-%s'            % (dbse, '7'                     )] =    0.69
BIND_HEAVY280['%s-%s'            % (dbse, '8'                     )] =    1.32
BIND_HEAVY280['%s-%s'            % (dbse, '9'                     )] =    0.68
BIND_HEAVY280['%s-%s'            % (dbse, '10'                    )] =    0.44
BIND_HEAVY280['%s-%s'            % (dbse, '11'                    )] =    0.80
BIND_HEAVY280['%s-%s'            % (dbse, '12'                    )] =    1.04
BIND_HEAVY280['%s-%s'            % (dbse, '13'                    )] =    1.29
BIND_HEAVY280['%s-%s'            % (dbse, '14'                    )] =    0.70
BIND_HEAVY280['%s-%s'            % (dbse, '15'                    )] =    1.30
BIND_HEAVY280['%s-%s'            % (dbse, '16'                    )] =    1.70
BIND_HEAVY280['%s-%s'            % (dbse, '17'                    )] =    1.14
BIND_HEAVY280['%s-%s'            % (dbse, '18'                    )] =    2.20
BIND_HEAVY280['%s-%s'            % (dbse, '19'                    )] =    2.07
BIND_HEAVY280['%s-%s'            % (dbse, '20'                    )] =    1.64
BIND_HEAVY280['%s-%s'            % (dbse, '21'                    )] =    2.80
BIND_HEAVY280['%s-%s'            % (dbse, '22'                    )] =    0.58
BIND_HEAVY280['%s-%s'            % (dbse, '23'                    )] =    0.68
BIND_HEAVY280['%s-%s'            % (dbse, '24'                    )] =    0.50
BIND_HEAVY280['%s-%s'            % (dbse, '25'                    )] =    1.24
BIND_HEAVY280['%s-%s'            % (dbse, '26'                    )] =    1.24
BIND_HEAVY280['%s-%s'            % (dbse, '27'                    )] =    0.84
BIND_HEAVY280['%s-%s'            % (dbse, '28'                    )] =    3.29
# Set default
BIND = BIND_HEAVY280
# Reference information
BINDINFO_HEAVY280 = {}
for rxn in HRXN:
    BINDINFO_HEAVY280['%s-%s' % (dbse, rxn)] = {'citation': 'dftd3', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'atqz???'}  # CC bas?, PP?, core?

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, '1'                     )] = """Bismuthine Dimer"""
TAGL['%s-%s'            % (dbse, '2'                     )] = """Bismuthine-Water"""
TAGL['%s-%s'            % (dbse, '3'                     )] = """Bismuthine-Hydrogen Sulfide"""
TAGL['%s-%s'            % (dbse, '4'                     )] = """Bismuthine-Hydrogen Chloride"""
TAGL['%s-%s'            % (dbse, '5'                     )] = """Bismuthine-Hydrogen Bromide"""
TAGL['%s-%s'            % (dbse, '6'                     )] = """Bismuthine-Hydrogen Iodide"""
TAGL['%s-%s'            % (dbse, '7'                     )] = """Bismuthine-Ammonia"""
TAGL['%s-%s'            % (dbse, '8'                     )] = """Plumbane Dimer"""
TAGL['%s-%s'            % (dbse, '9'                     )] = """Plumbane-Bismuthine"""
TAGL['%s-%s'            % (dbse, '10'                    )] = """Plumbane-Water"""
TAGL['%s-%s'            % (dbse, '11'                    )] = """Plumbane-Hydrogen Chloride"""
TAGL['%s-%s'            % (dbse, '12'                    )] = """Plumbane-Hydrogen Bromide"""
TAGL['%s-%s'            % (dbse, '13'                    )] = """Plumbane-Hydrogen Iodide"""
TAGL['%s-%s'            % (dbse, '14'                    )] = """Plumbane-Tellane"""
TAGL['%s-%s'            % (dbse, '15'                    )] = """Stilbine Dimer"""
TAGL['%s-%s'            % (dbse, '16'                    )] = """Stilbine-Water"""
TAGL['%s-%s'            % (dbse, '17'                    )] = """Stilbine-Hydrogen Sulfide"""
TAGL['%s-%s'            % (dbse, '18'                    )] = """Stilbine-Hydrogen Chloride"""
TAGL['%s-%s'            % (dbse, '19'                    )] = """Stilbine-Hydrogen Bromide"""
TAGL['%s-%s'            % (dbse, '20'                    )] = """Stilbine-Hydrogen Iodide"""
TAGL['%s-%s'            % (dbse, '21'                    )] = """Stilbine-Ammonia"""
TAGL['%s-%s'            % (dbse, '22'                    )] = """Tellane Dimer"""
TAGL['%s-%s'            % (dbse, '23'                    )] = """Tellane-Water"""
TAGL['%s-%s'            % (dbse, '24'                    )] = """Tellane-Hydrogen Sulfide"""
TAGL['%s-%s'            % (dbse, '25'                    )] = """Tellane-Hydrogen Chloride"""
TAGL['%s-%s'            % (dbse, '26'                    )] = """Tellane-Hydrogen Bromide"""
TAGL['%s-%s'            % (dbse, '27'                    )] = """Tellane-Hydrogen Iodide"""
TAGL['%s-%s'            % (dbse, '28'                    )] = """Tellane-Ammonia"""
TAGL['%s-%s-reagent'    % (dbse, 'BiH3'                  )] = """Bismuthine"""
TAGL['%s-%s-reagent'    % (dbse, 'BiH3_2'                )] = """Bismuthine dimer"""
TAGL['%s-%s-reagent'    % (dbse, 'BiH3_H2O'              )] = """Bismuthine-Water"""
TAGL['%s-%s-reagent'    % (dbse, 'BiH3_H2S'              )] = """Bismuthine-Hydrogen Sulfide"""
TAGL['%s-%s-reagent'    % (dbse, 'BiH3_HBr'              )] = """Bismuthine-Hydrogen Bromide"""
TAGL['%s-%s-reagent'    % (dbse, 'BiH3_HCl'              )] = """Bismuthine-Hydrogen Chloride"""
TAGL['%s-%s-reagent'    % (dbse, 'BiH3_HI'               )] = """Bismuthine-Hydrogen Iodide"""
TAGL['%s-%s-reagent'    % (dbse, 'BiH3_NH3'              )] = """Bismuthine-Ammonia"""
TAGL['%s-%s-reagent'    % (dbse, 'H2O'                   )] = """Water"""
TAGL['%s-%s-reagent'    % (dbse, 'H2S'                   )] = """Hydrogen Sulfide"""
TAGL['%s-%s-reagent'    % (dbse, 'HBr'                   )] = """Hydrogen Bromide"""
TAGL['%s-%s-reagent'    % (dbse, 'HCl'                   )] = """Hydrogen Chloride"""
TAGL['%s-%s-reagent'    % (dbse, 'HI'                    )] = """Hydrogen Iodide"""
TAGL['%s-%s-reagent'    % (dbse, 'NH3'                   )] = """Ammonia"""
TAGL['%s-%s-reagent'    % (dbse, 'PbH4'                  )] = """Plumbane"""
TAGL['%s-%s-reagent'    % (dbse, 'PbH4_2'                )] = """Plumbane dimer"""
TAGL['%s-%s-reagent'    % (dbse, 'PbH4_BiH3'             )] = """Plumbane-Bismuthine"""
TAGL['%s-%s-reagent'    % (dbse, 'PbH4_H2O'              )] = """Plumbane-Water"""
TAGL['%s-%s-reagent'    % (dbse, 'PbH4_HBr'              )] = """Plumbane-Hydrogen Bromide"""
TAGL['%s-%s-reagent'    % (dbse, 'PbH4_HCl'              )] = """Plumbane-Hydrogen Chloride"""
TAGL['%s-%s-reagent'    % (dbse, 'PbH4_HI'               )] = """Plumbane-Hydrogen Iodide"""
TAGL['%s-%s-reagent'    % (dbse, 'PbH4_TeH2'             )] = """Plumbane-Tellane"""
TAGL['%s-%s-reagent'    % (dbse, 'SbH3'                  )] = """Stilbine"""
TAGL['%s-%s-reagent'    % (dbse, 'SbH3_2'                )] = """Stilbine dimer"""
TAGL['%s-%s-reagent'    % (dbse, 'SbH3_H2O'              )] = """Stilbine-Water"""
TAGL['%s-%s-reagent'    % (dbse, 'SbH3_H2S'              )] = """Stilbine-Hydrogen Sulfide"""
TAGL['%s-%s-reagent'    % (dbse, 'SbH3_HBr'              )] = """Stilbine-Hydrogen Bromide"""
TAGL['%s-%s-reagent'    % (dbse, 'SbH3_HCl'              )] = """Stilbine-Hydrogen Chloride"""
TAGL['%s-%s-reagent'    % (dbse, 'SbH3_HI'               )] = """Stilbine-Hydrogen Iodide"""
TAGL['%s-%s-reagent'    % (dbse, 'SbH3_NH3'              )] = """Stilbine-Ammonia"""
TAGL['%s-%s-reagent'    % (dbse, 'TeH2'                  )] = """Tellane"""
TAGL['%s-%s-reagent'    % (dbse, 'TeH2_2'                )] = """Tellane dimer"""
TAGL['%s-%s-reagent'    % (dbse, 'TeH2_H2O'              )] = """Tellane-Water"""
TAGL['%s-%s-reagent'    % (dbse, 'TeH2_H2S'              )] = """Tellane-Hydrogen Sulfide"""
TAGL['%s-%s-reagent'    % (dbse, 'TeH2_HBr'              )] = """Tellane-Hydrogen Bromide"""
TAGL['%s-%s-reagent'    % (dbse, 'TeH2_HCl'              )] = """Tellane-Hydrogen Chloride"""
TAGL['%s-%s-reagent'    % (dbse, 'TeH2_HI'               )] = """Tellane-Hydrogen Iodide"""
TAGL['%s-%s-reagent'    % (dbse, 'TeH2_NH3'              )] = """Tellane-Ammonia"""

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'BiH3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BI               0.000000000000     0.000000000000     1.446025999101
    H               -1.357253687789    -2.350832346011    -0.482033890343
    H               -1.357253687789     2.350832346011    -0.482033890343
    H                2.714507375578     0.000000000000    -0.482033890343

""")

GEOS['%s-%s-%s' % (dbse, 'BiH3_2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BI               4.211509508066     0.659360628192     0.000000000000
    BI              -4.211509508066    -0.659360628192     0.000000000000
    H                6.112864634103    -0.746223135362     2.349061510047
    H                2.245098220765    -2.034108244024     0.000000000000
    H                6.112864634103    -0.746223135362    -2.349061510047
    H               -6.112864634103     0.746223135362    -2.349061510047
    H               -2.245098220765     2.034108244024     0.000000000000
    H               -6.112864634103     0.746223135362     2.349061510047

""")

GEOS['%s-%s-%s' % (dbse, 'BiH3_H2O', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BI               1.666237975583    -2.293355317636    -0.545084655425
    H                1.252981568219    -5.572442910008     0.012637473816
    H                0.220361244348    -1.644759073188     2.387366947801
    H               -1.312802268379    -2.206096061407    -2.033135861768
    O                0.054761063022     3.580482763555     0.020433757349
    H               -1.249901556760     3.298292003259     1.246911009471
    H               -0.631638026033     4.837878595424    -1.089128671243

""")

GEOS['%s-%s-%s' % (dbse, 'BiH3_H2S', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BI               1.534851422370     2.987385671895     0.000000000000
    H                0.465106527035     5.097203081183     2.348791379542
    H                0.465106527035     5.097203081183    -2.348791379542
    H               -1.450928128371     1.503164481318     0.000000000000
    H                0.242434374414    -6.716544432687     0.000000000000
    S                0.607502154980    -4.223694374074     0.000000000000
    H               -1.864072877463    -3.744717508818     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'BiH3_HBr', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BI               0.000000000000     0.000000000000    -1.950033917573
    H               -1.358187536950    -2.352449820205    -3.874329732623
    H               -1.358187536950     2.352449820205    -3.874329732623
    H                2.716375073901     0.000000000000    -3.874329732623
    BR               0.000000000000     0.000000000000     5.449328838825
    H                0.000000000000     0.000000000000     8.124731202173

""")

GEOS['%s-%s-%s' % (dbse, 'BiH3_HCl', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BI               0.000000000000     0.000000000000    -1.889594671890
    H               -1.357165532631    -2.350679656799    -3.819476086560
    H               -1.357165532631     2.350679656799    -3.819476086560
    H                2.714331065263     0.000000000000    -3.819476086560
    CL               0.000000000000     0.000000000000     5.471867442690
    H                0.000000000000     0.000000000000     7.877057583277

""")

GEOS['%s-%s-%s' % (dbse, 'BiH3_HI', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BI               0.000000000000     0.000000000000    -1.996592564632
    H               -1.360047695160    -2.355671708734    -3.912278832850
    H               -1.360047695160     2.355671708734    -3.912278832850
    H                2.720095390320     0.000000000000    -3.912278832850
    I                0.000000000000     0.000000000000     5.361059497320
    H                0.000000000000     0.000000000000     8.372809121123

""")

GEOS['%s-%s-%s' % (dbse, 'BiH3_NH3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BI               0.000000000000     0.000000000000    -3.379546144161
    H               -1.353013030975    -2.343487312951    -5.337461057795
    H               -1.353013030975     2.343487312951    -5.337461057795
    H                2.706026061949     0.000000000000    -5.337461057795
    N                0.000000000000     0.000000000000     4.309083302812
    H               -1.768181081307     0.000000000000     5.027915930040
    H                0.884090540653    -1.531289734903     5.027915930040
    H                0.884090540653     1.531289734903     5.027915930040

""")

GEOS['%s-%s-%s' % (dbse, 'H2O', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    O                0.000000000000     0.000000000000    -0.741141166511
    H                1.428822067888     0.000000000000     0.370570583255
    H               -1.428822067888    -0.000000000000     0.370570583255

""")

GEOS['%s-%s-%s' % (dbse, 'H2S', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    H               -1.814333975616    -0.000000000000    -0.581751760461
    S                0.000000000000     0.000000000000     1.163503520922
    H                1.814333975616     0.000000000000    -0.581751760461

""")

GEOS['%s-%s-%s' % (dbse, 'HBr', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    BR               0.000000000000     0.000000000000    -1.335750414149
    H                0.000000000000     0.000000000000     1.335750414149

""")

GEOS['%s-%s-%s' % (dbse, 'HCl', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    H                0.000000000000     0.000000000000     1.201321041337
    CL               0.000000000000     0.000000000000    -1.201321041337

""")

GEOS['%s-%s-%s' % (dbse, 'HI', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    I                0.000000000000     0.000000000000    -1.502195686012
    H                0.000000000000     0.000000000000     1.502195686012

""")

GEOS['%s-%s-%s' % (dbse, 'NH3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    N                0.000000000000     0.000000000000    -0.535255509545
    H               -0.884611497284     1.532192058255     0.178422455880
    H               -0.884611497284    -1.532192058255     0.178422455880
    H                1.769222994567     0.000000000000     0.178422455880

""")

GEOS['%s-%s-%s' % (dbse, 'PbH4', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    PB               0.000000000000     0.000000000000     0.000000000000
    H               -1.879233755259     1.879233755259     1.879233755259
    H                1.879233755259    -1.879233755259     1.879233755259
    H               -1.879233755259    -1.879233755259    -1.879233755259
    H                1.879233755259     1.879233755259    -1.879233755259

""")

GEOS['%s-%s-%s' % (dbse, 'PbH4_2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    PB               0.000000000000     0.000000000000    -4.219190605506
    H               -1.875810079479     1.875810079479    -2.331468390907
    H                1.875810079479    -1.875810079479    -2.331468390907
    H               -1.879312656802    -1.879312656802    -6.100581621626
    H                1.879312656802     1.879312656802    -6.100581621626
    PB               0.000000000000     0.000000000000     4.219190605506
    H               -1.879312656802     1.879312656802     6.100581621626
    H                1.879312656802    -1.879312656802     6.100581621626
    H               -1.875810079479    -1.875810079479     2.331468390907
    H                1.875810079479     1.875810079479     2.331468390907

""")

GEOS['%s-%s-%s' % (dbse, 'PbH4_BiH3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    PB               0.000000000000     0.000000000000    -4.980117851559
    H                1.534819974993    -2.658386177160    -6.062927849272
    H                1.534819974993     2.658386177160    -6.062927849272
    H               -3.069639949987     0.000000000000    -6.062927849272
    H                0.000000000000     0.000000000000    -1.721649436303
    BI               0.000000000000     0.000000000000     4.773659391931
    H                2.714004562416     0.000000000000     6.706580266371
    H               -1.357002281208     2.350396897040     6.706580266371
    H               -1.357002281208    -2.350396897040     6.706580266371

""")

GEOS['%s-%s-%s' % (dbse, 'PbH4_H2O', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    PB               2.749179618555     2.095565262318     0.000000000000
    H                2.794709383586     3.998545849974     2.646343638075
    H                5.438025061662     0.252969851540     0.000000000000
    H                0.100331323009     0.221871897192     0.000000000000
    H                2.794709383586     3.998545849974    -2.646343638075
    O               -4.018752895484    -3.096172910804     0.000000000000
    H               -5.750537926904    -2.566263315844     0.000000000000
    H               -4.107663948010    -4.905062484350     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'PbH4_HBr', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    PB               0.000000000000     0.000000000000    -2.616263811339
    H                0.000000000000    -2.650373840410    -0.724677721044
    H               -0.000000000000     2.650373840410    -0.724677721044
    H                2.658211312748     0.000000000000    -4.493796034380
    H               -2.658211312748    -0.000000000000    -4.493796034380
    H                0.000000000000     0.000000000000     7.863364534653
    BR               0.000000000000     0.000000000000     5.189846787535

""")

GEOS['%s-%s-%s' % (dbse, 'PbH4_HCl', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    PB               0.000000000000     0.000000000000    -2.535810103461
    H                0.000000000000    -2.654409912807    -0.652568220482
    H               -0.000000000000     2.654409912807    -0.652568220482
    H                2.655632061444     0.000000000000    -4.419186598615
    H               -2.655632061444    -0.000000000000    -4.419186598615
    H                0.000000000000     0.000000000000     7.541700128930
    CL               0.000000000000     0.000000000000     5.137619612725

""")

GEOS['%s-%s-%s' % (dbse, 'PbH4_HI', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    PB               0.000000000000     0.000000000000    -2.670208714536
    H                0.000000000000    -2.644139186255    -0.765594961507
    H               -0.000000000000     2.644139186255    -0.765594961507
    H                2.661970069377     0.000000000000    -4.539454529413
    H               -2.661970069377    -0.000000000000    -4.539454529413
    H                0.000000000000     0.000000000000     8.144199432670
    I                0.000000000000     0.000000000000     5.136108263706

""")

GEOS['%s-%s-%s' % (dbse, 'PbH4_TeH2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    PB              -0.005275109644    -4.065699399460     0.000000000000
    H               -0.105302851335    -0.810664635179     0.000000000000
    H               -1.504362443110    -5.197069098290     2.658573403788
    H               -1.504362443110    -5.197069098290    -2.658573403788
    H                3.097632876507    -5.050179129079     0.000000000000
    TE               0.020913656055     5.316180464234     0.000000000000
    H                2.181984086953     7.522856552549     0.000000000000
    H               -2.181227772316     7.481644343514     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'SbH3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    SB               0.000000000000     0.000000000000    -1.340998601775
    H               -1.313777404557    -2.275529214529     0.447009232630
    H               -1.313777404557     2.275529214529     0.447009232630
    H                2.627554809114     0.000000000000     0.447009232630

""")

GEOS['%s-%s-%s' % (dbse, 'SbH3_2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    SB              -1.364558987878     3.561776302126     0.000000000000
    H               -0.417272714330     5.571920881072     2.274948848935
    H               -0.417272714330     5.571920881072    -2.274948848935
    H                1.474232723624     2.118627711079     0.000000000000
    SB               1.492837263549    -3.923320526331     0.000000000000
    H               -0.537399664790    -3.012390498308     2.271690293004
    H                0.306833758945    -6.876144252402     0.000000000000
    H               -0.537399664790    -3.012390498308    -2.271690293004

""")

GEOS['%s-%s-%s' % (dbse, 'SbH3_H2O', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    SB              -0.082097780411    -2.209640942089     0.000000000000
    H                1.350328455332    -3.871020595738     2.289164968418
    H                1.350328455332    -3.871020595738    -2.289164968418
    H               -2.601169712214    -4.138079211543     0.000000000000
    O                0.676727697321     4.986950451609     0.000000000000
    H                0.209355789921     3.232006199821     0.000000000000
    H               -0.903472905281     5.870804693677     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'SbH3_H2S', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    SB              -0.075036696794    -2.541861903274     0.000000000000
    H                1.319574823344    -4.251669553228     2.282755636653
    H                1.319574823344    -4.251669553228    -2.282755636653
    H               -2.630381094063    -4.425512151949     0.000000000000
    H                1.512829073137     6.346504431386     0.000000000000
    S               -0.945732664278     5.803733647965     0.000000000000
    H               -0.500828264690     3.320475082329     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'SbH3_HBr', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    SB               0.000000000000     0.000000000000     1.350604250123
    H                0.000000000000     0.000000000000    -3.979770293287
    BR               0.000000000000     0.000000000000    -6.668577435731
    H               -1.322481372194    -2.290604928704     3.099122262118
    H               -1.322481372194     2.290604928704     3.099122262118
    H                2.644962744389     0.000000000000     3.099122262118

""")

GEOS['%s-%s-%s' % (dbse, 'SbH3_HCl', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    SB               0.000000000000     0.000000000000     1.287353837352
    H                0.000000000000     0.000000000000    -3.980905160866
    CL               0.000000000000     0.000000000000    -6.400307136296
    H               -1.323331028084    -2.292076575874     3.031110566508
    H               -1.323331028084     2.292076575874     3.031110566508
    H                2.646662056168     0.000000000000     3.031110566508

""")

GEOS['%s-%s-%s' % (dbse, 'SbH3_HI', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    SB               0.000000000000     0.000000000000     1.425479655830
    H                0.000000000000     0.000000000000    -3.974965525955
    I                0.000000000000     0.000000000000    -6.996409930809
    H               -1.321095958472    -2.288205321748     3.181566920154
    H               -1.321095958472     2.288205321748     3.181566920154
    H                2.642191916945     0.000000000000     3.181566920154

""")

#GEOS['%s-%s-%s' % (dbse, 'SbH3_NH3', 'reagent')] = qcdb.Molecule("""
#    units Bohr
#    no_com
#    no_reorient
#    0 1
#    SB               1.506573249011    -2.895042300647     0.000000000000
#    H               -0.653704873916    -2.493414313714     2.291629220690
#    H                0.998670757836    -6.057884036610     0.000000000000
#    H               -0.653704873916    -2.493414313714    -2.291629220690
#    N               -0.114977594517     2.982067358652     0.000000000000
#    H               -1.192770701027     3.348368106439     1.532995411457
#    H                1.302684737556     4.260951393156     0.000000000000
#    H               -1.192770701027     3.348368106439    -1.532995411457
#
#""")

GEOS['%s-%s-%s' % (dbse, 'TeH2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    TE               0.000000000000     0.000000000000    -1.454237107850
    H                2.183324845891     0.000000000000     0.727118553925
    H               -2.183324845891    -0.000000000000     0.727118553925

""")

GEOS['%s-%s-%s' % (dbse, 'TeH2_2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    TE               1.495489962102     5.569327594892     0.000000000000
    H                1.000109972583     2.521152865619     0.000000000000
    H               -1.551984876488     6.059635221257     0.000000000000
    TE              -0.213029679843    -3.268403571538     0.000000000000
    H                1.816584097634    -5.593798014651     0.000000000000
    H               -2.547169475988    -5.287914095579     0.000000000000

""")

#GEOS['%s-%s-%s' % (dbse, 'TeH2_H2O', 'reagent')] = qcdb.Molecule("""
#    units Bohr
#    no_com
#    no_reorient
#    0 1
#    TE               0.000000000000     0.000000000000    -3.061779351671
#    H               -2.168067233331    -0.000000000000    -5.266220756460
#    H                2.168067233331     0.000000000000    -5.266220756460
#    O                0.000000000000     0.000000000000     3.789426010565
#    H               -1.428826352623    -0.000000000000     4.902397427013
#    H                1.428826352623     0.000000000000     4.902397427013
#
#""")

GEOS['%s-%s-%s' % (dbse, 'TeH2_H2S', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    TE               0.000000000000     0.000000000000    -3.732528394571
    H               -2.179271331125    -0.000000000000    -5.920948715430
    H                2.179271331125     0.000000000000    -5.920948715430
    H                1.814617962832     0.000000000000     5.773433065052
    S                0.000000000000     0.000000000000     4.027559695326
    H               -1.814617962832    -0.000000000000     5.773433065052

""")

GEOS['%s-%s-%s' % (dbse, 'TeH2_HBr', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    TE              -2.109038714517    -2.029285951871     0.000000000000
    H               -2.102871026753     1.055673447010     0.000000000000
    H               -5.200158547634    -2.029028348194     0.000000000000
    H                4.565096419790     2.830262568906     0.000000000000
    BR               4.846971869113     0.172378284148     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'TeH2_HCl', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    TE              -2.090461576798    -1.119118367285     0.000000000000
    H               -2.277592820169     1.958495013397     0.000000000000
    H               -5.172897231464    -1.329660082877     0.000000000000
    H                5.069294946982    -0.919234127581     0.000000000000
    CL               4.471656681449     1.409517564346     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'TeH2_HI', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    TE              -3.343383976209    -2.349782560742     0.000000000000
    H               -0.945576528617    -0.401132072242     0.000000000000
    H               -5.293938270020     0.042040620884     0.000000000000
    H                5.463064727802     0.010291770458     0.000000000000
    I                4.119834047045     2.698582241642     0.000000000000

""")

#GEOS['%s-%s-%s' % (dbse, 'TeH2_NH3', 'reagent')] = qcdb.Molecule("""
#    units Bohr
#    no_com
#    no_reorient
#    0 1
#    TE               1.610398111735    -3.087256836699     0.000000000000
#    H                1.680949249606    -6.192827650173     0.000000000000
#    H               -1.467688716130    -3.164649100651     0.000000000000
#    N               -0.330236844539     2.590621920677     0.000000000000
#    H               -1.357340685976     3.075216031208    -1.534948021110
#    H               -1.357340685976     3.075216031208     1.534948021110
#    H                1.221259571279     3.703679604431     0.000000000000
#
#""")

