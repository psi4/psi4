import pytest
from utils import *

import collections

import qcdb

data = collections.defaultdict(dict)

data['C2H2']['ref'] = \
 [[    0.000000000000,     0.000000000000,     0.650000000000],
    [    0.000000000000,     0.000000000000,    -0.650000000000],
    [    0.000000000000,     0.000000000000,     1.750000000000],
    [    0.000000000000,     0.000000000000,    -1.750000000000]]


data['N2']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -0.550000000000],
    [    0.000000000000,     0.000000000000,     0.550000000000]]


data['CN']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -0.753922540200],
    [    0.000000000000,     0.000000000000,     0.646077459800]]


data['HCCCl']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -2.719188122301],
    [    0.000000000000,     0.000000000000,    -1.719188122301],
    [    0.000000000000,     0.000000000000,    -0.619188122301],
    [    0.000000000000,     0.000000000000,     0.880811877699]]

data['CHFClBr']['ref'] = \
 [[    0.251995331546,    -0.251226142375,    -1.228834841943],
    [    0.251995331546,    -0.251226142375,    -0.228834841943],
    [    0.251995331546,    -1.217151968664,     0.029984203160],
    [    1.088511635284,     0.231736770769,     0.029984203160],
    [   -0.584520972191,     0.231736770769,     0.029984203160]]


data['CH2ClBr']['ref'] = \
 [[   -0.588382367674,     0.890373072017,     0.000000000000],
    [   -0.588382367674,    -0.109626927983,     0.000000000000],
    [    0.377543458615,    -0.368445973085,     0.000000000000],
    [   -1.071345280819,    -0.368445973085,     0.836516303738],
    [   -1.071345280819,    -0.368445973085,    -0.836516303738]]


data['HOCl']['ref'] = \
 [[   -1.074855537230,     1.371823577282,     0.000000000000],
    [   -1.074855537230,     0.371823577282,     0.000000000000],
    [    0.522621918106,    -0.209610666372,     0.000000000000]]


data['C4H4Cl2F2']['ref'] = \
 [[    0.432781050498,     1.898774028282,     0.810337938486],
    [   -1.658744642774,     0.805191018766,    -0.984829058337],
    [    1.658744642774,    -0.805191018766,     0.984829058337],
    [   -0.432781050498,    -1.898774028282,    -0.810337938486],
    [   -0.317971784026,     2.532165941971,     2.640915161238],
    [   -1.615729990528,     1.614062700629,    -2.881498569657],
    [    1.615729990528,    -1.614062700629,     2.881498569657],
    [    0.317971784026,    -2.532165941971,    -2.640915161238],
    [   -4.852178875691,     1.024620478757,     0.190249941464],
    [    4.852178875691,    -1.024620478757,    -0.190249941464],
    [   -1.913713787211,    -3.739567959534,     0.258534542158],
    [    1.913713787211,     3.739567959534,    -0.258534542158]]


data['HOOH_dimer']['ref'] = \
 [[    0.991126228500,    -1.797922633300,     0.146518251500],
    [    2.769109309500,    -1.348521864900,    -0.007155768400],
    [    2.517803031100,     1.380837492300,    -0.115405801400],
    [    3.288320045300,     1.830859509500,     1.475770682500],
    [   -3.288320045300,    -1.830859509500,    -1.475770682500],
    [   -2.517803031100,    -1.380837492300,     0.115405801400],
    [   -2.769109309500,     1.348521864900,     0.007155768400],
    [   -0.991126228500,     1.797922633300,    -0.146518251500]]


data['HOOH']['ref'] = \
 [[   -0.657774170895,     0.990250821217,    -0.765560415408],
    [   -0.134283313545,     0.737880743551,     0.048237265941],
    [    0.134283313545,    -0.737880743551,     0.048237265941],
    [    0.657774170895,    -0.990250821217,    -0.765560415408]]


data['NOHOHOH']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -1.219526050935],  # Gh
    [    0.000000000000,     0.000000000000,    -0.219526050935],
    [    0.000000000000,    -1.477211629518,     0.040946215565],
    [    1.279302797929,     0.738605814759,     0.040946215565],
    [   -1.279302797929,     0.738605814759,     0.040946215565],
    [    0.899871989900,    -1.767037453836,     0.366877793280],
    [    1.080363329511,     1.662830730326,     0.366877793280],
    [   -1.980235319411,     0.104206723511,     0.366877793280]]


data['H2O']['ref'] = \
 [[    0.000000000000,     0.816641555162,    -0.512554059234],
    [    0.000000000000,     0.000000000000,     0.064591130803],
    [    0.000000000000,    -0.816641555162,    -0.512554059234]]


data['CH2F2']['ref'] = \
 [[    0.000000000000,     0.000000000000,     1.089095845660],
    [    0.000000000000,    -2.122315581200,    -0.459816147540],
    [   -0.000000000000,     2.122315581200,    -0.459816147540],
    [    1.708413985000,     0.000000000000,     2.184106800160],
    [   -1.708413985000,    -0.000000000000,     2.184106800160]]


data['NH3']['ref'] = \
 [[    0.000000000000,    -1.071293777318,     0.000000000000],  # Gh
    [    0.000000000000,    -0.071293777318,     0.000000000000],
    [   -0.430496198842,     0.330193571336,    -0.745641288860],
    [   -0.430496198842,     0.330193571336,     0.745641288860],
    [    0.860992397685,     0.330193571336,     0.000000000000]]


data['BrF5']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -1.514287735691],
    [    0.000000000000,     0.000000000000,     0.185712264309],
    [    0.000000000000,    -1.700000000000,     0.185712264309],
    [   -1.700000000000,     0.000000000000,     0.185712264309],
    [    1.700000000000,     0.000000000000,     0.185712264309],
    [    0.000000000000,     1.700000000000,     0.185712264309]]


data['N2H2']['ref'] = \
 [[    0.000000000000,     0.700000000000,     0.000000000000],
    [    0.000000000000,    -0.700000000000,     0.000000000000],
    [   -0.642787609687,     1.466044443119,     0.000000000000],
    [    0.642787609687,    -1.466044443119,     0.000000000000]]


data['NOHOHOH']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -1.000000000000],  # Gh
    [    0.000000000000,     0.000000000000,     0.000000000000],
    [    0.000000000000,    -1.500000000000,     0.000000000000],
    [    1.299038105677,     0.750000000000,     0.000000000000],
    [   -1.299038105677,     0.750000000000,     0.000000000000],
    [    0.939692620786,    -1.842020143326,     0.000000000000],
    [    1.125389928010,     1.734807753012,     0.000000000000],
    [   -2.065082548796,     0.107212390313,     0.000000000000]]


data['TFCOT']['ref'] = \
 [[   -1.618188000000,    -0.437140000000,    -0.409373000000],
    [   -1.394411000000,     0.896360000000,    -0.429596000000],
    [   -0.896360000000,    -1.394411000000,     0.429596000000],
    [   -0.437140000000,     1.618188000000,     0.409373000000],
    [    0.437140000000,    -1.618188000000,     0.409373000000],
    [    0.896360000000,     1.394411000000,     0.429596000000],
    [    1.394411000000,    -0.896360000000,    -0.429596000000],
    [    1.618188000000,     0.437140000000,    -0.409373000000],
    [    2.147277000000,    -1.690111000000,    -1.235043000000],
    [    1.690111000000,     2.147277000000,     1.235043000000],
    [   -2.147277000000,     1.690111000000,    -1.235043000000],
    [   -1.690111000000,    -2.147277000000,     1.235043000000],
    [    0.878010000000,    -2.418132000000,     1.029595000000],
    [   -2.418132000000,    -0.878010000000,    -1.029595000000],
    [   -0.878010000000,     2.418132000000,     1.029595000000],
    [    2.418132000000,     0.878010000000,    -1.029595000000]]


data['Li_H2O_4_p']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -1.000000000000],  # Gh
    [    0.000000000000,     0.000000000000,     0.000000000000],
    [    0.000000000000,    -1.000000000000,     0.000000000000],  # Gh
    [    0.000000000000,     0.000000000000,     1.000000000000],  # Gh
    [   -1.497220431853,     0.000000000000,    -1.169756803119],
    [    1.497220431853,     0.000000000000,    -1.169756803119],
    [    0.000000000000,    -1.497220431853,     1.169756803119],
    [    0.000000000000,     1.497220431853,     1.169756803119],
    [   -1.565808146965,     0.498804253130,    -1.977713511362],
    [   -2.264066335924,    -0.551608165437,    -1.068881676192],
    [    1.565808146965,    -0.498804253130,    -1.977713511362],
    [    2.264066335924,     0.551608165437,    -1.068881676192],
    [   -0.498804253130,    -1.565808146965,     1.977713511362],
    [    0.551608165437,    -2.264066335924,     1.068881676192],
    [    0.498804253130,     1.565808146965,     1.977713511362],
    [   -0.551608165437,     2.264066335924,     1.068881676192]]


data['ethylene_cation']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -0.705000000000],
    [    0.000000000000,     0.000000000000,     0.705000000000],
    [    0.353742012310,     0.854008763700,    -1.282611998014],
    [   -0.353742012310,    -0.854008763700,    -1.282611998014],
    [   -0.353742012310,     0.854008763700,     1.282611998014],
    [    0.353742012310,    -0.854008763700,     1.282611998014]]


data['ethane_gauche']['ref'] = \
 [[    1.092020143326,     0.163175911167,    -0.925416578398],
    [    0.750000000000,     0.000000000000,     0.000000000000],
    [   -0.750000000000,     0.000000000000,     0.000000000000],
    [   -1.092020143326,    -0.163175911167,    -0.925416578398],
    [   -1.092020143326,    -0.719846310393,     0.604022773555],
    [   -1.092020143326,     0.883022221559,     0.321393804843],
    [    1.092020143326,     0.719846310393,     0.604022773555],
    [    1.092020143326,    -0.883022221559,     0.321393804843]]


data['triplet_ethylene']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -0.705000000000],
    [    0.000000000000,     0.000000000000,     0.705000000000],
    [    0.000000000000,     0.924372424811,    -1.282611998014],
    [    0.000000000000,    -0.924372424811,    -1.282611998014],
    [   -0.924372424811,     0.000000000000,     1.282611998014],
    [    0.924372424811,     0.000000000000,     1.282611998014]]


data['allene']['ref'] = \
 [[   -2.000000000000,    -0.707106781187,     0.707106781187],
    [   -2.000000000000,     0.707106781187,    -0.707106781187],
    [   -1.500000000000,    -0.000000000000,     0.000000000000],
    [    0.000000000000,    -0.000000000000,     0.000000000000],
    [    1.500000000000,     0.000000000000,     0.000000000000],
    [    2.000000000000,     0.707106781187,     0.707106781187],
    [    2.000000000000,    -0.707106781187,    -0.707106781187]]


data['ethane_staggered']['ref'] = \
 [[   -0.657774170895,     0.990250821217,    -0.813797681349],
    [   -0.134283313545,     0.737880743551,     0.000000000000],
    [    0.134283313545,    -0.737880743551,     0.000000000000],
    [    0.657774170895,    -0.990250821217,    -0.813797681349],
    [   -0.728988008574,    -1.242620898884,     0.000000000000],
    [    0.657774170895,    -0.990250821217,     0.813797681349],
    [   -0.657774170895,     0.990250821217,     0.813797681349],
    [    0.728988008574,     1.242620898884,     0.000000000000]]


data['singlet_ethylene']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -0.705000000000],
    [    0.000000000000,     0.000000000000,     0.705000000000],
    [    0.000000000000,     0.924372424811,    -1.282611998014],
    [    0.000000000000,    -0.924372424811,    -1.282611998014],
    [    0.000000000000,     0.924372424811,     1.282611998014],
    [    0.000000000000,    -0.924372424811,     1.282611998014]]


data['ethane_eclipsed']['ref'] = \
 [[    0.000000000000,     1.092020143326,    -0.939692620786],
    [    0.000000000000,     0.750000000000,     0.000000000000],
    [    0.000000000000,    -0.750000000000,     0.000000000000],
    [    0.000000000000,    -1.092020143326,    -0.939692620786],
    [    0.813797681349,    -1.092020143326,     0.469846310393],
    [   -0.813797681349,    -1.092020143326,     0.469846310393],
    [   -0.813797681349,     1.092020143326,     0.469846310393],
    [    0.813797681349,     1.092020143326,     0.469846310393]]


data['BH4p']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -1.000000000000],  # Gh
    [    0.000000000000,     0.000000000000,     0.000000000000],
    [    0.000000000000,    -1.000000000000,     0.000000000000],
    [    1.000000000000,     0.000000000000,     0.000000000000],
    [    0.000000000000,     1.000000000000,     0.000000000000],
    [   -1.000000000000,     0.000000000000,     0.000000000000]]


data['CH4']['ref'] = \
 [[    0.000000000000,     0.000000000000,     0.000000000000],
    [    0.000000000000,     0.889981273211,    -0.629311793417],
    [    0.889981273211,     0.000000000000,     0.629311793417],
    [   -0.889981273211,     0.000000000000,     0.629311793417],
    [    0.000000000000,    -0.889981273211,    -0.629311793417]]


data['SF6']['ref'] = \
 [[    0.000000000000,     0.000000000000,    -1.800000000000],
    [    0.000000000000,     0.000000000000,     0.000000000000],
    [    0.000000000000,    -1.800000000000,     0.000000000000],
    [    1.800000000000,     0.000000000000,     0.000000000000],
    [    0.000000000000,     1.800000000000,     0.000000000000],
    [   -1.800000000000,     0.000000000000,     0.000000000000],
    [    0.000000000000,     0.000000000000,     1.800000000000]]


data['Ih']['ref'] = \
 [[   -0.000000000000,     1.000000000000,     1.618033988750],
    [    0.000000000000,    -1.000000000000,     1.618033988750],
    [   -0.000000000000,     1.000000000000,    -1.618033988750],
    [    0.000000000000,    -1.000000000000,    -1.618033988750],
    [    1.000000000000,     1.618033988750,     0.000000000000],
    [   -1.000000000000,     1.618033988750,     0.000000000000],
    [    1.000000000000,    -1.618033988750,     0.000000000000],
    [   -1.000000000000,    -1.618033988750,     0.000000000000],
    [    1.618033988750,     0.000000000000,     1.000000000000],
    [    1.618033988750,     0.000000000000,    -1.000000000000],
    [   -1.618033988750,    -0.000000000000,     1.000000000000],
    [   -1.618033988750,    -0.000000000000,    -1.000000000000]]

#! Tests to determine full point group symmetry.  Currently, these only matter
#! for the rotational symmetry number in thermodynamic computations.

data['C2H2']['pg'] = 'D_inf_h'
data['C2H2']['rsn'] = 2
data['C2H2']['mol'] = """
  C 0 0  r1
  C 0 0 -r1
  H 0 0  r2
  H 0 0 -r2
  r1 = 0.65
  r2 = 1.75
"""
data['isoC2H2']['idx'] = 0
data['isoC2H2']['pg'] = 'C_inf_v'
data['isoC2H2']['rsn'] = 1


def test_C2H2():
    sys_tester('C2H2')


def test_isoC2H2():
    sys_tester('isoC2H2')


data['N2']['pg'] = 'D_inf_h'
data['N2']['rsn'] = 2
data['N2']['mol'] = """
  N 0.0 0.0 0.0
  N 0.0 0.0 r
  r = 1.1
"""
data['isoN2']['idx'] = 0
data['isoN2']['pg'] = 'C_inf_v'
data['isoN2']['rsn'] = 1


def test_N2():
    sys_tester('N2')


def test_isoN2():
    sys_tester('isoN2')


data['CN']['pg'] = 'C_inf_v'
data['CN']['rsn'] = 1
data['CN']['mol'] = """
  0 2
  C 0.0 0.0 0.0
  N 0.0 0.0 r
  r = 1.4
"""
data['isoCN']['idx'] = 0
data['isoCN']['pg'] = 'C_inf_v'
data['isoCN']['rsn'] = 1


def test_CN():
    sys_tester('CN')


def test_isoCN():
    sys_tester('isoCN')


data['HCCCl']['pg'] = 'C_inf_v'
data['HCCCl']['rsn'] = 1
data['HCCCl']['mol'] = """
  H  0 0 -1.0
  C  0 0  0.0
  C  0 0  1.1
  Cl 0 0  2.6
"""
data['isoHCCCl']['idx'] = 1
data['isoHCCCl']['pg'] = 'C_inf_v'
data['isoHCCCl']['rsn'] = 1


def test_HCCCl():
    sys_tester('HCCCl')


def test_isoHCCCl():
    sys_tester('isoHCCCl')


data['CHFClBr']['pg'] = "C1"
data['CHFClBr']['rsn'] = 1
data["CHFClBr"]['mol'] = """
  H
  C  1 1.0
  F  2 1.0 1 105.0
  Cl 2 1.0 1 105.0 3  120.0
  Br 2 1.0 1 105.0 3 -120.0
"""
data['isoCHFClBr']['idx'] = 0
data['isoCHFClBr']['pg'] = 'C1'
data['isoCHFClBr']['rsn'] = 1


def test_CHFClBr():
    sys_tester('CHFClBr')


def test_isoCHFClBr():
    sys_tester('isoCHFClBr')


data['CH2ClBr']['pg'] = "Cs"
data['CH2ClBr']['rsn'] = 1
data['CH2ClBr']['mol'] = """
  Cl
  C  1 1.0
  Br 2 1.0 1 105.0
  H  2 1.0 1 105.0 3  120.0
  H  2 1.0 1 105.0 3 -120.0
"""
data['isoCH2ClBr']['idx'] = 4
data['isoCH2ClBr']['pg'] = 'C1'
data['isoCH2ClBr']['rsn'] = 1


def test_CH2ClBr():
    sys_tester('CH2ClBr')


def test_isoCH2ClBr():
    sys_tester('isoCH2ClBr')


data['HOCl']['pg'] = "Cs"
data['HOCl']['rsn'] = 1
data['HOCl']['mol'] = """
  H
  O 1 1.0
  Cl 2 1.7 1 110.0
"""
data['isoHOCl']['idx'] = 0
data['isoHOCl']['pg'] = 'Cs'
data['isoHOCl']['rsn'] = 1


def test_HOCl():
    sys_tester('HOCl')


def test_isoHOCl():
    sys_tester('isoHOCl')


data['C4H4Cl2F2']['au'] = True
data['C4H4Cl2F2']['pg'] = "Ci"
data['C4H4Cl2F2']['rsn'] = 1
data['C4H4Cl2F2']['mol'] = """
  units    bohr
  C     0.432781050498     1.898774028282     0.810337938486
  C    -1.658744642774     0.805191018766    -0.984829058337
  C     1.658744642774    -0.805191018766     0.984829058337
  C    -0.432781050498    -1.898774028282    -0.810337938486
  H    -0.317971784026     2.532165941971     2.640915161238
  H    -1.615729990528     1.614062700629    -2.881498569657
  H     1.615729990528    -1.614062700629     2.881498569657
  H     0.317971784026    -2.532165941971    -2.640915161238
  Cl   -4.852178875691     1.024620478757     0.190249941464
  Cl    4.852178875691    -1.024620478757    -0.190249941464
  F    -1.913713787211    -3.739567959534     0.258534542158
  F     1.913713787211     3.739567959534    -0.258534542158
"""
data['isoC4H4Cl2F2']['idx'] = 2
data['isoC4H4Cl2F2']['pg'] = 'C1'
data['isoC4H4Cl2F2']['rsn'] = 1


def test_C4H4Cl2F2():
    sys_tester('C4H4Cl2F2')


def test_isoC4H4Cl2F2():
    sys_tester('isoC4H4Cl2F2')


data['HOOH_dimer']['pg'] = "Ci"
data['HOOH_dimer']['rsn'] = 1
data['HOOH_dimer']['mol'] = """
  H   0.9911262285  -1.7979226333   0.1465182515
  O   2.7691093095  -1.3485218649  -0.0071557684
  O   2.5178030311   1.3808374923  -0.1154058014
  H   3.2883200453   1.8308595095   1.4757706825
  H  -3.2883200453  -1.8308595095  -1.4757706825
  O  -2.5178030311  -1.3808374923   0.1154058014
  O  -2.7691093095   1.3485218649   0.0071557684
  H  -0.9911262285   1.7979226333  -0.1465182515
"""
data['isoHOOH_dimer']['idx'] = 7
data['isoHOOH_dimer']['pg'] = 'C1'
data['isoHOOH_dimer']['rsn'] = 1


def test_HOOH_dimer():
    sys_tester('HOOH_dimer')


def test_isoHOOH_dimer():
    sys_tester('isoHOOH_dimer')


data['HOOH']['pg'] = "C2"
data['HOOH']['rsn'] = 2
data['HOOH']['mol'] = """
  H
  O 1 1.0
  O 2 1.5 1 110.0
  H 3 1.0 2 110.0 1 60.0
"""
data['isoHOOH']['idx'] = 2
data['isoHOOH']['pg'] = 'C1'
data['isoHOOH']['rsn'] = 1


def test_HOOH():
    sys_tester('HOOH')


def test_isoHOOH():
    sys_tester('isoHOOH')


data['NOHOHOH']['pg'] = 'C3'
data['NOHOHOH']['rsn'] = 3
data['NOHOHOH']['mol'] = """
  X
  N 1 1.0
  O 2 1.5 1 100.0
  O 2 1.5 1 100.0  3  120.0
  O 2 1.5 1 100.0  3 -120.0
  H 3 1.0 2 110.0 4 0.0
  H 4 1.0 2 110.0 5 0.0
  H 5 1.0 2 110.0 3 0.0
"""
data['isoNOHOHOH']['idx'] = 2
data['isoNOHOHOH']['pg'] = 'Cs'
data['isoNOHOHOH']['rsn'] = 1


def test_NOHOHOH():
    sys_tester('NOHOHOH')


def test_isoNOHOHOH():
    sys_tester('isoNOHOHOH')


data['H2O']['pg'] = "C2v"
data['H2O']['rsn'] = 2
data['H2O']['mol'] = """
  H
  O 1 1.0
  H 2 1.0 1 109.5
"""
data['isoH2O']['idx'] = 2
data['isoH2O']['pg'] = 'Cs'
data['isoH2O']['rsn'] = 1


def test_H2O():
    sys_tester('H2O')


def test_isoH2O():
    sys_tester('isoH2O')


data['CH2F2']['au'] = True
data['CH2F2']['pg'] = "C2v"
data['CH2F2']['rsn'] = 2
data['CH2F2']['mol'] = """
  units au
  C     0.0000000000  -0.0000000000   1.0890958457
  F     0.0000000000  -2.1223155812  -0.4598161475
  F    -0.0000000000   2.1223155812  -0.4598161475
  H     1.7084139850   0.0000000000   2.1841068002
  H    -1.7084139850  -0.0000000000   2.1841068002
"""
data['isoCH2F2']['idx'] = 3
data['isoCH2F2']['pg'] = 'Cs'
data['isoCH2F2']['rsn'] = 1


def test_CH2F2():
    sys_tester('CH2F2')


def test_isoCH2F2():
    sys_tester('isoCH2F2')


data['NH3']['pg'] = "C3v"
data['NH3']['rsn'] = 3
data['NH3']['mol'] = """
  X
  N 1 1.0
  H 2 rNH 1 aXNH
  H 2 rNH 1 aXNH 3 120.0
  H 2 rNH 1 aXNH 4 120.0

  rNH = 0.95
  aXNH = 115.0
"""
data['isoNH3']['idx'] = 3
data['isoNH3']['pg'] = 'Cs'
data['isoNH3']['rsn'] = 1


def test_NH3():
    sys_tester('NH3')


def test_isoNH3():
    sys_tester('isoNH3')


data['BrF5']['pg'] = "C4v"
data['BrF5']['rsn'] = 4
data['BrF5']['mol'] = """
 F
 Br 1 r
 F  2 r 1 90.0
 F  2 r 3 90.0 1  90.0
 F  2 r 3 90.0 1 -90.0
 F  2 r 1 90.0 3 180.0
 r = 1.7
"""
data['isoBrF5']['idx'] = 3
data['isoBrF5']['pg'] = 'Cs'
data['isoBrF5']['rsn'] = 1


def test_BrF5():
    sys_tester('BrF5')


def test_isoBrF5():
    sys_tester('isoBrF5')


data['N2H2']['pg'] = "C2h"
data['N2H2']['rsn'] = 2
data['N2H2']['mol'] = """
  N
  N 1 rNN
  H 1 rNH 2 aHNN
  H 2 rNH 1 aHNN 3 180.0
  rNH  = 1.0
  rNN  = 1.4
  aHNN = 140.0
"""


def test_N2H2():
    sys_tester('N2H2')


data['NOHOHOH']['pg'] = "C3h"
data['NOHOHOH']['rsn'] = 3
data['NOHOHOH']['mol'] = """
  X
  N 1 1.0
  O 2 1.5 1 90.0
  O 2 1.5 1 90.0  3  120.0
  O 2 1.5 1 90.0  3 -120.0
  H 3 1.0 2 110.0 4 0.0
  H 4 1.0 2 110.0 5 0.0
  H 5 1.0 2 110.0 3 0.0
"""
data['isoNOHOHOH']['idx'] = 5
data['isoNOHOHOH']['pg'] = 'Cs'
data['isoNOHOHOH']['rsn'] = 1


def test_NOHOHOH():
    sys_tester('NOHOHOH')


def test_isoNOHOHOH():
    sys_tester('isoNOHOHOH')


# 1,3,5,7-tetrafluorocyclooctatetraene
data['TFCOT']['pg'] = "S4"
data['TFCOT']['rsn'] = 2
data['TFCOT']['mol'] = """
  C       -1.618188     -0.437140     -0.409373
  C       -1.394411      0.896360     -0.429596
  C       -0.896360     -1.394411      0.429596
  C       -0.437140      1.618188      0.409373
  C        0.437140     -1.618188      0.409373
  C        0.896360      1.394411      0.429596
  C        1.394411     -0.896360     -0.429596
  C        1.618188      0.437140     -0.409373
  F        2.147277     -1.690111     -1.235043
  F        1.690111      2.147277      1.235043
  F       -2.147277      1.690111     -1.235043
  F       -1.690111     -2.147277      1.235043
  H        0.878010     -2.418132      1.029595
  H       -2.418132     -0.878010     -1.029595
  H       -0.878010      2.418132      1.029595
  H        2.418132      0.878010     -1.029595
"""


def test_TFCOT():
    sys_tester('TFCOT')


data['Li_H2O_4_p']['pg'] = "S4"
data['Li_H2O_4_p']['rsn'] = 2
data['Li_H2O_4_p']['mol'] = """
   1 1
   X
   Li 1 1.0
   X 2 1.0 1 90.0
   X 2 1.0 3 90.0 1 180.0
   O 2 oli 1 olix 3 -90.0
   O 2 oli 1 olix 3 90.0
   O 2 oli 4 olix 3 0.0
   O 2 oli 4 olix 3 180.0
   H 5 oh1 2 lioh1 1 xlioh1
   H 5 oh2 2 lioh2 1 xlioh2
   H 6 oh1 2 lioh1 1 xlioh1
   H 6 oh2 2 lioh2 1 xlioh2
   H 7 oh1 2 lioh1 4 -xlioh1
   H 7 oh2 2 lioh2 4 -xlioh2
   H 8 oh1 2 lioh1 4 -xlioh1
   H 8 oh2 2 lioh2 4 -xlioh2
   olix=52.0
   oli=1.9
   oh1=0.952
   oh2=0.950
   lioh1=125.4
   lioh2=124.8
   xlioh1=-40.0
   xlioh2=135.0
"""


def test_Li_H2O_4_p():
    sys_tester('Li_H2O_4_p')


data['ethylene_cation']['pg'] = "D2"
data['ethylene_cation']['rsn'] = 4
data['ethylene_cation']['mol'] = """
  C1
  C2 C1 rCC
  H1 C1 rCH C2 aHCC
  H2 C1 rCH C2 aHCC H1 180.0
  H3 C2 rCH C1 aHCC H1 D
  H4 C2 rCH C1 aHCC H3 180.0
  rCC  = 1.41
  rCH  = 1.09
  aHCC = 122.0
  D    = 45.0
"""


def test_ethylene_cation():
    sys_tester('ethylene_cation')


data['ethane_gauche']['pg'] = "D3"
data['ethane_gauche']['rsn'] = 6
data['ethane_gauche']['mol'] = """
  H
  C 1 1.0
  C 2 1.5 1 110.0
  H 3 1.0 2 110.0 1   20.0
  H 3 1.0 2 110.0 1  140.0
  H 3 1.0 2 110.0 1 -100.0
  H 2 1.0 3 110.0 1  120.0
  H 2 1.0 3 110.0 1 -120.0
"""


def test_ethane_gauche():
    sys_tester('ethane_gauche')


data['triplet_ethylene']['pg'] = "D2d"
data['triplet_ethylene']['rsn'] = 4
data['triplet_ethylene']['mol'] = """
  C1
  C2 C1 rCC
  H1 C1 rCH C2 aHCC
  H2 C1 rCH C2 aHCC H1 180.0
  H3 C2 rCH C1 aHCC H1 D
  H4 C2 rCH C1 aHCC H3 180.0
  rCC  = 1.41
  rCH  = 1.09
  aHCC = 122.0
  D    = 90.0
"""


def test_triplet_ethylene():
    sys_tester('triplet_ethylene')


data['allene']['pg'] = "D2d"
data['allene']['rsn'] = 4
data['allene']['mol'] = """
  H -2.0  0.0  1.0
  H -2.0  0.0 -1.0
  C -1.5  0.0  0.0
  C  0.0  0.0  0.0
  C  1.5  0.0  0.0
  H  2.0  1.0  0.0
  H  2.0 -1.0  0.0
"""
data['isoallene']['idx'] = [0, 6]
data['isoallene']['pg'] = 'C2'
data['isoallene']['rsn'] = 2


def test_allene():
    sys_tester('allene')


def test_isoallene():
    sys_tester('isoallene')


data['ethane_staggered']['pg'] = "D3d"
data['ethane_staggered']['rsn'] = 6
data['ethane_staggered']['mol'] = """
  H
  C 1 1.0
  C 2 1.5 1 110.0
  H 3 1.0 2 110.0 1   60.0
  H 3 1.0 2 110.0 1  -60.0
  H 3 1.0 2 110.0 1  180.0
  H 2 1.0 3 110.0 1  120.0
  H 2 1.0 3 110.0 1 -120.0
"""


def test_ethane_staggered():
    sys_tester('ethane_staggered')


data['singlet_ethylene']['pg'] = "D2h"
data['singlet_ethylene']['rsn'] = 4
data['singlet_ethylene']['mol'] = """
    C1
    C2 C1 rCC
    H1 C1 rCH C2 aHCC
    H2 C1 rCH C2 aHCC H1 180.0
    H3 C2 rCH C1 aHCC H1 D
    H4 C2 rCH C1 aHCC H3 180.0
    rCC  = 1.41
    rCH  = 1.09
    aHCC = 122.0
    D    = 0.0
"""
data['isosinglet_ethylene']['idx'] = [2, 4]
data['isosinglet_ethylene']['pg'] = 'C2v'
data['isosinglet_ethylene']['rsn'] = 2


def test_singlet_ethylene():
    sys_tester('singlet_ethylene')


def test_isosinglet_ethylene():
    sys_tester('isosinglet_ethylene')


data['ethane_eclipsed']['pg'] = "D3h"
data['ethane_eclipsed']['rsn'] = 6
data['ethane_eclipsed']['mol'] = """
  H
  C 1 1.0
  C 2 1.5 1 110.0
  H 3 1.0 2 110.0 1   00.0
  H 3 1.0 2 110.0 1  120.0
  H 3 1.0 2 110.0 1 -120.0
  H 2 1.0 3 110.0 1  120.0
  H 2 1.0 3 110.0 1 -120.0
"""


def test_ethane_eclipsed():
    sys_tester('ethane_eclipsed')


data['BH4p']['pg'] = "D4h"
data['BH4p']['rsn'] = 8
data['BH4p']['mol'] = """
 1 1
 X
 B 1 1.0
 H 2 1.0 1 90.0
 H 2 1.0 1 90.0 3  90.0
 H 2 1.0 1 90.0 3 180.0
 H 2 1.0 1 90.0 3 -90.0
"""
data['isoBH4p']['idx'] = [2, 4]
data['isoBH4p']['pg'] = 'D2h'
data['isoBH4p']['rsn'] = 4


def test_BH4p():
    sys_tester('BH4p')


def test_isoBH4p():
    sys_tester('isoBH4p')


data['CH4']['pg'] = "Td"
data['CH4']['rsn'] = 12
data['CH4']['mol'] = """
   C
   H 1 r
   H 1 r 2 TDA
   H 1 r 2 TDA 3 120
   H 1 r 2 TDA 4 120
   r = 1.09
"""
data['isoCH4']['idx'] = 1
data['isoCH4']['pg'] = 'C3v'
data['isoCH4']['rsn'] = 3


def test_CH4():
    sys_tester('CH4')


def test_isoCH4():
    sys_tester('isoCH4')


data['SF6']['pg'] = "Oh"
data['SF6']['rsn'] = 24
data['SF6']['mol'] = """
  F
  S 1 r
  F 2 r 1 90.0
  F 2 r 1 90.0 3  90.0
  F 2 r 1 90.0 3 180.0
  F 2 r 1 90.0 3 -90.0
  F 2 r 5 90.0 1 180.0
  r = 1.8
"""
data['isoSF6']['idx'] = [0, 2, 4, 6]
data['isoSF6']['pg'] = 'D4h'
data['isoSF6']['rsn'] = 8


def test_SF6():
    sys_tester('SF6')


def test_isoSF6():
    sys_tester('isoSF6')


data['Ih']['au'] = True
data['Ih']['pg'] = 'Ih'
data['Ih']['rsn'] = 60
data['Ih']['mol'] = """
  unit = au
  0 1
  H   0   1   x
  H   0  -1   x
  H   0   1  -x
  H   0  -1  -x
  H   1   x   0
  H  -1   x   0
  H   1  -x   0
  H  -1  -x   0
  H   x   0   1
  H   x   0  -1
  H  -x   0   1
  H  -x   0  -1
  x = 1.618033988749894848
"""
data['isoIh']['idx'] = [1, 4, 5, 8, 10]
data['isoIh']['pg'] = 'C5v'
data['isoIh']['rsn'] = 5


def test_Ih():
    sys_tester('Ih')


def test_isoIh():
    sys_tester('isoIh')


def mol_tester(lbl, molstr, pg, sigma, refgeomang, isbohr=False, iso=False):
    symmol = qcdb.Molecule(molstr)
    if iso is not False:
        if isinstance(iso, int):
            iso = [iso]
        for at in iso:
            # mass needn't make sense for element, just breaking the symmetry
            symmol.set_mass(at, 2.014)
    symmol.update_geometry()
    symmol.axis_representation()
    assert compare_strings(pg, symmol.get_full_point_group(), pg + " point group: " + lbl)
    assert compare_integers(sigma, symmol.rotational_symmetry_number(), pg + " sigma")
    if isbohr:
        geom_now = symmol.full_geometry()
    else:
        geom_now = qcdb.mscale(symmol.full_geometry(), qcdb.constants.bohr2angstroms)
    if refgeomang:
        assert compare_matrices(refgeomang, geom_now, 6, pg + " orientation")


def sys_tester(sys):
    if sys.startswith('iso'):
        base = sys[3:]
        isbohr = data[base].get('au', False)
        mol_tester(sys, data[base]['mol'], data[sys]['pg'], data[sys]['rsn'], None, isbohr, iso=data[sys]['idx'])
    else:
        isbohr = data[sys].get('au', False)
        mol_tester(sys, data[sys]['mol'], data[sys]['pg'], data[sys]['rsn'], data[sys]['ref'], isbohr)
