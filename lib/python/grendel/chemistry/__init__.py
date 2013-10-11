from copy import copy
import sys
from grendel.util.containers import AttributeLookupTable

__submodules__ = [
    'element_data',
    'atom',
    'molecule',
    'element',
    'molecular_properties',
    'derivative_properties',
]

__all__ = [
    "SampleMolecules",
    "ElementData"
]

for name in __submodules__:
    __import__(__name__ + "." + name)
    m = sys.modules[__name__ + "." + name]
    globals()[name] = m
    if hasattr(m, '__all__'):
        attrlist = copy(m.__all__)
    else:
        attrlist = list(filter(lambda x: x[0]!='_', dir(m)))
    for attr in attrlist:
        globals()[attr] = getattr(m, attr)
    if hasattr(m, '__not_parent_all__'):
        for item in m.__not_parent_all__:
            attrlist.remove(item)
    __all__.extend(attrlist)

__all__.extend(__submodules__)

#import atom; from atom import *
#__all__.extend(atom.__all__)
#
#import molecule; from molecule import *
#__all__.extend(molecule.__all__)


SampleMolecules = dict()
from molecule import Molecule
from element_data import Elements


ElementData = AttributeLookupTable(
    attribute_names=["symbol", "atomic_number", "atomic_weight"],
    initial_values=set(Elements.values())
)


def init_sample_molecules():
    SampleMolecules = globals()['SampleMolecules']
    SampleMolecules["quantum water"] = Molecule.from_z_matrix(
        """
        O
        H 1 1.0
        H 2 1.0 1 90.0
        """
    )

    SampleMolecules["H2O"] = \
    SampleMolecules["water"] =\
    SampleMolecules["Water"] = Molecule(
        description = "CCSD(T)/aug-cc-pVTZ Water",
        xyz_string = """
        3

        O  0.000000  0.000000  0.118154
        H  0.000000  0.758734 -0.472614
        H  0.000000 -0.758734 -0.472614
    """
    )

    SampleMolecules["H2"] = Molecule(
        description = "CCSD(T)/aug-cc-pVTZ H2",
        xyz_string = """
        2

        H  0.000000 0.000000 0.000000
        H  0.000000 0.000000 0.743000
        """
    )

    SampleMolecules["CO2"] = \
    SampleMolecules["Carbon Dioxide"] = Molecule(
        description="Carbon Dioxide from NIH CIR server.",
        xyz_string = """
        3

        O    -1.20800000  -0.00000000  -0.00000000
        O     1.20800000  -0.00000000  -0.00000000
        C    -0.00000000   0.00000000   0.00000000
        """
    )

    SampleMolecules["benzene"] = \
    SampleMolecules["Benzene"] = Molecule(
        description="B3LYP/cc-pVTZ Benzene",
        xyz_string = """
        12

        C    0.000000    1.390732    0.000000
        C    1.204409    0.695366    0.000000
        C    1.204409   -0.695366    0.000000
        C    0.000000   -1.390732    0.000000
        C   -1.204409   -0.695366    0.000000
        C   -1.204409    0.695366    0.000000
        H    0.000000    2.472794    0.000000
        H    2.141502    1.236397    0.000000
        H    2.141502   -1.236397    0.000000
        H    0.000000   -2.472794    0.000000
        H   -2.141502   -1.236397    0.000000
        H   -2.141502    1.236397    0.000000
        """
    )

    SampleMolecules['ethylene'] = \
    SampleMolecules['C2H4'] = Molecule(
        description="Ethylene from NIH CIR server.",
        xyz_string = """
        6

        C     0.65500000   0.00000000   0.00000000
        C    -0.65500000   0.00000000  -0.00000000
        H     1.19500000  -0.93530000   0.00000000
        H     1.19500000   0.93530000   0.00000000
        H    -1.19500000   0.93530000  -0.00000000
        H    -1.19500000  -0.93530000  -0.00000000
        """
    )

    SampleMolecules['methane'] = \
    SampleMolecules['CH4'] = Molecule(
        description="CCSD(T) aug-cc-pVTZ methane from CCCBDB",
        xyz_string = """
        5

        C    0.000000    0.000000    0.000000
        H    0.629271    0.629271    0.629271
        H   -0.629271   -0.629271    0.629271
        H   -0.629271    0.629271   -0.629271
        H    0.629271   -0.629271   -0.629271
        """
    )