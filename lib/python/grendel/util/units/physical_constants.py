""" Simple structure for retrieval of physical constants with their units
"""

import math
from grendel.util.units.unit import Angstrom, Second, Kilogram, Coulomb, Meter, Second, Joule

__all__ = [
    "PhysicalConstants",
]


PhysicalConstants = dict()

def define_physical_constant(name, value):
    globals()['PhysicalConstants'][name] = globals()[name] = value
    globals()['__all__'].append(name)

# NIST value
define_physical_constant('BohrRadius', 0.52917721092 * Angstrom)
# Old grendel value
#define_physical_constant('BohrRadius', 0.52917720859 * Angstrom)
define_physical_constant('PlanckConstant', 6.62606957e-34 * Joule / Second)
define_physical_constant('ReducedPlanckConstant', PlanckConstant / (2.0 * math.pi))
define_physical_constant('ElectronMass', 9.10938291e-31 * Kilogram)
define_physical_constant('ElementaryCharge', 1.602176565e-19 * Coulomb)
define_physical_constant('SpeedOfLight', 299792458.0 * Meter / Second )
define_physical_constant('AvogadrosNumber', 6.02214129e23 )

