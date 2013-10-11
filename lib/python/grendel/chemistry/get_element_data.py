#!/usr/bin/env python
from __future__ import print_function

# Note:  Not a grendel module.  This is a script that generates the element_data.py file from
#   the NIST data site:
#   http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=all

from collections import defaultdict
from os.path import dirname, abspath, join as path_join
import re
import math
from grendel.util.strings import strip_quotes

header = """
from grendel.chemistry.element import Element, Isotope
from grendel.util.units import AtomicMassUnit, ElectronVolt, Picometers, KilojoulesPerMole

MASS_UNIT = AtomicMassUnit


Elements = dict()

# Ghost atoms
Gh = X = Element(
    symbol="X",
    atomic_number=0,
    atomic_weight=0,
    isotopes=[
        Isotope(
            mass=0.0,
            mass_uncertainty=0,
            abundance=1.0,
            abundance_uncertainty=None,
            mass_number=0,
            special_symbol=None
        )
    ]
)
Elements['X'] = X
Elements['Gh'] = X

"""

def get_val_uncert(string):
    if len(string) > 0:
        if "(" in string:
            val, rest = string.split("(")
            val_uncert = rest.split(")")[0]
            digits = len(val.split('.')[-1])
            val_uncert = "0." + "0" * (digits-len(val_uncert)) + val_uncert
            return val, val_uncert
        else:
            return string, None
    else:
        return None, None

class MathematicaElementDatum(object):

    NotApplicable = "NotApplicable"
    NotAvailable = "NotAvailable"
    MissingOther = "Missing"

    def __init__(self, atom_number, prop_name, value_string, units_string, interval, notes):
        self.atom_number = atom_number
        self.property_name = strip_quotes(prop_name)
        #----------------------------------------#
        self.value = self._get_missing_status(value_string)
        if self.value is None:
            self.value = self._parse_value_string(value_string)
        #----------------------------------------#
        self.units = self._get_missing_status(units_string)
        if self.units is None:
            self.units = strip_quotes(units_string)
        #----------------------------------------#
        self.interval = self._get_missing_status(interval)
        if self.interval is None:
            self.interval = self._parse_interval(interval)
        #----------------------------------------#
        self.notes = self._get_missing_status(notes)
        if self.notes is None:
            self.notes = strip_quotes(notes)

    def value_not_missing(self):
        return self.value not in [
            MathematicaElementDatum.NotAvailable,
            MathematicaElementDatum.NotApplicable,
            MathematicaElementDatum.MissingOther
        ]

    @staticmethod
    def _parse_interval(interval):
        m = re.match(r"Interval\[\{(.+), (.+)\}\]", interval)
        if m:
            num1, num2 = m.groups()
            return "({}, {})".format(num1, num2)
        else:
            raise ValueError("Interval '{}' does not match expected interval pattern".format(interval))

    @classmethod
    def _get_missing_status(cls, value_string):
        if value_string == 'Missing["NotApplicable"]':
            return cls.NotApplicable
        elif value_string == 'Missing["NotAvailable"]':
            return cls.NotAvailable
        elif 'Missing[' in value_string:
            return cls.NotAvailable
        elif value_string == "MissingOther":
            return cls.MissingOther
        else:
            return None

    @staticmethod
    def _parse_value_string(value_string):
        m = re.match(r'(-?\d+\.?\d*)(?:`(-?\d+?\.?\d*))?(?:\*^(-?\d+\.?\d*))?', value_string)
        if m:
            num = m.group(1)
            if m.groups()[2] is not None:
                float_expon = float(m.group(3))
                expon = int(round(float_expon))
                value_string = num + "e{}".format(expon)
            else:
                value_string = num
            return value_string
        elif value_string[0] == "{" and value_string[-1] == "}":
            if len(value_string) == 2:
                return []
            else:
                return [MathematicaElementDatum._parse_value_string(v) for v in re.split("\s*,\s*", value_string[1:-1])]
        else:
            return strip_quotes(value_string)

if __name__ == "__main__":

    mydir = abspath(dirname(__file__))

    # Grab some other stuff from the Mathematica element data
    atominfo = path_join(mydir, "mathematica_data.txt")
    math_data = defaultdict(lambda: dict())
    # remove the first line that gives the catagories
    lines = open(atominfo).readlines()[1:]
    for line in lines:
        line = line.strip()
        line = re.sub(r'ElementData\[[^\]]+\]', 'MissingOther', line)
        m = re.match(r"\{([^,]+),\s*([^,]+),\s*([^,]+),\s*(\{(?:[^,}]+,)*[^,}]*\}|[^,]*|Row\[.+\]*]),\s*([^,]+),\s*([^,]+),\s*(Interval\[.*\]|[^,]+),\s*(.+)\}", line)
        if m:
            if m.group(3) in [
                '"ElectronConfiguration"',
                '"QuantumNumbers"',
                '"SpaceGroupName"',
                '"ElectronConfigurationString"',
                '"IsotopeAbundances"'
            ]:
                continue
            atom_number, name, prop, value, units, __, interval, notes = m.groups()
        else:
            raise ValueError("Unrecognized mathematica line {}, {}, {}".format(line, line[0], line[-1]))
        try:
            datum = MathematicaElementDatum(
                atom_number, prop, value, units, interval, notes
            )
        except:
            print(line)
            raise
        math_data[int(atom_number)][datum.property_name] = datum

    # Grab stuff from the isotopes file that comes from NIST
    isofile = path_join(mydir, "isotopes.txt")

    data = open(isofile).read()

    section_re = r'Atomic Number = (.*)\nAtomic Symbol = (.*)\nMass Number = (.*)\nRelative Atomic Mass = (.*)\nIsotopic Composition = (.*)\nStandard Atomic Weight = (.*)\nNotes = (.*)\n'
    section_re = re.compile(section_re)

    isotopes = defaultdict(lambda: [])
    element_data = dict()
    symbols = []
    other_keywords = dict()

    for m in section_re.finditer(data):
        atomic_number = int(m.group(1))
        atomic_symbol = m.group(2).strip()
        mass_number = m.group(3).strip()
        atomic_mass = m.group(4).strip()
        isotopic_composition = m.group(5).strip()
        standard_atomic_weight = m.group(6).strip()
        notes = m.group(7).strip()

        if len(atomic_mass) == 0:
            # this entry is useless to me
            continue

        if len(symbols) < atomic_number:
            symbols.append(atomic_symbol)
        standard_symbol = symbols[atomic_number-1]

        mass, mass_uncert = get_val_uncert(atomic_mass)
        abund, abund_uncert = get_val_uncert(isotopic_composition)

        mass_number = int(mass_number)

        if standard_atomic_weight[0] == '[':
            is_synthetic = True
            atomic_weight = standard_atomic_weight[1:-1]
            atomic_weight_uncertainty = None
        else:
            atomic_weight, atomic_weight_uncertainty = get_val_uncert(standard_atomic_weight)
            is_synthetic = False

        if standard_symbol in element_data:
            if element_data[standard_symbol] != (atomic_number, atomic_weight, atomic_weight_uncertainty, is_synthetic):
                raise ValueError("Element data mismatch for {}{}".format(mass_number, standard_symbol))
        else:
            element_data[standard_symbol] = (atomic_number, atomic_weight, atomic_weight_uncertainty, is_synthetic)
            # This is our first time on this element.  Let's grab some of the mathematica data
            other_keywords[standard_symbol] = []
            mdat = math_data[atomic_number]
            #----------------------------------------#
            if "VanDerWaalsRadius" in mdat:
                vdw = mdat["VanDerWaalsRadius"]
                if vdw.value_not_missing():
                    other_keywords[standard_symbol].append("vdw_radius=" + vdw.value + " * " + vdw.units)
            else:
                print(atomic_number)
            #----------------------------------------#
            ie = math_data[atomic_number]["IonizationEnergies"]
            if ie.value_not_missing():
                vals = [v + " * " + ie.units for v in ie.value]
                other_keywords[standard_symbol].append("ionization_energies=[" + ", ".join(vals) + "]")
            #----------------------------------------#
            eneg = math_data[atomic_number]["Electronegativity"]
            if eneg.value_not_missing():
                other_keywords[standard_symbol].append("electronegativity=" + eneg.value)
            #----------------------------------------#







        if atomic_symbol != standard_symbol:
            special_symbol = atomic_symbol
        else:
            special_symbol = None

        isotopes[standard_symbol].append("""Isotope(
            mass={},
            mass_uncertainty={},
            abundance={},
            abundance_uncertainty={},
            mass_number={},
            special_symbol={}
        )""".format(
            mass, mass_uncert, abund, abund_uncert, mass_number,
            ('"' + special_symbol + '"') if special_symbol is not None else None
        ))
        print("Parsed data for {}{}".format(mass_number, standard_symbol))

    out = open("element_data.py", "w+")
    out.write(header)
    for symbol in symbols:
        eldat = element_data[symbol]
        if eldat[-1]:
            # synthetic
            out.write("""
{symbol} = Element(
    symbol="{symbol}",
    atomic_number={atomic_number},
    atomic_weight={atomic_weight},
    is_synthetic=True,{other_keywords}
    isotopes=[
        {isotopes}
    ]
)
Elements['{symbol}'] = {symbol}

""".format(
                symbol=symbol,
                atomic_number=eldat[0],
                atomic_weight=eldat[1],
                isotopes=",\n        ".join(isotopes[symbol]),
                other_keywords="" if len(other_keywords[symbol]) == 0 else ("\n    " + ",\n    ".join(other_keywords[symbol]) + ",")
            )
        )
        else:
            out.write("""
{symbol} = Element(
    symbol="{symbol}",
    atomic_number={atomic_number},
    atomic_weight={atomic_weight},
    atomic_weight_uncertainty={aw_uncert},{other_keywords}
    isotopes=[
        {isotopes}
    ]
)
Elements['{symbol}'] = {symbol}

""".format(
                symbol=symbol,
                atomic_number=eldat[0],
                atomic_weight=eldat[1],
                aw_uncert=eldat[2],
                isotopes=",\n        ".join(isotopes[symbol]),
                other_keywords="" if len(other_keywords[symbol]) == 0 else ("\n    " + ",\n    ".join(other_keywords[symbol]) + ",")
            ))










