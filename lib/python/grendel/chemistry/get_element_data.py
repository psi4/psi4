#!/usr/bin/env python
from __future__ import print_function

# Note:  Not a grendel module.  This is a script that generates the element_data.py file from
#   the NIST data site:
#   http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=all

from collections import defaultdict
from os.path import dirname, abspath, join as path_join
import re

header = """
from grendel.chemistry.element import Element, Isotope
from grendel.util.units import AtomicMassUnit

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


if __name__ == "__main__":

    mydir = abspath(dirname(__file__))
    isofile = path_join(mydir, "iso.txt")

    data = open(isofile).read()

    section_re = r'Atomic Number = (.*)\nAtomic Symbol = (.*)\nMass Number = (.*)\nRelative Atomic Mass = (.*)\nIsotopic Composition = (.*)\nStandard Atomic Weight = (.*)\nNotes = (.*)\n'
    section_re = re.compile(section_re)

    isotopes = defaultdict(lambda: [])
    element_data = dict()
    symbols = []

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
    is_synthetic=True,
    isotopes=[
        {isotopes}
    ]
)
Elements['{symbol}'] = {symbol}

""".format(
                symbol=symbol,
                atomic_number=eldat[0],
                atomic_weight=eldat[1],
                isotopes=",\n        ".join(isotopes[symbol])
            )
        )
        else:
            out.write("""
{symbol} = Element(
    symbol="{symbol}",
    atomic_number={atomic_number},
    atomic_weight={atomic_weight},
    atomic_weight_uncertainty={aw_uncert},
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
                isotopes=",\n        ".join(isotopes[symbol])
            ))










