from functools import partial
import re
from grendel.coordinates.cartesian_coordinate import CartesianCoordinate
from grendel.gmath.matrix import Matrix
from grendel.gmath.vector import Vector

from grendel.interface.computation_details import ComputationDetails, Methods
from grendel.interface.output_parser import PropertyParser, OutputParser, OutputParsingError, RepresentationDependentPropertyParser, OptimizedGeometryParser
from grendel.chemistry.molecular_properties import Energy
from grendel.chemistry.derivative_properties import PropertyDerivative
from grendel.representations.cartesian_representation import CartesianRepresentation, X, Y, Z
from grendel.util.iteration import flattened
from grendel.util.parsing import RegexSequence, MatrixRegexSequence
from grendel.util.strings import indented
from grendel.util.units import Hartrees, Angstrom, Bohr

###########
# Getters #
###########


def hf_dft_sp_energy_parser(file_contents):
    sequence = RegexSequence(
        start = [
            r' ---------------------------',
            r'  Cycle       Energy        ',
            r' ---------------------------'
        ],
        main = r'\s*\d+\s+(-\d+\.\d+)\s+\d+\.\d+',
        stop = r' ---------------------------------------',
        last_only=True
    )
    for line in file_contents.splitlines():
        sequence.parse_line(line)
    if len(sequence.groups()) == 0:
        raise OutputParsingError
    return float(sequence.groups()[0]) * Hartrees

def parse_qchem_matrix(data):
    ret_val = [[]]
    def parse_line(cols, m):
        row = int(m.group(1)) - 1
        data = [float(f) for f in m.group(2).strip().split()]
        while len(ret_val) < row + 1:
            ret_val.append([])
        while len(ret_val[row]) < cols[-1] + 1:
            ret_val[row].append(None)
        if len(cols) != len(data):
            raise OutputParsingError("ragged matrix in output:\n{}".format(indented(data)))
        for d, col in zip(data, cols):
            ret_val[row][col] = d

    def parse_group_of_lines(m):
        columns = [int(s.strip()) - 1 for s in m.group(1).split()]
        re.sub(r'(\d)((?:\s+-?\d+\.\d+)+)', partial(parse_line, columns), m.group(2))

    re.sub(r'((?:\d+[ \t]*)+\n)[ \t]+((?:\d+(?:\s+-?\d+\.\d+)+\s+)+)', parse_group_of_lines, data)

    if None in flattened(ret_val):
        raise OutputParsingError("output matrix missing value:\n{}".format(indented(data)))

    return Matrix(ret_val)

#TODO copy parts of this to some library function
def hf_dft_gradient_parser(file_contents):
    section_re = re.search(r'''
        Gradient.of.SCF.Energy.*\s+                      # The title
        ((?:(?:\d+\s+)+(?:\d+(?:\s+-?\d+\.\d+)+\s+)+)+)  # The data
        Max.gradient.component                           # The value signaling the end
        ''',
        file_contents,
        re.MULTILINE|re.VERBOSE
    )
    if section_re is None:
        raise OutputParsingError("could not find gradient section")
    data = section_re.group(1)
    ret_val = parse_qchem_matrix(data)
    ret_val = Matrix(ret_val, units=Hartrees/Bohr).transpose().flatten(order='C').view(Vector)
    return ret_val

def hf_dft_hessian_parser(file_contents):
    section_re = re.search(r'''
        Hessian.of.the.SCF.Energy.*\s+                   # The title
        ((?:(?:\d+\s+)+(?:\d+(?:\s+-?\d+\.\d+)+\s+)+)+)  # The data
        Gradient.time:                                   # The value signaling the end
        ''',
        file_contents,
        re.MULTILINE|re.VERBOSE
    )
    if section_re is None:
        raise OutputParsingError("could not find hessian section")
    data = section_re.group(1)
    ret_val = parse_qchem_matrix(data)
    ret_val = Matrix(ret_val, units=Hartrees/Bohr**2)
    return ret_val


def standard_orientation_getter(file_contents, molecule):
    m = re.search(r"""
        -+.*\s+                                                  # A divider line
        Standard\sNuclear\sOrientation\ \((Angstroms|Bohr)\)\s+  # Title line
        .*\s+                                                    # A line I don't care about
        -+.*\s*                                                  # A divider line
        ((?:\d\s+\w{1,3}\s+(?:-?\d+\.\d+\s+){3})+)               # Geometry lines
        -+                                                       # Ending divider line
        """,
        file_contents,
        re.VERBOSE|re.MULTILINE
    )
    units = Angstrom if m.group(1) == 'Angstroms' else Bohr
    if not m:
        raise OutputParsingError
    ret_val = CartesianRepresentation(molecule)
    atom_num = [0]
    def parse_coord(match):
        # TODO Sanity checking for atom symbol
        x, y, z = tuple(
            CartesianCoordinate(
                atom = molecule.atoms[atom_num[0]],
                direction = direct,
                parent = ret_val,
                units = molecule.cartesian_units,
                freeze_value=True,
                value = float(match.groups()[direct])*units)
                    for direct in [X, Y, Z])
        ret_val.add_coordinate_copy(x)
        ret_val.add_coordinate_copy(y)
        ret_val.add_coordinate_copy(z)
        atom_num[0] += 1
        return ''
    re.sub(r'\d\s+\w{1,3}\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', parse_coord, m.group(2))
    return ret_val

# UNTESTED!!!
def opt_geom_getter(file_contents):
    m = re.search(r"""
        Optimization.Cycle:\s+\d+\s+                    # Optimization cycle line
        Coordinates.\((Angstroms|Bohr)\)\s+             # Coordinates, Bohr or Angstroms
        .*\s+                                           # XYZ head line
        ((\d+\s+\s+\w{1,3}\s+(?:-?\d+\.\d+\s+){3})+)    # Geometry lines
        Point.Group                                     # Signals end of optimization
        """,
        file_contents,
        re.VERBOSE|re.MULTILINE
    )
    if not m:
        raise OutputParsingError
    atoms = []
    def parse_coord(match):
        atoms.append(
            Atom(match.group(1), Vector(map(float, match.groups()[1:3])))
        )

    re.sub(r'\d\s+(\w{1,3})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', parse_coord, m.group(2))
    return Molecule(atoms)


def opt_energy_getter(file_contents):
    raise NotImplementedError



class QChemOutputParser(OutputParser):
    """
    """

    # Put more specific parsers first; the requested MoleculerPropertyDetails
    property_parsers = [
        PropertyParser(
            property=Energy,
            details=ComputationDetails(
                method=(
                    # The "dotted keyword" syntax is on it's way out, since it's too confusing
                    Methods.HF, 'HF',
                    Methods.DFT, 'DFT'
                )
            ),
            units=Hartrees,
            getter=hf_dft_sp_energy_parser
        ),
        RepresentationDependentPropertyParser(
            property = PropertyDerivative(Energy, CartesianRepresentation),
            details = ComputationDetails(
                method=(
                    # The "dotted keyword" syntax is on it's way out, since it's too confusing
                    Methods.HF, 'HF',
                    Methods.DFT, 'DFT'
                )
            ),
            units=Hartrees / Bohr,
            getter=hf_dft_gradient_parser,
            representation_getter=standard_orientation_getter
        ),
        RepresentationDependentPropertyParser(
            property=PropertyDerivative(Energy, CartesianRepresentation, order=2),
            details=ComputationDetails(
                method = (
                    # The "dotted keyword" syntax is on it's way out, since it's too confusing
                    Methods.HF, 'HF',
                    Methods.DFT, 'DFT'
                )
            ),
            units=Hartrees / Bohr ** 2,
            getter=hf_dft_hessian_parser,
            representation_getter=standard_orientation_getter
        ),
        # UNTESTED!!!
        #OptimizedGeometryParser(
        #    property=Energy,
        #    details=ComputationDetails(
        #        method=(
        #            # The "dotted keyword" syntax is on it's way out, since it's too confusing
        #            Methods.HF, 'HF',
        #            Methods.DFT, 'DFT'
        #        )
        #    ),
        #    geometry_getter=opt_geom_getter,
        #    getter=opt_energy_getter
        #)
    ]


#####################
# Dependent Imports #
#####################

from grendel.chemistry.atom import Atom
from grendel.chemistry.molecule import Molecule
