"""
Simplified version of the interfacing part of PyGrendel that only parses geometries out
of a given file.  Some subclasses may support parsing of other properties like charge and
multiplicity, but the bare minimum simply parses out a geometry.  This is intended to
make it easier to support a wider variety of files than can be fully implemented in the
general input_generator/output_parser motif.
"""
from __future__ import print_function
from abc import ABCMeta, abstractmethod, abstractproperty
import re
from grendel.chemistry.molecule import Molecule
from grendel.util.parsing import re_float_or_int_or_scientific
from tempfile import TemporaryFile
from os.path import exists
from grendel.chemistry import ElementData
from grendel.util.strings import indent

known_parsers = []

def parse_geometry(filename_or_string, filetype=None):
    if exists(filename_or_string):
        return parse_geometry_file(filename_or_string, filetype=filetype)
    else:
        try:
            return parse_geometry_string(filename_or_string, filetype=filetype)
        except:
            print("Couldn't parse string, and no file by that name exists:\n{}\nReraising error...".format(
                filename_or_string
            ))
            raise

def parse_geometry_file(filename, filetype=None):
    # TODO check file extension
    possible_parsers = []
    for parser in known_parsers:
        if filetype is None:
            if parser.is_this_filetype(filename):
                possible_parsers.append(parser)
        else:
            if parser._matches_filetype_name(filetype):
                possible_parsers.append(parser)
    if len(possible_parsers) == 0:
        if filetype is None:
            raise ValueError("Can't automatically determine file type of file {}".format(
                filename
            ))
        else:
            raise ValueError("No known parsers for filetype '{}'".format(
                filetype
            ))
    elif len(possible_parsers) > 1:
        raise ValueError("File type of {} is ambiguous.  File could be interpreted as:"
                         " include:\n{}".format(
            filename,
            indent("\n".join(p.format_description for p in possible_parsers))
        ))
    else:
        return possible_parsers[0].parse_file(filename)

def parse_geometry_string(string, filetype=None):
    possible_parsers = []
    for parser in known_parsers:
        if filetype is None:
            if parser.string_is_filetype(string):
                possible_parsers.append(parser)
        else:
            if parser._matches_filetype_name(filetype):
                possible_parsers.append(parser)
    if len(possible_parsers) == 0:
        if filetype is None:
            raise ValueError("Can't automatically determine file type of file with data:\n{}".format(
                string
            ))
        else:
            raise ValueError("No known parsers for filetype '{}'".format(
                filetype
            ))
    elif len(possible_parsers) > 1:
        raise ValueError("File type is ambiguous.  File could be interpreted as:"
                         " include:\n{}\nfor file with data:\n{}".format(
            indent("\n".join(p.format_description for p in possible_parsers)), string
        ))
    else:
        return possible_parsers[0].parse_string(string)

#region | Abstract Base Classes                                                     {{{1 |

class GeometryParser(object):
    __metaclass__ = ABCMeta

    def __init__(self, file_path=None, string=None):
        if (file_path is not None and string is not None) or (file_path is None and string is None):
            raise TypeError("Geometry parsers require a file_path or a string, but not both.")
        self.file_path = file_path
        self.string = string

    @abstractproperty
    def format_description(self):
        """
        A short, human-readable description of the format that the parser handles.  For
        instance, "Q-Chem input file" or "MPQC output file"
        """
        return NotImplemented

    @abstractproperty
    def filetype_strings(self):
        """
        A list of strings that the user can specify for the "filetype" argument to parse_geometry()
        and mean that this parser should be used.  Strings are case insensative, and whitespace as
        well as all other non-alphanumeric input is ignored.  This should, of course, not contain
        any strings describing other classes
        """

    def parse(self):
        """
        Parse out the geometry.  Handles string version and file version and then
        calls the subclass method that implements the actual parsing
        @rtype: Molecule
        @return: Returns a Molecule instance with the geometry from the given file.  Some
            subclasses may also set the charge, multiplicity, description, and other attributes
        """
        if self.file_path is not None:
            return self._parse_file(self.file_path)
        else:
            return self._parse_string(self.string)

    def _parse_file(self, filename):
        # TODO User friendly errors
        with open(filename) as f:
            rv = self._parse_stream(f)
        return rv

    def _parse_string(self, string):
        with TemporaryFile("w+") as f:
            f.write(string)
            f.flush()
            f.seek(0)
            rv = self._parse_stream(f)
        return rv

    @abstractmethod
    def _parse_stream(self, f):
        """
        Parse out the geometry from file-like object f.  Must be implemented by subclasses
        @rtype: Molecule
        @return: Returns a Molecule instance with the geometry from the given file.  Some
            subclasses may also set the charge, multiplicity, description, and other attributes
        """
        return NotImplemented

    @classmethod
    def _matches_filetype_name(cls, string):
        string = str(string)
        def _strip(s):
            return re.sub(r'\W', '', s.lower())
        string = _strip(string)
        for tystr in cls.filetype_strings:
            if string == _strip(tystr):
                return True
        #----------------------------------------#
        return False

    @classmethod
    def parse_file(cls, file_path):
        return cls(file_path=file_path).parse()

    @classmethod
    def parse_string(cls, string):
        return cls(string=string).parse()

    @classmethod
    def is_this_filetype(cls, filename):
        """
        Subclasses should override this so that parse_geometry() can figure out whether or
        not the file `filename` is the kind of thing that can be parsed by the given parser.
        It is *not* a validation; it only needs to be thorough enough to determine unambiguously
        if the file is useable by the parser.  The default implementation calls string_is_filetype,
        which in turn returns False by default, meaning that the parser will always be left
        out of parse_geometry().

        @see GeometryParser.string_is_filetype()
        """
        with open(filename) as f:
            return cls.string_is_filetype(f.read())

    @classmethod
    def string_is_filetype(cls, string):
        """
        Subclasses should override this so that parse_geometry() can figure out whether or
        not the data in `string` is from a file of the type that the given parser can parse.

        @see GeometryParser.is_this_filetype()
        """
        return False

class RegexBasedGeometryParser(GeometryParser):

    @abstractproperty
    def regex(self):
        """
        Should be a compiled regex for which group(1) is the xyz-like geometry
        """
        return NotImplemented

    def _parse_stream(self, f):
        data = f.read()
        matches = list(self.regex.finditer(data))
        if len(matches) == 0:
            raise ValueError("No geometry found in file {}".format(self.file_path))
        elif len(matches) > 1:
            raise ValueError("More than 1 geometry found in file {} (total of {} matches found)".format(
                self.file_path,
                len(matches)
            ))
        else:
            return Molecule(
                xyz_string=matches[0].group(1)
            )

class DoubleRegexBasedGeometryParser(GeometryParser):
    """
    Base class for parsing by first getting the molecule section using an outer_regex
     and then finding atoms within that section using an inner regex
    """

    @abstractproperty
    def outer_regex(self):
        """
        Should be a compiled regex for which group(1) is a geometry section
        """
        return NotImplemented

    @abstractproperty
    def inner_regex(self):
        """
        Should be a compiled regex for which group(1) the atom's symbol,
        group(2) is the x coordinate,  group(3) is the y coordinate, and
        group(4) is the x coordinate.  Should match exactly once for each atom
        in a molecule section
        """
        return NotImplemented

    def _parse_stream(self, f):
        if self.inner_regex.groups < 4:
            raise ValueError("Not enough capturing groups in inner regex.")
        #----------------------------------------#
        data = f.read()
        matches = list(self.outer_regex.finditer(data))
        if len(matches) == 0:
            raise ValueError("No molecule section found in file {}".format(self.file_path))
        elif len(matches) > 1:
            raise ValueError("More than one molecule section in file {} (total of {} matches found)".format(
                self.file_path,
                len(matches)
            ))
        else:
            xyz_lines = []
            natoms = 0
            mol_section = matches[0].group(1)
            for m in self.inner_regex.finditer(mol_section):
                xyz_lines.append("{0} {1} {2} {3}".format(*m.groups()))
                natoms += 1
            if natoms == 0:
                raise ValueError("No atoms found in geometry section of {}".format(self.file_path))
            return Molecule(
                xyz_string="\n".join(xyz_lines)
            )


#endregion }}}1

#================================================================================#

#region | Parsers                                                                   {{{1 |

class MPQCInputGeometryParser(DoubleRegexBasedGeometryParser):
    """
    Parse a geometry out of an MPQC input file.  If the file does not contain exactly one geometry,
    throw up our hands and complain.
    """

    inner_regex = re.compile(r'(\w{1,3})\s*\[\s*('
        + ")\s+(".join([re_float_or_int_or_scientific]*3)
        + ")\s*\]"
    )
    outer_regex = re.compile(r'\{\s*atoms\s+geometry\s*\}\s*=\s*\{((?:\s*[^}])+)\}', re.MULTILINE)
    format_description = "MPQC input file"
    filetype_strings = ["mpqci", "mpqcin", "MPQC input", "MPQC input file", "keyval", "MPQC keyval", "MPQC keyval input"]

    @classmethod
    def string_is_filetype(cls, string):
        # Make sure there's a mpqc section
        m = re.search(r'mpqc:\(.*?mole(?:\s*=|\s*<\s*\w+\s*>\s*:\s*\().*?\)', string, re.MULTILINE | re.DOTALL)
        if m:
            # Make sure there's a Molecule section
            m2 = re.search(r'<\s*Molecule\s*>\s*:\s*\(.+?\)', string, re.MULTILINE | re.DOTALL)
            if m2:
                return True
        #----------------------------------------#
        return False


known_parsers.append(MPQCInputGeometryParser)

#endregion }}}1

#================================================================================#

#region | OpenBabel interface                                                       {{{1 |

try:
    import openbabel
    have_openbabel = has_openbabel = True

    class OpenBabelGeometryParser(GeometryParser):
        """
        Abstract base class for parsers that use OpenBabel to do their dirty work.
        All Molecule objects output by subclasses of this will have an obmol attribute
        that is set to the openbabel.OBMol object from which it is constructed.
        Idealy, though, this attribute won't be used since all relevent data will
        already be parsed out of it.
        """

        @abstractproperty
        def filename_extensions(self):
            """
            A list of filename extensions that openbabel recognizes as identifying
            the type of the file passed in.
            """


        def __init__(self, *args, **kwargs):
            self._obc = openbabel.OBConversion()
            self._obc.SetInFormat(self.main_filename_extension)
            super(OpenBabelGeometryParser, self).__init__(*args, **kwargs)

        @property
        def main_filename_extension(self):
            return self.filename_extensions[0]

        def _parse_file(self, filename):
            obmol = openbabel.OBMol()
            success = self._obc.ReadFile(obmol, filename)
            if not success:
                raise RuntimeError("OpenBabel parser failed when parsing"
                                   " file {}.  (ReadFile() returned False)".format(filename))
            return self._convert_obmol(obmol)

        def _parse_string(self, string):
            obmol = openbabel.OBMol()
            success = self._obc.ReadString(obmol, string)
            if not success:
                raise RuntimeError("OpenBabel parser failed"
                                   " (ReadString() returned False) when parsing data:\n{}".format(
                    string
                ))
            return self._convert_obmol(obmol)

        def _parse_stream(self, f):
            return self._parse_string(f.data())

        def _convert_obmol(self, obmol):
            symbols = []
            positions = []
            for atom in openbabel.OBMolAtomIter(obmol):
                symbols.append(ElementData.get(atomic_number=atom.GetAtomicNum()).symbol)
                positions.append((atom.x(), atom.y(), atom.z()))
            charge = self._get_charge(obmol)
            multiplicity = self._get_multiplicity(obmol)
            description = self._get_description(obmol)
            rv = Molecule(
                atom_names=symbols,
                cart_mat=positions,
                charge=charge,
                multiplicity=multiplicity,
                description=description
            )
            rv.obmol = obmol
            return rv

        def _get_charge(self, obmol):
            return obmol.GetTotalCharge()

        def _get_multiplicity(self, obmol):
            if obmol.HasSpinMultiplicityAssigned():
                return obmol.GetTotalSpinMultiplicity()
            else:
                return None

        def _get_description(self, obmol):
            return obmol.GetTitle()

    class MPQCOutputGeometryParser(OpenBabelGeometryParser):
        filename_extensions = ["mpqc"]
        format_description = "MPQC output file"
        filetype_strings = ["mpqco", "mpqcout", "MPQC output", "MPQC output file"]

        @classmethod
        def string_is_filetype(cls, string):
            lines = string.splitlines()
            # see if the title is in the first 50 lines
            first = "\n".join(lines[:50])
            m = re.search(r'MPQC: Massively Parallel Quantum Chemistry', first)
            if m:
                return True
            #----------------------------------------#
            return False

    known_parsers.append(MPQCOutputGeometryParser)

    class MoldenGeometryParser(OpenBabelGeometryParser):
        filename_extensions = ["molden"]
        format_description = "Molden format"
        specification_url = "http://www.cmbi.ru.nl/molden/molden_format.html"
        filetype_strings = ["molden"]

        @classmethod
        def string_is_filetype(cls, string):
            lines = string.splitlines()
            # see if the title is in the first 3 lines (really needs to be the
            #   first one, but we'll let it slide a bit
            first = "\n".join(lines[:3])
            m = re.search(r'\[Molden Format\]', first, re.IGNORECASE)
            if m:
                return True
            #----------------------------------------#
            return False
    known_parsers.append(MoldenGeometryParser)

    class QChemGeometryParser(OpenBabelGeometryParser):
        filename_extensions = ["qcout"]
        format_description = "Q-Chem output file"
        specification_url = "http://www.q-chem.com/"
        filetype_strings = ["qchem", "qchem output" , "qcout", "qc", "qchem output file", "qchemout"]
    known_parsers.append(QChemGeometryParser)


except ImportError:
    have_openbabel = has_openbabel = False


#endregion }}}1
