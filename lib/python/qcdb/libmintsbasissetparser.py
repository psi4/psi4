import os
import re
from exceptions import *
from libmintsgshell import *


class Gaussian94BasisSetParser(object):
    """Class for parsing basis sets from a text file in Gaussian 94
    format. Translated directly from the Psi4 libmints class written
    by Justin M. Turney and Andrew C. Simmonett.

    """

    def __init__(self, forced_puream=None):
        """Constructor"""
        # If the parser needs to force spherical or cartesian (e.g., loading old guess)
        self.force_puream_or_cartesian = False if forced_puream is None else True
        # Is the forced value to use puream?  (Otherwise force Cartesian).
        self.forced_is_puream = False if forced_puream is None else forced_puream
        # string filename
        self.filename = None

    def load_file(self, filename, basisname=None):
        """Load and return the file to be used by parse.  Return only
        portion of *filename* pertaining to *basisname* if specified (for
        multi-basisset files) otherwise entire file as list of strings.

        """
        # string filename
        self.filename = filename

        given_basisname = False if basisname is None else True
        found_basisname = False
        basis_separator = re.compile(r'^\s*\[\s*(.*?)\s*\]\s*$')

        # Loads an entire file.
        try:
            infile = open(filename, 'r')
        except IOError:
            raise BasisSetFileNotFound("""BasisSetParser::parse: Unable to open basis set file: %s""" % (filename))
        if os.stat(filename).st_size == 0:
            raise ValidationError("""BasisSetParser::parse: given filename '%s' is blank.""" % (filename))
        contents = infile.readlines()

        lines = []
        for text in contents:
            text = text.strip()
            # If no basisname was given always save the line.
            if given_basisname is False:
                lines.append(text)

            if found_basisname:
                # If we find another [*] we're done.
                if basis_separator.match(text):
                    what = basis_separator.match(text).group(1)
                    break
                lines.append(text)
                continue

            # If the user gave a basisname AND text matches the basisname we want to trigger to retain
            if given_basisname and basis_separator.match(text):
                if basisname == basis_separator.match(text).group(1):
                    found_basisname = True

        return lines

    def parse(self, symbol, dataset):
        """Given a string, parse for the basis set needed for atom.
        * @param symbol atom symbol to look for in dataset
        * @param dataset data set to look through
        dataset can be list of lines or a single string which will be converted to list of lines

        """
        if isinstance(dataset, basestring):
            lines = dataset.split('\n')
        else:
            lines = dataset

        # Regular expressions that we'll be checking for.
        cartesian = re.compile(r'^\s*cartesian\s*', re.IGNORECASE)
        spherical = re.compile(r'^\s*spherical\s*', re.IGNORECASE)
        comment = re.compile(r'^\s*\!.*')  # line starts with !
        separator = re.compile(r'^\s*\*\*\*\*')  # line starts with ****
        ATOM = '(([A-Z]{1,3}\d*)|([A-Z]{1,3}_\w+))'  # match 'C 0', 'Al c 0', 'P p88 p_pass 0' not 'Ofail 0', 'h99_text 0'
        atom_array = re.compile(r'^\s*((' + ATOM + '\s+)+)0\s*$', re.IGNORECASE)  # array of atomic symbols terminated by 0
        shell = re.compile(r'^\s*(\w+)\s*(\d+)\s*(-?\d+\.\d+)')  # Match beginning of contraction
        blank = re.compile(r'^\s*$')
        NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
        primitives1 = re.compile(r'^\s*' + NUMBER + '\s+' + NUMBER + '.*')  # Match s, p, d, f, g, ... functions
        primitives2 = re.compile(r'^\s*' + NUMBER + '\s+' + NUMBER + '\s+' + NUMBER + '.*')  # match sp functions

        # s, p and s, p, d can be grouped together in Pople-style basis sets
        sp = 'SP'
        spd = 'SPD'

        #                a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
        #shell_to_am = [-1,-1,-1, 2,-1, 3, 4, 5, 6,-1, 7, 8, 9,10,11, 1,12,13, 0,14,15,16,17,18,19,20]
        alpha = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
            'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        angmo = [-1, -1, -1, 2, -1, 3, 4, 5, 6, -1, 7, 8,
            9, 10, 11, 1, 12, 13, 0, 14, 15, 16, 17, 18, 19, 20]
        shell_to_am = dict(zip(alpha, angmo))

        # Basis type.
        gaussian_type = 'Pure'

        if self.force_puream_or_cartesian:
            if self.forced_is_puream == False:
                gaussian_type = 'Cartesian'

        # Need a dummy center for the shell.
        center = [0.0, 0.0, 0.0]

        shell_list = []
        lineno = 0
        found = False

        while lineno < len(lines):
            line = lines[lineno]
            lineno += 1

            # Ignore blank lines
            if blank.match(line):
                continue

            # Look for Cartesian or Spherical
            if not self.force_puream_or_cartesian:
                if cartesian.match(line):
                    gaussian_type = 'Cartesian'
#TODO                    if psi4.get_global_option('PUREAM').has_changed():
#TODO                        gaussian_type = 'Pure' if int(psi4.get_global('PUREAM')) else 'Cartesian'
                    continue
                elif spherical.match(line):
                    gaussian_type = 'Pure'
#TODO                    if psi4.get_global_option('PUREAM').has_changed():
#TODO                        gaussian_type = 'Pure' if int(psi4.get_global('PUREAM')) else 'Cartesian'
                    continue
                #end case where puream setting wasn't forced by caller

            # Do some matches
            if comment.match(line):
                continue
            if separator.match(line):
                continue

            # Match: H    0
            # or:    H    O...     0
            if atom_array.match(line):
                what = atom_array.match(line).group(1).split()
                # Check the captures and see if this basis set is for the atom we need.
                found = False
                if symbol in [x.upper() for x in what]:
                    found = True
                    msg = """line %5d""" % (lineno)

                    # Read in the next line
                    line = lines[lineno]
                    lineno += 1

                    # Need to do the following until we match a "****" which is the end of the basis set
                    while not separator.match(line):
                        # Match shell information
                        if shell.match(line):
                            what = shell.match(line)
                            shell_type = str(what.group(1)).upper()
                            nprimitive = int(what.group(2))
                            scale = float(what.group(3))

                            if len(shell_type) == 1:
                                am = shell_to_am[shell_type[0]]

                                exponents = [0.0] * nprimitive
                                contractions = [0.0] * nprimitive

                                for p in range(nprimitive):
                                    line = lines[lineno]
                                    lineno += 1
                                    line = line.replace('D', 'e', 1)
                                    line = line.replace('d', 'e', 1)

                                    what = primitives1.match(line)
                                    # Must match primitives1; will work on the others later
                                    if not what:
                                        raise ValidationError("""Gaussian94BasisSetParser::parse: Unable to match an exponent with one contraction: line %d: %s""" % (lineno, line))
                                    exponent = float(what.group(1))
                                    contraction = float(what.group(2))

                                    # Scale the contraction and save the information
                                    contraction *= scale
                                    exponents[p] = exponent
                                    contractions[p] = contraction

                                # We have a full shell, push it to the basis set
                                shell_list.append(ShellInfo(am, contractions, exponents,
                                    gaussian_type, 0, center, 0, 'Unnormalized'))

                            elif len(shell_type) == 2:
                                # This is to handle instances of SP, PD, DF, FG, ...
                                am1 = shell_to_am[shell_type[0]]
                                am2 = shell_to_am[shell_type[1]]

                                exponents = [0.0] * nprimitive
                                contractions1 = [0.0] * nprimitive
                                contractions2 = [0.0] * nprimitive

                                for p in range(nprimitive):
                                    line = lines[lineno]
                                    lineno += 1
                                    line = line.replace('D', 'e', 1)
                                    line = line.replace('d', 'e', 1)

                                    what = primitives2.match(line)
                                    # Must match primitivies2
                                    if not what:
                                        raise ValidationError("Gaussian94BasisSetParser::parse: Unable to match an exponent with two contractions: line %d: %s" % (lineno, line))
                                    exponent = float(what.group(1))
                                    contraction = float(what.group(2))

                                    # Scale the contraction and save the information
                                    contraction *= scale
                                    exponents[p] = exponent
                                    contractions1[p] = contraction

                                    # Do the other contraction
                                    contraction = float(what.group(3))

                                    # Scale the contraction and save the information
                                    contraction *= scale
                                    contractions2[p] = contraction

                                shell_list.append(ShellInfo(am1, contractions1, exponents,
                                    gaussian_type, 0, center, 0, 'Unnormalized'))
                                shell_list.append(ShellInfo(am2, contractions2, exponents,
                                    gaussian_type, 0, center, 0, 'Unnormalized'))
                            else:
                                raise ValidationError("""Gaussian94BasisSetParser::parse: Unable to parse basis sets with spd, or higher grouping""")
                        else:
                            raise ValidationError("""Gaussian94BasisSetParser::parse: Expected shell information, but got: line %d: %s""" % (lineno, line))
                        line = lines[lineno]
                        lineno += 1

                    break

        if not found:
            #raise BasisSetNotFound("Gaussian94BasisSetParser::parser: Unable to find the basis set for %s in %s" % \
            #   (symbol, self.filename), silent=True)
            return None, None

        return shell_list, msg
