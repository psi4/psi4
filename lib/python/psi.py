import PsiMod

import signal
import sys

def signal_handler(signal, frame):
        print 'PSI4 exiting due to user thrown Ctrl+C'
        sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)

class PsiException: pass
class ValueNotSet (PsiException): pass
class RowAlignmentError(PsiException):

    def __init__(self, l1name, l1, l2name, l2):
        msg = "Rows %s and %s not aligned. Length %d != %d" % (l1name, l2name, l1, l2)
        PsiException.__init__(self, msg)

molecule = None

def new_set_attr(self, name, value):
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)
    if isvar:
        fxn = object.__getattribute__(self, "set_variable")
        fxn(name, value)
        return

    object.__setattr__(self, name, value)

def new_get_attr(self, name):
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)

    if isvar:
        fxn = object.__getattribute__(self, "get_variable")
        return fxn(name)

    return object.__getattribute__(self, name)


def dynamic_variable_bind(cls):
    #class specific
    cls.__setattr__ = new_set_attr
    cls.__getattr__ = new_get_attr


dynamic_variable_bind(PsiMod.Molecule) #pass class type, not class instance

#
# Define geometry to be used by PSI4.
# The molecule created by this will be set in options.
#
# geometry("
#   O  1.0 0.0 0.0
#   H  0.0 1.0 0.0
#   H  0.0 0.0 0.0
#
def geometry(geom, reorient = True, shiftToCOM = True):
    # Create a Molecule object
    global molecule
    molecule = PsiMod.Molecule.create_molecule_from_string(geom)

    # If requested, shift molecule to center of mass
    if shiftToCOM == True:
        molecule.move_to_com()

    # If requested, reorient the molecule.
    if reorient == True:
        molecule.reorient()

    PsiMod.set_active_molecule(molecule)

    return molecule

def dummify():
    if not molecule:
        raise ValueNotSet("no default molecule found")

    #molecule.set_dummy_atom

#
# Set options
#
# options({
#    'WFN': 'SCF',
#    'BASIS': 'STO-3G'
# })
def options(dict, module = ''):
    m = module.strip().upper()
    if len(m) > 0:
        PsiMod.set_default_options_for_module(m)

    for key in dict.keys():
        PsiMod.set_option(key.upper(), dict[key]);

#
# Define geometry to be used by PSI4.
# This geometry will be saved to the checkpoint file.
#
# geometry( [[ "O", 0.0, 0.0, 0.0 ],
#            [ "H", 1.0, 0.0, 0.0 ]] )
#
def geometry_old(geom, reorient = True, prefix = "", chkpt = None, shiftToCOM = True):
    # Make sure the user passed in what we expect
    if isinstance(geom, list) == False:
        raise TypeError("geometry must be a list")

    # Create a Molecule object from PSI4 C++
    molecule = PsiMod.Molecule()

    # For each atom in the geometry given add it to the molecule
    for atom in geom:
        # Make sure the atom is a list
        if isinstance(atom, list) == False:
            raise TypeError("an atom entry is not a list")

        # Make sure the length of the atom entry is long enough
        if len(atom) < 4:
            print "Insufficient information for atom entry:"
            print atom
            raise AttributeError("atom entries must have at least 4 elements")

        molecule.add_atom(1, atom[1], atom[2], atom[3], atom[0], 1.0, 0, 0.0)

    # Print the molecule:
    molecule.print_to_output()

    # If requested, shift molecule to center of mass
    if shiftToCOM == True:
        molecule.move_to_com()

    # If requested, reorient the molecule.
    if reorient == True:
        molecule.reorient()

        molecule.print_to_output()

    # Save the molecule to the checkpoint file using prefix
    if chkpt == None:
        # If the user didn't provide a checkpoint object create one
        # using the shared psio object
        psio = PsiMod.IO.shared_object()
        chkpt = PsiMod.Checkpoint(psio, 1)

    molecule.save_to_checkpoint(chkpt, prefix)

    return molecule

class Table:

    def __init__(self, rows=(),
                 row_label_width=10,
                 row_label_precision=4,
                 cols=(),
                 width=16, precision=10):
        self.row_label_width = row_label_width
        self.row_label_precision = row_label_precision
        self.width = width
        self.precision = precision
        self.rows = rows

        if isinstance(cols, str):
            self.cols = (cols,)
        else:
            self.cols = cols

        self.labels = []
        self.data = []

    def format_label(self):
        str = lambda x: (('%%%d.%df' % (self.row_label_width, self.row_label_precision)) % x)
        return " ".join(map(str, self.labels))

    def format_values(self, values):
        str = lambda x: (('%%%d.%df' % (self.width, self.precision)) % x)
        return " ".join(map(str, values))

    def __getitem__(self, value):
        self.labels.append(value)
        return self

    def __setitem__(self, name, value):
        self.labels.append(name)
        label = self.format_label()
        self.labels = []

        #if hasattr(value, "__iter__"):
        self.data.append( (label, list(value) ) )
        #else:
            #self.data.append( (label, list(value) ) )

    def save(self, file):
        import pickle
        pickle_str = pickle.dumps(self)
        fileobj = open(file, "w")
        fileobj.write(str(self))
        fileobj.close()

    def __str__(self):
        rowstr = lambda x: '%%%ds' % self.row_label_width % x
        colstr = lambda x: '%%%ds' % self.width % x

        lines = []

        row_header = " ".join(map(rowstr, self.rows))
        row_header += " ".join(map(colstr, self.cols))

        lines.append(row_header)

        for datarow in self.data:
            row_data = datarow[0]
            row_data += self.format_values(datarow[1])
            lines.append(row_data)

        return "\n".join(lines)

    def copy(self):
        import copy
        return copy.deepcopy(self)

    def absolute_to_relative(self, Factor = 627.51):
        import copy

        if len(self.data) == 0:
            return

        current_min = list(copy.deepcopy(self.data[0][1]))
        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                if current_min[col] > datarow[1][col]:
                    current_min[col] = datarow[1][col]

        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                datarow[1][col] = (datarow[1][col] - current_min[col]) * Factor

