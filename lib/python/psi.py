import PsiMod

class ValueNotSet (Exception): pass

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
