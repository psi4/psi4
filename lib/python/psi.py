import PsiMod

def psidatadir(key):
    return PsiMod.Process.environment[key]

#
# Define geometry to be used by PSI4.
# This geometry will be saved to the checkpoint file.
#
# geometry( [[ "O", 0.0, 0.0, 0.0 ],
#            [ "H", 1.0, 0.0, 0.0 ]] )
#
def geometry(geom, reorient = True, prefix = "", chkpt = None, shiftToCOM = True):
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
