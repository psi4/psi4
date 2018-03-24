from __future__ import division
from psi4.driver.p4util.exceptions import ValidationError
from psi4 import core
from itertools import combinations
import math


def displace_cart(mol, geom, salc_list, i_m, step_size):
    """
    Given a molecule, a geometry, a list of salcs, an iterable of (SALC displacement index, # steps to displace) tuples,
    and a step size, name the geometry according to the salc displacements and displace the geometry accordingly.
    """
    name = ""
    # This for loop and tuple unpacking is why the function can handle an arbitrary number of SALCs.
    # This is the sole motivator for the use of tuples in geom_generator.
    for salc_index, disp_steps in i_m:
        # We get rid of the initial ", " later.
        name += ", Coord: {:d}, Disp: {:d}".format(salc_index, disp_steps)
        for component in salc_list[salc_index]:
            geom.add(0, component.atom, component.xyz,
                     disp_steps * step_size * component.coef / math.sqrt(mol.mass(component.atom)))
    # Now let's get rid of that initial ", ".
    geom.name = name[2:]


def geom_generator(mol, freq_irrep_only, mode):
    """
    Generate geometries for a finite difference calculation and print information about the geometries generated.
    You probably want to use the fd_geoms_1_0, fd_geoms_freq_2_0, or fd_geoms_freq_2_1 convenience functions. 

    mol: The molecule to generate displaced geometries for
    freq_irrep_only: An int specifying the Cotton ordered irrep to take frequencies for.
    mode: A string that specifies the quantities determined at and from the displacements. See msg_dist.
    """

    print_lvl = core.get_option("FINDIF", "PRINT")
    num_pts = core.get_option("FINDIF", "POINTS")
    disp_size = core.get_option("FINDIF", "DISP_SIZE")

    msg_dict = {
        "1_0":
        "energies to determine gradients",
        "2_1":
        "gradients to determine vibrational frequencies and \n"
        "  normal modes.  Resulting frequencies are only valid at stationary points.",
        "2_0":
        "gradients to determine vibrational frequencies and \n"
        "  normal modes.  Resulting frequencies are only valid at stationary points."
    }

    try:
        print_msg = msg_dict[mode]
    except KeyError:
        raise ValidationError("FINDIF: Mode {} not recognized.".format(mode))

    if print_lvl:
        core.print_out("\n-------------------------------------------------------------\n\n"
                       "  Using finite-differences of {:s}\n"
                       "\tGenerating geometries for use with {:d}-point formula.\n"
                       "\tDisplacement size will be {:6.2e}.\n".format(print_msg, num_pts, disp_size))

    # TODO: Replace with a generator from a stencil to a set of points.
    pts_dict = {
        3: {
            "sym": ((-1, ), (1, )),
            "asym": ((-1, ), ),
            "off": ((1, 1), (-1, -1))
        },
        5: {
            "sym": ((-2, ), (-1, ), (1, ), (2, )),
            "asym": ((-2, ), (-1, )),
            "off": ((-1, -2), (-2, -1), (-1, -1), (1, -1), (-1, 1), (1, 1), (2, 1), (1, 2))
        }
    }

    # We won't need this variable for a while, but let's validate num_pts early.
    try:
        disps = pts_dict[num_pts]
    except KeyError:
        raise ValidationError("FINDIF: Invalid number of points!")

    # Get settings for CdSalcList, then get it.
    method_allowed_irreps = 0x1 if mode == "1_0" else 0xFF
    t_project = not core.get_global_option("EXTERN") and not core.get_global_option("PERTURB_H")
    # core.get_local_option returns an int, but CdSalcList expect a bool... Re-cast as bool.
    r_project = t_project and bool(core.get_local_option("FINDIF", "FD_PROJECT"))
    salc_list = core.CdSalcList(mol, method_allowed_irreps, t_project, r_project)

    # Convention: A capital N at the start of variable means "number of"
    Natom = mol.natom()
    Nirrep = salc_list.nirrep()
    Nsalc = salc_list.ncd()

    if print_lvl:
        core.print_out("\tNumber of atoms is {:d}.\n".format(Natom))
        if method_allowed_irreps != 0x1:
            core.print_out("\tNumber of irreps is {:d}.\n".format(Nirrep))
        core.print_out("\tNumber of {!s}SALCs is {:d}.\n".format("" if method_allowed_irreps != 0x1 else "symmetric ",
                                                                 Nsalc))
        core.print_out("\tTranslations projected? {:d}. Rotations projected? {:d}.\n".format(t_project, r_project))

    # Convention: x_pi means x_per_irrep. The ith element is x for irrep i, with Cotton ordering.
    # Here, x is a list of indices of CdSALCs belonging to the specified irrep.
    salc_indices_pi = [[] for irrep in range(Nirrep)]

    # Validate that we have an irrep matching the user-specified irrep, if any.
    try:
        salc_indices_pi[freq_irrep_only]
    except (TypeError, IndexError):
        if freq_irrep_only != -1:
            raise ValidationError("FINDIF: Irrep value not in valid range.")

    # Populate salc_indices_pi for all irreps.
    for i, salc in enumerate(salc_list):
        salc_indices_pi[salc.irrep_index()].append(i)

    # If the method allows more than one irrep, print how the irreps partition the SALCS.
    if print_lvl and method_allowed_irreps != 0x1:
        core.print_out("\tIndex of SALCs per irrep:\n")
        for i in range(Nirrep):
            if print_lvl > 1 or freq_irrep_only in {i, -1}:
                tmp = (" {:d} " * len(salc_indices_pi[i])).format(*salc_indices_pi[i])
                core.print_out("\t {:d} : ".format(i + 1) + tmp + "\n")
        core.print_out("\tNumber of SALCs per irrep:\n")
        for i in range(Nirrep):
            if print_lvl > 1 or freq_irrep_only in {i, -1}:
                core.print_out("\t Irrep {:d}: {:d}\n".format(i + 1, len(salc_indices_pi[i])))

    # Now that we've printed the SALCs, clear any that are not of user-specified symmetry.
    if freq_irrep_only != -1:
        for i in range(Nirrep):
            if i != freq_irrep_only:
                salc_indices_pi[i].clear()

    # Define the number of displacements per irrep.
    Ndisp_pi = []

    for irrep, indices in enumerate(salc_indices_pi):
        Ndisp = len(indices) * len(disps["asym" if irrep != 0 else "sym"])
        if mode == "2_0":
            # Either len(indices) or len(indices)-1 is even, so dividing by two is safe.
            Ndisp += len(indices) * (len(indices) - 1) // 2 * len(disps["off"])
        Ndisp_pi.append(Ndisp)

    # Let's print out the number of geometries, the displacement multiplicity, and the CdSALCs!
    if print_lvl:
        core.print_out("\tNumber of geometries (including reference) is {:d}.\n".format(sum(Ndisp_pi) + 1))
        if method_allowed_irreps != 0x1:
            core.print_out("\tNumber of displacements per irrep:\n")
            for i, ndisp in enumerate(Ndisp_pi, start=1):
                core.print_out("\t  Irrep {:d}: {:d}\n".format(i, ndisp))

    if print_lvl > 1:
        for salc in salc_list:
            salc.print()

    # We can finally start generating displacements.
    ref_geom_temp = mol.geometry()
    ref_geom = ref_geom_temp.clone()
    ref_geom.name = "Reference geometry"

    disp_geoms = []

    def append_geoms(indices, steps):
        """Given a list of indices and a list of steps to displace each, append the corresponding geometry to the list."""
        new_geom = ref_geom.clone()
        # Next, to make this salc/magnitude composite.
        index_steps = zip(indices, steps)
        displace_cart(mol, new_geom, salc_list, index_steps, disp_size)
        disp_geoms.append(new_geom)

    for h in range(Nirrep):
        active_indices = salc_indices_pi[h]

        for index in active_indices:
            # Displace along the diagonal.
            # Remember that the totally symmetric irrep has special displacements.
            for val in disps["sym" if h == 0 else "asym"]:
                append_geoms((index,), val)

        # Hessian from energies? We have off-diagonal displacements to worry about.
        if mode == "2_0":
            # i indexes SALC indices of the current irrep.
            for i, index in enumerate(active_indices):
                for index2 in active_indices[:i]:
                    for val in disps["off"]:
                        append_geoms((index, index2), val)

    disp_geoms.append(ref_geom)

    if print_lvl > 2:
        for disp in disp_geoms:
            disp.print_out()

    if print_lvl > 1:
        core.print_out("\n-------------------------------------------------------------\n")

    return disp_geoms


# The gradient only depends on totally symmetric irreps, so the user shouldn't specify irreps further.
def fd_geoms_1_0(molecule):
    """Generate geometries for gradient by finite difference of energies."""
    return geom_generator(molecule, -1, "1_0")


def fd_geoms_freq_1(molecule, irrep):
    """Generate geometries for hessian by finite difference of gradients."""
    return geom_generator(molecule, irrep, "2_1")


def fd_geoms_freq_0(molecule, irrep):
    """Generate geometries for a hessian by finite difference of energies."""
    return geom_generator(molecule, irrep, "2_0")
