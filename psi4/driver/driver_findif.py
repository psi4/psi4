from __future__ import division
from psi4.driver.p4util.exceptions import ValidationError
from psi4.driver.qcdb import libmintsmolecule
from psi4 import core
import numpy as np
#import scipy.linalg as la
from psi4.driver.p4util import block_diag

def _displace_cart(mol, geom, salc_list, i_m, step_size):
    """
    Return a geometry corresponding to the specified displacement along SALCs.

    Parameters
    ----------
    mol : qcdb.molecule or psi4.core.Molecule
        The molecule to displace
    geom : psi4.core.Matrix
        The geometry of the molecule.
    salc_list : psi4.core.CdSalcList
        A list of Cartesian displacement SALCs
    i_m : iterator of tuples
        An iterator containing tuples. Each tuple has the index of a salc in
        salc_list and the number of steps (positive or negative) to displace
        the salc at that index.
    step_size : float
        The size of a single "step," i.e., the stencil size.

    Returns
    ------
    None
    """
    name = ""
    # This for loop and tuple unpacking is why the function can handle
    # an arbitrary number of SALCs.
    for salc_index, disp_steps in i_m:
        # We get rid of the initial ", " later.
        name += ", Coord: {:d}, Disp: {:d}".format(salc_index, disp_steps)
        for component in salc_list[salc_index]:
            geom.add(0, component.atom, component.xyz,
                     disp_steps * step_size * component.coef / np.sqrt(mol.mass(component.atom)))
    # Now let's get rid of that initial ", ".
    geom.name = name[2:]


def _initialize_findif(mol, freq_irrep_only, mode, initialize_string, verbose=0):
    """
    Perform several initialization tasks needed by all primary functions.

    Parameters
    ----------
    mol : qcdb.molecule or psi4.core.Molecule
        The molecule to displace
    freq_irrep_only : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.
    mode : {"1_0", "2_0", "2_1"}
         The first number specifies the derivative level determined from
         displacements, and the second number is the level determined at.
    initialize_string : function
         A function that returns the string to print to show the caller was entered.
         Necessary because the strings depends on values found in this function.
    verbose : int
         Set to 0 to silence extra print information, regardless of the print level.
         Used so the information is printed only during geometry generation, and not
         during the derivative computation as well.

    :returns: *dict* 
    """

    print_lvl = core.get_option("FINDIF", "PRINT")
    num_pts = core.get_option("FINDIF", "POINTS")
    disp_size = core.get_option("FINDIF", "DISP_SIZE")

    data = {"print_lvl": print_lvl, "num_pts": num_pts, "disp_size": disp_size}

    if print_lvl:
        core.print_out(initialize_string(data))

    # Get settings for CdSalcList, then get it.
    method_allowed_irreps = 0x1 if mode == "1_0" else 0xFF
    t_project = not core.get_global_option("EXTERN") and (not core.get_global_option("PERTURB_H"))
    # core.get_option returns an int, but CdSalcList expect a bool, so re-cast
    r_project = t_project and bool(core.get_option("FINDIF", "FD_PROJECT"))
    salc_list = core.CdSalcList(mol, method_allowed_irreps, t_project, r_project)

    # Convention: A capital N at the start of variable means "number of"
    Natom = mol.natom()
    Nirrep = salc_list.nirrep()
    Nsalc = salc_list.ncd()

    if print_lvl and verbose:
        core.print_out("\tNumber of atoms is {:d}.\n".format(Natom))
        if method_allowed_irreps != 0x1:
            core.print_out("\tNumber of irreps is {:d}.\n".format(Nirrep))
        core.print_out("\tNumber of {!s}SALCs is {:d}.\n".format("" if method_allowed_irreps != 0x1 else "symmetric ",
                                                                 Nsalc))
        core.print_out("\tTranslations projected? {:d}. Rotations projected? {:d}.\n".format(t_project, r_project))

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
    if print_lvl and method_allowed_irreps != 0x1 and verbose:
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
    if print_lvl and verbose:
        core.print_out("\tNumber of geometries (including reference) is {:d}.\n".format(sum(Ndisp_pi) + 1))
        if method_allowed_irreps != 0x1:
            core.print_out("\tNumber of displacements per irrep:\n")
            for i, ndisp in enumerate(Ndisp_pi, start=1):
                core.print_out("\t  Irrep {:d}: {:d}\n".format(i, ndisp))

    if print_lvl > 1 and verbose:
        for salc in salc_list:
            salc.print()

    data.update({
        "Ndisp_pi": Ndisp_pi,
        "Nirrep": Nirrep,
        "Nsalc": Nsalc,
        "Natom": Natom,
        "salc_list": salc_list,
        "salc_indices_pi": salc_indices_pi,
        "disps": disps
    })

    return data


def _geom_generator(mol, freq_irrep_only, mode):
    """
    Generate geometries for the specified molecule and derivative levels.
    You probably want to use the geoms_grad_from_energy, fd_geoms_freq_2_0,
    or fd_geoms_freq_2_1 convenience functions. 

    Parameters
    ----------
    mol : qcdb.molecule or psi4.core.Molecule
        The molecule to perform finite difference calculations on.
    freq_irrep_only : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.
    mode : {"1_0", "2_0", "2_1"}
         The first number specifies the derivative level determined from
         displacements, and the second number is the level determined at.

    Returns
    -------
    disp_geoms : list of psi4.core.Matrix
        A list of displaced geometries to compute quantities at.
    """

    # TODO: Switch disp_geoms to list of ndarray?

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

    def init_string(data):
        return ("\n-------------------------------------------------------------\n\n"
                "  Using finite-differences of {:s}.\n"
                "\tGenerating geometries for use with {:d}-point formula.\n"
                "\tDisplacement size will be {:6.2e}.\n".format(print_msg, data["num_pts"], data["disp_size"]))

    data = _initialize_findif(mol, freq_irrep_only, mode, init_string, 1)

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
        _displace_cart(mol, new_geom, data["salc_list"], index_steps, data["disp_size"])
        disp_geoms.append(new_geom)

    for h in range(data["Nirrep"]):
        active_indices = data["salc_indices_pi"][h]

        for index in active_indices:
            # Displace along the diagonal.
            # Remember that the totally symmetric irrep has special displacements.
            for val in data["disps"]["sym" if h == 0 else "asym"]:
                append_geoms((index, ), val)

        # Hessian from energies? We have off-diagonal displacements to worry about.
        if mode == "2_0":
            # i indexes SALC indices of the current irrep.
            for i, index in enumerate(active_indices):
                for index2 in active_indices[:i]:
                    for val in data["disps"]["off"]:
                        append_geoms((index, index2), val)

    disp_geoms.append(ref_geom)

    if data["print_lvl"] > 2:
        for disp in disp_geoms:
            disp.print_out()

    if data["print_lvl"] > 1:
        core.print_out("\n-------------------------------------------------------------\n")

    return disp_geoms


def comp_grad_from_energy(mol, E):
    """Compute the gradient by finite difference of energies.

    Parameters
    ----------
    mol : qcdb.molecule or psi4.core.Molecule
        The molecule to compute the gradient of.
    E : list of floats
        A list of energies of the molecule at displaced geometries.

    Returns
    -------
    gradient : psi4.core.Matrix
        The gradient in Cartesians, as a matrix with dimensions
        number-of-atoms by 3. """

    # TODO: Do we really want to return a Matrix here, instead of a Numpy array?

    def init_string(data):
        return ("  Computing gradient from energies.\n"
                "\tUsing {:d}-point formula.\n"
                "\tEnergy without displacement: {:15.10f}\n"
                "\tCheck energies below for precision!\n"
                "\tForces are for mass-weighted, symmetry-adapted cartesians (in au).\n".format(
                    data["num_pts"], E[-1]))

    data = _initialize_findif(mol, -1, "1_0", init_string)

    Ndisp = sum(data["Ndisp_pi"]) + 1  # +1 for the reference
    if len(E) != Ndisp:
        raise ValidationError("FINDIF: Received {} energies for {} geometries.".format(len(E), Ndisp))

    salc_to_grad = {
        3:
        lambda i: (E[2 * i + 1] - E[2 * i]) / (2.0 * data["disp_size"]),
        5:
        lambda i: (E[4 * i] - 8.0 * E[4 * i + 1] + 8.0 * E[4 * i + 2] - E[4 * i + 3]) / (12.0 * data["disp_size"])
    }

    try:
        g_q = [salc_to_grad[data["num_pts"]](i) for i in range(data["Nsalc"])]
    except KeyError:
        raise ValidationError("FINDIF: {} is an invalid number of points.".format(data["num_pts"]))
    g_q = np.asarray(g_q)

    if data["print_lvl"]:
        max_disp = (data["num_pts"] - 1) // 2  # The numerator had better be divisible by two.
        e_per_salc = 2 * max_disp
        energy_string = ""
        for i in range(1, max_disp + 1):
            energy_string = "Energy(-{})        ".format(i) + energy_string + "Energy(+{})        ".format(i)
        core.print_out("\n\t Coord      " + energy_string + "    Force\n")
        for salc in range(data["Nsalc"]):
            print_str = "\t{:5d}" + " {:17.10f}" * (e_per_salc + 1) + "\n"
            core.print_out(print_str.format(salc, *E[e_per_salc*salc:e_per_salc*(salc+1)], g_q[salc]))
        core.print_out("\n")

    Bmat = data["salc_list"].matrix()
    g_cart = g_q @ Bmat

    gradient_matrix = core.Matrix("F-D gradient", data["Natom"], 3)

    for atom in range(data["Natom"]):
        for xyz in range(3):
            gradient_matrix.set(atom, xyz, g_cart[3 * atom + xyz] * np.sqrt(mol.mass(atom)))

    if core.get_option("FINDIF", "GRADIENT_WRITE"):
        grad = core.GradientWriter(mol, gradient_matrix)
        gradfile = core.get_writer_file_prefix(mol.name()) + ".grad"
        grad.write(gradfile)
        core.print_out("\t Gradient written.\n")

    sgradient = gradient_matrix.clone()

    if data["print_lvl"]:
        core.print_out("\n-------------------------------------------------------------\n")

    return sgradient


def comp_hess_from_grad(mol, G, freq_irrep_only):
    """Compute the Hessian by finite difference of gradients.

    Parameters
    ----------
    mol : qcdb.molecule or psi4.core.Molecule
        The molecule to compute the Hessian of.
    G : list of psi4.core.Matrix
        A list of gradients of the molecule at displaced geometries
    freq_irrep_only
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.

    Returns
    -------
    hessian : np.array
        The hessian in Cartesians, as a matrix with dimensions of
        3*number-of-atoms by 3*number-of-atoms. """

    def init_string(data):
        return ("  Computing second-derivative from gradients using projected, \n"
                "  symmetry-adapted, cartesian coordinates.\n\n"
                "  {:d} gradients passed in, including the reference geometry.\n".format(len(G)))

    data = _initialize_findif(mol, freq_irrep_only, "2_1", init_string)

    Ndisp = sum(data["Ndisp_pi"]) + 1  # +1 for the reference
    if len(G) != Ndisp:
        raise ValidationError("FINDIF: Received {} gradients for {} geometries.".format(len(G), Ndisp))

    # For non-totally symmetric CdSALCs, a symmetry operation can convert + and - displacements.
    # Good News: By taking advantage of that, we (potentially) ran less computations.
    # Bad News: We need to find the - displacements from the + computations now.
    # The next ~200 lines of code are dedicated to that task.
    if data["print_lvl"]:
        core.print_out("  Generating complete list of displacements from unique ones.\n\n")

    pg = mol.point_group()
    ct = pg.char_table()
    order = pg.order()

    # Determine what atoms map to what other atoms under the point group operations.
    # The py-side compute_atom_map will work whether mol is a Py-side or C-side object.
    atom_map = libmintsmolecule.compute_atom_map(mol)
    if data["print_lvl"] >= 3:
        core.print_out("\tThe atom map:\n")
        for atom, sym_image_list in enumerate(atom_map):
            core.print_out("\t {:d} : ".format(atom + 1))
            for image_atom in sym_image_list:
                core.print_out("{:4d}".format(image_atom + 1))
            core.print_out("\n")
        core.print_out("\n")

    # A list of lists of gradients, per irrep
    gradients_pi = []
    # Extract and print the symmetric gradients. This need no additional processing.
    gradients_pi.append([np.array(grad) for grad in G[0:data["Ndisp_pi"][0]]])
    # Asymmetric gradients need additional processing. For future convenience, we discard the symmetric ones.
    G = G[data["Ndisp_pi"][0]:]

    if data["print_lvl"] >= 3:
        core.print_out("\tSymmetric gradients\n")
        for gradient in gradients_pi[0]:
            core.Matrix.from_array(gradient).print_out()
            #gradient.print_out()

    # Asymmetric gradient. There's always SOME operation that transforms a positive
    # into a negative displacement.By doing extra things here, we can find the
    # gradients at the positive displacements.
    for h in range(1, data["Nirrep"]):

        # If there are no CdSALCs in this irrep, let's skip it.
        if not data["Ndisp_pi"][h]:
            gradients_pi.append([])
            continue

        gamma = ct.gamma(h)
        if data["print_lvl"] >= 3:
            core.print_out("Characters for irrep {}\n".format(h))
            for group_op in range(order):
                core.print_out(" {:5.1f}".format(gamma.character(group_op)))
            core.print_out("\n")

        # Find the group operation that converts + to - displacements.
        for group_op in range(order):
            if gamma.character(group_op) == -1:
                break
        else:
            raise ValidationError("A symmetric gradient passed for a non-symmetric one.")
        if data["print_lvl"]:
            core.print_out("\tOperation {} takes plus displacements of irrep {} to minus ones.\n".format(
                group_op + 1, gamma.symbol()))

        so = np.array(ct.symm_operation(group_op).matrix())
        gradients = []

        def recursive_gradients(n):
            """Populate gradients, with step -n, -n+1, ... -1, 1, ... n.
               Positive displacements are computed."""
            gradients.append(np.array(G.pop(0)))
            new_grad = np.zeros((data["Natom"], 3))
            for atom, image in enumerate(atom_map):
                atom2 = image[group_op]
                new_grad[atom2] = np.einsum("xy,y->x", so, gradients[-1][atom])
            if n > 1:
                recursive_gradients(n - 1)
            gradients.append(new_grad)

        max_disp = (data["num_pts"] - 1) // 2  # The numerator had better be divisible by two.
        # i is just for counting the number of times to run recursive_gradients.
        # recursive_gradients "knows" due to the structure of G the number of computed gradients per CdSALC.
        for i in data["salc_indices_pi"][h]:
            recursive_gradients(max_disp)
        gradients_pi.append(gradients)

    # Massweight all gradients.
    # Remember, the atom currently corresponds to our 0 index, hence these transpose tricks.
    massweighter = np.asarray([mol.mass(a) for a in range(data["Natom"])])**(-0.5)
    gradients_pi = [[(grad.T * massweighter).T for grad in gradients] for gradients in gradients_pi]

    if data["print_lvl"] >= 3:
        core.print_out("\tAll mass-weighted gradients\n")
        for gradients in gradients_pi:
            for grad in gradients:
                core.Matrix.from_array(grad).print_out()

    # We have all our gradients generated now!
    # Next, time to get our Hessian.

    H_pi = []
    B_pi = []
    irrep_lbls = mol.irrep_labels()
    massweighter = np.repeat(massweighter, 3)

    for h in range(data["Nirrep"]):
        Ndisp = data["Ndisp_pi"][h]
        Nindices = len(data["salc_indices_pi"][h])
        gradients = gradients_pi[h]

        if not Nindices:
            continue

        # Flatten each gradient, and turn it into a COLUMN of the matrix.
        gradient_matrix = np.array([grad.flatten() for grad in gradients]).T
        # Transform disps from Cartesian to CdSalc coordinates.
        # For future convenience, we transpose.
        # Rows are gradients and columns are coordinates with respect to a particulr CdSALC.
        B_pi.append(data["salc_list"].matrix_irrep(h))
        grads_adapted = (B_pi[-1] @ gradient_matrix).T

        if data["print_lvl"] >= 3:
            core.print_out("Gradients in B-matrix coordinates\n")
            for disp in range(Ndisp):
                core.print_out(" disp {:d}: ".format(disp))
                for salc in grads_adapted[disp]:
                    core.print_out("{:15.10f}".format(salc))
                core.print_out("\n")

        H_pi.append(np.empty([Nindices, Nindices]))

        if data["num_pts"] == 3:
            H_pi[-1] = (grads_adapted[1::2] - grads_adapted[::2]) / (2.0 * data["disp_size"])
        elif data["num_pts"] == 5:
            H_pi[-1] = (grads_adapted[::4] - 8 * grads_adapted[1::4] + 8 * grads_adapted[2::4] - grads_adapted[3::4]
                        ) / (12.0 * data["disp_size"])

        # Symmetrize our Hessian block.
        # The symmetric structure is lost due to errors in the gradients
        temp_hess = core.Matrix.from_array(H_pi[-1])
        temp_hess.hermitize()
        H_pi[-1] = np.array(temp_hess)

        if data["print_lvl"] >= 3:
            core.print_out("\n\tForce Constants for irrep {} in mass-weighted, ".format(irrep_lbls[h]))
            core.print_out("symmetry-adapted cartesian coordinates.\n")
            core.Matrix.from_array(H_pi[-1]).print_out()

        evals, evects = np.linalg.eigh(H_pi[-1])
        # Get our eigenvalues and eigenvectors in descending order.
        idx = evals.argsort()[::-1]
        evals = evals[idx]
        evects = evects[:, idx]

        normal_irr = (B_pi[-1] * massweighter).T @ evects

        if data["print_lvl"] >= 2:
            core.print_out("\n\tNormal coordinates (non-mass-weighted) for irrep {}:\n".format(irrep_lbls[h]))
            core.Matrix.from_array(normal_irr).print_out()

    # We have the Hessian in each irrep! The final task is to perform coordinate transforms.
    #H = la.block_diag(*H_pi)
    H = block_diag(*H_pi)
    B = np.vstack(B_pi)

    if data["print_lvl"] >= 3:
        core.print_out("\n\tForce constant matrix for all computed irreps in mass-weighted SALCS.\n")
        core.Matrix.from_array(H).print_out()

    # Transform the massweighted Hessian from the CdSalc basis to Cartesians.
    # The Hessian is the matrix not of a linear transformation, but of a (symmetric) bilinear form
    # As such, the change of basis is formula A' = Xt A X, no inverses!
    Hx = B.T @ H @ B
    if data["print_lvl"] >= 3:
        core.print_out("\n\tForce constants in mass-weighted Cartesian coordinates.\n")
        core.Matrix.from_array(Hx).print_out()

    # Un-massweight the Hessian.
    Hx = np.transpose(Hx / massweighter) / massweighter

    if data["print_lvl"] >= 3:
        core.print_out("\n\tForce constants in Cartesian coordinates.\n")
        core.Matrix.from_array(Hx).print_out()

    if data["print_lvl"]:
        core.print_out("\n-------------------------------------------------------------\n")

    return Hx


def geoms_grad_from_energy(molecule):
    """
    Generate geometries for a gradient by finite difference of energies.
    
    Parameters
    ----------
    molecule : qcdb.molecule or psi4.core.Molecule
        The molecule to compute the gradient of.

    Returns
    -------
    disp_geoms : list of psi4.core.Matrix
        A list of displaced geometries to compute energies at.

    Notes
    -----
    Only symmetric displacements are necessary, so user specification of
    symmetry is disabled.
    """
    return _geom_generator(molecule, -1, "1_0")


def geoms_hess_from_grad(molecule, irrep):
    """
    Generate geometries for a hessian by finite difference of energies.
    
    Parameters
    ----------
    molecule : qcdb.molecule or psi4.core.Molecule
        The molecule to compute the frequencies of.
    irrep : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.

    Returns
    -------
    disp_geoms : list of psi4.core.Matrix
        A list of displaced geometries to compute gradients at.
    """
    return _geom_generator(molecule, irrep, "2_1")


def fd_geoms_freq_0(molecule, irrep):
    """
    Generate geometries for a hessian by finite difference of energies.
    
    Parameters
    ----------
    molecule : qcdb.molecule or psi4.core.Molecule
        The molecule to compute the frequencies of.
    irrep : int
        The Cotton ordered irrep to get frequencies for. Choose -1 for all
        irreps.

    Returns
    -------
    disp_geoms : list of psi4.core.Matrix
        A list of displaced geometries to compute energies at.
    """
    return _geom_generator(molecule, irrep, "2_0")
