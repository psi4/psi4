# FAILED ATTEMPT AT REFORMING MOLECULESTUB
# Perhaps I'll come back to this later...
#from grendel import type_checking_enabled, sanity_checking_enabled
#from grendel.chemistry.atom import Atom
#from grendel.gmath import magnitude, angle_between_vectors
#from grendel.gmath.matrix import Matrix
#from grendel.util.decorators import with_flexible_arguments, typechecked, IterableOf
#from grendel.util.exceptions import ChemistryError
#from grendel.util.overloading import overloaded, OverloadedFunctionCallError
#from grendel.util.strings import indented
#from grendel.util.units import strip_units, DistanceUnit, AngularUnit, Radians, Degrees, Angstroms, isunit


## Immutable "friend class" of Molecule
#class MoleculeStub(object):
#    """
#    Immutable "friend" class of `Molecule`, used for hashing.
#    """
#
#    ####################
#    # Class Attributes #
#    ####################
#
#    eq_precision = 8
#    same_internal_tol = {AngularUnit: 0.0001*Degrees, DistanceUnit: 1e-7*Angstroms}
#
#    ##############
#    # Attributes #
#    ##############
#
#    multiplicity = None
#    """ The multiplicity of the electronic state of the molecule.
#      (i.e. 2S+1 where S is the total spin).  Defaults to singlet. """
#
#    charge = None
#    """ The charge on the molecule.  Defaults to neutral. """
#
#    reoriented_matrix = None
#
#    ######################
#    # Private Attributes #
#    ######################
#
#    _hash = None
#    _cartesian_representation = None
#    _cartesian_units = None
#    _internal_representation = None
#    _from_molecule = None
#    _xyz = None
#    _element_list = None
#
#    ##################
#    # Initialization #
#    ##################
#
#    @overloaded
#    def __init__(self, *args, **kwargs):
#        raise OverloadedFunctionCallError
#
#    @__init__.overload_with(
#        atoms=IterableOf('Atom'),
#    )
#    def __init__(self,
#            atoms,
#            **kwargs):
#        self.__init__(
#            [(atom.element, atom.isotope) for atom in atoms],
#            Matrix([atom.position for atom in atoms]),
#            **kwargs)
#
#    @__init__.overload_with(
#        cartesian_units=isunit,
#        charge=(int, None),
#        multiplicity=(int, None)
#    )
#    def __init__(self,
#            elements_and_isotopes,
#            xyz,
#            cartesian_units=DistanceUnit.default,
#            charge=None,
#            multiplicity=None):
#        self._cartesian_units = cartesian_units
#        self.charge = charge if charge is not None else Molecule.default_charge
#        self.multiplicity = multiplicity if multiplicity is not None else Molecule.default_multiplicity
#        self._xyz = xyz
#        self._element_list = elements_and_isotopes
#        # TODO strip units
#        tmpmol = Molecule(
#            [Atom(el, iso, pos) for (el, iso), pos in zip(self._element_list, self._xyz.iter_rows)],
#            charge=self.charge,
#            multiplicity=self.multiplicity
#        )
#        self.reoriented_matrix = tmpmol.reoriented().xyz
#
#    ###################
#    # Special Methods #
#    ###################
#
#    def __hash__(self):
#        if self._hash is not None:
#            return self._hash
#        self._hash = MoleculeDict.hash_for(self)
#        return self._hash
#
#    def __eq__(self, other):
#        if isinstance(other, MoleculeStub):
#            if [a.isotope for a in self] != [a.isotope for a in other]:
#                return False
#            elif (self.multiplicity, self.charge) != (other.multiplicity, other.charge):
#                return False
#            else:
#                reoriented = self.reoriented_matrix * self._cartesian_units.to(DistanceUnit.default)
#                rounded = [round(v, MoleculeStub.eq_precision) for v in reoriented.ravel()]
#                other_oriented = other.reoriented_matrix * other._cartesian_units.to(DistanceUnit.default)
#                other_rounded = [round(v, MoleculeStub.eq_precision) for v in other_oriented.ravel()]
#                return rounded == other_rounded
#        else:
#            return NotImplemented
#
#
#    ###########
#    # Methods #
#    ###########
#
#    # TODO document this!
#    # TODO class variables for default tolerances
#    def is_valid_stub_for(self, other, cart_tol=None, internal_tol=None, ang_tol=None):
#        """
#        """
#        # cart_tol is used for comparison between cartesian positions
#        cart_tol = cart_tol or 1e-8*Angstroms
#        # internal_tol should be unitless, since the difference between internal coordinates could have multiple
#        #   units, and we're taking the magnitude across these units
#        internal_tol = internal_tol or MoleculeStub.same_internal_tol
#        # ang_tol is used when comparing cartesian geometries.  If the angle between two corresponding atoms
#        #   in self and other differs from the angle between the first two corresponding atoms in self and other
#        #   by more than ang_tol, we assume they do not have the same geometry and thus return False
#        ang_tol = ang_tol or 1e-5*Degrees
#        #--------------------------------------------------------------------------------#
#        if type_checking_enabled:
#            if not isinstance(other, Molecule):
#                raise TypeError
#            if isinstance(other, MoleculeStub):
#                raise TypeError
#        #--------------------------------------------------------------------------------#
#        if self.multiplicity != other.multiplicity:
#            return False
#        elif self.charge != other.charge:
#            return False
#        else:
#            if len(self._element_list) != other.natoms:
#                for num, (element, isotope) in enumerate(self._element_list):
#                    if (other[num].element, other[num].isotope) != (element, isotope):
#                        return False
#        if other.natoms <= 1:
#            # if we have 1 or 0 atoms and we've gotten this far, we have a match
#            return True
#        #--------------------------------------------------------------------------------#
#        # no failures yet, so we have to compare geometries
#        # if self has an internal_representation, use it
#        if self._internal_representation is not None:
#            diff = self._internal_representation.values - self._internal_representation.values_for_molecule(other)
#            if any(abs(d) > internal_tol[c.units.genre].in_units(c.units) for d, c in zip(diff, self._internal_representation)):
#                return False
#            return True
#        # if mol has an internal representation, use it:
#        elif other.internal_representation is not None:
#            diff = other.internal_representation.values - other.internal_representation.values_for_molecule(self)
#            if any(abs(d) > internal_tol[c.units.genre].in_units(c.units) for d, c in zip(diff, other.internal_representation)):
#                return False
#            return True
#        else:
#            # They're both fully cartesian.  This could take a while...
#            # We should first try to short-circuit as many ways as possible
#            #----------------------------------------#
#            # first strip units and store stripped versions to speed up the rest of the work
#            # strip the units off of ang_tol
#            ang_tol = strip_units(ang_tol, Radians)
#            # strip the units off of cart_tol
#            cart_tol = strip_units(cart_tol, self._cartesian_units)
#            # make a list of positions with stripped units, since we'll use it up to three times
#            stripped = [strip_units(atom, self._cartesian_units) for atom in (self if self.is_centered() else self.recentered()) ]
#            other_stripped = [strip_units(atom, self._cartesian_units) for atom in (other if other.is_centered() else other.recentered())]
#            #----------------------------------------#
#            # Try to short-circuit negatively by looking for an inconsistancy in the angles
#            #   between pairs of corresponding atoms
#            # If the first atom is at the origin, use the second one
#            offset = 1
#            if stripped[0].is_zero():
#                if magnitude(stripped[0] - other_stripped[0]) > cart_tol:
#                    return False
#                else:
#                    if sanity_checking_enabled and stripped[1].is_zero():
#                        raise ChemistryError, "FrozenMolecule:\n{}\nhas two atoms on top of each other.".format(indented(str(self)))
#                    if other_stripped[1].is_zero():
#                        return False
#                    else:
#                        offset = 2
#                        first_ang = angle_between_vectors(self.atoms[1].pos, other.atoms[1].pos)
#            else:
#                if other_stripped[0].is_zero():
#                    return False
#                else:
#                    first_ang = angle_between_vectors(self.atoms[0].pos, other.atoms[0].pos)
#            for apos, opos in zip(stripped[offset:], other_stripped[offset:]):
#                if apos.is_zero():
#                    if magnitude(apos - opos) > cart_tol:
#                        return False
#                elif opos.is_zero():
#                    # Try again, since Tensor.zero_cutoff could smaller than cart_tol, causing a false zero
#                    if magnitude(apos - opos) > cart_tol:
#                        return False
#                else:
#                    ang = angle_between_vectors(apos, opos)
#                    if abs(ang - first_ang) > ang_tol:
#                        return False
#                    # Also, the magnitude of the distance from the center of mass should be the same:
#                    if abs(apos.magnitude() - opos.magnitude()) > cart_tol:
#                        return False
#            #----------------------------------------#
#            # Try to short-circuit positively:
#            exact_match = True
#            for apos, opos in zip(stripped, other_stripped):
#                if magnitude(apos - opos) > cart_tol:
#                    exact_match = False
#                    break
#            if exact_match:
#                return True
#            exact_match = True
#            # Check negative version
#            for apos, opos in zip(stripped, other_stripped):
#                if magnitude(apos + opos) > cart_tol:
#                    exact_match = False
#                    break
#            if exact_match:
#                return True
#            #----------------------------------------#
#            # We can't short-circuit, so this is the only means we have left
#            # It's far more expensive than the rest, but it always works.
#            return self.has_same_geometry(other, cart_tol)
#
#
######################
## Dependent Imports #
######################
#
#from grendel.chemistry.molecule_dict import MoleculeDict
#from grendel.chemistry.molecule import Molecule

