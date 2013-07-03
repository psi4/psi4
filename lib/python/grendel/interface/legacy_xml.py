""" Legacy support for parsing results from Grendel++ data.xml files
"""
from __future__ import print_function
from collections import defaultdict
import sys

from xml.etree import ElementTree

from grendel import sanity_checking_enabled
from grendel.chemistry.molecular_properties import Energy, MolecularProperty
from grendel.chemistry.molecule_dict import MoleculeDict
from grendel.chemistry.molecule import MoleculeStub
from grendel.interface.computation_details import ComputationDetails
from grendel.interface.result_getter import ResultGetter
from grendel.util.metaprogramming import ReadOnlyAttribute
from grendel.util.overloading import listify_args
from grendel.util.strings import indented
from grendel.util.units import *
from grendel.util.units.unit import DistanceUnit

class InvalidLegacyXMLFileError(ValueError): pass

class PropertyUnavailableError(ValueError): pass

class XMLResultGetterError(ValueError): pass

class LegacyXMLResultGetter(ResultGetter):

    ##############
    # Attributes #
    ##############

    files = ReadOnlyAttribute('files')
    properties = None
    properties_for_molecules = None

    ##################
    # Initialization #
    ##################

    def __init__(self, comparison_representation, *files):
        self.started = True
        self._files = listify_args(*files)
        self.properties = MoleculeDict(comparison_representation, default=lambda: [])
        self.properties_for_molecules = defaultdict(lambda: [])
        for file in files:
            self._parse_file(file)

    ###################
    # Private Methods #
    ###################

    def _parse_file(self, file):
        def _get_at_least_one(parent, tag, dispnum):
            ret_val = parent.findall(tag)
            if sanity_checking_enabled:
                if len(ret_val) == 0:
                    raise InvalidLegacyXMLFileError("missing {} section "
                                                    "for displacement number {}".format(tag, dispnum))
            return ret_val
        def _get_exactly_one(parent, tag, dispnum):
            ret_val = _get_at_least_one(parent, tag, dispnum)
            if sanity_checking_enabled:
                if len(ret_val) > 1:
                    raise InvalidLegacyXMLFileError("multiple {} sections "
                                                    "for displacement number {}".format(tag, dispnum))
            return ret_val[0]
        #========================================#
        etr = ElementTree.parse(file)
        for disp in etr.iter('displacement'):
            disp_number = disp.get('number', '<unnumbered>')
            # Get the molecule part
            mol_sect = _get_exactly_one(disp, 'molecule', disp_number)
            # Get the XYZ section
            xyz_sect = _get_exactly_one(mol_sect, 'xyz', disp_number)
            if 'units' in xyz_sect.keys():
                unitstr = xyz_sect.get('units')
                units = eval(unitstr.title(), globals())
            else:
                units = DistanceUnit.default
            energy_el = _get_at_least_one(disp, 'energy', disp_number)
            # for now just use the "molecular" energy
            energy_el = [e for e in energy_el if e.get('type', '') == 'molecular']
            if len(energy_el) == 0:
                raise InvalidLegacyXMLFileError("missing energy with type='molecular' "
                                                "for displacement number {}".format(disp_number))
            elif len(energy_el) > 1:
                raise InvalidLegacyXMLFileError("multiple energy elements with type='molecular' "
                                                "for displacement number {}".format(disp_number))
            energy_el = energy_el[0]
            if 'units' in energy_el.keys():
                unitstr = energy_el.get('units')
                energy_units = eval(unitstr.title(), globals())
            else:
                energy_units = Hartrees
            energy_val = float(energy_el.get('value')) * energy_units
            mol_stub = MoleculeStub(xyz_sect.text, units=units)
            energy = Energy(mol_stub, units=energy_units, details=ComputationDetails(type='molecular'))
            energy.value =  energy_val
            self.properties[mol_stub].append(energy)

    def can_get_property_for_molecule(self, molecule, property, details=None):
        return self.has_property_for_molecule(molecule, property, details)

    def has_property_for_molecule(self, molecule, property, details=None, verbose=True):
        if molecule in self.properties:
            props = self.properties[molecule]
            for p in props:
                pcopy = copy(p)
                pcopy.molecule = molecule
                self.properties_for_molecules[molecule].append(pcopy)
                if MolecularProperty.is_same_property(property, p):
                    if ComputationDetails.is_compatible_details(details, p.details):
                        return True
        return False

    def get_property_for_molecule(self, molecule, property, details=None):
        for p in self.properties_for_molecules[molecule]:
            if MolecularProperty.is_same_property(property, p):
                if ComputationDetails.is_compatible_details(details, p.details):
                    return p
        props = self.properties[molecule]
        for p in props:
            pcopy = copy(p)
            pcopy.molecule = molecule
            self.properties_for_molecules[molecule].append(pcopy)
            if MolecularProperty.is_same_property(property, p):
                if ComputationDetails.is_compatible_details(details, p.details):
                    return pcopy
        raise PropertyUnavailableError

    ###################
    # Private Methods #
    ###################

