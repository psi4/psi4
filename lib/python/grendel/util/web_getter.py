from copy import copy
import cirpy

# TODO CCCBDB data retrieval:
#  Use curl --cookie-jar session --location -o results.html http://cccbdb.nist.gov/guide.asp to get the cookie (replacing session with a logical filename or - for stdout
# Then use the cookie with:  curl -d "formula=<formula>&submit1=Submit" -e http://cccbdb.nist.gov/geom1.asp --location --cookie <cookie contents> -o <outputfile> http://cccbdb.nist.gov/getform.asp
#       to get a list of available results
# Then find the method/basis pair and get the geometry



__all__ = [
    "MoleculeNotFoundError"
]

valid_identifiers = copy(cirpy.available_properties)
input_identifiers = copy(cirpy.available_properties)
retrieve_only_identifiers = [
    'sdf',
    'formula',
    'h_bond_donor_count',
    'h_bond_acceptor_count',
    'h_bond_center_count',
    'rule_of_5_violation_count',
    'rotor_count',
    'effective_rotor_count',
    'ring_count',
    'ringsys_count'
]
input_only_identifiers = [
]
name_types = [
    'opsin_name',
    'cir_name',
    'chemspider_name'
]
input_only_identifiers.extend(name_types)
for i in cirpy.file_formats: input_identifiers.remove(i)
for i in retrieve_only_identifiers: input_identifiers.remove(i)
input_identifiers.extend(input_only_identifiers)

nice_names = {
    'stdinchikey'   : ["Standard InChIKey"],
    'formula'       : ["Chemical Formula"],
    'mw'            : ["Molecular Weight"],
}



# TODO Link in openbabel using Pybel?

class MoleculeNotFoundError(Exception):
    """ Raised when a function within the `web_getter` module cannot find a molecule with the given identifier.
    """
    # TODO: Check for general internet connection (e.g. ping google.com) and raise a different error if this is the problem
    def __init__(self, identifier, value, url):
        if isinstance(identifier, list):
            super(MoleculeNotFoundError, self).__init__("Could not find molecule with a value of " + str(value)
                                                        + " for identifiers [" + ", ".join(identifier) + "]:  The url "
                                                        + str(url) + " returned a 404 error.")
        elif identifier == "__any__":
            super(MoleculeNotFoundError, self).__init__("Could not find molecule with a value of " + str(value)
                                                        + " for any input identifier:  The url "
                                                        + str(url) + " returned a 404 error.")
        else:
            super(MoleculeNotFoundError, self).__init__("Could not find molecule with a value of " + str(value)
                                                    + " for identifier " + str(identifier) + ":  The url "
                                                    + str(url) + " returned a 404 error.")










