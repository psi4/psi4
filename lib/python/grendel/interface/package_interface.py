from grendel.interface.computation import Computation
from grendel.interface.computation_details import ComputationDetails

class PackageInterface(object):
    """ Container for the various custom classes needed to run a calculation in a given software package.
    """

    ##############
    # Attributes #
    ##############

    # Things that have no logical defaults that could be used for all packages...
    input_generator = None
    runner = None
    output_parser = None

    # Things that _can_ be overridden, but need not be in most cases...
    computation_class = Computation
    details_class = ComputationDetails


    ##################
    # Initialization #
    ##################

    def __init__(self, input_generator, runner, output_parser):
        self.input_generator = input_generator
        self.runner = runner
        self.output_parser = output_parser

    @classmethod
    def load_package(cls, name,
                     input_generator=None,
                     output_parser=None,
                     runner=None,
                     details_class=None,
                     computation_class=None):
        """ Load a package with the given name to be used for computations.  Optionally, specify specific subclasses
        to use for certain things.
        """
        if name.upper() == "Q-CHEM" or name.upper() == "QCHEM":
            import qchem
            ret_val = PackageInterface(
                qchem.QChemInputGenerator,
                qchem.LocalQChemRunner,
                qchem.QChemOutputParser
            )
        else:
            raise ValueError("Unknown computational package '{0}'".format(name))
        #--------------------------------------------------------------------------------#
        if input_generator is not None:
            ret_val.input_generator = input_generator
        if runner is not None:
            ret_val.runner = runner
        if output_parser is not None:
            ret_val.output_parser = runner
        if details_class is not None:
            ret_val.details_class = details_class
        if computation_class is not None:
            ret_val.computation_class = computation_class
        #--------------------------------------------------------------------------------#
        return ret_val





