import os
from grendel.chemistry.molecular_properties import Energy
from grendel.chemistry.derivative_properties import PropertyDerivative, Hessian, Gradient
from grendel.interface.computation_details import Methods, ComputationDetails, Basis

from grendel.interface.input_generator import TemplateInputGenerator, PropertyGenerator
from grendel.representations.cartesian_representation import CartesianRepresentation


template_dir = os.path.dirname(__file__) + os.sep + "templates" + os.sep

class QChemInputGenerator(TemplateInputGenerator):
    """
    """

    property_generators = [
        PropertyGenerator(
            filename = template_dir + "basic.mako",
            property = [Energy, Gradient],
            details = ComputationDetails(
                method = (
                    # The "dotted keyword" syntax is on it's way out, since it's too confusing
                    Methods.HF, 'HF',
                    Methods.DFT, 'DFT'
                ),
                required_details = [
                    'basis'
                ],
            ),
            unsupported_details = {
                'basis' : Basis.custom_basis,
                'optimization' : True
            }
        ),
        PropertyGenerator(
            filename = template_dir + "basic.mako",
            property = [Energy, Hessian],
            details = ComputationDetails(
                method = (
                    # The "dotted keyword" syntax is on it's way out, since it's too confusing
                    Methods.HF, 'HF',
                    Methods.DFT, 'DFT'
                    ),
                required_details = [
                    'basis'
                ],
            ),
            unsupported_details = {
                'basis' : Basis.custom_basis,
                'optimization' : True
            }
        ),
        PropertyGenerator(
            filename = template_dir + "qcopt.mako",
            property = [Energy],
            details= ComputationDetails(
                method = (
                    'HF', 'DFT'
                ),
                optimization = ( True, ),
                required_details = [
                    'basis'
                ]
            ),
            unsupported_details = {
                'basis' : Basis.custom_basis,
                'optimization' : False
            }
        )
    ]