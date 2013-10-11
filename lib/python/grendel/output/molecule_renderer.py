"""
Functions for rendering molecules in general, with out a predetermined backend
"""
from __future__ import print_function, division
from grendel.output.renderer_base import MoleculeRenderer, MoleculeSeriesRenderer
from grendel.output.jmol_interface import JmolRenderer, JmolSeriesRenderer
from grendel.chemistry.molecule import Molecule

default_renderer = JmolRenderer
default_series_renderer = JmolSeriesRenderer
default_width = 500
default_height = 500

def _get_renderer(molecule, renderer_type):
    if renderer_type is None:
        renderer_type = default_renderer
    #----------------------------------------#
    if not (isinstance(renderer_type, type) and issubclass(renderer_type, MoleculeRenderer)):
        renderer_type = str(renderer_type)
        if default_renderer.lower() in ["jmol", "jmolrenderer"]:
            renderer_type = JmolRenderer
        else:
            raise ValueError("Unknown renderer '{}'".format(renderer_type))
    #----------------------------------------#
    return renderer_type(molecule)

def _get_series_renderer(molecules, renderer_type):
    if renderer_type is None:
        renderer_type = default_series_renderer
    #----------------------------------------#
    if not (isinstance(renderer_type, type) and issubclass(renderer_type, MoleculeSeriesRenderer)):
        renderer_type = str(renderer_type)
        if default_renderer.lower() in ["jmol", "jmolrenderer", "jmolseriesrenderer"]:
            renderer_type = JmolSeriesRenderer
        else:
            raise ValueError("Unknown renderer '{}'".format(renderer_type))
    #----------------------------------------#
    return renderer_type(molecules)

def write_png(molecule, filename,
        width=None,
        height=None,
        renderer_type=None,
        **kwargs
):
    if width is None:
        width = default_width
    if height is None:
        height = default_height
    #----------------------------------------#
    renderer = _get_renderer(molecule, renderer_type)
    #----------------------------------------#
    return renderer.write_png(
        filename=filename,
        width=width,
        height=height,
        **kwargs
    )

def write_jpg(molecule, filename,
        width=None,
        height=None,
        renderer_type=None,
        **kwargs
):
    if width is None:
        width = default_width
    if height is None:
        height = default_height
    #----------------------------------------#
    renderer = _get_renderer(molecule, renderer_type)
    #----------------------------------------#
    return renderer.write_jpg(
        filename=filename,
        width=width,
        height=height,
        **kwargs
    )

def write_3d_pdf(molecule, filename,
        renderer_type=None,
        quiet=False,
        **kwargs
):
    if isinstance(molecule, Molecule):
        renderer = _get_renderer(molecule, renderer_type)
        return renderer.write_3d_pdf(filename, quiet=quiet, **kwargs)
    else:
        renderer = _get_series_renderer(molecule, renderer_type)
        return renderer.write_3d_pdf(filename, quiet=quiet, **kwargs)

