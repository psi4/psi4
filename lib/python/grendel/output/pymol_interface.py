#!/usr/bin/env python

import __main__

__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()

from grendel.output.renderer_base import MoleculeRenderer

class PymolRenderer(MoleculeRenderer):


    def write_png(self, filename, width, height, **kwargs):
        self._render_molecule(**kwargs)
        # TODO finish this
        raise NotImplementedError()

