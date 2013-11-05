"""
Contains the MoleculeRenderer and MoleculeSeriesRenderer base class
"""
from tempfile import mkdtemp
from shutil import rmtree
import atexit
from os.path import join as path_join

_keep_tmp_files = __debug__


class MoleculeRenderer(object):
    """
    Abstract base class for interfaces that render pictures and such of molecules
    """

    def __init__(self, mol):
        """
        @param mol: Molecule to be rendered
        @type mol: Molecule
        """
        self.molecule = mol

    def write_png(self, filename, width, height, **kwargs):
        raise NotImplementedError("MoleculeRenderer subclass '{}' has not implemented write_png()".format(
            self.__class__.__name__
        ))

    @classmethod
    def render_png(cls, mol, *args, **kwargs):
        tmp = cls(mol)
        return tmp.write_png(*args, **kwargs)

    def write_pdf(self, filename, width, height, **kwargs):
        raise NotImplementedError("MoleculeRenderer subclass '{}' has not implemented write_pdf()".format(
            self.__class__.__name__
        ))

    @classmethod
    def render_pdf(cls, mol, *args, **kwargs):
        tmp = cls(mol)
        return tmp.write_pdf(*args, **kwargs)

    def write_jpg(self, filename, width, height, **kwargs):
        raise NotImplementedError("MoleculeRenderer subclass '{}' has not implemented write_jpg()".format(
            self.__class__.__name__
        ))

    @classmethod
    def render_jpg(cls, mol, *args, **kwargs):
        tmp = cls(mol)
        return tmp.write_jpg(*args, **kwargs)

class _TempDirMoleculeRenderer(MoleculeRenderer):

    def __init__(self, *args, **kwargs):
        # Make a temporary directory for feeding things into jmol
        self.tempdir = mkdtemp()
        if not _keep_tmp_files:
            @atexit.register
            def rmtmpdir(tmpdir = self.tempdir):
                rmtree(tmpdir)
        self._geom_path = path_join(self.tempdir, "geom.xyz")
        #----------------------------------------#
        super(_TempDirMoleculeRenderer, self).__init__(*args, **kwargs)


class MoleculeSeriesRenderer(object):
    """
    Abstract base class for rendering movies of molecules and such.
    """

    def __init__(self, mols):
        self.molecules = list(mols)


