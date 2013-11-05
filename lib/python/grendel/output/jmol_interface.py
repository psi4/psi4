"""
Do stuff with jmol
"""
from __future__ import print_function, division
import atexit
import shlex
from shutil import rmtree, move
import sys
from grendel.output.renderer_base import _TempDirMoleculeRenderer, MoleculeSeriesRenderer
from os.path import join as path_join, exists
from os import remove
from tempfile import mkdtemp, NamedTemporaryFile
from subprocess import Popen, PIPE
from grendel.util.strings import indent
from grendel.util.context_managers import working_directory
from grendel.util.misc import ProgressBar

try:
    from matplotlib.colors import rgb2hex, colorConverter
    _have_mpl_colors = True
except ImportError:
    _have_mpl_colors = False

# TODO Uh... I can do a lot better at portability
_idtfconverter_sh_path = '/Applications/meshlab.app/Contents/plugins/U3D_OSX/'
_have_meshlab = exists(_idtfconverter_sh_path)
_idtf_to_u3d_template = _idtfconverter_sh_path + 'IDTFConverter.sh "' + _idtfconverter_sh_path + '" -en 1 -rzf 0 -pq 500 -input "{input}" -output "{output}"'

_keep_tmp_files = __debug__

# TODO use setuptools or something to make this search for the jmol path and version etc.
jmol_path = "jmol"
pdflatex_path = "pdflatex"
pdflatex_args = ['-interaction', 'nonstopmode']

#TODO Decouple the IDTF to U3D conversion from this
class JmolRenderer(_TempDirMoleculeRenderer):

    default_keywords = dict(
        # black background
        background="[x000000]"
    )

    default_3d_keywords = dict(
        # black background
        background="[xFFFFFF]",
        show_axes=False,
    )

    def __init__(self, *args, **kwargs):
        super(JmolRenderer, self).__init__(*args, **kwargs)


    def _get_kwarg_script_lines(self, **kwargs):
        no_kwarg = "__no_kwarg__"
        rv_lines = []
        #----------------------------------------#
        background_color = kwargs.pop("background", no_kwarg)
        if background_color is not no_kwarg:
            bg = JmolRenderer.Jmol_color(background_color)
            rv_lines.append("background " + bg)
        #----------------------------------------#
        show_axes = kwargs.pop("show_axes", no_kwarg)
        if show_axes is not no_kwarg:
            if show_axes:
                rv_lines.append("axes ON")
            else:
                rv_lines.append("axes OFF")
        #----------------------------------------#
        if len(kwargs) > 0:
            raise TypeError("Unknown render keywords: (\"{}\")".format('", "'.join(kwargs.keys())))
        #- - - - - - - - - - - - - - - - - - - - #
        return "\n".join(rv_lines)


    @staticmethod
    def Jmol_color(in_color):
        if isinstance(in_color, str) and in_color[:2] == "[x" and in_color[-1] == ']':
            # I guess the user knows what they're doing
            return in_color
        #----------------------------------------#
        if _have_mpl_colors:
            color = colorConverter.to_rgb(in_color)
            hex_str = rgb2hex(color)
            hex_str = "[x" + hex_str[1:] + "]"
            return hex_str
        else:
            # I can't help you, sorry
            return in_color


    def write_png(self, filename, width, height,
        compression_amount=0,
        transparent_background=True,
        **kwargs
    ):
        new_kw = JmolRenderer.default_keywords.copy()
        if transparent_background:
            # make the background something unusual to ensure only the background is made transparent, per the jmol manual
            new_kw.update(background="[x101010]")
        # user preferences always have veto power
        new_kw.update(kwargs)
        #----------------------------------------#
        if not 0.0 <= float(compression_amount) <= 10.0:
            raise ValueError("compression_amount must be a number between 0 and 10 (got value of {})".format(
                float(compression_amount)
            ))
        #----------------------------------------#
        self._write_image(
            image_type="PNG" if not transparent_background else "PNGT",
            filename=filename,
            width=width,
            height=height,
            number_param=compression_amount,
            extra_lines=self._get_kwarg_script_lines(**new_kw)
        )


    def write_jpg(self, filename, width, height,
            quality=100,
            **kwargs
    ):
        new_kw = JmolRenderer.default_keywords.copy()
        # user preferences always have veto power
        new_kw.update(kwargs)
        #----------------------------------------#
        if not 1 <= int(quality) <= 100:
            raise ValueError("quality must be an integer between 1 and 100 (got {})".format(
                int(quality)
            ))
        #----------------------------------------#
        self._write_image(
            image_type="JPG",
            filename=filename,
            width=width,
            height=height,
            number_param=quality,
            extra_lines=self._get_kwarg_script_lines(**new_kw)
        )

    def write_3d_pdf(self, filename,
            verbose=False,
            extra_latex_before="",
            extra_latex_after="",
            latex_file_handler=None,
            toolbar=False,
            quiet=False,
            **kwargs
    ):
        self._ensure_tmp_xyz()
        #----------------------------------------#
        if not _have_meshlab:
            raise RuntimeError("You need to have meshlab.app installed (currently it must be installed in /Applications/ to work...) to do make 3D pdfs.  Get it from Sourceforge")
        #----------------------------------------#
        new_kw = JmolRenderer.default_keywords.copy()
        # 3d keywords are higher priority
        new_kw.update(JmolRenderer.default_3d_keywords)
        # user preferences always have veto power
        new_kw.update(kwargs)
        extra_lines = self._get_kwarg_script_lines(**new_kw)
        #----------------------------------------#
        if verbose: print("Generating idtf file with Jmol...")
        jscript = """
            load "{molfile}"
            {extra_lines}
            write "{tempdir}/output.idtf"
        """.format(
            molfile=self._geom_path,
            tempdir=self.tempdir,
            extra_lines=extra_lines
        )
        self._run_jmol_script(jscript)
        idtf_path = path_join(self.tempdir, "output.idtf")
        if not exists(idtf_path):
            raise RuntimeError("Jmol command failed to produce "
                               "file '{}'.  Script was:\n{}".format(
                idtf_path,
                indent(jscript)
            ))
        #----------------------------------------#
        # now convert to u3d format
        # TODO make this portable
        if verbose: print("Converting idtf file to u3d file...")
        idtf_cmd = _idtf_to_u3d_template.format(
            input=path_join(self.tempdir, "output.idtf"),
            output=path_join(self.tempdir, "output.u3d")
        )
        JmolRenderer._safe_run(shlex.split(idtf_cmd), descr="IDTF to u3d conversion")
        #shlex.os.system(idtf_cmd)
        #----------------------------------------#
        # Now embed in a latex document and run PDFLatex
        if verbose: print("Constructing latex document...".format(pdflatex_path, self.tempdir))
        with open(path_join(self.tempdir, "output.idtf.tex")) as f:
            jmol_latex_data = f.read()
        stuff3d = jmol_latex_data.split("toolbar=false,")[1].split("inline=true,")[0].strip()
        if latex_file_handler is None:
            f = open(path_join(self.tempdir, "output.tex"), "w+")
            run_pdflatex = True
        else:
            f = latex_file_handler
            run_pdflatex = False
        if run_pdflatex:
            f.write(pdf_3d_header)
        pdf3d_section = pdf_3d_entry_template\
            .replace("@@@LABEL@@@", "molecule")\
            .replace("@@@3DSTUFF@@@", stuff3d)\
            .replace("@@@FNAME@@@", path_join(self.tempdir, "output.u3d"))
        if toolbar:
            pdf3d_section = pdf3d_section.replace("toolbar=false", "toolbar=true")
        f.write(extra_latex_before)
        f.write(pdf3d_section)
        f.write(extra_latex_after)
        if run_pdflatex:
            f.write(pdf_3d_footer)
        #- - - - - - - - - - - - - - - - - - - - #
        if run_pdflatex:
            if verbose: print("Running '{} {} output' in directory {}...".format(pdflatex_path, " ".join(pdflatex_args), self.tempdir))
            pdflatex_cmd_list = [pdflatex_path] + pdflatex_args + ["output"]
            with working_directory(self.tempdir):
                # Run 3 times to get references right
                JmolRenderer._safe_run(pdflatex_cmd_list, "pdflatex",
                    require_empty_stderr=True, require_zero_exit=False)
                if verbose: print("  ...run 1 complete")
                JmolRenderer._safe_run(pdflatex_cmd_list, "output",
                    require_empty_stderr=True, require_zero_exit=False)
                if verbose: print("  ...run 2 complete")
                JmolRenderer._safe_run(pdflatex_cmd_list, "output",
                    require_empty_stderr=True, require_zero_exit=False)
                if verbose: print("  ...run 3 complete")
            #----------------------------------------#
            # Finally, copy back the temporary output to the requested file path
            #   but only if we wrote the file
            move(path_join(self.tempdir, "output.pdf"), filename)


    def _write_image(self, image_type, filename, width, height, number_param, extra_lines=''):
        if image_type not in ["PNG", "PNGT", "JPG"]:
            raise ValueError("Unsupported image type '{}'".format(image_type))
        #----------------------------------------#
        self._ensure_tmp_xyz()
        #----------------------------------------#
        self._run_jmol_script("""
            load "{molfile}"
            {extra_lines}
            write IMAGE {width} {height} {image_type} {n} "{filename}"
        """.format(
            molfile=self._geom_path,
            width=width, height=height,
            n=number_param,
            filename=filename,
            image_type=image_type,
            extra_lines=extra_lines
        ))
        #----------------------------------------#


    def _run_jmol_script(self, script,
            nodisplay=True,
            background_transparent=True,
            silent_startup=True,
            noconsole=True
    ):
        """
        @return: Returns a tuple of stdout (as a string), stderror (as a string), and returncode
        """
        with NamedTemporaryFile("w+", dir=self.tempdir, delete=False) as f:
            f.write(script)
            name = f.name
        #----------------------------------------#
        cmd_args = [jmol_path]
        cmd_args.extend(['--script', name])
        if nodisplay:
            cmd_args.append("--nodisplay")
        if background_transparent:
            cmd_args.append("--backgroundtransparent")
        if silent_startup:
            cmd_args.append("--silent")
        if noconsole:
            cmd_args.append("--noconsole")
        #----------------------------------------#
        try:
            stdoutstr, stderrstr = JmolRenderer._safe_run(cmd_args)
        except RuntimeError:
            print("Jmol script failed.  Script was:\n{}".format(
                script
            ), file=sys.stderr)
            raise
        #----------------------------------------#
        if not _keep_tmp_files:
            remove(name)
        #----------------------------------------#
        return stdoutstr, stderrstr

    # TODO move this somewhere more appropriate
    @staticmethod
    def _safe_run(cmd_arg_list, descr="External command",
            require_empty_stderr=False, require_zero_exit=True
    ):
        process = Popen(cmd_arg_list,
            stdout=PIPE,
            stderr=PIPE
        )
        stdoutstr, stderrstr = process.communicate()
        return_code = process.returncode
        if (require_zero_exit and return_code != 0) or (require_empty_stderr and len(stderrstr) > 0):
            raise RuntimeError("{} run failed.  Command was:\n{}\nstdout was:\n{}\nstderr was:\n{}\nexit code was: {}".format(
                descr,
                indent(" ".join(cmd_arg_list)),
                indent(stdoutstr),
                indent(stderrstr),
                return_code
            ))
        return stdoutstr, stderrstr


    def _ensure_tmp_xyz(self):
        if not exists(self._geom_path):
            self.molecule.write_xyz(self._geom_path)


class JmolSeriesRenderer(MoleculeSeriesRenderer):

    def __init__(self, *args, **kwargs):
        # Make a temporary directory for feeding things into jmol
        self.tempdir = mkdtemp()
        if not _keep_tmp_files:
            @atexit.register
            def rmtmpdir(tmpdir = self.tempdir):
                rmtree(tmpdir)
        self._geom_path = path_join(self.tempdir, "geom.xyz")
        #----------------------------------------#
        super(JmolSeriesRenderer, self).__init__(*args, **kwargs)

    def write_3d_pdf(self, filename,
            verbose=False,
            quiet=False,
            **kwargs
    ):
        """
        Write interactive 3D pictures of the molecules to a pdf file, one per page.
        @param filename: File to write
        @param kwargs: Keywords to pass to JMolRenderer.write_3d_pdf
        """
        latex_out = path_join(self.tempdir, "output.tex")
        #----------------------------------------#
        pbar = ProgressBar(len(self.molecules))
        if not quiet:
            pbar.update(0)
        with open(latex_out, "w+") as f:
            f.write(pdf_3d_header)
            for mol in self.molecules:
                renderer = JmolRenderer(mol)
                if mol.description is not None:
                    f.write(pdf_description_template.replace("@@@LABEL@@@", mol.description))
                renderer.write_3d_pdf(
                    filename=None, verbose=verbose,
                    latex_file_handler=f,
                    **kwargs
                )
                if mol is not self.molecules[-1]:
                    f.write(r"\clearpage")
                if not quiet:
                    pbar += 1
            f.write(pdf_3d_footer)
        #----------------------------------------#
        pdflatex_cmd_list = [pdflatex_path] + pdflatex_args + ["output"]
        with working_directory(self.tempdir):
            # Run 3 times to get references right
            JmolRenderer._safe_run(pdflatex_cmd_list, "pdflatex",
                require_empty_stderr=True, require_zero_exit=False)
            if verbose: print("  ...run 1 complete")
            JmolRenderer._safe_run(pdflatex_cmd_list, "output",
                require_empty_stderr=True, require_zero_exit=False)
            if verbose: print("  ...run 2 complete")
            JmolRenderer._safe_run(pdflatex_cmd_list, "output",
                require_empty_stderr=True, require_zero_exit=False)
            if verbose: print("  ...run 3 complete")
        #----------------------------------------#
        move(path_join(self.tempdir, "output.pdf"), filename)


pdf_description_template = r'''
\begin{center}{\huge @@@LABEL@@@}\end{center}
'''

pdf_3d_header = r"""
\documentclass[12pt,letter]{article}
\usepackage{hyperref}
\usepackage[3D]{movie15}
\usepackage{verbatim}
\usepackage{fullpage}
\pagestyle{empty}
\begin{document}
"""

pdf_3d_entry_template = r"""
 \begin{center}
  \includemovie[
   label=@@@LABEL@@@,
    autoplay,
    repeat=1,
    toolbar=false,
    @@@3DSTUFF@@@
    inline=true,
  ]{0.95\textwidth}{0.95\textwidth}{@@@FNAME@@@}
\end{center}
"""

pdf_3d_footer = r"""
\end{document}
"""
