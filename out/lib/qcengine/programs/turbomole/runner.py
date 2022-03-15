"""
Calls the Turbomole executable.
"""
import os
import re
from decimal import Decimal
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import safe_version, which

from ...util import execute, temporary_directory
from ..model import ProgramHarness
from ..qcvar_identities_resources import build_atomicproperties, build_out
from .define import execute_define, prepare_stdin
from .harvester import harvest
from .methods import KEYWORDS, METHODS


class TurbomoleHarness(ProgramHarness):

    _defaults = {
        "name": "Turbomole",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": False,
        "node_parallel": True,
        "managed_memory": True,
    }

    version_cache: Dict[str, str] = {}

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "define",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via http://www.cosmologic.de/turbomole/home.html",
        )

    def get_version(self) -> str:
        which_prog = which("define")
        if which_prog not in self.version_cache:
            # We use basically a dummy stdin as we dont want to pipe any real
            # input into define. We only want to parse the version number from
            # the string.
            with temporary_directory(suffix="_define_scratch") as tmpdir:
                tmpdir = Path(tmpdir)
                stdout = execute_define("\n", cwd=tmpdir)
            # Tested with V7.3 and V7.4.0
            version_re = re.compile("TURBOMOLE (?:rev\. )?(V.+?)\s+")
            mobj = version_re.search(stdout)
            version = mobj[1]
            self.version_cache[which_prog] = safe_version(version)
        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        # TODO: handle input errors?! But then define probably already crashed...
        # if 'There is an error in the input file' in dexe["stdout"]:
        # raise InputError(dexe["stdout"])

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)

    def sub_control(self, control, pattern, repl, **kwargs):
        control_subbed = re.sub(pattern, repl, control, **kwargs)
        return control_subbed

    def append_control(self, control, to_append):
        return self.sub_control(control, "\$end", f"{to_append}\n$end")

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        turbomolerec = {
            "infiles": {},
            "outfiles": {"control": "control"},
            "scratch_directory": config.scratch_directory,
        }

        # Handle molecule
        # TODO: what's up with moldata? Do I need it?
        coord_str, moldata = input_model.molecule.to_string(dtype="turbomole", return_data=True)

        # Prepare stdin for define call
        model = input_model.model
        # geeopt will hold the for which to calculate the gradient.
        # 'x' corresponds to the ground state, 'a 1' would be the GS too.
        # 'a1 2' would be the 1st excited state of the irreducible group A1.
        # Right now only GS are supported, so this is hardcoded as 'x'.
        geoopt = "x" if input_model.driver.derivative_int() > 0 else ""
        stdin, subs = prepare_stdin(
            model.method,
            model.basis,
            input_model.keywords,
            input_model.molecule.molecular_charge,
            input_model.molecule.molecular_multiplicity,
            geoopt,
        )
        with temporary_directory(suffix="_define_scratch") as tmpdir:
            tmpdir = Path(tmpdir)
            with open(tmpdir / "coord", "w") as handle:
                handle.write(coord_str)
            stdout = execute_define(stdin, cwd=tmpdir)
            # The define scratch will be populated by some files that we want to keep
            to_keep = "basis auxbasis coord control alpha beta mos".split()

            for fn in to_keep:
                full_fn = tmpdir / fn
                if not full_fn.exists():
                    continue
                with open(full_fn) as handle:
                    turbomolerec["infiles"][fn] = handle.read()

        env = os.environ.copy()
        env["PARA_ARCH"] = "SMP"
        env["PARNODES"] = str(config.ncores)
        env["SMPCPUS"] = str(config.ncores)
        turbomolerec["environment"] = env
        # Memory is set in the control file

        keywords = input_model.keywords

        ########################
        # DETERMINE SOME FLAGS #
        ########################

        ri_calculation = any([keywords.get(ri_kw, False) for ri_kw in KEYWORDS["ri"]])
        ricc2_calculation = model.method in METHODS["ricc2"]

        ###################
        # MEMORY HANDLING #
        ###################

        # Central file that controls Turbomole. We assign it here to the "control"
        # variable as we may need to modify it, e.g. for a Hessian calculation.
        control = turbomolerec["infiles"]["control"]

        # Calculate total available memory in MB
        mem_mb = config.memory * (1024 ** 3) / 1e6
        ri_fraction = 0.25
        # Total amount of memory allocated to ricore
        ricore = 0
        if ri_calculation:
            # This is the default given by Turbomole
            ricore = mem_mb * ri_fraction
            ri_per_core = int(ricore / config.ncores)
            # Update $ricore entry in the control file
            control = self.sub_control(control, "\$ricore\s+(\d+)", f"$ricore {ri_per_core} MiB per_core")
        # Calculate remaining memory
        maxcor = mem_mb - ricore
        assert maxcor > 0, "Not enough memory for maxcor! Need {-maxcor} MB more!"

        # maxcore per_core
        per_core = int(maxcor / config.ncores)
        # Update $maxcor entry in the control file
        control = self.sub_control(control, "\$maxcor\s+(\d+)\s+MiB\s+per_core", f"$maxcor {per_core} MiB per_core")

        ############################
        # DETERMINE SHELL COMMANDS #
        ############################

        # ----------------------#
        # | Energy calculations |
        # ----------------------#

        # Set appropriate commands. We always need a reference wavefunction
        # so the first command will be dscf or ridft to converge the SCF.
        commands = ["ridft"] if ri_calculation else ["dscf"]

        # ------------------------#
        # | Gradient calculations |
        # ------------------------#

        # Keep the gradient file for parsing
        if input_model.driver.derivative_int() == 1:
            turbomolerec["outfiles"]["gradient"] = "gradient"

        # ricc2 will also calculate the gradient. But this requires setting
        # 'geoopt (state)' in the control file. This is currently handled in the
        # 'define' call.
        if ricc2_calculation:
            commands.append("ricc2")
        # Gradient calculation for DFT/HF
        elif input_model.driver.derivative_int() == 1:
            grad_command = "rdgrad" if ri_calculation else "grad"
            commands.append(grad_command)

        # -----------------------#
        # | Hessian calculations |
        # -----------------------#

        if input_model.driver.derivative_int() == 2:
            freq_command = "NumForce -level cc2" if ricc2_calculation else "aoforce"
            # NumForce seems to ignore the nprhessian command and will always
            # write to hessian
            hessian_outfile = "hessian" if ricc2_calculation else "nprhessian"
            commands.append(freq_command)
            # Add some keywords to the control file
            #   noproj: Don't project out translation and rotation
            #   nprhessian: Set filename of un-projected hessian
            control = self.append_control(control, "$noproj\n$nprhessian file=nprhessian")
            turbomolerec["outfiles"][hessian_outfile] = None

        # Build the full shell command and set it
        command = ["; ".join(commands)]
        turbomolerec["command"] = command
        # Re-assign the potentially modified control file, e.g. for a Hessian calculation
        turbomolerec["infiles"]["control"] = control

        # TODO: check if the chosen commands are available with which()?

        return turbomolerec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            inputs["outfiles"],
            shell=True,
            # TODO: scratch_messy?
            # scratch_messy=False,
        )
        return success, dexe

    def parse_output(
        self, outfiles: Dict[str, str], input_model: "AtomicInput"
    ) -> "AtomicResult":  # lgtm: [py/similar-function]

        stdout = outfiles.pop("stdout")

        qcvars, gradient, hessian = harvest(input_model.molecule, stdout, **outfiles)

        if gradient is not None:
            qcvars["CURRENT GRADIENT"] = gradient

        if hessian is not None:
            qcvars["CURRENT HESSIAN"] = hessian

        retres = qcvars[f"CURRENT {input_model.driver.upper()}"]
        if isinstance(retres, Decimal):
            retres = float(retres)

        build_out(qcvars)
        atprop = build_atomicproperties(qcvars)

        output_data = input_model.dict()
        output_data["extras"]["outfiles"] = outfiles
        output_data["properties"] = atprop
        output_data["provenance"] = Provenance(creator="Turbomole", version=self.get_version(), routine="turbomole")
        output_data["return_result"] = retres
        output_data["stdout"] = stdout
        output_data["success"] = True

        return AtomicResult(**output_data)
