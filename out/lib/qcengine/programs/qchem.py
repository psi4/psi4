"""
Calls the Q-Chem executable.
"""

import os
import re
import tempfile
import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from qcelemental import constants
from qcelemental.models import AtomicInput, AtomicResult, Molecule, Provenance
from qcelemental.molparse import regex
from qcelemental.util import parse_version, safe_version, which

from qcengine.config import TaskConfig, get_config

from ..exceptions import InputError, UnknownError
from ..util import disk_files, execute, temporary_directory
from .model import ProgramHarness

NUMBER = r"(?x:" + regex.NUMBER + ")"


class QChemHarness(ProgramHarness):
    _defaults: Dict[str, Any] = {
        "name": "QChem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        return which(
            "qchem",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install by visiting the Q-Chem website. Check it's in your PATH with `which qchem'.",
        )

    def _get_qc_path(self, config: Optional["TaskConfig"] = None):
        paths = os.environ.copy()
        paths["QCSCRATCH"] = tempfile.gettempdir()
        if config and config.scratch_directory:
            paths["QCSCRATCH"] = config.scratch_directory

        # Nothing else to pass in
        if "QC" in os.environ:
            return paths

        # Assume QC path is
        qchem_path = which("qchem")
        if qchem_path:
            paths["QC"] = os.path.dirname(os.path.dirname(qchem_path))

        return paths

    def get_version(self) -> str:
        # excuse missing Q-Chem for QCEngineRecords tests
        found = self.found()
        if not found:
            return safe_version("0.0.0")

        self.found(raise_error=True)

        # Get the node configuration
        config = get_config()

        which_prog = which("qchem")
        if which_prog not in self.version_cache:
            success, exc = execute(
                [which_prog, "v.in"],
                {"v.in": "$rem\n"},
                scratch_directory=config.scratch_directory,
                environment=self._get_qc_path(),
                timeout=15,
            )

            mobj = re.search(r"Q-Chem\s+([\d.]+)\s+for", exc["stdout"])
            if not mobj:
                mobj = re.search(r"Q-Chem version:\s+([\d.]+)\s+", exc["stdout"])

            if mobj:
                self.version_cache[which_prog] = safe_version(mobj.group(1))

            # if "QC not defined" in exc["stdout"]:
            else:
                return safe_version("0.0.0")

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: TaskConfig) -> "AtomicResult":
        """
        Run qchem
        """
        # Check if qchem executable is found
        self.found(raise_error=True)

        # Check qchem version
        qceng_ver = "5.1"
        if parse_version(self.get_version()) < parse_version(qceng_ver):
            raise TypeError(f"Q-Chem version <{qceng_ver} not supported (found version {self.get_version()})")

        # Setup the job
        job_inputs = self.build_input(input_model, config)

        # Run qchem
        exe_success, proc = self.execute(job_inputs)

        # Determine whether the calculation succeeded
        if exe_success:
            # If execution succeeded, collect results
            result = self.parse_output(proc["outfiles"], input_model)
            return result
        else:
            outfile = proc["outfiles"]["dispatch.out"]
            if "fatal error occurred in module qparser" in outfile:
                raise InputError(proc["outfiles"]["dispatch.out"])
            else:
                # Return UnknownError for error propagation
                raise UnknownError(proc["outfiles"]["dispatch.out"])

    def execute(
        self,
        inputs: Dict[str, Any],
        *,
        extra_infiles: Optional[Dict[str, str]] = None,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        scratch_messy: bool = False,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:
        """
        For option documentation go look at qcengine/util.execute
        """

        # Collect all input files and update with extra_infiles
        infiles = inputs["infiles"]
        if extra_infiles is not None:
            infiles.update(extra_infiles)

        binary_files = [os.path.join("savepath", x) for x in ["99.0", "131.0", "132.0"]]

        # Collect all output files and extend with with extra_outfiles
        outfiles = ["dispatch.out"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        # Replace commands with extra_commands if present
        commands = inputs["commands"] + ["savepath"]
        if extra_commands is not None:
            commands = extra_commands

        envs = self._get_qc_path()

        with temporary_directory(parent=inputs["scratch_directory"], suffix="_qchem_scratch") as tmpdir:
            envs["QCSCRATCH"] = tmpdir
            bdict = {x: None for x in binary_files}

            with disk_files({}, bdict, cwd=tmpdir, as_binary=binary_files):
                exe_success, proc = execute(
                    commands,
                    infiles=infiles,
                    outfiles=outfiles,
                    scratch_name=scratch_name,
                    scratch_directory=tmpdir,
                    scratch_messy=scratch_messy,
                    timeout=timeout,
                    environment=envs,
                )

            proc["outfiles"].update({os.path.split(k)[-1]: v for k, v in bdict.items()})

        if (proc["outfiles"]["dispatch.out"] is None) or (
            "Thank you very much for using Q-Chem" not in proc["outfiles"]["dispatch.out"]
        ):
            exe_success = False

        # QChem does not create an output file and only prints to stdout
        return exe_success, proc

    def build_input(
        self, input_model: "AtomicInput", config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:

        # Check some bounds on what cannot be parsed
        if "ccsd" in input_model.model.method.lower() or "ccd" in input_model.model.method.lower():
            raise InputError("Cannot handle CC* methods currently.")

        # Build keywords
        keywords = {k.upper(): v for k, v in input_model.keywords.items()}
        keywords["INPUT_BOHR"] = "TRUE"
        keywords["MEM_TOTAL"] = str(int(config.memory * 1024))  # In MB

        if input_model.driver == "energy":
            keywords["JOBTYPE"] = "sp"
        elif input_model.driver == "gradient":
            keywords["JOBTYPE"] = "force"
        elif input_model.driver == "hessian":
            keywords["JOBTYPE"] = "freq"
        else:
            raise InputError(f"Driver {input_model.driver} not implemented for Q-Chem.")

        if input_model.molecule.fix_com or input_model.molecule.fix_orientation:
            keywords["SYM_IGNORE"] = "TRUE"

        keywords["METHOD"] = input_model.model.method
        if input_model.model.basis:
            keywords["BASIS"] = input_model.model.basis

        # Begin the input file
        input_file = []
        input_file.append(
            f"""$comment
Automatically generated Q-Chem input file by QCEngine
$end
            """
        )

        # Add Molecule, TODO: Add to QCElemental
        mol = input_model.molecule
        input_file.append("$molecule")
        input_file.append(f"""{int(mol.molecular_charge)} {mol.molecular_multiplicity}""")

        for real, sym, geom in zip(mol.real, mol.symbols, mol.geometry):
            if real is False:
                raise InputError("Cannot handle ghost atoms yet.")
            input_file.append(f"{sym} {geom[0]:14.8f} {geom[1]:14.8f} {geom[2]:14.8f}")

        input_file.append("$end\n")

        # Write out the keywords
        input_file.append("$rem")
        for k, v in keywords.items():
            input_file.append(f"{k:20s}  {v}")
        input_file.append("$end\n")

        ret = {
            "infiles": {"dispatch.in": "\n".join(input_file)},
            "commands": [which("qchem"), "-nt", str(config.ncores), "dispatch.in", "dispatch.out"],
            "scratch_directory": config.scratch_directory,
        }

        return ret

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> AtomicResult:

        output_data = {}

        bdata = {}
        outtext = ""
        for k, v in outfiles.items():
            if k == "dispatch.out":
                outtext = v
                continue
            if v is None:
                continue
            bdata[k] = np.frombuffer(v)

        if input_model.driver == "energy":
            output_data["return_result"] = bdata["99.0"][-1]
        elif input_model.driver == "gradient":
            output_data["return_result"] = bdata["131.0"]
        elif input_model.driver == "hessian":
            output_data["return_result"] = bdata["132.0"]
        else:
            raise ValueError(f"Could not parse driver of type {input_model.driver}.")

        properties = {
            "nuclear_repulsion_energy": bdata["99.0"][0],
            "scf_total_energy": bdata["99.0"][1],
            "return_energy": bdata["99.0"][-1],
        }

        qcvars = {}
        _mp2_methods = {"mp2", "rimp2"}
        if input_model.model.method.lower() in _mp2_methods:
            emp2 = bdata["99.0"][-1]
            properties["mp2_total_energy"] = emp2
            qcvars["MP2 TOTAL ENERGY"] = emp2
            qcvars["CURRENT ENERGY"] = emp2
            ehf = bdata["99.0"][1]
            qcvars["HF TOTAL ENERGY"] = ehf
            qcvars["SCF TOTAL ENERGY"] = ehf
            qcvars["CURRENT REFERENCE ENERGY"] = ehf
            qcvars["MP2 CORRELATION ENERGY"] = emp2 - ehf
            qcvars["CURRENT CORRELATION ENERGY"] = emp2 - ehf
            properties["mp2_correlation_energy"] = emp2 - ehf

        # Correct CCSD because its odd?
        # if input_model.model.method.lower() == "ccsd":
        #     m1 = re.findall(" CCSD correlation energy.+=.+\d+\.\d+", outfiles["dispatch.out"])
        #     m2 = re.findall(" CCSD total energy.+=.+\d+\.\d+", outfiles["dispatch.out"])

        props, prov = self._parse_logfile_common(outtext, input_model.dict())
        output_data["provenance"] = prov
        output_data["properties"] = properties
        output_data["properties"].update(props)
        output_data["stdout"] = outfiles["dispatch.out"]
        output_data["success"] = True

        merged_data = {**input_model.dict(), **output_data}
        merged_data["extras"]["qcvars"] = qcvars

        return AtomicResult(**merged_data)

    def _parse_logfile_common(self, outtext: str, input_dict: Dict[str, Any]):
        """
        Parses fields from log file that are not parsed from QCSCRATCH in parse_output
        """

        properties = {}
        provenance = Provenance(creator="QChem", version=self.get_version(), routine="qchem").dict()
        mobj = re.search(r"This is a multi-thread run using ([0-9]+) threads", outtext)
        if mobj:
            provenance["nthreads"] = int(mobj.group(1))

        mobj = re.search(r"Total job time:\s*" + NUMBER + r"s\(wall\)", outtext)
        if mobj:
            provenance["wall_time"] = float(mobj.group(1))

        mobj = re.search(r"Archival summary:\s*\n[0-9]+\\[0-9+]\\([\w\.\-]+)\\", outtext)
        if mobj:
            provenance["hostname"] = mobj.group(1)

        mobj = re.search(r"\n\s*There are\s+(\d+) alpha and\s+(\d+) beta electrons\s*\n", outtext)
        if mobj:
            properties["calcinfo_nalpha"] = int(mobj.group(1))
            properties["calcinfo_nbeta"] = int(mobj.group(2))

        mobj = re.search(r"\n\s*There are\s+\d+ shells and\s+(\d+) basis functions\s*\n", outtext)
        if mobj:
            properties["calcinfo_nbasis"] = int(mobj.group(1))

        mobj = re.search(r"\n\s*RI-MP2 CORRELATION ENERGY\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj:
            properties["mp2_correlation_energy"] = float(mobj.group(1))

        mobj = re.search(r"\n\s*RI-MP2 SINGLES ENERGY\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj:
            properties["mp2_singles_energy"] = float(mobj.group(1))

        mobj_aaaa = re.search(r"\n\s*RI-MP2 ENERGY \(aa\|aa\)\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        mobj_bbbb = re.search(r"\n\s*RI-MP2 ENERGY \(bb\|bb\)\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj_aaaa and mobj_bbbb:
            properties["mp2_same_spin_correlation_energy"] = float(mobj_aaaa.group(1)) + float(mobj_bbbb.group(1))

        mobj_aabb = re.search(r"\n\s*RI-MP2 ENERGY \(aa\|bb\)\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        mobj_bbaa = re.search(r"\n\s*RI-MP2 ENERGY \(bb\|aa\)\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj_aaaa and mobj_bbbb:
            properties["mp2_opposite_spin_correlation_energy"] = float(mobj_aabb.group(1)) + float(mobj_bbaa.group(1))

        properties["calcinfo_natom"] = len(input_dict["molecule"]["symbols"])

        mobj = re.search(r"\n\s*(\d+)\s+" + NUMBER + "\s+" + NUMBER + r"\s+Convergence criterion met\s*\n", outtext)
        if mobj:
            properties["scf_iterations"] = int(mobj.group(1))

        mobj = re.search(
            r"\n\s+Dipole Moment \(Debye\)\s*\n\s+X\s+" + NUMBER + r"\s+Y\s+" + NUMBER + r"\s+Z\s+" + NUMBER + r"\s*\n",
            outtext,
        )
        if mobj:
            cf = constants.conversion_factor("debye", "e * bohr")
            properties["scf_dipole_moment"] = [float(mobj.group(i)) * cf for i in range(1, 4)]

        return properties, provenance

    def parse_logfile(self, outfiles: Dict[str, str]) -> AtomicResult:
        """
        Parses a log file.
        """
        warnings.warn(
            "parse_logfile will result in precision loss for some fields due to trunctation in " "Q-Chem output files."
        )

        outtext = outfiles["dispatch.out"]

        mobj = re.search(r"(?:User input\:|Running Job?)\s+\d+\s+of\s+(\d+)", outtext)
        if mobj:
            if int(mobj.group(1)) > 1:
                raise InputError("Multi-job Q-Chem log files not supported.")

        input_dict = {}
        mobj = re.search(r"\n-{20,}\nUser input:\n-{20,}\n(.+)\n-{20,}", outtext, re.DOTALL)
        if mobj:
            inputtext = mobj.group(1)

            rem_match = re.search(r"\$rem\s*\n([^\$]+)\n\s*\$end", inputtext, re.DOTALL | re.IGNORECASE)
            if rem_match:
                rem_text = rem_match.group(1)
                lines = rem_text.split("\n")
                keywords = {}
                for line in lines:
                    s = re.sub(r"(^|[^\\])!.*", "", line).split()
                    if len(s) == 0:
                        continue
                    keywords[s[0].lower()] = s[1].lower()
                input_dict["model"] = {}
                input_dict["model"]["method"] = keywords.pop("method").lower()
                input_dict["model"]["basis"] = keywords.pop("basis").lower()
                if "jobtype" in keywords:
                    jobtype = keywords.pop("jobtype")
                else:
                    jobtype = "sp"
                _jobtype_to_driver = {
                    "sp": "energy",
                    "force": "gradient",
                    "freq": "hessian",
                }  # optimization intentionally not supported
                try:
                    input_dict["driver"] = _jobtype_to_driver[jobtype]
                except KeyError:
                    raise KeyError(f"Jobtype {jobtype} not supported in qchem log file parser.")

                for key in keywords:
                    if keywords[key] == "false":
                        keywords[key] = False
                    if keywords[key] == "true":
                        keywords[key] = True
                input_dict["keywords"] = keywords

            molecule_match = re.search(r"\$molecule\s*\n([^\$]+)\n\s*\$end", inputtext, re.DOTALL | re.IGNORECASE)
            if molecule_match:
                molecule_text = molecule_match.group(1)
                if keywords.get("input_bohr", False):
                    molecule_text += "\nunits au"
                molecule = Molecule.from_data(molecule_text, dtype="psi4")
                input_dict["molecule"] = molecule.dict()

            _excluded_rem = {
                "input_bohr",
                "mem_total",
                "mem_static",
            }  # parts of the rem normally written by build_input, and which do not affect results
            for item in _excluded_rem:
                if item in input_dict["keywords"]:
                    input_dict["keywords"].pop(item)

        try:
            qcscr_result = self.parse_output(outfiles, AtomicInput(**input_dict)).dict()
        except KeyError:
            props, prov = self._parse_logfile_common(outtext, input_dict)
            qcscr_result = {"properties": props, "provenance": prov, **input_dict}

        mobj = re.search(r"\n\s*Total\s+energy in the final basis set =\s+" + NUMBER + r"\s*\n", outtext)
        if mobj and qcscr_result["properties"].get("scf_total_energy", None) is None:
            qcscr_result["properties"]["scf_total_energy"] = float(mobj.group(1))

        mobj = re.search(r"\n\s*Nuclear Repulsion Energy =\s+" + NUMBER + r"\s+hartrees\s*\n", outtext)
        if mobj and qcscr_result["properties"].get("nuclear_repulsion_energy", None) is None:
            qcscr_result["properties"]["nuclear_repulsion_energy"] = float(mobj.group(1))

        mobj = re.search(r"\n\s*RI-MP2 TOTAL ENERGY\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj and qcscr_result["properties"].get("mp2_total_energy", None) is None:
            qcscr_result["properties"]["mp2_total_energy"] = float(mobj.group(1))

        _scf_methods = {
            "hf",
            "spw92",
            "lda",
            "svwn5",
            "b97-d3(0)",
            "b97-d",
            "pbe",
            "blyp",
            "revpbe",
            "beef-vdw",
            "bop",
            "bp86",
            "bp86vwn",
            "bpbe",
            "edf1",
            "edf2",
            "gam",
            "hcth93",
            "hcth120",
            "hcth147",
            "hcth407",
            "hle16",
            "kt1",
            "kt2",
            "kt3",
            "mpw91",
            "n12",
            "olyp",
            "pbeop",
            "pbesol",
            "pw91",
            "rpbe",
            "rvv10",
            "sogga",
            "sogga11",
            "vv10",
            "b97m-v",
            "b97m-rv",
            "m06-l",
            "tpss",
            "revtpss",
            "bloc",
            "m11-l",
            "mbeef",
            "mgga_ms0",
            "mgga_ms1",
            "mgga_ms2",
            "mgga_mvs",
            "mn12-l",
            "mn15-l",
            "otpss",
            "pkzb",
            "scan",
            "t-hcth",
            "tm",
            "vsxc",
            "b3lyp",
            "pbe0",
            "revpbe0",
            "b97",
            "b1lyp",
            "b1pw91",
            "b3lyp5",
            "b3p86",
            "b1lyp",
            "b1pw91",
            "b3lyp5",
            "b3p86",
            "b3pw91",
            "b5050lyp",
            "b97-1",
            "b97-2",
            "b97-3",
            "b97-k",
            "bhhlyp",
            "hflyp",
            "mpw1k",
            "mpw1lyp",
            "mpw1pbe",
            "mpw1pw91",
            "o3lyp",
            "pbeh-3c",
            "pbe50",
            "sogga11-x",
            "wc04",
            "wp04",
            "x3lyp",
            "m06-2x",
            "m08-hx",
            "tpssh",
            "revtpssh",
            "b1b95",
            "b3tlap",
            "bb1k",
            "bmk",
            "dldf",
            "m05",
            "m05-2x",
            "m06",
            "m06-hf",
            "m08-so",
            "mgga_ms2h",
            "mgga_mvsh",
            "mn15",
            "mpw1b95",
            "mpwb1k",
            "pw6b95",
            "pwb6k",
            "scan0",
            "t-hcthh",
            "tpss0",
            "wb97x-v",
            "wb97x-d3",
            "wb97x-d",
            "cam-b3lyp",
            "cam-qtp00",
            "cam-qtp01",
            "hse-hjs",
            "lc-rvv10",
            "lc-vv10",
            "lc-wpbe08",
            "lrc-bop",
            "lrc-wpbe",
            "lrc-wpbeh",
            "n12-sx",
            "rcam-b3lyp",
            "wb97",
            "wb97x",
            "wb97x-rv",
            "wb97m-v",
            "m11",
            "mn12-sx",
            "wb97m-rv",
            "wm05-d",
            "wm06-d3",
            "dsd-pbepbe-d3",
            "wb97x-2(lp)",
            "wb97x-2(tqz)",
            "xyg3",
            "xygj-os",
            "b2plyp",
            "b2gpplyp",
            "dsd-pbep86-d3",
            "ls1dh-pbe",
            "pbe-qidh",
            "pbe0-2",
            "pbe0-dh",
            "ptpss-d3",
            "dsd-pbeb95-d3",
            "pwpb95-d3",
        }
        _mp2_methods = {"rimp2"}

        method = input_dict["model"]["method"].lower()
        if qcscr_result["properties"].get("return_energy", None) is None:
            if method in _mp2_methods:
                qcscr_result["properties"]["return_energy"] = qcscr_result["properties"]["mp2_total_energy"]
            elif method in _scf_methods:
                qcscr_result["properties"]["return_energy"] = qcscr_result["properties"]["scf_total_energy"]
            else:
                raise NotImplementedError(f"Method {method} not supported by logfile parser for energy driver.")

        if input_dict["driver"] == "gradient" and qcscr_result.get("return_result", None) is None:

            def read_matrix(text):
                lines = text.split("\n")
                i = 0
                mdict = defaultdict(dict)
                maxcol = 0
                maxrow = 0
                while i < len(lines):
                    cols = [int(idx) for idx in lines[i].split()]
                    maxcol = max(maxcol, *cols)
                    i += 1
                    while i < len(lines):
                        s = lines[i].split()
                        if len(s) <= len(cols):
                            break
                        assert len(s) == len(cols) + 1, s
                        row = int(s[0])
                        maxrow = max(maxrow, row)
                        data = [float(field) for field in s[1:]]
                        for col_idx, col in enumerate(cols):
                            mdict[row - 1][col - 1] = data[col_idx]
                        i += 1

                ret = np.zeros((maxrow, maxcol))
                for row in mdict:
                    for col in mdict[row]:
                        ret[row, col] = mdict[row][col]
                return ret

            if method in _mp2_methods:
                mobj = re.search(
                    r"\n\s+Full Analytical Gradient of MP2 Energy \(in au.\)\s*\n"
                    r"([\s\d\.EDed\+\-]+)\n"
                    r"\s*Gradient time:",
                    outtext,
                )

            elif method in _scf_methods:
                mobj = re.search(
                    r"\n\s+Gradient of SCF Energy\s*\n([\s\d\.EDed\+\-]+)\n\s*Max gradient component =", outtext
                )

            else:
                raise NotImplementedError(f"Method {method} not supported by the logfile parser for gradient driver.")

            if mobj:
                qcscr_result["return_result"] = read_matrix(mobj.group(1)).T

        qcscr_result["success"] = True  # XXX: have a nice day?
        qcscr_result["stdout"] = outtext
        if input_dict["driver"] == "energy" and qcscr_result.get("return_result", None) is None:
            qcscr_result["return_result"] = qcscr_result["properties"]["return_energy"]

        return AtomicResult(**qcscr_result)
