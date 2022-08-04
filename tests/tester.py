import re
from pathlib import Path

# Goals for Psithon Tests
# -----------------------
# [x] Notice if test undocumented
# [x] Notice if test unregistered to CMake/CTest or registered with wrong name
# [x] Notice if test unregistered to pytest or registered with wrong name
# [x] If test unregistered to pytest, fill in runner file (completely for most, but add'l files are manual)
# [x] Check ctest/pytest marks consistency
# [ ] Write Sphinx files
# [ ] Translate docs comment from LaTeX to reST
# [ ] Update samples/
# [ ] Replace document_tests.pl


complaints = []

internal_dirs = [
    "psithon2/psiaux1/myplugin1",
    "dftd3/psithon2/psiaux1/myplugin1",
]

skip_nested_dirs = [
    "brianqc",  # can't test locally
    "cfour",  # imminent deletion
    "erd",  # imminent deletion
    "pasture-ccsorttransqt2",  # retired
]

model_testinputpy = """from addons import *

@ctest_labeler("{markstr}")
def test_{pytest_name}():
    ctest_runner(__file__{copyfilesstr})

"""

def filter_marks(marks, nested):
    markstr = []
    for m in marks:
        if m == "psi":
            pass
        elif m == "quicktests":
            markstr.append("quick")
        elif m == "longtests":
            markstr.append("long")
        elif m == "smoketests":
            markstr.append("smoke")
        elif m == "addon":
            pass
        elif m == nested:
            if m in ["json"]:
                markstr.append(m)
            else:
                pass
        else:
            markstr.append(m)
    return markstr


with open("CMakeLists.txt") as fp:
    cm = fp.read()
cm = cm.replace("add_subdirectory(", "").replace(")", "").split()

tests = Path(".")
for fl in sorted(tests.rglob("*")):
    if fl.name in ["input.dat", "input.py"]:
        testdir = fl.parent

        legit = True
        ctest_name = str(testdir).replace("/", "-")
        pytest_name = ctest_name.replace("-", "_").replace("+", "plus")

        nested = "/" in str(testdir)
        if nested:
            nested = str(testdir.parent)

        if str(testdir) in internal_dirs:
            continue
        if nested in skip_nested_dirs:
            continue

        if nested:
            fnested_cm = testdir.parent / "CMakeLists.txt"
            with open(fnested_cm) as fp:
                nested_cm = fp.read()
            nested_cm = nested_cm.replace("add_subdirectory(", "").replace(")", "").split()
            if str(testdir.name) not in nested_cm:
                complaints.append(f"{testdir}: missing cmake directory registration. `vi {fnested_cm}`")
                legit = False
        else:
            if ctest_name not in cm:
                complaints.append(f"{testdir}: missing cmake directory registration. `vi CMakeLists.txt`")
                legit = False

        with open(fl) as fp:
            inputdat = fp.readlines()

        description = ""
        for ln in inputdat:
            if ln.startswith("#!"):
                description += ln[2:]
        description = description.strip()
        if not description:
            complaints.append(f"{testdir}: missing docs comment. `vi {fl}`")

        marks = []
        copyfiles = False
        cml = fl.with_name("CMakeLists.txt")
        if cml.is_file():
            with open(cml) as fp:
                cmakeliststxt = fp.read()

            mobj = re.search(r"^\s*" + r"add_regression_test\(" + r"(?P<name>([a-zA-Z0-9-+_]+))" + r'\s+\"' + r"(?P<marks>([a-z0-9-_;]+))" + r'\"\)', cmakeliststxt, re.MULTILINE)
            if mobj:
                if mobj.group("name") != ctest_name:
                    complaints.append(f"{testdir}: mismatched directory ({ctest_name}) and ctest registration name ({mobj.group('name')}). `vi {cml}`")
                marks = mobj.group("marks").split(";")
            else:
                complaints.append(f"{testdir}: missing ctest registration. `vi {cml}`")
                legit = False

            mobj = re.search(f"^\s*" + r"file\(COPY", cmakeliststxt, re.MULTILINE)
            if mobj:
                copyfiles = True

        else:
            complaints.append(f"{testdir}: missing CMakeLists. `vi {cml}`")
            legit = False

        marks = filter_marks(marks, nested=nested)
        markstr = ";".join(marks)

        tipy = fl.with_name("test_input.py")
        if tipy.is_file():
            with open(tipy) as fp:
                testinputpy = fp.read()

            pymarks = []
            pyplugins = []
            mobj = re.findall(r'^@uusing\("' + r"([a-z0-9-_;]+)" + r'\"\)', testinputpy, re.MULTILINE)
            if mobj:
                pyplugins = mobj

            mobj = re.search(r'^@ctest_labeler\("' + r"(?P<pymarks>([a-z0-9-_;]+))" + r'\"\)', testinputpy, re.MULTILINE)
            if mobj:
                pymarks = mobj.group("pymarks").split(";")

            if (set(marks) - set(pyplugins)) != set(pymarks):
                complaints.append(f"{testdir}: mismatched marks ctest ({markstr}) and pytest ({';'.join(pymarks)}). `vi {testdir}/CMakeLists.txt {testdir}/test_input.py`")

            mobj = re.search(r"^def test_" + r"(?P<name>([a-zA-Z0-9_]+))" + r"\(", testinputpy, re.MULTILINE)
            if mobj:
                if mobj.group("name") != pytest_name:
                    complaints.append(f"{testdir}: mismatched directory ({pytest_name}) and pytest registration ({mobj.group('name')}). `vi {tipy}`")

        else:
            if legit:
                complaints.append(f"{testdir}: missing pytest input generated. check it! `vi {tipy}`")

                copyfilesstr = ", [\n    ]" if copyfiles else ""
                with open(tipy, "w") as fp:
                    fp.write(model_testinputpy.format(markstr=markstr, pytest_name=pytest_name, copyfilesstr=copyfilesstr))


if complaints:
    print("Complaints\n----------")
    for idx, item in enumerate(complaints):
        print(f"- [ ] {idx + 1:3}. {item}")
