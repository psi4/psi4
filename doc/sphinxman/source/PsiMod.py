"""PsiMod module
The Sphinx documentation system scans all Python modules for docstrings
and class structure (good) but in doing so imports modules it finds in
input statements and follows any exposed code like the procedures
dictionary (bad for Psi4 b/c of PsiMod). So, this fake PsiMod module
exists to appease Sphinx when it looks for PsiMod.py to import. Any PsiMod
commands that aren't protected by functions need to have skeleton versions
here.

"""

class Molecule:
    pass

def sldkfjksl():
    return 4.4

def plugin_load(sofile):
    return 4


class SuperFunctional():

    def __init__(self):
        pass

    @staticmethod
    def blank():
        return SuperFunctional()

    @staticmethod
    def build(sstr, iint, iint2):
        return SuperFunctional()

    def add_c_functional(self, Functional):
        pass

    def add_x_functional(self, Functional):
        pass

    def allocate(self):
        pass

    def ansatz(self):
        pass

    def c_alpha(self):
        pass

    def c_functional(self, sstr):
        pass

    def c_omega(self):
        pass

    def citation(self):
        pass

    def deriv(self):
        pass

    def description(self):
        pass

    def dispersion(self):
        pass

    def is_c_hybrid(self):
        pass

    def is_c_lrc(self):
        pass

    def is_gga(self):
        pass

    def is_meta(self):
        pass

    def is_x_hybrid(self):
        pass

    def is_x_lrc(self):
        pass

    def max_points(self):
        pass

    def name(self):
        return 'SuperFunctionalName'

    def print_detail(self, iint):
        pass

    def print_out(self):
        pass

    def set_c_alpha(self, ffloat):
        pass

    def set_c_omega(self, ffloat):
        pass

    def set_citation(self, sstr):
        pass

    def set_deriv(self, iint):
        pass

    def set_description(self, sstr):
        pass

    def set_dispersion(self, Dispersion):
        pass

    def set_max_points(self, iint):
        pass

    def set_name(self, sstr):
        pass

    def set_x_alpha(self, ffloat):
        pass

    def set_x_omega(self, ffloat):
        pass

    #def test_functional(self, ):
    #    pass

    def value(self, sstr):
        pass

    def x_alpha(self):
        pass

    def x_functional(self, sstr):
        pass

    def x_omega(self):
        pass


class Functional():

    def __init__(self):
        pass

    @staticmethod
    def build_base(sstr):
        return Functional()
  
    def alpha(self):
        pass
    
    def citation(self):
        pass

    def description(self):
        pass

    def is_gga(self):
        pass

    def is_lrc(self):
        pass

    def is_meta(self):
        pass

    def lsda_cutoff(self):
        pass

    def meta_cutoff(self):
        pass

    def name(self):
        pass

    def omega(self):
        pass

    def print_detail(SuperFunctional, iint):
        pass

    def print_out(self):
        pass

    def set_alpha(self, ffloat):
        pass

    def set_citation(self, sstr):
        pass

    def set_description(self, sstr):
        pass

    def set_gga(self, bbool):
        pass

    def set_lsda_cutoff(self, ffloat):
        pass

    def set_meta(self, bbool):
        pass

    def set_meta_cutoff(self, ffloat):
        pass

    def set_name(self, sstr):
        pass

    def set_omega(self, ffloat):
        pass

    def set_parameter(self, sstr, ffloat):
        pass


class Dispersion():

    def __init__(self):
        pass

    @staticmethod
    def build(sstr, ffloat):
        return Dispersion()
  
