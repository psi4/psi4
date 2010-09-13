import PsiMod
import psiopt

class PsiException: pass
class ValueNotSet (PsiException): pass
class RowAlignmentError(PsiException):

    def __init__(self, l1name, l1, l2name, l2):
        msg = "Rows %s and %s not aligned. Length %d != %d" % (l1name, l2name, l1, l2)
        PsiException.__init__(self, msg)

molecule = None

def new_set_attr(self, name, value):
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)
    if isvar:
        fxn = object.__getattribute__(self, "set_variable")
        fxn(name, value)
        return

    object.__setattr__(self, name, value)

def new_get_attr(self, name):
    fxn = object.__getattribute__(self, "is_variable")
    isvar = fxn(name)

    if isvar:
        fxn = object.__getattribute__(self, "get_variable")
        return fxn(name)

    return object.__getattribute__(self, name)


def dynamic_variable_bind(cls):
    #class specific
    cls.__setattr__ = new_set_attr
    cls.__getattr__ = new_get_attr


dynamic_variable_bind(PsiMod.Molecule) #pass class type, not class instance

#
# Define geometry to be used by PSI4.
# The molecule created by this will be set in options.
#
# geometry("
#   O  1.0 0.0 0.0
#   H  0.0 1.0 0.0
#   H  0.0 0.0 0.0
#
def geometry(geom, name = "default"):
    # Create a Molecule object
    global molecule
    molecule = PsiMod.Molecule.create_molecule_from_string(geom)
    molecule.set_name(name)

    PsiMod.set_active_molecule(molecule)
    PsiMod.IO.set_default_namespace(molecule.get_name())

    return molecule

def dummify():
    if not molecule:
        raise ValueNotSet("no default molecule found")
    #molecule.set_dummy_atom

# extract_clusters(shared_ptr<Molecule> molecule, bool ghost, int cluster_size)
# returns all subclusters of the molecule of real size cluster_size
# and all other atoms ghosted if ghost equals true, all other atoms 
# discarded if ghost equals false
#
# if cluster_size = 0, returns all possible combinations of cluster size
# 
# defaults: 
#  ghost = True 
#  cluster_size = 0 
#
def extract_clusters(mol, ghost = True, cluster_size = 0):
   
    if not mol:
        raise ValueNotSet("no molecule found")
    
    #How many levels of clusters are possible?
    nfrag = mol.nfragments()
    
    #Initialize the cluster array
    clusters = []

    #scope the arrays
    reals = []
    ghosts = []    

    #counter
    counter = 0

    #loop over all possible cluster sizes
    for nreal in range(nfrag,0,-1):

        #if a specific cluster size size is requested, only do that 
        if (nreal != cluster_size and cluster_size > 0):
            continue
        
        # initialize the reals list
        reals = []
        
        #setup first combination [3,2,1] lexical ordering
        #fragments indexing is 1's based, bloody hell
        for index in range(nreal,0,-1):
            reals.append(index)
       
        #start loop through lexical promotion
        while True:

            counter = counter + 1

            #Generate cluster from last iteration
            if (ghost):
                ghosts = []
                for g in range(nfrag,0,-1):
                    if (g not in reals): 
                        ghosts.append(g)
                #print "Cluster #%d: %s reals, %s ghosts" % (counter,str(reals), str(ghosts))
                clusters.append(mol.extract_subsets(reals,ghosts)) 
            else:
                #print "Cluster #%d: %s reals" % (counter,str(reals))
                clusters.append(mol.extract_subsets(reals)) 

            #reset rank
            rank = 0;
            
            #look for lexical promotion opportunity
            #i.e.: [4 2 1] has a promotion opportunity at 
            # index 1 to produce [4 3 1] 
            for k in range(nreal - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1;
                    break

            #do the promotion
            reals[rank] = reals[rank] + 1;

            #demote the right portion of the register
            val = 1
            for k in range(nreal-1,rank,-1):
                reals[k] = val
                val = val + 1
            
            #boundary condition is promotion into
            #[nfrag+1 nfrag-1 ...]
            if (reals[0] > nfrag):
                break
            
    return clusters

# wxtract_cluster_indexing(shared_ptr<Molecule> molecule, bool ghost, int cluster_size)
# returns a LIST of all subclusters of the molecule of real size cluster_size
# and all other atoms ghosted if ghost equals true, all other atoms 
# discarded if ghost equals false
#
# if cluster_size = 0, returns all possible combinations of cluster size
# 
# defaults: 
#  ghost = True 
#  cluster_size = 0 
#
def extract_cluster_indexing(mol,cluster_size = 0):

    import copy   

    if not mol:
        raise ValueNotSet("no molecule found")
    
    #How many levels of clusters are possible?
    nfrag = mol.nfragments()
 
    #Initialize the cluster array
    clusters = []

    #scope the arrays
    reals = []

    #counter
    counter = 0

    #loop over all possible cluster sizes
    for nreal in range(nfrag,0,-1):

        #if a specific cluster size size is requested, only do that 
        if (nreal != cluster_size and cluster_size > 0):
            continue
        
        # initialize the reals list
        reals = []
        
        #setup first combination [3,2,1] lexical ordering
        #fragments indexing is 1's based, bloody hell
        for index in range(nreal,0,-1):
            reals.append(index)
       
        #start loop through lexical promotion
        while True:

            counter = counter + 1

            #Generate cluster from last iteration
            clusters.append(copy.deepcopy(reals)) 

            #reset rank
            rank = 0;
            
            #look for lexical promotion opportunity
            #i.e.: [4 2 1] has a promotion opportunity at 
            # index 1 to produce [4 3 1] 
            for k in range(nreal - 2, -1, -1):
                if (reals[k] != reals[k + 1] + 1):
                    rank = k + 1;
                    break

            #do the promotion
            reals[rank] = reals[rank] + 1;

            #demote the right portion of the register
            val = 1
            for k in range(nreal-1,rank,-1):
                reals[k] = val
                val = val + 1
            
            #boundary condition is promotion into
            #[nfrag+1 nfrag-1 ...]
            if (reals[0] > nfrag):
                break
            
    return clusters

#
# Set options
#
# options({
#    'WFN': 'SCF',
#    'BASIS': 'STO-3G'
# })
def options(dict, module = ''):
    m = module.strip().upper()
    if len(m) > 0:
        PsiMod.set_default_options_for_module(m)

    for key in dict.keys():
        PsiMod.set_option(key.upper(), dict[key]);

#
# Define geometry to be used by PSI4.
# This geometry will be saved to the checkpoint file.
#
# geometry( [[ "O", 0.0, 0.0, 0.0 ],
#            [ "H", 1.0, 0.0, 0.0 ]] )
#
def geometry_old(geom, reorient = True, prefix = "", chkpt = None, shiftToCOM = True):
    # Make sure the user passed in what we expect
    if isinstance(geom, list) == False:
        raise TypeError("geometry must be a list")

    # Create a Molecule object from PSI4 C++
    molecule = PsiMod.Molecule()

    # For each atom in the geometry given add it to the molecule
    for atom in geom:
        # Make sure the atom is a list
        if isinstance(atom, list) == False:
            raise TypeError("an atom entry is not a list")

        # Make sure the length of the atom entry is long enough
        if len(atom) < 4:
            print "Insufficient information for atom entry:"
            print atom
            raise AttributeError("atom entries must have at least 4 elements")

        molecule.add_atom(1, atom[1], atom[2], atom[3], atom[0], 1.0, 0, 0.0)

    # Print the molecule:
    molecule.print_to_output()

    # If requested, shift molecule to center of mass
    if shiftToCOM == True:
        molecule.move_to_com()

    # If requested, reorient the molecule.
    if reorient == True:
        molecule.reorient()

        molecule.print_to_output()

    # Save the molecule to the checkpoint file using prefix
    if chkpt == None:
        # If the user didn't provide a checkpoint object create one
        # using the shared psio object
        psio = PsiMod.IO.shared_object()
        chkpt = PsiMod.Checkpoint(psio, 1)

    molecule.save_to_checkpoint(chkpt, prefix)

    return molecule

def activate(mol):
    PsiMod.set_active_molecule(mol)
    PsiMod.IO.set_default_namespace(mol.get_name())

class Table:

    def __init__(self, rows=(),
                 row_label_width=10,
                 row_label_precision=4,
                 cols=(),
                 width=16, precision=10):
        self.row_label_width = row_label_width
        self.row_label_precision = row_label_precision
        self.width = width
        self.precision = precision
        self.rows = rows

        if isinstance(cols, str):
            self.cols = (cols,)
        else:
            self.cols = cols

        self.labels = []
        self.data = []

    def format_label(self):
        #str = lambda x: (('%%%d.%df' % (self.row_label_width, self.row_label_precision)) % x)
        str = lambda x: (('%%%ds' % (self.row_label_width)) % x)
        return " ".join(map(str, self.labels))

    def format_values(self, values):
        str = lambda x: (('%%%d.%df' % (self.width, self.precision)) % x)
        return " ".join(map(str, values))

    def __getitem__(self, value):
        self.labels.append(value)
        return self

    def __setitem__(self, name, value):
        self.labels.append(name)
        label = self.format_label()
        self.labels = []

        if isinstance(value, list):
            self.data.append( (label, value ) )
        else:
            self.data.append( (label, [value] ) )

    def save(self, file):
        import pickle
        pickle_str = pickle.dumps(self)
        fileobj = open(file, "w")
        fileobj.write(str(self))
        fileobj.close()

    def __str__(self):
        rowstr = lambda x: '%%%ds' % self.row_label_width % x
        colstr = lambda x: '%%%ds' % self.width % x

        lines = []

        table_header = ""
        if isinstance(self.rows, str):
            table_header += "%%%ds" % self.row_label_width % self.rows
        else:
            table_header += " ".join(map(rowstr, self.rows))
        table_header += " ".join(map(colstr, self.cols))

        lines.append(table_header)

        for datarow in self.data:
            #print datarow
            row_data = datarow[0]
            row_data += self.format_values(datarow[1])
            lines.append(row_data)

        return "\n".join(lines)

    def copy(self):
        import copy
        return copy.deepcopy(self)

    def absolute_to_relative(self, Factor = 627.51):
        import copy

        if len(self.data) == 0:
            return

        current_min = list(copy.deepcopy(self.data[0][1]))
        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                if current_min[col] > datarow[1][col]:
                    current_min[col] = datarow[1][col]

        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                #print datarow[1][col]
                datarow[1][col] = (datarow[1][col] - current_min[col]) * Factor

def banner(text, type = 1, width = 35):
    lines = text.split('\n')
    max_length = 0
    for line in lines:
        if (len(line) > max_length):
            max_length = len(line)
    
    max_length = max([width, max_length]);
   
    null = '' 
    if type == 1:
        banner  = '  //' + null.center(max_length,'>') + '//\n'
        for line in lines:
            banner += '  //' + line.center(max_length) + '//\n'
        banner += '  //' + null.center(max_length,'<') + '//\n'
    
    if type == 2:
        banner = ''
        for line in lines:
            banner += (' ' + line + ' ').center(max_length,'=') 

    PsiMod.print_out(banner)

def energy(name):
    
    handle = energy_procedures(name)
    handle(name)

def energy_procedures(name):
    if (name.lower() == 'sapt0'):
        return lambda name : run_sapt('sapt0') 
    elif (name.lower() == 'sapt2'):
        return lambda name : run_sapt('sapt2')
    elif (name.lower() == 'sapt2+'):
        return lambda name : run_sapt('sapt2+')
    elif (name.lower() == 'sapt2+(3)'):
        return lambda name : run_sapt('sapt2+(3)')
    else:
        raise 'Undefined Energy Procedure' 

def run_sapt(name):
    if not molecule:
        raise ValueNotSet("no molecule found")
     
    molecule.update_geometry()
    monomerA = molecule.extract_subsets(1,2)
    monomerA.set_name("monomerA")   
    monomerB = molecule.extract_subsets(2,1)    
    monomerB.set_name("monomerB")   

    mol = molecule

    PsiMod.set_active_molecule(mol)
    PsiMod.IO.set_default_namespace("dimer")
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("SAPT","2-dimer")
    PsiMod.print_out("\n")
    banner('Dimer HF')
    PsiMod.print_out("\n")
    e_dimer = PsiMod.scf()

    PsiMod.set_active_molecule(monomerA)
    PsiMod.IO.set_default_namespace("monomerA")
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("SAPT","2-monomer_A")
    PsiMod.print_out("\n")
    banner('Monomer A HF')
    PsiMod.print_out("\n")
    e_monomerA = PsiMod.scf()
    
    PsiMod.set_active_molecule(monomerB)
    PsiMod.IO.set_default_namespace("monomerB")
    PsiMod.set_default_options_for_module("SCF")
    PsiMod.set_option("NO_INPUT",True)
    PsiMod.set_option("SAPT","2-monomer_B")
    PsiMod.print_out("\n")
    banner('Monomer B HF')
    PsiMod.print_out("\n")
    e_monomerB = PsiMod.scf()

    PsiMod.IO.change_file_namespace(121,"monomerA","dimer")
    PsiMod.IO.change_file_namespace(122,"monomerB","dimer")

    PsiMod.set_active_molecule(mol)
    PsiMod.IO.set_default_namespace("dimer")
    PsiMod.set_default_options_for_module("SAPT")
    PsiMod.set_option("NO_INPUT",True)
    if (name.lower() == 'sapt0'):
        PsiMod.set_option("SAPT_LEVEL","SAPT0")
    elif (name.lower() == 'sapt2'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2")
    elif (name.lower() == 'sapt2+'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2+")
    elif (name.lower() == 'sapt2+(3)'):
        PsiMod.set_option("SAPT_LEVEL","SAPT2+(3)")
    PsiMod.print_out("\n")
    banner(name.upper())
    PsiMod.print_out("\n")
    e_sapt = PsiMod.sapt()
