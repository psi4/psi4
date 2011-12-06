import PsiMod

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
    molecule = PsiMod.Molecule.create_molecule_from_string(geom)
    molecule.set_name(name)

    activate(molecule)

    return molecule

def activate(mol):
    PsiMod.set_active_molecule(mol)
    #PsiMod.IO.set_default_namespace(mol.get_name())

def activate_potential(pot):
    PsiMod.set_active_potential(pot)


