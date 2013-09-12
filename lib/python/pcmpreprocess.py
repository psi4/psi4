import pcmgetkw as getkw
from copy import deepcopy

isAngstrom = False

toAngstrom =  0.52917721092; # CODATA 2010 recommended data
toAtomicUnits = 1.0/toAngstrom;

def	setup_keywords():
	top=getkw.Section('toplevel', callback=verify_top)
	top.set_status(True)
	top.add_kw('Units',     'STR', 'AU')

	cavity=getkw.Section('Cavity', callback=verify_cavity)
	cavity.add_kw('Type','STR')
	cavity.add_kw('PatchLevel', 'INT', 2)
	cavity.add_kw('Coarsity', 'DBL', 0.5)
	cavity.add_kw('Area','DBL', 0.3)
	cavity.add_kw('Scaling', 'BOOL', True)
	cavity.add_kw('AddSpheres', 'BOOL', True)
        cavity.add_kw('Mode','STR','Explicit')
        cavity.add_kw('Atoms','INT_ARRAY')
        cavity.add_kw('Radii','DBL_ARRAY')
	cavity.add_kw('RadiiSet', 'STR', 'Bondi')
	cavity.add_kw('Spheres','DBL_ARRAY', callback=verify_spheres)
	top.add_sect(cavity)
        
	medium=getkw.Section('Medium', callback=verify_medium)
	medium.add_kw('Solvent',     'STR', 'Water')
	medium.add_kw('SolverType',  'STR', 'IEFPCM')
	medium.add_kw('EquationType','STR', 'SecondKind')
	medium.add_kw('Correction', 'DBL', 0.0)
	medium.add_kw('ProbeRadius', 'DBL', 1.0)
	top.add_sect(medium)
	
	green=getkw.Section('Green', callback=verify_green)
	green.add_kw('Type',           'STR', 'Vacuum')
	green.add_kw('Der',            'STR', 'Derivative')
	green.add_kw('Eps',            'DBL', 1.0)
	green.add_kw('EpsRe',          'DBL', 1.0)
	green.add_kw('EpsImg',         'DBL', 1.0)
	green.add_kw('SphereRadius',   'DBL', 1.0)
	green.add_kw('SpherePosition', 'DBL_ARRAY')
	medium.add_sect(green)

	green_part = deepcopy(green)
	green.add_sect(green_part)

	return top

def verify_top(section):
	global isAngstrom
	allowed_units = ('AU', 'Angstrom')
	key = section.get('Units')
	val = key.get()
	if (val not in allowed_units):
	    	print "Allowed units are: ", allowed_units
		sys.exit(1)
	if (val == 'Angstrom'):
		isAngstrom = True

def verify_cavity(section):
	allowed = ('GePol', 'Wavelet')
        type = section.get('Type')
        if (type.get() not in allowed):
        	print "Allowed types are: ", allowed
        	sys.exit(1)

	if section['Area'].is_set(): convert_area_scalar(section['Area'])
	if (type.get() == 'GePol'):
        	area=section.get('Area')
        	a=area.get()
        	if (a < 0.01):
        		print "Area value is too small"
        		print "Minimum value: 0.01"
        		sys.exit(1)
        elif (type.get() == 'Wavelet'):
        	key = section.get('PatchLevel')
        	if (key.get() < 1):
        		print "Patch level must be > 0"
        		sys.exit(1)
        	key = section.get('Coarsity')
        	if (key.get() < 0.0 or key.get() >= 1.0):
        		print "Coarsity has to be within ]0,1["
        		sys.exit(1)
        radiiSet = section.get('RadiiSet')
        allowed_sets = ('Bondi', 'UFF')
        if (radiiSet.get() not in allowed_sets):
                print "Allowed radii sets are: ", allowed_sets
                sys.exit(1)
	allowed_modes = ("Explicit", "Atoms", "Implicit")
	mode = section.get('Mode')
	if (mode.get() not in allowed_modes):
		print "Allowed modes are: ", allowed_modes
		sys.exit(1)                         

        atoms=section.get('Atoms')
        at=atoms.get()
        radii = section.get('Radii')
	convert_length_array(radii);
        r = radii.get()
	
        if (mode.get() == 'Atoms'):
		if (len(r) != len(at) or len(at) == 0):
			print "Incoherent input for Atoms keyword."
			print "Check that Atoms and Radii are consistent."
			sys.exit(1)
		else:
			for i, v in enumerate(at):
				if (at.count(v) > 1):
					print "Incoherent input for Atoms keyword."
					print "Too many spheres on the same atom(s)."
					sys.exit(1)
	
	

def verify_medium(section):
	allowedSolvents = {'Water':                ('Water',                'water',                'H2O'),
			   'Methanol':             ('Methanol',             'methanol',             'CH3OH'),
			   'Ethanol':              ('Ethanol',              'ethanol',              'CH3CH2OH'),
			   'Chloroform':           ('Chloroform',           'chloroform',           'CHCl3'),
			   'Methylenechloride':    ('Methylenechloride',    'methylenechloride',    'CH2Cl2'),
			   '1,2-Dichloroethane':   ('1,2-Dichloroethane',   '1,2-dichloroethane',   'C2H4Cl2'),
			   'Carbon Tetrachloride': ('Carbon Tetrachloride', 'carbon tetrachloride', 'CCl4'),
			   'Benzene':              ('Benzene',              'benzene',              'C6H6'),
			   'Toluene':              ('Toluene',              'toluene',              'C6H5CH3'),
			   'Chlorobenzene':        ('Chlorobenzene',        'chlorobenzene',        'C6H5Cl'),
			   'Nitromethane':         ('Nitromethane',         'nitromethane',         'CH3NO2'),
			   'N-heptane':            ('N-heptane',            'n-heptane',            'C7H16'),
			   'Cyclohexane':          ('Cyclohexane',          'cyclohexane',          'C6H12'),
			   'Aniline':              ('Aniline',              'aniline',              'C6H5NH2'),
			   'Acetone':              ('Acetone',              'acetone',              'C2H6CO'),
			   'Tetrahydrofurane':     ('Tetrahydrofurane',     'tetrahydrofurane',     'THF'),
			   'Dimethylsulfoxide':    ('Dimethylsulfoxide',    'dimethylsulfoxide',    'DMSO'),
			   'Acetonitrile':         ('Acetonitrile',         'acetonitrile',         'CH3CN'),
			   'Explicit':             ('Explicit',             'explicit')}
	solvent = section.get('Solvent')
	explicitSolvent = solvent.get() in allowedSolvents['Explicit']
	if(explicitSolvent):
		PRF = section.is_set('ProbeRadius')
		GIF = section.is_set('Green<inside>')
		GOF = section.is_set('Green<outside>')
		if (not PRF):
			print "Error: Explicit solvent chosen but ProbeRadius not specified"
		if (not GIF):
			print "Error: Explicit solvent chosen but Green<inside> not specified"
		if (not GOF):
			print "Error: Explicit solvent chosen but Green<outside> not specified"
		if (not GIF or not GOF or not PRF):
			sys.exit(1)
	solventFound = False
	for i, v in allowedSolvents.iteritems():
		if (solvent.get() in v):
			solventName = i
			solventFound = True
			break
	if (not solventFound):
		print "Unknown solvent"
		print "Choose a solvent from the following list: "
		print allowedSolvents.keys()
		print "or specify the solvent data explicitly."
		sys.exit(1)

        correction = section.get('Correction')
	if (correction.get() < 0.0):
		print "Correction for CPCM solver must be greater than 0.0"
		sys.exit(1)

	convert_length_scalar(section.get('ProbeRadius'))
	radius = section.get('ProbeRadius')
	if (radius.get() < 0.1 or radius.get() > 100):
		print "Probe radius has to be within [0.1,100] Atomic Units"
		sys.exit(1)

	allowed_types = ('IEFPCM', 'CPCM', 'Wavelet', 'Linear')
        key = section.get('SolverType')
        val = key.get()
        if (val not in allowed_types):  
                print "Allowed types are: ", allowed_types
        	sys.exit(1)
	allowed_equations = ('FirstKind', 'SecondKind', 'Full')
        key = section.get('EquationType')
        val = key.get()
        if (val not in allowed_equations):  
                print "Allowed equations are: ", allowed_equations
        	sys.exit(1)



def verify_green(section):
	required = ('Type',)
	allowed = ('Vacuum', 'UniformDielectric', 'MetalSphere', 'GreensFunctionSum')
	allowed_der = ('Numerical', 'Derivative', 'Gradient', 'Hessian')

	green1 = section.fetch_sect('Green<one>')
	green2 = section.fetch_sect('Green<two>')
	eps = section.get('Eps')
	epsimg = section.get('EpsImg')
	epsre = section.get('EpsRe')

	convert_length_array(section.get('SpherePosition'))
        position = section.get('SpherePosition')
	convert_length_scalar(section.get('SphereRadius'))
        radius = section.get('SphereRadius')

	type=section.get('Type')
	if (type.get() not in allowed):
		print "Allowed Green's functions are:", allowed
		sys.exit(1)

	type=section.get('Der')
	if (type.get() not in allowed_der):
		print "Allowed Derivatives are:", allowed
		sys.exit(1)

	if (type.get() == 'UniformDielectric'):
		if not eps.is_set():
			print "Eps not defined for UniformDielectric"
			sys.exit(1)

	if (type.get() == 'MetalSphere'):
		if not (eps.is_set() and epsre.is_set and epsimg.is_set()):
			print "Eps and/or EpsImg not defined for MetalSphere"
			sys.exit(1)
		if not (position.is_set() and radius.is_set()):
			print "SpherePosition and/or SphereRadius not defined for MetalSphere"
			sys.exit(1)
		if (len(position.get()) != 3):
			print "SpherePosition error"
			sys.exit(1)
		if (radius.get()  < 0.1):
			print "Minimum value allowed for Radius is 0.1"
			sys.exit(1)

	if (type.get() == 'GreensFunctionSum'):
		if not (green1.is_set() and green2.is_set()):
			print "One or both components not defined for GreensFunctionSum"
			sys.exit(1)

def verify_spheres(keyword):
	length=len(keyword.get())
	if (length % 4 != 0):
                print "Empty or incoherent Spheres list."
		sys.exit(1)
	convert_length_array(keyword)

def convert_length_array(keyword):
	length=len(keyword.get())
	if (isAngstrom):
		for i in range(length):
			keyword[i] *= toAtomicUnits

def convert_length_scalar(keyword):
	if (isAngstrom):
		keyword[0] *= toAtomicUnits

def convert_area_scalar(keyword):
	if (isAngstrom):
		keyword[0] *= toAtomicUnits * toAtomicUnits


def preprocess():
    """ Takes the PCM input file in @pcmsolver.inp, and preprocesses to make it machine-readable."""
    valid_keywords = setup_keywords()
    input=getkw.GetkwParser()
    inkw=input.parseFile('@pcmsolver.inp')
    inkw.sanitize(valid_keywords)
    topsect=inkw.get_topsect()
    inkw.run_callbacks(valid_keywords)

    xfile='@pcmsolver.inp'
    fd=open(xfile,'w')
    fd.write(str(inkw.top))
    fd.close()
