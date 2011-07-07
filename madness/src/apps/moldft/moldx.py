import math, sys

# Element info taken mostly from Rasmol but the effective covalent
# radii are from NWChem with those probably coming from Hondo.

cpk = {0: "light grey",
       1: "sky blue",
       2: "red",
       3: "yellow",
       4: "white",
       5: "pink",
       6: "gold",
       7: "blue",
       8: "orange",
       9: "dark grey",
       10: "brown",
       11: "purple",
       12: "deep pink",
       13: "green",
       14: "fire brick",
       15: "mid green"}

# Colors from X11/rgb.txt
colors = {
    "light grey": (211,211,211),
    "sky blue": (135,206,235),
    "red": (255,0,0),
    "yellow": (255,255,0),
    "white": (255,255,255),
    "pink": (255,105,180),
    "gold": (255,215,0),
    "blue": (0,0,255),
    "orange": (255,165,0),
    "dark grey": (169,169,169),
    "brown": (165,42,42),
    "purple": (160,32,240),
    "deep pink": (205,16,118),
    "green": (0,255,0),
    "fire brick": (178,34,34),
    "mid green": (60,179,113)}

for color in colors.keys():
    r,g,b = colors[color]
    colors[color]=r/255.0,g/255.0,b/255.0
    
class Element:
    def __init__(self,symbol,covrad,vdwrad,col,name,atno):
        self.symbol = symbol
        self.covrad = covrad
        self.vdwrad = vdwrad
        self.col    = colors[cpk[col]]
        self.name   = name
        self.atno   = atno

    def __str__(self):
        return "%3d %2s %16s %.1f %.1f (%.2f %.2f %.2f) " % (
            self.atno,self.symbol,self.name,self.covrad,self.vdwrad,
            self.col[0],self.col[1],self.col[2])

elements = range(104)
elements[  0] = Element('Bq' ,0.00, 0.0,  0, "Dummy"    ,  0 )
elements[  1] = Element('H' ,0.30, 2.75,  4, "HYDROGEN"    ,  1 )
elements[  2] = Element('He',1.22, 5.50,  5, "HELIUM"      ,  2 )
elements[  3] = Element('Li',1.23, 3.05, 14, "LITHIUM"     ,  3 )
elements[  4] = Element('Be',0.89, 1.57, 12, "BERYLLIUM"   ,  4 )
elements[  5] = Element('B' ,0.88, 3.87, 13, "BORON"       ,  5 )
elements[  6] = Element('C' ,0.77, 3.87,  9, "CARBON"      ,  6 ) # was color 0
elements[  7] = Element('N' ,0.70, 3.50,  1, "NITROGEN"    ,  7 )
elements[  8] = Element('O' ,0.66, 3.37,  2, "OXYGEN"      ,  8 )
elements[  9] = Element('F' ,0.58, 3.25,  6, "FLUORINE"    ,  9 )
elements[ 10] = Element('Ne',1.60, 5.05, 12, "NEON"        , 10 )
elements[ 11] = Element('Na',1.66, 5.50,  7, "SODIUM"      , 11 )
elements[ 12] = Element('Mg',1.36, 3.75, 15, "MAGNESIUM"   , 12 )
elements[ 13] = Element('Al',1.25, 3.75,  9, "ALUMINIUM"   , 13 )
elements[ 14] = Element('Si',1.17, 5.50,  6, "SILICON"     , 14 )
elements[ 15] = Element('P' ,1.10, 4.70,  8, "PHOSPHORUS"  , 15 )
elements[ 16] = Element('S' ,1.04, 4.52,  3, "SULPHUR"     , 16 )
elements[ 17] = Element('Cl',0.99, 4.37, 13, "CHLORINE"    , 17 )
elements[ 18] = Element('Ar',1.91, 6.92, 12, "ARGON"       , 18 )
elements[ 19] = Element('K' ,2.03, 5.97, 12, "POTASSIUM"   , 19 )
elements[ 20] = Element('Ca',1.74, 4.87,  9, "CALCIUM"     , 20 )
elements[ 21] = Element('Sc',1.44, 3.30, 12, "SCANDIUM"    , 21 )
elements[ 22] = Element('Ti',1.32, 4.87,  9, "TITANIUM"    , 22 )
elements[ 23] = Element('V' ,1.22, 2.65, 12, "VANADIUM"    , 23 )
elements[ 24] = Element('Cr',1.19, 2.82,  9, "CHROMIUM"    , 24 )
elements[ 25] = Element('Mn',1.17, 2.97,  9, "MANGANESE"   , 25 )
elements[ 26] = Element('Fe',1.165,4.87,  8, "IRON"        , 26 )
elements[ 27] = Element('Co',1.16, 2.82, 12, "COBALT"      , 27 )
elements[ 28] = Element('Ni',1.15, 3.10, 10, "NICKEL"      , 28 )
elements[ 29] = Element('Cu',1.17, 2.87, 10, "COPPER"      , 29 )
elements[ 30] = Element('Zn',1.25, 2.87, 10, "ZINC"        , 30 )
elements[ 31] = Element('Ga',1.25, 3.87, 12, "GALLIUM"     , 31 )
elements[ 32] = Element('Ge',1.22, 9.99, 12, "GERMANIUM"   , 32 )
elements[ 33] = Element('As',1.21, 2.07, 12, "ARSENIC"     , 33 )
elements[ 34] = Element('Se',1.17, 2.25, 12, "SELENIUM"    , 34 )
elements[ 35] = Element('Br',1.14, 4.37, 10, "BROMINE"     , 35 )
elements[ 36] = Element('Kr',1.98, 4.75, 12, "KRYPTON"     , 36 )
elements[ 37] = Element('Rb',2.22, 6.62, 12, "RUBIDIUM"    , 37 )
elements[ 38] = Element('Sr',1.92, 5.05, 12, "STRONTIUM"   , 38 )
elements[ 39] = Element('Y' ,1.62, 4.02, 12, "YTTRIUM"     , 39 )
elements[ 40] = Element('Zr',1.45, 3.55, 12, "ZIRCONIUM"   , 40 )
elements[ 41] = Element('Nb',1.34, 3.32, 12, "NIOBIUM"     , 41 )
elements[ 42] = Element('Mo',1.29, 4.37, 12, "MOLYBDENUM"  , 42 )
elements[ 43] = Element('Tc',1.27, 4.50, 12, "TECHNETIUM"  , 43 )
elements[ 44] = Element('Ru',1.24, 3.00, 12, "RUTHENIUM"   , 44 )
elements[ 45] = Element('Rh',1.25, 3.05, 12, "RHODIUM"     , 45 )
elements[ 46] = Element('Pd',1.28, 3.60, 12, "PALLADIUM"   , 46 )
elements[ 47] = Element('Ag',1.34, 3.87,  9, "SILVER"      , 47 )
elements[ 48] = Element('Cd',1.41, 4.37, 12, "CADMIUM"     , 48 )
elements[ 49] = Element('In',1.50, 3.62, 12, "INDIUM"      , 49 )
elements[ 50] = Element('Sn',1.40, 4.17, 12, "TIN"         , 50 )
elements[ 51] = Element('Sb',1.41, 2.80, 12, "ANTIMONY"    , 51 )
elements[ 52] = Element('Te',1.37, 3.15, 12, "TELLURIUM"   , 52 )
elements[ 53] = Element('I' ,1.33, 4.37, 11, "IODINE"      , 53 )
elements[ 54] = Element('Xe',2.09, 5.25, 12, "XENON"       , 54 )
elements[ 55] = Element('Cs',2.35, 7.52, 12, "CAESIUM"     , 55 )
elements[ 56] = Element('Ba',1.98, 6.02,  8, "BARIUM"      , 56 )
elements[ 57] = Element('La',1.69, 4.57, 12, "LANTHANUM"   , 57 )
elements[ 58] = Element('Ce',1.65, 4.65, 12, "CERIUM"      , 58 )
elements[ 59] = Element('Pr',1.65, 4.05, 12, "PRASEODYMIUM", 59 )
elements[ 60] = Element('Nd',1.64, 4.47, 12, "NEODYMIUM"   , 60 )
elements[ 61] = Element('Pm',1.65, 4.40, 12, "PROMETHIUM"  , 61 )
elements[ 62] = Element('Sm',1.66, 4.35, 12, "SAMARIUM"    , 62 )
elements[ 63] = Element('Eu',1.65, 4.90, 12, "EUROPIUM"    , 63 )
elements[ 64] = Element('Gd',1.61, 4.22, 12, "GADOLINIUM"  , 64 )
elements[ 65] = Element('Tb',1.59, 4.15, 12, "TERBIUM"     , 65 )
elements[ 66] = Element('Dy',1.59, 4.07, 12, "DYSPROSIUM"  , 66 )
elements[ 67] = Element('Ho',1.58, 4.02, 12, "HOLMIUM"     , 67 )
elements[ 68] = Element('Er',1.57, 3.97, 12, "ERBIUM"      , 68 )
elements[ 69] = Element('Tm',1.56, 3.92, 12, "THULIUM"     , 69 )
elements[ 70] = Element('Yb',1.56, 3.85, 12, "YTTERBIUM"   , 70 )
elements[ 71] = Element('Lu',1.56, 3.82, 12, "LUTETIUM"    , 71 )
elements[ 72] = Element('Hf',1.44, 3.50, 12, "HAFNIUM"     , 72 )
elements[ 73] = Element('Ta',1.34, 3.05, 12, "TANTALUM"    , 73 )
elements[ 74] = Element('W' ,1.30, 3.15, 12, "TUNGSTEN"    , 74 )
elements[ 75] = Element('Re',1.28, 3.25, 12, "RHENIUM"     , 75 )
elements[ 76] = Element('Os',1.26, 3.95, 12, "OSMIUM"      , 76 )
elements[ 77] = Element('Ir',1.26, 3.05, 12, "IRIDIUM"     , 77 )
elements[ 78] = Element('Pt',1.29, 3.87, 12, "PLATINUM"    , 78 )
elements[ 79] = Element('Au',1.34, 3.62,  6, "GOLD"        , 79 )
elements[ 80] = Element('Hg',1.44, 4.95, 12, "MERCURY"     , 80 )
elements[ 81] = Element('Tl',1.55, 4.27, 12, "THALLIUM"    , 81 )
elements[ 82] = Element('Pb',1.54, 5.40, 12, "LEAD"        , 82 )
elements[ 83] = Element('Bi',1.52, 4.32, 12, "BISMUTH"     , 83 )
elements[ 84] = Element('Po',1.53, 3.02, 12, "POLONIUM"    , 84 )
elements[ 85] = Element('At',1.50, 2.80, 12, "ASTATINE"    , 85 )
elements[ 86] = Element('Rn',2.20, 5.75, 12, "RADON"       , 86 )
elements[ 87] = Element('Fr',3.24, 8.10, 12, "FRANCIUM"    , 87 )
elements[ 88] = Element('Ra',2.68, 6.42, 12, "RADIUM"      , 88 )
elements[ 89] = Element('Ac',2.25, 5.30, 12, "ACTINIUM"    , 89 )
elements[ 90] = Element('Th',2.16, 4.60, 12, "THORIUM"     , 90 )
elements[ 91] = Element('Pa',1.93, 4.00, 12, "PROTACTINIUM", 91 )
elements[ 92] = Element('U' ,3.00, 4.37, 12, "URANIUM"     , 92 )
elements[ 93] = Element('Np',1.57, 4.27, 12, "NEPTUNIUM"   , 93 )
elements[ 94] = Element('Pu',1.81, 4.17, 12, "PLUTONIUM"   , 94 )
elements[ 95] = Element('Am',2.21, 4.15, 12, "AMERICIUM"   , 95 )
elements[ 96] = Element('Cm',1.43, 4.12, 12, "CURIUM"      , 96 )
elements[ 97] = Element('Bk',1.42, 4.10, 12, "BERKELIUM"   , 97 )
elements[ 98] = Element('Cf',1.40, 4.07, 12, "CALIFORNIUM" , 98 )
elements[ 99] = Element('Es',1.39, 4.05, 12, "EINSTEINIUM" , 99 )
elements[100] = Element('Fm',1.38, 4.02, 12, "FERMIUM"     ,100 )
elements[101] = Element('Md',1.37, 4.00, 12, "MENDELEVIUM" ,101 )
elements[102] = Element('No',1.36, 3.97, 12, "NOBELIUM"    ,102 )
elements[103] = Element('Lr',1.34, 3.95, 12, "LAWRENCIUM"  ,103 )

default_radii = range(104)
default_colors = range(104)
for i in range(1,len(elements)):
    default_radii[i]  = elements[i].covrad
    default_colors[i] = elements[i].col
    
def scale_coords(coords, fac):
    natom = len(coords)
    for atom in range(natom):
        for i in range(3):
            coords[atom][i] = coords[atom][i]*fac

def dist(i,j):
    return math.sqrt((i[0]-j[0])**2+
                     (i[1]-j[1])**2+
                     (i[2]-j[2])**2)

def make_bonds(coords, z, radii):
    bonds = []
    for i in range(len(coords)):
        for j in range(i):
            rij = dist(coords[i],coords[j])
            if rij < 1.1*(radii[z[i]]+radii[z[j]]):
                #print "bond",i,j,rij
                bonds.append([i,j])
    return bonds

def moldx(coords,z,bonds=None,units="angstrom",fname="mol.dx",
          radii=default_radii,colors=default_colors,
          scale=None,shift=None):
    '''

    Make an Opendx format file with a stick representation of the
    molecule.  Additional representations will be easy to add.

    coords = coords[natom][3] input coords (unchanged on output)

    z = z[natom] input atomic number of each atom (unchanged on output)

    bonds = optional list of precomputed connections.  If None, these
    .       are computed from the covalent radii.

    units = the input units.  If the input units are atomic, they will be
    .       converted to angstrom.  The output units are always angstrom
    .       unless additional scaling is specified.

    fname = the name of the output file

    radii = covalent radii for forming bonds.  Defaults will be used
    .       if not speficied.

    colors = colors for the atoms (3-vector float in [0,1]).  Defaults
    .        will be used if not specified.

    scale = if scale (scalar float) is specified, the *output*
    .       coordinates are scaled by that factor.

    shift = if shift (3-vector float) is specified, the *output*
    .       coordinates are shifted *after* scaling
    .       (e.g., x*scale+shift[0])

    '''

    au_to_angs = 0.529177249
    if len(z) != len(coords): raise ValueError

    if not bonds:
        # Need angstroms for making bonds
        if units in ["au","a.u.","atomic","bohr"]:
            scale_coords(coords,au_to_angs)

        bonds = make_bonds(coords,z,radii)

        if units in ["au","a.u.","atomic","bohr"]:
            scale_coords(coords,1.0/au_to_angs)

    # We have the atomic coordinates, numbers, radii, and
    # connectivity.  The opendx stick representation is formed
    # by splitting the bonds at the midpoint (weighted by
    # the ratio of the radii) with appropriate coloring.

    if not units in ["au","a.u.","atomic","bohr"]:
        # Need atomic units for the output
        scale_coords(coords,1.0/au_to_angs)

    pos = []
    col = []
    con = []
    
    for i,j in bonds:
        npos = len(pos)
        rij = dist(coords[i],coords[j])
        si = radii[z[i]]/(radii[z[i]]+radii[z[j]])
        sj = 1.0-si
        mid = [coords[i][0]*sj+coords[j][0]*si,
               coords[i][1]*sj+coords[j][1]*si,
               coords[i][2]*sj+coords[j][2]*si]
        
        pos.append(coords[i]); col.append(colors[z[i]])
        pos.append(mid); col.append(colors[z[i]])

        pos.append(mid); col.append(colors[z[j]])
        pos.append(coords[j]); col.append(colors[z[j]])

        con.append([npos,npos+1])
        con.append([npos+2,npos+3])

    f = open(fname,'w+')

    f.write("object 1 class array type float rank 1 shape 3 items %d data follows\n" % len(pos))
    for x,y,z in pos:
        if scale: x, y, z = x*scale, y*scale, z*scale
        if shift: x, y, z = x+shift[0], y+shift[1], z+shift[2]
        f.write("%f %f %f\n" % (x,y,z))
    f.write('attribute "dep" string "positions"\n')

    f.write("object 2 class array type int rank 1 shape 2 items %d data follows\n" % len(con))
    for i,j in con: f.write("%d %d \n" % (i,j))
    f.write('attribute "element type" string "lines"\n')
    f.write('attribute "ref" string "positions"\n')
    f.write('attribute "dep" string "connections"\n')

    f.write("object 3 class array type float rank 0 items %d data follows\n" % len(pos))
    for p in pos: f.write("0 ")
    f.write("\n")
    f.write('attribute "dep" string "positions"\n')

    f.write("object 4 class array type float rank 1 shape 3 items %d data follows\n" % len(col))
    for r,g,b in col: f.write("%f %f %f\n" % (r,g,b))
    f.write('attribute "dep" string "positions"\n')

    f.write('object "irregular positions irregular connections" class field\n')
    f.write('component "positions" value 1\n')
    f.write('component "connections" value 2\n')
    f.write('component "data" value 3\n')
    f.write('component "colors" value 4\n')
    f.write('attribute "name" string "irregular positions irregular connections"\n')
    f.write('end\n')

    f.close()

    # Undo-Force a.u. while making the bonds (cannot do until after
    # writing file since pos contains references to these values).
    if not units in ["au","a.u.","atomic","bohr"]:
        scale_coords(coords,au_to_angs)
    

def symtoatn(sym):
    sym = sym.lower()
    for el in elements:
        if el.symbol.lower() == sym:
            return el.atno
    raise "Uh?"

if __name__ == "__main__":
    # Read first molecule from the madness output file
    f = sys.stdin
    while 1:
        line = f.readline().split()
        if len(line) and line[0] == 'geometry':
            coords = []
            z = []
            units = 'au'
            while 1:
                line = f.readline().split()
                if len(line) == 0: continue
                if line[0] == 'end':
                    break
                elif line[0] == 'units':
                    units = line[1]
                elif not (line[0] == 'units' or line[0] == 'eprec'):
                    z.append(symtoatn(line[0]))
                    coords.append([float(line[1]), float(line[2]), float(line[3])])
            moldx(coords, z, units=units, fname='molecule.dx')
            sys.exit(0)

    
