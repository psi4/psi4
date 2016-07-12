#!/usr/bin/env python

# Francesco Evangelista
# Emory University

from __future__ import print_function
import argparse
import sys
import re
import subprocess
import os
import datetime

from os import listdir, environ
from os.path import isfile, join

vmd_cube_help = """vmd_cube is a script to render cube files with vmd.
To generate cube files with Psi4, add the command cubeprop(wfn) at the end
of your input file, where *wfn* is a Wavefunction object that may be
retrieved from any calculation and used following the pattern "ene, wave =
energy('pbe', return_wfn=True)\\n cubeprop(wave)"."""

vmd_exe = ""

vmd_script_name = "vmd_mo_script.vmd"

vmd_template = """#
# VMD script to plot MOs from cube files
#

# Load the molecule and change the atom style
mol load cube PARAM_CUBEFILE.cube
mol modcolor 0 PARAM_CUBENUM Element
mol modstyle 0 PARAM_CUBENUM CPK 0.400000 0.40000 30.000000 16.000000

# Define the material
material change ambient Opaque 0.310000
material change diffuse Opaque 0.720000
material change specular Opaque 0.500000
material change shininess Opaque 0.480000
material change opacity Opaque 1.000000
material change outline Opaque 0.000000
material change outlinewidth Opaque 0.000000
material change transmode Opaque 0.000000
material change specular Opaque 0.750000

material change ambient EdgyShiny 0.310000
material change diffuse EdgyShiny 0.720000
material change shininess EdgyShiny 1.0000
material change opacity EdgyShiny PARAM_OPACITY

# Customize atom colors
color Element C silver
color Element H white

# Rotate and translate the molecule
rotate x by PARAM_RX
rotate y by PARAM_RY
rotate z by PARAM_RZ
translate by PARAM_TX PARAM_TY PARAM_TZ
scale by PARAM_SCALE

# Eliminate the axis and perfect the view
axes location Off
display projection Orthographic
display depthcue off
color Display Background white"""


vmd_template_surface = """#
# Add the surfaces
mol color ColorID PARAM_SURF1ID
mol representation Isosurface PARAM_ISOVALUE1 0 0 0 1 1
mol selection all
mol material EdgyShiny
mol addrep PARAM_CUBENUM
mol color ColorID PARAM_SURF2ID
mol representation Isosurface PARAM_ISOVALUE2 0 0 0 1 1
mol selection all
mol material EdgyShiny
mol addrep PARAM_CUBENUM

# Render
render TachyonInternal PARAM_CUBEFILE.tga
mol delete PARAM_CUBENUM
"""

vmd_template_rotate = """
light 1 off
light 0 rot y  30.0
light 0 rot x -30.0
"""

default_path = os.getcwd()

# Default parameters
options = {"SURF1ID"    : [None,"Surface1 Color Id"],
           "SURF2ID"    : [None,"Surface2 Color Id"],
           "ISOVALUE1"  : [None,"Isosurface1 Value"],
           "ISOVALUE2"  : [None,"Isosurface2 Value"],
           "RX"         : [None,"X-axis Rotation"],
           "RY"         : [None,"Y-axis Rotation"],
           "RZ"         : [None,"Z-axis Rotation"],
           "TX"         : [None,"X-axis Translation"],
           "TY"         : [None,"Y-axis Translation"],
           "TZ"         : [None,"Z-axis Translation"],
           "OPACITY"    : [None,"Opacity"],
           "CUBEDIR"    : [None,"Cubefile Directory"],
           "SCALE"      : [None,"Scaling Factor"],
           "MONTAGE"    : [None,"Montage"],
           "FONTSIZE"   : [None,"Font size"],
           "IMAGESIZE"  : [None,"Image size"],
           "VMDPATH"    : [None,"VMD Path"]}


def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def multigsub(subs,str):
    for k,v in subs.items():
        str = re.sub(k,v,str)
    return str


def find_vmd(options):
    if environ['VMDPATH']:
        vmdpath = environ['VMDPATH']
        vmdpath = multigsub({" " : r"\ "},vmdpath)
        options["VMDPATH"][0] = vmdpath
    else:
        print("Please set the VMDPATH environmental variable to the path of VMD.")
        exit(1)


def save_setup_command(argv):
    file_name = join(default_path, 'vmd_cube_command')
    f = open(file_name, 'w')
    f.write('# setup command was executed '+datetime.datetime.now().strftime("%d-%B-%Y %H:%M:%S"+"\n"))
    f.write(" ".join(argv[:])+"\n")
    f.close()


def read_options(options):
    parser = argparse.ArgumentParser(description=vmd_cube_help)
    parser.add_argument('data', metavar='<cubefile dir>', type=str, nargs='?',default=".",
                   help='The directory containing the cube files.')
                   
    parser.add_argument('--color1', metavar='<integer>', type=int, nargs='?',default=3,
                   help='the color ID of surface 1 (integer, default = 3)')
    parser.add_argument('--color2', metavar='<integer>', type=int, nargs='?',default=23,
                   help='the color ID of surface 2 (integer, default = 23)')
                   
    parser.add_argument('--iso', metavar='<isovalue>', type=float, nargs='?',default=0.05,
                   help='the isosurface value (float, default = 0.05)')
                   
    parser.add_argument('--rx', metavar='<angle>', type=float, nargs='?',default=30.0,
                   help='the x-axis rotation angle (float, default = 30.0)')
    parser.add_argument('--ry', metavar='<angle>', type=float, nargs='?',default=40.0,
                   help='the y-axis rotation angle (float, default = 40.0)')
    parser.add_argument('--rz', metavar='<angle>', type=float, nargs='?',default=15.0,
                   help='the z-axis rotation angle (float, default = 15.0)')
                   
    parser.add_argument('--tx', metavar='<angle>', type=float, nargs='?',default=0.0,
                   help='the x-axis translation (float, default = 0.0)')
    parser.add_argument('--ty', metavar='<angle>', type=float, nargs='?',default=0.0,
                   help='the y-axis translation (float, default = 0.0)')
    parser.add_argument('--tz', metavar='<angle>', type=float, nargs='?',default=0.0,
                   help='the z-axis translation (float, default = 0.0)')

    parser.add_argument('--opacity', metavar='<opacity>', type=float, nargs='?',default=1.0,
                   help='opacity of the isosurface (float, default = 1.0)')

    parser.add_argument('--scale', metavar='<factor>', type=float, nargs='?',default=1.0,
                   help='the scaling factor (float, default = 1.0)')
    parser.add_argument('--montage', const=True, default=False, nargs='?',
                   help='call montage to combine images. (string, default = false)')

    parser.add_argument('--imagesize', metavar='<integer>', type=int, nargs='?',default=250,
                   help='the size of each image (integer, default = 250)')
    parser.add_argument('--fontsize', metavar='<integer>', type=int, nargs='?',default=20,
                   help='the font size (integer, default = 20)')


    args = parser.parse_args()

    options["CUBEDIR"][0] = str(args.data)
    options["SURF1ID"][0] = str(args.color1)
    options["SURF2ID"][0] = str(args.color2)
    options["ISOVALUE1"][0] = str(args.iso)
    options["ISOVALUE2"][0] = str(-args.iso)
    options["RX"][0] = str(args.rx)
    options["RY"][0] = str(args.ry)
    options["RZ"][0] = str(args.rz)
    options["TX"][0] = str(args.tx)
    options["TY"][0] = str(args.ty)
    options["TZ"][0] = str(args.tz)
    options["OPACITY"][0] = str(args.opacity)
    options["SCALE"][0] = str(args.scale)
    options["MONTAGE"][0] = str(args.montage)
    options["FONTSIZE"][0] = str(args.fontsize)
    options["IMAGESIZE"][0] = str(args.imagesize)

    print("Parameters:")
    for k,v in options.items():
        print("  %-20s %s" % (v[1],v[0]))


def find_cubes(options):
    # Find all the cube files in a given directory
    dir = options["CUBEDIR"][0]
    sorted_files = []
    for f in listdir(options["CUBEDIR"][0]):
        if ".cube" in f:
            sorted_files.append(f)

    return sorted(sorted_files)


def write_and_run_vmd_script(options,cube_files):
    vmd_script = open(vmd_script_name,"w+")
    vmd_script.write(vmd_template_rotate)

    # Define a map that contains all the values of the VMD parameters
    replacement_map = {}
    for k,v in options.items():
        key = "PARAM_" + k.upper()
        replacement_map[key] = v[0]

    for n,f in enumerate(cube_files):
        replacement_map["PARAM_CUBENUM"] = "%03d" % n
        replacement_map["PARAM_CUBEFILE"] = options["CUBEDIR"][0] + "/" + f[:-5]

        vmd_script_surface = multigsub(replacement_map,vmd_template_surface)
        vmd_script_head = multigsub(replacement_map,vmd_template)
        vmd_script.write(vmd_script_head + "\n" + vmd_script_surface)

    vmd_script.write("quit")
    vmd_script.close()

    # Call VMD
    FNULL = open(os.devnull, 'w')
    subprocess.call(("%s -dispdev text -e %s" % (options["VMDPATH"][0],vmd_script_name)),stdout=FNULL, shell=True)


def call_montage(options,cube_files):
    if options["MONTAGE"][0] == 'True':
        # Optionally, combine all figures into one image using montage
        montage_exe = which("montage")
        if montage_exe:
            alpha_mos = []
            beta_mos = []
            densities = []
            basis_functions = []
            for f in cube_files:
                tga_file = f[:-5] + ".tga"
                if "Psi_a" in f:
                    alpha_mos.append(tga_file)
                if "Psi_b" in f:
                    beta_mos.append(tga_file)
                if "D" in f:
                    densities.append(tga_file)
                if "Phi" in f:
                    basis_functions.append(tga_file)


            # Sort the MOs
            sorted_mos = []
            for set in [alpha_mos,beta_mos]:
                sorted_set = []
                for s in set:
                    s_split = s.split('_')
                    sorted_set.append((int(s_split[2]),"Psi_a_%s_%s" % (s_split[2],s_split[3])))
                sorted_set = sorted(sorted_set)
                sorted_mos.append([s[1] for s in sorted_set])
                    
            os.chdir(options["CUBEDIR"][0])

            for f in sorted_mos[0]:
                f_split = f.split('_')
                label = "%s\ \(%s\)" % (f_split[3][:-4],f_split[2])
                subprocess.call(("montage -pointsize %s -label %s %s -geometry '%sx%s+0+0>' %s" %
                    (options["FONTSIZE"][0],label,f,options["IMAGESIZE"][0],options["IMAGESIZE"][0],f)), shell=True)
            
            if len(alpha_mos) > 0:
                subprocess.call(("%s %s -geometry +2+2 AlphaMOs.tga" % (montage_exe," ".join(sorted_mos[0]))), shell=True)
            if len(beta_mos) > 0:
                subprocess.call(("%s %s -geometry +2+2 BetaMOs.tga" % (montage_exe," ".join(sorted_mos[1]))), shell=True)
            if len(densities) > 0:
                subprocess.call(("%s %s -geometry +2+2 Densities.tga" % (montage_exe," ".join(densities))), shell=True)
            if len(basis_functions) > 0:
                subprocess.call(("%s %s -geometry +2+2 BasisFunctions.tga" % (montage_exe," ".join(basis_functions))), shell=True)


def main(argv):
    read_options(options)
    find_vmd(options)
    save_setup_command(argv)
    cube_files = find_cubes(options)
    write_and_run_vmd_script(options,cube_files)
    call_montage(options,cube_files)


if __name__ == '__main__':
    main(sys.argv)
