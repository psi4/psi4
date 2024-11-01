#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2017 Francesco Evangelista
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import print_function

import argparse
import datetime
import os
import re
import subprocess
import sys
from os import environ, listdir
from os.path import isfile, join

vmd_cube_help = """vmd_cube is a script to render cube files with vmd.
To generate cube files with Psi4 add the command cubeprop() at the end of your input file."""

vmd_exe = ""

vmd_script_name = "vmd_mo_script.vmd"

vmd_template = """#
# VMD script to plot MOs from cube files
#

# Load the molecule and change the atom style
mol load cube PARAM_CUBEFILE.cube
mol modcolor 0 PARAM_CUBENUM Element
mol modstyle 0 PARAM_CUBENUM Licorice 0.110000 10.000000 10.000000
#mol modstyle 0 PARAM_CUBENUM CPK 0.400000 0.40000 30.000000 16.000000

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

material change ambient   EdgyShiny 0.310000
material change diffuse   EdgyShiny 0.720000
material change shininess EdgyShiny 1.0000
material change opacity   EdgyShiny PARAM_OPACITY

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
display resize PARAM_IMAGEW PARAM_IMAGEH
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
"""

vmd_template_interactive = """#
# Disable rendering
mol off PARAM_CUBENUM
"""

vmd_template_render = """
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
           "LABEL_MOS"  : [None,"Label MOs"],
           "FONTSIZE"   : [None,"Font size"],
           "IMAGEW"     : [None,"Image width"],
           "IMAGEH"     : [None,"Image height"],
           "VMDPATH"    : [None,"VMD Path"],
           "INTERACTIVE": [None,"Interactive Mode"],
           "GZIP"       : [None,"Gzip Cube Files"]}


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
                   
    parser.add_argument('--tx', metavar='<length>', type=float, nargs='?',default=0.0,
                   help='the x-axis translation (float, default = 0.0)')
    parser.add_argument('--ty', metavar='<length>', type=float, nargs='?',default=0.0,
                   help='the y-axis translation (float, default = 0.0)')
    parser.add_argument('--tz', metavar='<length>', type=float, nargs='?',default=0.0,
                   help='the z-axis translation (float, default = 0.0)')

    parser.add_argument('--opacity', metavar='<opacity>', type=float, nargs='?',default=1.0,
                   help='opacity of the isosurface (float, default = 1.0)')

    parser.add_argument('--scale', metavar='<factor>', type=float, nargs='?',default=1.0,
                   help='the scaling factor (float, default = 1.0)')
    parser.add_argument('--no-montage', action="store_true",
                   help='call montage to combine images. (string, default = false)')
    parser.add_argument('--no-labels', action="store_true",
                   help='do not add labels to images. (string, default = false)')

    parser.add_argument('--imagesize', metavar='<integer>', type=int, nargs='?',default=250,
                   help='the size of each image (integer, default = 250)')
    parser.add_argument('--imagew', metavar='<integer>', type=int, nargs='?',default=250,
                   help='the width of images (integer, default = 250)')
    parser.add_argument('--imageh', metavar='<integer>', type=int, nargs='?',default=250,
                   help='the height of images (integer, default = 250)')
    parser.add_argument('--fontsize', metavar='<integer>', type=int, nargs='?',default=20,
                   help='the font size (integer, default = 20)')

    parser.add_argument('--interactive', action="store_true",
                   help='run in interactive mode (default = false)')

    parser.add_argument('--gzip', action="store_true",
                   help='gzip cube files (default = false)')

    parser.add_argument('--national_scheme', action="store_true",
                   help='use a red/blue color scheme. (string, default = false)')
    parser.add_argument('--silver_scheme', action="store_true",
                   help='use a gray/white color scheme. (string, default = false)')
    parser.add_argument('--bright_scheme', action="store_true",
                   help='use a soft yellow/blue color scheme. (string, default = false)')
    parser.add_argument('--electron_scheme', action="store_true",
                   help='use a purple/green color scheme. (string, default = false)')

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
    options["LABEL_MOS"][0] = str(not args.no_labels)
    options["MONTAGE"][0] = str(not args.no_montage)
    options["FONTSIZE"][0] = str(args.fontsize)
    options["IMAGEW"][0] = str(args.imagew)
    options["IMAGEH"][0] = str(args.imageh)
    options["INTERACTIVE"][0] = str(args.interactive)
    options["GZIP"][0] = str(args.gzip)

    if args.national_scheme:
        options["SURF1ID"][0] = '23'
        options["SURF2ID"][0] = '30'

    if args.silver_scheme:
        options["SURF1ID"][0] = '2'
        options["SURF2ID"][0] = '8'

    if args.electron_scheme:
        options["SURF1ID"][0] = '13'
        options["SURF2ID"][0] = '12'

    if args.bright_scheme:
        options["SURF1ID"][0] = '32'
        options["SURF2ID"][0] = '22'

    print("Parameters:")
    sorted_parameters = sorted(options.keys())
    for k in sorted_parameters:
        print("  %-20s %s" % (options[k][1],options[k][0]))

def find_cubes(options):
    # Find all the cube files in a given directory
    dir = options["CUBEDIR"][0]
    sorted_files = []
    zipped_files = []

    for f in listdir(options["CUBEDIR"][0]):
        if "\'" in f:
            nf = f.replace("\'", "p")
            os.rename(f,nf)
            f = nf
        if "\"" in f:
            nf = f.replace("\"", "pp")
            os.rename(f,nf)
            f = nf
        if f[-5:] == '.cube':
            sorted_files.append(f)
        elif f[-8:] == '.cube.gz':
            found_zipped = True
            # unzip file
            sorted_files.append(f[:-3])
            zipped_files.append(f)

    if len(zipped_files) > 0:
        print("\nDecompressing gzipped cube files")
        FNULL = open(os.devnull, 'w')
        subprocess.call(("gzip -d %s" % " ".join(zipped_files)),stdout=FNULL, shell=True)
        options["GZIP"][0] = 'True'

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
        
        if options["INTERACTIVE"][0] == 'True':
            vmd_script_render = multigsub(replacement_map,vmd_template_interactive)
        else:
            vmd_script_render = multigsub(replacement_map,vmd_template_render)

        vmd_script.write(vmd_script_head + "\n" + vmd_script_surface + "\n" + vmd_script_render)

    if options["INTERACTIVE"][0] == 'False':
        vmd_script.write("quit")
        vmd_script.close()
        # Call VMD in text mode
        FNULL = open(os.devnull, 'w')
        subprocess.call(("%s -dispdev text -e %s" % (options["VMDPATH"][0],vmd_script_name)),stdout=FNULL, shell=True)
    else:
        vmd_script.close()
        # Call VMD in graphic mode
        FNULL = open(os.devnull, 'w')
        subprocess.call(("%s -e %s" % (options["VMDPATH"][0],vmd_script_name)),stdout=FNULL, shell=True)


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
                    
            # Add labels
            if options["LABEL_MOS"][0] == 'True':
                for f in sorted_mos[0]:
                    f_split = f.split('_')
                    label = r'%s\ \(%s\)' % (f_split[3][:-4],f_split[2])
                    subprocess.call(("montage -pointsize %s -label %s %s -geometry '%sx%s+0+0>' %s" %
                        (options["FONTSIZE"][0],label,f,options["IMAGEW"][0],options["IMAGEH"][0],f)), shell=True)

            # Combine together in one image
            if len(alpha_mos) > 0:
                subprocess.call(("%s %s -geometry +2+2 AlphaMOs.tga" % (montage_exe," ".join(sorted_mos[0]))), shell=True)
            if len(beta_mos) > 0:
                subprocess.call(("%s %s -geometry +2+2 BetaMOs.tga" % (montage_exe," ".join(sorted_mos[1]))), shell=True)
            if len(densities) > 0:
                subprocess.call(("%s %s -geometry +2+2 Densities.tga" % (montage_exe," ".join(densities))), shell=True)
            if len(basis_functions) > 0:
                subprocess.call(("%s %s -geometry +2+2 BasisFunctions.tga" % (montage_exe," ".join(basis_functions))), shell=True)


def zip_files(cube_files,options):
    """Gzip cube files if requested or necessary."""
    if options["GZIP"][0] == 'True':
        print("\nCompressing cube files")
        FNULL = open(os.devnull, 'w')
        subprocess.call(("gzip %s" % " ".join(cube_files)),stdout=FNULL, shell=True)


def main(argv):
    find_vmd(options)
    read_options(options)
    save_setup_command(argv)
    cube_files = find_cubes(options)
    write_and_run_vmd_script(options,cube_files)
    call_montage(options,cube_files)
    zip_files(cube_files,options)

if __name__ == '__main__':
    main(sys.argv)
