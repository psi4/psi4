#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from __future__ import absolute_import

import sys, os
import math
from math import *
from datetime import date

class InPsight:

    # POV-Ray defines
    defines = {}
    defines['Shadows'] = 'false'
    defines['Background_Color'] = '<0.6,0.6,0.6>'
    defines['Output_File_Type'] = 'N'
    defines['Output_Alpha'] = 'true'
    defines['Light_Color'] = '<1,1,1>'
    defines['Filename'] = 'inpsight'
    defines['Filepath'] = os.getcwd()
    defines['Antialias'] = 'true'
    defines['Antialias_Threshold'] = '0.1'

    # Molecule geometry
    atoms = [] # (Z,x,y,z,R,r,g,b,t) in bohr
    bonds = [] # (x1,y1,z1,R1,x2,y2,z2,R2,r,g,b,t)

    # Molecular geometry defines
    colors = []
    radii = []
    radial_scale = 0.25
    bond_width = 0.2 # bohr
    bohr_per_ang = 1.8897161646320724
    bonding_alpha = 0.65 # Used to select/reject bonds via sum of vDW radii

    # View defines (high-level)
    azimuth = 0.0
    elevation = 0.0
    zoom = 0.5
    height = 900
    width = 1200

    # Camera positions (low-level)
    location = [1.0,0.0,0.0]
    up = [0.0,0.75,0.0]
    right = [1.0,0.0,0.0]
    sky = [0.0,-1.0,0.0]
    look_at = [0.0,0.0,0.0]
    light = [1.0,0.0,0.0]
    light_color = [0.6,0.6,0.6]

    # Standard Jmol colors, 256-based
    colors.append([0,0,0])
    colors.append([255,255,255])
    colors.append([217,255,255])
    colors.append([204,128,255])
    colors.append([194,255,0])
    colors.append([255,181,181])
    colors.append([144,144,144])
    colors.append([48,80,248])
    colors.append([255,13,13])
    colors.append([144,224,80])
    colors.append([179,227,245])
    colors.append([171,92,242])
    colors.append([138,255,0])
    colors.append([191,166,166])
    colors.append([240,200,160])
    colors.append([255,128,0])
    colors.append([255,255,48])
    colors.append([31,240,31])
    colors.append([128,209,227])
    colors.append([143,64,212])
    colors.append([61,255,0])
    colors.append([230,230,230])
    colors.append([191,194,199])
    colors.append([166,166,171])
    colors.append([138,153,199])
    colors.append([156,122,199])
    colors.append([224,102,51])
    colors.append([240,144,160])
    colors.append([80,208,80])
    colors.append([200,128,51])
    colors.append([125,128,176])
    colors.append([194,143,143])
    colors.append([102,143,143])
    colors.append([189,128,227])
    colors.append([255,161,0])
    colors.append([166,41,41])
    colors.append([92,184,209])
    colors.append([112,46,176])
    colors.append([0,255,0])
    colors.append([148,255,255])
    colors.append([148,224,224])
    colors.append([115,194,201])
    colors.append([84,181,181])
    colors.append([59,158,158])
    colors.append([36,143,143])
    colors.append([10,125,140])
    colors.append([0,105,133])
    colors.append([192,192,192])
    colors.append([255,217,143])
    colors.append([166,117,115])
    colors.append([102,128,128])
    colors.append([158,99,181])
    colors.append([212,122,0])
    colors.append([148,0,148])
    colors.append([66,158,176])
    colors.append([87,23,143])
    colors.append([0,201,0])
    colors.append([112,212,255])
    colors.append([255,255,199])
    colors.append([217,255,199])
    colors.append([199,255,199])
    colors.append([163,255,199])
    colors.append([143,255,199])
    colors.append([97,255,199])
    colors.append([69,255,199])
    colors.append([48,255,199])
    colors.append([31,255,199])
    colors.append([0,255,156])
    colors.append([0,230,117])
    colors.append([0,212,82])
    colors.append([0,191,56])
    colors.append([0,171,36])
    colors.append([77,194,255])
    colors.append([77,166,255])
    colors.append([33,148,214])
    colors.append([38,125,171])
    colors.append([38,102,150])
    colors.append([23,84,135])
    colors.append([208,208,224])
    colors.append([255,209,35])
    colors.append([184,184,208])
    colors.append([166,84,77])
    colors.append([87,89,97])
    colors.append([158,79,181])
    colors.append([171,92,0])
    colors.append([117,79,69])
    colors.append([66,130,150])
    colors.append([66,0,102])
    colors.append([0,125,0])
    colors.append([112,171,250])
    colors.append([0,186,255])
    colors.append([0,161,255])
    colors.append([0,143,255])
    colors.append([0,128,255])
    colors.append([0,107,255])
    colors.append([84,92,242])
    colors.append([120,92,227])
    colors.append([138,79,227])
    colors.append([161,54,212])
    colors.append([179,31,212])
    colors.append([179,31,186])
    colors.append([179,13,166])
    colors.append([189,13,135])
    colors.append([199,0,102])
    colors.append([204,0,89])
    colors.append([209,0,79])
    colors.append([217,0,69])
    colors.append([224,0,56])
    colors.append([230,0,46])
    colors.append([235,0,38])

    # Approximate vDW radii in angstrom
    radii.append(2.0)
    radii.append(1.001)
    radii.append(1.012)
    radii.append(0.825)
    radii.append(1.408)
    radii.append(1.485)
    radii.append(1.452)
    radii.append(1.397)
    radii.append(1.342)
    radii.append(1.287)
    radii.append(1.243)
    radii.append(1.144)
    radii.append(1.364)
    radii.append(1.639)
    radii.append(1.716)
    radii.append(1.705)
    radii.append(1.683)
    radii.append(1.639)
    radii.append(1.595)
    radii.append(1.485)
    radii.append(1.474)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.562)
    radii.append(1.650)
    radii.append(1.727)
    radii.append(1.760)
    radii.append(1.771)
    radii.append(1.749)
    radii.append(1.727)
    radii.append(1.628)
    radii.append(1.606)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.639)
    radii.append(1.672)
    radii.append(1.804)
    radii.append(1.881)
    radii.append(1.892)
    radii.append(1.892)
    radii.append(1.881)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)
    radii.append(2.0)

    def __init__(self,molecule):
        self.molecule = molecule
        self.molecule.update_geometry()
        self.update_geometry()

    def update_geometry(self):

        # Atoms
        natom = self.molecule.natom()
        self.atoms = []
        for k in range(0,natom):
            x = self.molecule.x(k)
            y = self.molecule.y(k)
            z = self.molecule.z(k)
            Z = self.molecule.Z(k)
            atom = Z, x, y, z, self.radial_scale * self.bohr_per_ang * self.radii[Z], self.colors[Z][0] / 256.0, \
                self.colors[Z][1] / 256.0, self.colors[Z][2] / 256.0, 0.0

            self.atoms.append(atom)

        # Bonds
        self.bonds = []
        for k in range(1,natom):
            for l in range (0, k):
                Z1 = self.atoms[k][0]
                Z2 = self.atoms[l][0]
                R1 = self.bohr_per_ang*self.radii[Z1]
                R2 = self.bohr_per_ang*self.radii[Z2]
                x1 = self.atoms[k][1]
                y1 = self.atoms[k][2]
                z1 = self.atoms[k][3]
                x2 = self.atoms[l][1]
                y2 = self.atoms[l][2]
                z2 = self.atoms[l][3]

                r1 = self.atoms[k][5];
                g1 = self.atoms[k][6];
                b1 = self.atoms[k][7];
                t1 = self.atoms[k][8];
                r2 = self.atoms[l][5];
                g2 = self.atoms[l][6];
                b2 = self.atoms[l][7];
                t2 = self.atoms[l][8];

                R = sqrt((x1-x2)*(x1-x2) + \
                         (y1-y2)*(y1-y2) + \
                         (z1-z2)*(z1-z2))

                if (R < self.bonding_alpha*(R1 + R2)):
                    omega = R2 / (R1 + R2)
                    xc = omega * (x1 - x2) + x2
                    yc = omega * (y1 - y2) + y2
                    zc = omega * (z1 - z2) + z2
                    bond1 = x1,y1,z1,self.bond_width, xc,yc,zc,self.bond_width,r1,g1,b1,t1
                    bond2 = x2,y2,z2,self.bond_width, xc,yc,zc,self.bond_width,r2,g2,b2,t2

                    self.bonds.append(bond1)
                    self.bonds.append(bond2)

    def set_define(self, key, value):
        self.defines[key] = value

    def set_color(self, Z, color):
        self.colors[Z] = color

    def set_radius(self, Z, radius):
        self.radii[Z] = radius

    def position_camera(self):
        xc = self.molecule.center_of_mass()
        self.look_at = [xc[0], xc[1], xc[2]]
        Rmax = 0.0
        natom = self.molecule.natom()
        for k in range(0,natom):
            x = [self.molecule.x(k), self.molecule.y(k), self.molecule.z(k)]
            R = sqrt((x[0] - xc[0])*(x[0] - xc[0]) + \
                     (x[1] - xc[1])*(x[1] - xc[1]) + \
                     (x[2] - xc[2])*(x[2] - xc[2]))
            if R > Rmax:
                Rmax = R
        Rmax = Rmax / self.zoom
        if (self.width < self.height):
            self.right = [Rmax, 0.0, 0.0]
            self.up = [0.0, self.right[0]*self.height/self.width, 0.0]
        else:
            self.up = [0.0, Rmax, 0.0]
            self.right = [self.up[1]*self.width/self.height, 0.0, 0.0]

        phi = math.pi*(-self.azimuth)/180.0
        theta = math.pi*(90.0 - self.elevation)/180.0
        delta = [Rmax*cos(phi)*sin(theta), Rmax*sin(phi)*sin(theta), Rmax*cos(theta)]
        self.location = [xc[0] + delta[0], xc[1] + delta[1], xc[2] + delta[2]]
        phi = math.pi*(-(self.azimuth + 30.0))/180.0
        theta = math.pi*(90.0 - (self.elevation + 30.0))/180.0
        delta = [Rmax*cos(phi)*sin(theta), Rmax*sin(phi)*sin(theta), Rmax*cos(theta)]
        self.light = [xc[0] + delta[0], xc[1] + delta[1], xc[2] + delta[2]]

    def set_view(self,azimuth, elevation, zoom = 0.7):
        self.azimuth = azimuth
        self.elevation = elevation
        self.zoom = zoom

        self.position_camera()

    def set_size(self, width,height):
        self.width = width
        self.height = height

    def set_camera(self, location, sky, up, right, look_at, light, light_color):
        self.location = location
        self.sky  = sky
        self.up = up
        self.right = right
        self.look_at = look_at
        self.light = light
        self.light_color = light_color

    def save_molecule(self, filename):
        if (filename != ''):
            self.defines['Filename'] = filename

        ini_filename = self.defines['Filepath'] + '/' + self.defines['Filename'] + '.pov.ini'
        pov_filename = self.defines['Filepath'] + '/' + self.defines['Filename'] + '.pov'
        png_filename = self.defines['Filepath'] + '/' + self.defines['Filename'] + '.png'
        pov_file = self.defines['Filename'] + '.pov'
        png_file = self.defines['Filename'] + '.png'

        # Write the pov.ini file
        fh = open(ini_filename,'w')
        fh.write('; InPsight: visualization in Psi4\n')
        fh.write(';  by Rob Parrish\n')
        fh.write('; .pov.ini file\n')
        fh.write('; Created %s\n' % str(date.today()))
        fh.write('\n')
        fh.write('Input_File_Name=%s\n' % pov_file)
        fh.write('Output_to_File=true\n')
        fh.write('Output_File_Type=%s\n' % self.defines['Output_File_Type'])
        fh.write('Output_File_Name=%s\n' % png_file)
        fh.write('Height=%s\n' % str(self.height))
        fh.write('Width=%s\n' % str(self.width))
        fh.write('Output_Alpha=%s\n' % self.defines['Output_Alpha'])
        fh.write('Antialias=%s\n' % self.defines['Antialias'])
        fh.write('Antialias_Threshold=%s\n' % self.defines['Antialias_Threshold'])
        fh.write('Display=true\n')
        fh.write('Warning_Level=5\n')
        fh.write('Verbose=false\n')

        fh.close()

        # Write the pov file
        fh = open(pov_filename, 'w')

        fh.write('// InPsight: visualization in Psi4\n')
        fh.write('//  by Rob Parrish\n')
        fh.write('// .pov file (adopted from Jmol)\n')
        fh.write('// Created %s\n' % str(date.today()))

        fh.write('#declare Width = %s;\n' % str(self.width))
        fh.write('#declare Height = %s;\n' % str(self.height))
        fh.write('#declare Shadows = %s; \n' % self.defines['Shadows'])
        fh.write('\n')
        fh.write('camera{\n')
        fh.write('  orthographic\n')
        fh.write('  location < %s, %s, %s>\n' %(str(self.location[0]),str(self.location[1]),str(self.location[2]) ))
        fh.write('  sky      < %s, %s, %s>\n' %(str(self.sky[0]),     str(self.sky[1]),     str(self.sky[2]) ))
        fh.write('  up       < %s, %s, %s>\n' %(str(self.up[0]),      str(self.up[1]),      str(self.up[2]) ))
        fh.write('  right    < %s, %s, %s>\n' %(str(self.right[0]),   str(self.right[1]),   str(self.right[2]) ))
        fh.write('  look_at  < %s, %s, %s>\n' %(str(self.look_at[0]), str(self.look_at[1]), str(self.look_at[2]) ))
        fh.write('}\n')
        fh.write('\n')
        fh.write('background { color rgb %s }\n' % self.defines['Background_Color'])
        fh.write('light_source { <%s,%s,%s>  rgb <%s,%s,%s> }\n' \
            %(str(self.light[0]),str(self.light[1]),str(self.light[2]),\
              str(self.light_color[0]),str(self.light_color[1]),str(self.light_color[2])))
        fh.write('\n')
        fh.write('// ***********************************************\n')
        fh.write('// macros for atom/bond shapes\n')
        fh.write('// ***********************************************\n')
        fh.write('#macro check_shadow()\n')
        fh.write(' #if (!Shadows)\n')
        fh.write('  no_shadow \n')
        fh.write(' #end\n')
        fh.write('#end\n')
        fh.write('\n')
        fh.write('#macro translucentFinish(T)\n')
        fh.write(' #local shineFactor = T;\n')
        fh.write(' #if (T <= 0.25)\n')
        fh.write('  #declare shineFactor = (1.0-4*T);\n')
        fh.write(' #end\n')
        fh.write(' #if (T > 0.25)\n')
        fh.write('  #declare shineFactor = 0;\n')
        fh.write(' #end\n')
        fh.write(' finish {\n')
        fh.write('  ambient 0.45\n')
        fh.write('  diffuse 0.84\n')
        fh.write('  specular 0.22\n')
        fh.write('  roughness .00001\n')
        fh.write('  metallic shineFactor\n')
        fh.write('  phong 0.9*shineFactor\n')
        fh.write('  phong_size 120*shineFactor\n')
        fh.write('}#end\n')
        fh.write('\n')
        fh.write('#macro a(X,Y,Z,RADIUS,R,G,B,T)\n')
        fh.write(' sphere{<X,Y,Z>,RADIUS\n')
        fh.write('  pigment{rgbt<R,G,B,T>}\n')
        fh.write('  translucentFinish(T)\n')
        fh.write('  check_shadow()}\n')
        fh.write('#end\n')
        fh.write('\n')
        fh.write('#macro b(X1,Y1,Z1,RADIUS1,X2,Y2,Z2,RADIUS2,R,G,B,T)\n')
        fh.write(' cone{<X1,Y1,Z1>,RADIUS1,<X2,Y2,Z2>,RADIUS2\n')
        fh.write('  pigment{rgbt<R,G,B,T>}\n')
        fh.write('  translucentFinish(T)\n')
        fh.write('  check_shadow()}\n')
        fh.write('#end \n')

        for bond in self.bonds:
            fh.write('b(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)\n' % \
                (str(bond[0]),str(bond[1]),str(bond[2]),str(bond[3]),\
                str(bond[4]),str(bond[5]),str(bond[6]),str(bond[7]),\
                str(bond[8]),str(bond[9]),str(bond[10]),str(bond[11])))
        for atom in self.atoms:
            fh.write('a(%s,%s,%s,%s,%s,%s,%s,%s)\n' % \
                (str(atom[1]),str(atom[2]),str(atom[3]),str(atom[4]),\
                str(atom[5]),str(atom[6]),str(atom[7]),str(atom[8])))

        fh.close()

    def save_density(self,filename='rho',overlap = 2.0,n = [40,40,40],caxis = [0.0,1.0]):
        if (filename != ''):
            self.defines['Filename'] = filename

        grid = GridProp()
        grid.set_n(n[0],n[1],n[2])
        grid.set_caxis(caxis[0],caxis[1])
        grid.set_filename(self.defines['Filename'])
        grid.add('RHO')
        grid.compute()

        df3_file = filename + '.RHO.df3'
        l = [grid.get_l(0),grid.get_l(1),grid.get_l(2)]
        o = [grid.get_o(0),grid.get_o(1),grid.get_o(2)]


        ini_filename = self.defines['Filepath'] + '/' + self.defines['Filename'] + '.pov.ini'
        pov_filename = self.defines['Filepath'] + '/' + self.defines['Filename'] + '.pov'
        png_filename = self.defines['Filepath'] + '/' + self.defines['Filename'] + '.png'
        pov_file = self.defines['Filename'] + '.pov'
        png_file = self.defines['Filename'] + '.png'

        # Write the pov.ini file
        fh = open(ini_filename,'w')
        fh.write('; InPsight: visualization in Psi4\n')
        fh.write(';  by Rob Parrish\n')
        fh.write('; .pov.ini file\n')
        fh.write('; Created %s\n' % str(date.today()))
        fh.write('\n')
        fh.write('Input_File_Name=%s\n' % pov_file)
        fh.write('Output_to_File=true\n')
        fh.write('Output_File_Type=%s\n' % self.defines['Output_File_Type'])
        fh.write('Output_File_Name=%s\n' % png_file)
        fh.write('Height=%s\n' % str(self.height))
        fh.write('Width=%s\n' % str(self.width))
        fh.write('Output_Alpha=%s\n' % self.defines['Output_Alpha'])
        fh.write('Antialias=%s\n' % self.defines['Antialias'])
        fh.write('Antialias_Threshold=%s\n' % self.defines['Antialias_Threshold'])
        fh.write('Display=true\n')
        fh.write('Warning_Level=5\n')
        fh.write('Verbose=false\n')

        fh.close()

        # Write the pov file
        fh = open(pov_filename, 'w')

        fh.write('// InPsight: visualization in Psi4\n')
        fh.write('//  by Rob Parrish\n')
        fh.write('// .pov file (adopted from Jmol)\n')
        fh.write('// Created %s\n' % str(date.today()))

        fh.write('#declare Shadows = %s; \n' % self.defines['Shadows'])
        fh.write('\n')
        fh.write('camera{\n')
        fh.write('  orthographic\n')
        fh.write('  location < %s, %s, %s>\n' %(str(self.location[0]),str(self.location[1]),str(self.location[2]) ))
        fh.write('  sky      < %s, %s, %s>\n' %(str(self.sky[0]),     str(self.sky[1]),     str(self.sky[2]) ))
        fh.write('  up       < %s, %s, %s>\n' %(str(self.up[0]),      str(self.up[1]),      str(self.up[2]) ))
        fh.write('  right    < %s, %s, %s>\n' %(str(self.right[0]),   str(self.right[1]),   str(self.right[2]) ))
        fh.write('  look_at  < %s, %s, %s>\n' %(str(self.look_at[0]), str(self.look_at[1]), str(self.look_at[2]) ))
        fh.write('}\n')
        fh.write('\n')
        fh.write('background { color rgb %s }\n' % self.defines['Background_Color'])
        fh.write('light_source { <%s,%s,%s>  rgb <%s,%s,%s> }\n' \
            %(str(self.light[0]),str(self.light[1]),str(self.light[2]),\
              str(self.light_color[0]),str(self.light_color[1]),str(self.light_color[2])))
        fh.write('\n')
        fh.write('// ***********************************************\n')
        fh.write('// macros for atom/bond shapes\n')
        fh.write('// ***********************************************\n')
        fh.write('#macro check_shadow()\n')
        fh.write(' #if (!Shadows)\n')
        fh.write('  no_shadow \n')
        fh.write(' #end\n')
        fh.write('#end\n')
        fh.write('\n')
        fh.write('#macro translucentFinish(T)\n')
        fh.write(' #local shineFactor = T;\n')
        fh.write(' #if (T <= 0.25)\n')
        fh.write('  #declare shineFactor = (1.0-4*T);\n')
        fh.write(' #end\n')
        fh.write(' #if (T > 0.25)\n')
        fh.write('  #declare shineFactor = 0;\n')
        fh.write(' #end\n')
        fh.write(' finish {\n')
        fh.write('  ambient 0.45\n')
        fh.write('  diffuse 0.84\n')
        fh.write('  specular 0.22\n')
        fh.write('  roughness .00001\n')
        fh.write('  metallic shineFactor\n')
        fh.write('  phong 0.9*shineFactor\n')
        fh.write('  phong_size 120*shineFactor\n')
        fh.write('}#end\n')
        fh.write('\n')
        fh.write('#macro a(X,Y,Z,RADIUS,R,G,B,T)\n')
        fh.write(' sphere{<X,Y,Z>,RADIUS\n')
        fh.write('  pigment{rgbt<R,G,B,T>}\n')
        fh.write('  translucentFinish(T)\n')
        fh.write('  check_shadow()}\n')
        fh.write('#end\n')
        fh.write('\n')
        fh.write('#macro b(X1,Y1,Z1,RADIUS1,X2,Y2,Z2,RADIUS2,R,G,B,T)\n')
        fh.write(' cone{<X1,Y1,Z1>,RADIUS1,<X2,Y2,Z2>,RADIUS2\n')
        fh.write('  pigment{rgbt<R,G,B,T>}\n')
        fh.write('  translucentFinish(T)\n')
        fh.write('  check_shadow()}\n')
        fh.write('#end \n')

        for bond in self.bonds:
            fh.write('b(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)\n' % \
                (str(bond[0]),str(bond[1]),str(bond[2]),str(bond[3]),\
                str(bond[4]),str(bond[5]),str(bond[6]),str(bond[7]),\
                str(bond[8]),str(bond[9]),str(bond[10]),str(bond[11])))
        for atom in self.atoms:
            fh.write('a(%s,%s,%s,%s,%s,%s,%s,%s)\n' % \
                (str(atom[1]),str(atom[2]),str(atom[3]),str(atom[4]),\
                str(atom[5]),str(atom[6]),str(atom[7]),str(atom[8])))

        fh.close()
