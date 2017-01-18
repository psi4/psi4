#!/usr/bin/perl

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

system ("ls -1 *cc-*.gbs | grep -v 'autogen' | grep -v 'tight' | grep -v 'polarization' | grep -v 'molpro' | grep -v 'diffuse' | grep -v 'basis' | grep -v 'corevalence' | grep -v 'hold' > realdunnings.txt");

system ("echo 'changed dunning basis sets\n\n' > basisdunningdiffer.txt");
open IN, "realdunnings.txt";
@text = <IN>;
close(IN);

foreach $file (@text) {

   chomp($file);
   print "$file\n";

   open BASFILE, $file;
   @bastext = <BASFILE>;
   close(BASFILE);

   if (!($bastext[0] =~ /spherical/)) {
      system ("./prepend_puream_onefile.sh $file");
   }

   #system("diff -bw -I '^! merged' $file ../../$file");  # less restrictive- shows changes in comment lines that arise when basis sets change but affected atoms not present
   system("diff -bw -I '^!' $file ../../$file");
   #system("diff -bw -I '^!' $file ../$file");
   system("diff -qbw -I '^!' $file ../../$file >> basisdunningdiffer.txt");
   #system("diff -qbw -I '^!' $file ../$file >> basisdunningdiffer.txt");
   print "\n\n\n";
}


