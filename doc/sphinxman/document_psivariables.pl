#!/usr/bin/env perl

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

use strict;
use warnings;
use File::Path qw(remove_tree);

# Then, for each module detected in read_options.cc, we collect the process environment
# variables set in src/bin/module and src/lib/libmodule_solver and place them in a TeX file,
# which is inlined into the manual.

#
# First, read the options for each module
#
my $DriverPath = "";
if ($#ARGV == 0) { $DriverPath = $ARGV[0] . "/"; }
my $DriverFile = $DriverPath . "../../psi4/src/read_options.cc";


my $CurrentModule;
my %Keywords;
my %ModuleDescriptions;
open(DRIVER, $DriverFile) or die "\nI can't read the PSI driver file\n";

# Start looping over the driver file
while(<DRIVER>){
    # Look for the name of the module for which options are being added
    if(/name\s*\=\=\s*\"(\w+)\"/){
        $CurrentModule = $1;
        $Keywords{$CurrentModule} = 1;
    }
    if(/\/\*-\s*MODULEDESCRIPTION/ and $CurrentModule){
        $ModuleDescriptions{$CurrentModule} = get_description($_);
    }
}
close DRIVER;


my @temp = ();
print_hash(\%Keywords);
my @PSIMODULES = @temp;
push(@PSIMODULES, "OEPROP");

sub print_hash
{
 my %hash     = %{$_[0]};
 my @RearrModules = sort {$a cmp $b} keys %hash;
 @RearrModules = grep { $_ ne "GLOBALS"} @RearrModules;
 unshift(@RearrModules, "GLOBALS");
 foreach my $Module (@RearrModules){
     $Module =~ s/_/-/g; # Things like plugin_module_name will screw things up...
     push(@temp, $Module);
 }  # module
}


sub get_description
{
 my $Line = shift;
 chomp $Line;
 my $String;
 #Process the inputted line
 if($Line =~ /(?:\/\*\-\s*MODULEDESCRIPTION)\s+(.*)/){
     $String = $1;
 }else{
     die "\nI don't know what to do with $Line\n";
 }
 if(!($String =~ s/\s*-\*\///g)){
     # 'twas more than a one-liner, let's keep searching
     while(<DRIVER>){
         # Add on the current line, sans the newline
         chomp;
         $String .= $_;
         # Attempt to nuke any -*/ patterns, if successful we're done
         last if $String =~ s/\s*-\*\///g;
     }
 }
 # Search and replace multiple spaces with a single space, not that LaTeX cares
 $String =~ s/ +/ /g;
 return $String;
}



#
# Scan the source for Process::Environment variables
#

my $SrcFolder = $DriverPath . "../../psi4/src/psi4/";
my $TexSummary = "variables_list.tex";
my $RstSummary = "source/autodoc_psivariables_bymodule.rst";
open(TEXOUT,">$TexSummary") or die "I can't write to $TexSummary\n";
print TEXOUT "{\n \\footnotesize\n";
open(VOUT,">$RstSummary") or die "I can't write to $RstSummary\n";
print VOUT "\n.. _`apdx:psivariables_module`:\n\n";
print VOUT "PSI Variables by Module\n=======================\n\n";
print VOUT ".. note:: Lowercase letters in PSI variable names represent variable portions of the name.\n";
print VOUT "   See :ref:`apdx:psivariables_alpha` for fuller description.\n\n";
print VOUT ".. toctree::\n   :maxdepth: 1\n\n";
# Grab psi modules and ordering from options parsing above
foreach my $Module (@PSIMODULES) {
    # Set path for each module of bin/module and lib/libmodule_solver
    #     Assign stray variables as for OEPROP below
    my @RelevantDirs = ($SrcFolder . lc($Module), $SrcFolder . "lib" . lc($Module) . "_solver");
    if ($Module eq "OEPROP") { push(@RelevantDirs, $SrcFolder . "libmints"); }
    if ($Module eq "GDMA") { push(@RelevantDirs, $SrcFolder . "libgdma"); }
    my @EnvVariables = ();
    my @EnvArrays = ();
    printf TEXOUT "\n\\subsection{%s}\n",$Module;
    # Search each line in each file in each module-relevant directory for environment variables
    foreach my $Dir (@RelevantDirs) {
        if (opendir(SRC, $Dir)) {
            while (my $file = readdir(SRC)) {
                if (open(CODE, "<$Dir/$file")) {
                    my @text = <CODE>;
                    foreach my $line (@text) {
                        if ($line =~ /\QProcess::environment.globals\E/) {
                            my @ltemp = split( /"/, $line);
                            if ($ltemp[0] =~ /\QProcess::environment.globals\E/) {
                                push(@EnvVariables, $ltemp[1]);
                            }
                        }
                        if ($line =~ /\QProcess::environment.arrays\E/) {
                            my @ltemp = split( /"/, $line);
                            if ($ltemp[0] =~ /\QProcess::environment.arrays\E/) {
                                push(@EnvArrays, $ltemp[1]);
                            }
                        }
                    }
                    close(CODE);
                }
            }
            closedir SRC;
        }
    }
    # Remove duplicate env variables, sort into alphabetical order, and print to tex file
    my %EnvHash;
    foreach my $EnvVar (@EnvVariables){
         $EnvHash{$EnvVar} = 1 if $EnvVar;
    }
    foreach my $EnvVar (@EnvArrays){
         $EnvHash{$EnvVar} = 2 if $EnvVar;
    }
    if (scalar keys %EnvHash > 0) {
       print VOUT "   autodir_psivariables/module__" . lc($Module) . "\n";
       my $vvout = "source/autodir_psivariables/module__" . lc($Module) . ".rst";
       open(VVOUT,">$vvout") or die "I can't write to $vvout\n";
       printf VVOUT ".. _`apdx:%s_psivar`:\n\n", lc($Module);
       my $Moddivider = "=" x length($Module);
       printf VVOUT "\n%s\n%s\n\n", uc($Module), $Moddivider;
       if (exists $ModuleDescriptions{$Module}) { printf VVOUT "$ModuleDescriptions{$Module}\n\n"; }
       print VVOUT ".. hlist::\n   :columns: 1\n\n";
       foreach my $Var (sort keys %EnvHash) {
           printf TEXOUT '\\begin{tabular*}{\\textwidth}[tb]{p{1.0\\textwidth}}';
           printf TEXOUT "\n\t %s \\\\ \n", $Var;
           print TEXOUT "\\end{tabular*}\n";
           my $squashedVar = $Var;
           $squashedVar =~ s/ //g;
           if ($EnvHash{$Var} == 2) {
               printf VVOUT "   * :psivar:`%s <%s>` (array)\n\n", $Var, $squashedVar;
           } else {
               printf VVOUT "   * :psivar:`%s <%s>`\n\n", $Var, $squashedVar;
           }
       }
       print VVOUT "\n";
       close VVOUT;
       printf "Auto-documenting psi variables in module %s\n", lc($Module);
    }
}
print TEXOUT "}\n";
unlink("variables_list.tex");
close TEXOUT;
print VOUT "   autodir_psivariables/module__cfour\n";
print VOUT "\n";
close VOUT;

