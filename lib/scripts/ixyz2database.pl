#!/usr/bin/env perl
use Cwd 'abs_path';

#  Utility: This script converts a set of geometry files in XYZ format into a perl resource file for Sherrill
#     group SETS scripts and a python database file for psi4.
#  Instructions: Detailed instructions may be found on the psi4 trac page. In short, move all XYZ files
#     intended for a database into a directory and run this script from that directory. Answer a few questions
#     about the intended database.
#     PSI4: Edit the resulting database.py file if necessary, then copy it into psi4/lib/databases/ . Its contents
#        can be accessed as normal through the db() wrapper with no further configuration or recompiling.
#     SETS: Edit the resulting resource_database.pl if necessary, then copy it into cdsgroup
#        cdsgroup/scripts/Lori_SETS_scripts/resources/ . Add two lines to resources/resources_sets.pl as 
#        described therein to register the new set. Its contents can be accessed as usual through imake-SETS.pl .
#  Created: Monday, December 21, 2009, LAB
#  Last Modified: Sunday, June 19, 2011, LAB



#   /////                \\\\\
#         CONFIGURE HERE
#   \\\\\                /////

$string2screen = <<END;

# Welcome to ixyz2database.
#    Just fill in the usual variables when prompted. Hit ENTER to accept default.
#    Strings should not be in quotes (all manner of ' " , [] () {} are forbidden). 
#    Elements in arrays should be space-delimited (notation (n .. m) is forbidden).
#    ENTER to continue.
END
print "\n$string2screen";
chomp($trash = <>);

$string2screen = <<END;
# name your database
#    Recommend 3-8 characters, all capitalized letters.
#    e.g., \$set = MINE
#    [SCALAR (required)]
END
print "\n$string2screen";
print "\$set = ";
chomp ($set = <>);

$string2screen = <<END;
# XYZ file extension
#    all files with this extension in the current directory will be processed
#    [SCALAR (default = xyz)]
END
print "\n$string2screen";
print "\$fext = ";
chomp ($fext = <>);
if ($fext eq "") { $fext = "xyz"; }

gobackA:
$string2screen = <<END;
# what should line two of the XYZ file be used for?
#    [comment]   treat content as text for the comment line
#    [cgmp]      treat first item in line as system charge, second as system multiplicity, remainder as comment line
#    [SCALAR (default = comment)]
END
print "\n$string2screen";
print "\$line2 = ";
chomp ($line2 = <>);
if ($line2 eq "") { $line2 = "comment"; }
if ( ($line2 ne "comment") && ($line2 ne "cgmp") ) { print "\nnot an option- try again.\n"; goto gobackA; }

gobackB:
$string2screen = <<END;
# open-shell or non-singlets are present among your systems (or subsystems in the case of dimers) (true/false)?
#    [SCALAR (default = false)]
END
print "\n$string2screen";
print "\$isOS = ";
chomp ($isOS = <>);
if ($isOS eq "") { $isOS = "false"; }
if ( ($isOS ne "true") && ($isOS ne "false") ) { print "\nnot an option- try again.\n"; goto gobackB; }

gobackC:
$string2screen = <<END;
# what is the nature of the systems in your incipient database?
#    [1]   I have a bunch of plain molecules (no need to act on any subsystems) that I want to be able to act upon in parallel.
#    [2]   I have a bunch of molecules that I want to form into a database whose reference quantity corresponds to various combinations thereof.
#    [3]   I have a bunch of dimers that I want to form into a database whose reference quantity is interaction energy.
#    Your final database (for psi4) may of course resemble any combination of these choices. This is but a humble script to get you started.
#    [SCALAR (required)]
END
print "\n$string2screen";
print "\$route = ";
chomp ($route = <>);
if ( ($route != 1) && ($route != 2) && ($route != 3) ) { print "\nnot an option- try again.\n"; goto gobackC; }

if ($route == 2) {
$string2screen = <<END;
# how many reactions (things that have a reference quantity, as opposed to reagents that have a geometry) are in the database?
#    [SCALAR (required)]
END
print "\n$string2screen";
print "\$Nrxn = ";
chomp ($Nrxn = <>);
}

#   /////                \\\\\
#         TERM CONFIGURE
#   \\\\\                /////



open(LAB_OUT, ">stats.txt");
open(SPL_OUT, ">temp1.txt");
open(RPL_OUT, ">temp5.txt");
open(GPL_OUT, ">temp2.txt");
open(SPY_OUT, ">temp3.txt");
open(GPY_OUT, ">temp4.txt");

if ($route == 3) {
   print SPL_OUT "\nsub load_indices_$set {\n\n   \$settype = \"complex-xyz-parted\";\n\n   \@HSYS = (";
}
else {
   print SPL_OUT "\nsub load_indices_$set {\n\n   \$settype = \"reaction\";\n\n   \@HRXN = (1 .. $Nrxn);\n\n";
   print SPL_OUT "   \$numberTT = $Nrxn;\n";
   print SPL_OUT "   \$numberHB = 0;\n";
   print SPL_OUT "   \$numberMX = $Nrxn;\n";
   print SPL_OUT "   \$numberDD = 0;\n";
   print SPL_OUT "\n   foreach \$reaction (\@HRXN) { push(\@setindex,\"\$set-\$reaction\"); }\n}";
   print SPL_OUT "\n\n\n\nsub load_reagents_$set {\n\n   \$settype = \"reaction\";\n\n   \@HSYS = (";
}

print GPL_OUT "\nsub load_xyz_$set {\n\n";

print GPL_OUT "   \$isHB = \"no\";\n";
print GPL_OUT "   \$isMX = \"yes\";\n";
print GPL_OUT "   \$isDD = \"no\";\n";
print GPL_OUT "   \$bas2 = \"no\";\n";
print GPL_OUT "   \$ABun = \"no\";\n";
print GPL_OUT "   \$ABcp = \"no\";\n";
print GPL_OUT "   \$crat = 1.0;\n";
print GPL_OUT "   \$EQyy = \"NaN\";\n";
print GPL_OUT "   \$bind = \"NaN\";\n\n";

@HRXN = (1 .. $Nrxn);
%BINDRXN = ();
%TAGLRXN = ();
foreach $rxn (@HRXN) { $BINDRXN{$rxn} = "nan"; }
foreach $rxn (@HRXN) { $TAGLRXN{$rxn} = "Reaction $rxn"; }

if ($route != 3) {
   print RPL_OUT "\nsub load_rxn_$set {\n\n   \%rxnm = ();\n\n";
   for ($i = 1; $i <= $Nrxn; $i++) {
   
      if ($i == 1) { print RPL_OUT "   if (\$reaction eq \"$i\") {\n\n"; }
      else         { print RPL_OUT "   elsif (\$reaction eq \"$i\") {\n\n"; }
   
      print RPL_OUT "      \$bind = NaN;\n\n";
      print RPL_OUT "      \@HACTIVE = (\"\",\"\",\"\");\n";
      print RPL_OUT "      \@rxnm{\@HACTIVE} = ();\n\n   }\n";
      #print RPL_OUT "      \%rxnm = (\n      \"\"            => ,\n      \"\"            =>      );\n\n   }\n";
   }
   
   print RPL_OUT "   else {\n\n      \$tagl = \"NONSENSE\";\n      \$bind = \"NaN\";\n      \%rxnm = ();\n      \@HACTIVE = ();\n\n   }\n\n";
   print RPL_OUT "foreach my \$component (\@HACTIVE) { \$component = \"\$component-reagent\"; }\n";
   print RPL_OUT "for my \$key (keys %rxnm) { \$rxnm{\"\$key-reagent\"} = delete \$rxnm{\$key}; }\n\n}\n\n\n";

}

# reagent section
print GPY_OUT "\n# <<< Molecule Specifications >>>\n";
if ($route == 3) {
   print GPY_OUT "monoA_unCP = 'monoA = dimer.extract_subsets(1)\\nmonoA.set_name(\"monoA\")\\nPsiMod.set_active_molecule(monoA)\\nPsiMod.IO.set_default_namespace(\"monoA\")\\n'\n";
   print GPY_OUT "monoB_unCP = 'monoB = dimer.extract_subsets(2)\\nmonoB.set_name(\"monoB\")\\nPsiMod.set_active_molecule(monoB)\\nPsiMod.IO.set_default_namespace(\"monoB\")\\n'\n";
   print GPY_OUT "monoA_CP   = 'monoA = dimer.extract_subsets(1,2)\\nmonoA.set_name(\"monoA\")\\nPsiMod.set_active_molecule(monoA)\\nPsiMod.IO.set_default_namespace(\"monoA\")\\n'\n";
   print GPY_OUT "monoB_CP   = 'monoB = dimer.extract_subsets(2,1)\\nmonoB.set_name(\"monoB\")\\nPsiMod.set_active_molecule(monoB)\\nPsiMod.IO.set_default_namespace(\"monoB\")\\n'\n\n";
}
printf "\n%-25s %6s %6s %6s %6s %6s\t\t%-50s %-50s\n", "system", "CHGsyst", "MLPsyst", "Natom", "Nmol1", "Nmol2", "Atoms1", "Atoms2";

$count = 0;
@HRGT = ();
%TAGLRGT = ();
%BINDRGT = ();

foreach $filename (<*.$fext>) {

   $NH = 0;  $NC = 0;  $NN = 0;  $NO = 0;  $NF = 0;  $NSi = 0;  $NP = 0;  $NS = 0;  $NCl = 0;

   # ascertain system name and open file
   ($system) = $filename =~ m/(.*).$fext$/;
   print SPL_OUT "\"$system\",";
   push (@HRGT, $system);

   open(IN, "$system.$fext");
   @text = <IN>;
   close(IN);
   chomp(@text);

   # process first line
   $Nsyst = $text[0];

   # process second line
   if ($line2 eq "cgmp") {

      $tagl = "";
      $text[1] =~ s/(^\s+)|(\s+$)//g;
      @ltemp = split( /\s+/, $text[1], 3);
      $CHGsyst = $ltemp[0];
      $MLPsyst = $ltemp[1];
      $tagl    = $ltemp[2];
   }
   else { 

      $CHGsyst = 0;
      $MLPsyst = 1;
      $tagl = $text[1];
   }
   $TAGLRGT{$system} = $tagl;
   $BINDRGT{$system} = "nan";

   if ($route == 3) {

      # print system introduction
      if ($count == 0) { print GPL_OUT "   if (\$system eq \"$system\") {\n\n"; }
      else             { print GPL_OUT "   elsif (\$system eq \"$system\") {\n\n"; }
   
      print GPL_OUT "      \$tagl = \"$tagl\";\n\n";
      print GPL_OUT "      \%cgmp = (\n";
      print GPL_OUT "      \"CHGsyst\"  => $CHGsyst,\n";
      print GPL_OUT "      \"CHGmol1\"  => $CHGsyst,\n";
      print GPL_OUT "      \"CHGmol2\"  => 0,\n";
      print GPL_OUT "      \"MLPsyst\"  => $MLPsyst,\n";
      print GPL_OUT "      \"MLPmol1\"  => $MLPsyst,\n";
      print GPL_OUT "      \"MLPmol2\"  => 1  );\n\n";
   
      print GPY_OUT "$set" . "_" . $system . " = input.process_input(\"\"\"\n";
      print GPY_OUT "molecule dimer {\n";
      print GPY_OUT "$CHGsyst $MLPsyst\n";
   
      # separate fragments with BFS script by SMS
      @atoms1 = ();
      @atoms2 = ();
      @elems1 = ();
      @elems2 = ();
      @fragment1X = ();
      @fragment2X = ();
      @fragment1Y = ();
      @fragment2Y = ();
      @fragment1Z = ();
      @fragment2Z = ();

      ($callingpath) = abs_path($0) =~ m/(.*)\/ixyz2database\.pl/;
      system("$callingpath/BFS4xyz2db.py $system.$fext > BFStemp.txt");
   
      open(IN, "BFStemp.txt");
      @textF = <IN>;
      close(IN);
   
      chomp($textF[0]);
      if ($textF[0] ne "Found 2 fragments") { print "ERROR: 2 fragments not detected for system $system\n"; }
   
      @ltemp = split( /\s+/, $textF[1]);
      $Nmol1 = $ltemp[3];
      if ( ($ltemp[0] ne "Fragment") || ($Nmol1 le 0) ) { print "ERROR: fragment 1 has $Nmol1 atoms\n"; }
   
      @ltemp = split( /\s+/, $textF[2]);
      foreach $atom (@ltemp) { push(@atoms1, $atom); }
   
      @ltemp = split( /\s+/, $textF[3]);
      $Nmol2 = $ltemp[3];
      if ( ($ltemp[0] ne "Fragment") || ($Nmol2 le 0) ) { print "ERROR: fragment 2 has $Nmol2 atoms\n"; }
   
      @ltemp = split( /\s+/, $textF[4]);
      foreach $atom (@ltemp) { push(@atoms2, $atom); }
   
      if ( ($Nmol1 + $Nmol2) != $Nsyst) { print "ERROR: Nsyst atoms ($Nsyst) unequal to Nmol1 ($Nmol1) + Nmol2 ($Nmol2)\n"; }
      printf "%-25s %6d %6d %6d %6d %6d", $system,$CHGsyst,$MLPsyst,$Nsyst,$Nmol1,$Nmol2;
      print "\t\t@atoms1\t\t@atoms2\n";
   
      print GPL_OUT "      \%stat = (\n";
      print GPL_OUT "      \"Nsyst\"  => $Nsyst,\n";
      print GPL_OUT "      \"Nmol1\"  => $Nmol1,\n";
      print GPL_OUT "      \"Nmol2\"  => $Nmol2  );\n\n";
   
      # separate fragments
      print GPL_OUT "      \@geom = ();\n";
   
      foreach $atom (@atoms1) {
   
         $text[$atom+1] =~ s/(^\s+)|(\s+$)//g;
         @ltemp = split( /\s+/, $text[$atom+1] );
         push(@elems1,$ltemp[0]);
         push(@fragment1X,$ltemp[1]);
         push(@fragment1Y,$ltemp[2]);
         push(@fragment1Z,$ltemp[3]);
      }
   
      foreach $atom (@atoms2) {
      
         $text[$atom+1] =~ s/(^\s+)|(\s+$)//g;
         @ltemp = split( /\s+/, $text[$atom+1] );
         push(@elems2,$ltemp[0]);
         push(@fragment2X,$ltemp[1]);
         push(@fragment2Y,$ltemp[2]);
         push(@fragment2Z,$ltemp[3]);
      }
   
      # print geometry
      for($i = 0; $i < $Nmol1; $i++) { printf GPL_OUT "      push(\@geom, \"%-3s  % 14.8f % 14.8f % 14.8f\");\n", $elems1[$i], $fragment1X[$i], $fragment1Y[$i], $fragment1Z[$i]; }
      for($i = 0; $i < $Nmol2; $i++) { printf GPL_OUT "      push(\@geom, \"%-3s  % 14.8f % 14.8f % 14.8f\");\n", $elems2[$i], $fragment2X[$i], $fragment2Y[$i], $fragment2Z[$i]; }
      print GPL_OUT "\n   }\n";
   
      for($i = 0; $i < $Nmol1; $i++) { printf GPY_OUT "%-3s  % 14.8f % 14.8f % 14.8f\n", $elems1[$i], $fragment1X[$i], $fragment1Y[$i], $fragment1Z[$i]; }
      print GPY_OUT "--\n0 1\n";
      for($i = 0; $i < $Nmol2; $i++) { printf GPY_OUT "%-3s  % 14.8f % 14.8f % 14.8f\n", $elems2[$i], $fragment2X[$i], $fragment2Y[$i], $fragment2Z[$i]; }
      print GPY_OUT "units angstrom\n}\n\"\"\", 0)\n\n";

   }
   else {

      # print system introduction
      if ($count == 0) { print GPL_OUT "   if (\$system eq \"$system\") {\n\n"; }
      else             { print GPL_OUT "   elsif (\$system eq \"$system\") {\n\n"; }
   
      print GPL_OUT "      \$tagl = \"$system $tagl\";\n\n";
      print GPL_OUT "      \%cgmp = (\n";
      print GPL_OUT "      \"CHGsyst\"  => $CHGsyst,\n";
      print GPL_OUT "      \"MLPsyst\"  => $MLPsyst  );\n\n";

      print GPY_OUT "$set" . "_" . $system . " = input.process_input(\"\"\"\n";
      print GPY_OUT "molecule dimer {\n";
      print GPY_OUT "$CHGsyst $MLPsyst\n";
   
      # separate fragments with BFS script by SMS
      $Nmol1 = $Nsyst;
      $Nmol2 = 0;
      @atoms1 = (1 .. $Nsyst);
      @elems1 = ();
      @fragment1X = ();
      @fragment1Y = ();
      @fragment1Z = ();
   
      if ( ($Nmol1 + $Nmol2) != $Nsyst) { print "ERROR: Nsyst atoms ($Nsyst) unequal to Nmol1 ($Nmol1) + Nmol2 ($Nmol2)\n"; }
      printf "%-25s %6d %6d %6d %6d %6d", $system,$CHGsyst,$MLPsyst,$Nsyst,$Nmol1,$Nmol2;
      print "\t\t@atoms1\t\t@atoms2\n";
   
      print GPL_OUT "      \%stat = (\n";
      print GPL_OUT "      \"Nsyst\"  => $Nsyst  );\n\n";
   
      # separate fragments
      print GPL_OUT "      \@geom = ();\n";
   
      foreach $atom (@atoms1) {
   
         $text[$atom+1] =~ s/(^\s+)|(\s+$)//g;
         @ltemp = split( /\s+/, $text[$atom+1] );
         push(@elems1,$ltemp[0]);
         push(@fragment1X,$ltemp[1]);
         push(@fragment1Y,$ltemp[2]);
         push(@fragment1Z,$ltemp[3]);

         if ($ltemp[0] eq "H")  { $NH++;  }
         if ($ltemp[0] eq "C")  { $NC++;  }
         if ($ltemp[0] eq "N")  { $NN++;  }
         if ($ltemp[0] eq "O")  { $NO++;  }
         if ($ltemp[0] eq "F")  { $NF++;  }
         if ($ltemp[0] eq "Si") { $NSi++; }
         if ($ltemp[0] eq "P")  { $NP++;  }
         if ($ltemp[0] eq "S")  { $NS++;  }
         if ($ltemp[0] eq "Cl") { $NCl++; }
      }
   
      # print geometry
      for($i = 0; $i < $Nmol1; $i++) { printf GPL_OUT "      push(\@geom, \"%-3s  % 14.8f % 14.8f % 14.8f\");\n", $elems1[$i], $fragment1X[$i], $fragment1Y[$i], $fragment1Z[$i]; }
      print GPL_OUT "\n   }\n";
   
      for($i = 0; $i < $Nmol1; $i++) { printf GPY_OUT "%-3s  % 14.8f % 14.8f % 14.8f\n", $elems1[$i], $fragment1X[$i], $fragment1Y[$i], $fragment1Z[$i]; }
      print GPY_OUT "units angstrom\n}\n\"\"\", 0)\n\n";

      printf LAB_OUT "%-25s\t%d\t%d\t%d\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", $system, $Nsyst, $CHGsyst, $MLPsyst, $NH, $NC, $NN, $NO, $NF, $NSi, $NP, $NS, $NCl, $tagl;

   }

   $count++;
}
$Nrgt = $#HRGT + 1;
if ($Nrgt != $count) { print "ERROR: discrepancy in counting systems $Nrgt vs $count!\n"; }


# perl resource file
print SPL_OUT ");\n\n";
print SPL_OUT "   \$numberTT = $Nrgt;\n";
print SPL_OUT "   \$numberHB = 0;\n   \$numberMX = 0;\n   \$numberDD = 0;\n\n";
print SPL_OUT "   foreach \$system (\@HSYS) { push(\@setindex,\"\$set-\$system\"); }\n}\n\n\n";

print GPL_OUT "   else {\n\n";
print GPL_OUT "      \$tagl = \"NONSENSE\";\n";
print GPL_OUT "      \%cgmp = ();\n";
print GPL_OUT "      \%stat = ();\n";
print GPL_OUT "      \@geom = ();\n";
print GPL_OUT "      \$crat = \"NaN\";\n\n   }\n\n";
      
if ($route != 3) {
   print GPL_OUT "   \$stat{\"Nmol1\"} = \$stat{\"Nsyst\"};\n";
   print GPL_OUT "   \$stat{\"Nmol2\"} = 0;\n";
}
print GPL_OUT "   \$stat{\"Msyst\"} = \$stat{\"Nsyst\"};\n";
print GPL_OUT "   \$stat{\"Mmol1\"} = \$stat{\"Nmol1\"};\n";
print GPL_OUT "   \$stat{\"Mmol2\"} = \$stat{\"Nmol2\"};\n";
print GPL_OUT "   \$stat{\"param\"} = 0;\n\n";
print GPL_OUT "   \$texl = \$tagl;\n\n";

if ($route == 3) {
   print GPL_OUT "   if (\$ABun eq \"yes\") { \@HACTIVEun = (\"\$system-dimer\",\"\$system-monoA-unCP\",\"\$system-monoA-unCP\"); }\n";
   print GPL_OUT "   else                { \@HACTIVEun = (\"\$system-dimer\",\"\$system-monoA-unCP\",\"\$system-monoB-unCP\"); }\n\n";

   print GPL_OUT "   \$rxnmUN{\"\$system-dimer\"}      = +1;\n";
   print GPL_OUT "   \$rxnmUN{\"\$system-monoA-unCP\"} = -1;\n";
   print GPL_OUT "   \$rxnmUN{\"\$system-monoB-unCP\"} = -1;\n\n";
   
   print GPL_OUT "   if (\$ABcp eq \"yes\") { \@HACTIVEcp = (\"\$system-dimer\",\"\$system-monoA-CP\",\"\$system-monoA-CP\"); }\n";
   print GPL_OUT "   else                { \@HACTIVEcp = (\"\$system-dimer\",\"\$system-monoA-CP\",\"\$system-monoB-CP\"); }\n\n";

   print GPL_OUT "   \$rxnmCP{\"\$system-dimer\"}      = +1;\n";
   print GPL_OUT "   \$rxnmCP{\"\$system-monoA-CP\"}   = -1;\n";
   print GPL_OUT "   \$rxnmCP{\"\$system-monoB-CP\"}   = -1;\n\n";
}

print GPL_OUT "}\n\nreturn 1;\n";


# python database file
print SPY_OUT "\"\"\"\n";
print SPY_OUT "**$set**\n\n";
print SPY_OUT "| Database of <description of members and reference energy type>.\n";
print SPY_OUT "| Geometries from <Reference>.\n";
print SPY_OUT "| Reference interaction energies from <Reference>.\n\n";
if ($route == 3) {
    print SPY_OUT "- **cp**  ``'off'`` <erase this comment and after unless on is a valid option> || ``'on'``\n\n";
    print SPY_OUT "- **rlxd** ``'off'`` <erase this comment and after unless on is valid option> || ``'on'``\n\n";
}
print SPY_OUT "- **benchmark**\n\n";
print SPY_OUT "  - ``'<benchmark_name>'`` <Reference>.\n";
print SPY_OUT "  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.\n\n";
print SPY_OUT "- **subset**\n\n";
print SPY_OUT "  - ``'small'``\n";
print SPY_OUT "  - ``'large'``\n";
print SPY_OUT "  - ``'<subset>'`` <members_description>\n\n";
print SPY_OUT "----\n\n";
print SPY_OUT "\"\"\"\n";
print SPY_OUT "import re\n";
print SPY_OUT "import input\n";

print SPY_OUT "\n# <<< $set Database Module >>>\n";
print SPY_OUT "dbse = '$set'\n";
if ($isOS eq "true") { print SPY_OUT "isOS = '$isOS'\n"; }

print SPY_OUT "\n# <<< Database Members >>>\n";
print SPY_OUT "HRXN = [";
if    ($route == 1) { foreach $rgt (@HRGT) { print SPY_OUT "'$rgt', "; } }
elsif ($route == 2) { foreach $rxn (@HRXN) { print SPY_OUT "'$rxn', "; } }
elsif ($route == 3) { foreach $rgt (@HRGT) { print SPY_OUT "'$rgt', "; } }
print SPY_OUT "]\n";
print SPY_OUT "HRXN_SM = []\n";
print SPY_OUT "HRXN_LG = []\n";

print SPY_OUT "\n# <<< Chemical Systems Involved >>>\n";
print SPY_OUT "RXNM = {}     # reaction matrix of reagent contributions per reaction\n";
print SPY_OUT "ACTV = {}     # order of active reagents per reaction\n";
if    ($route == 1) {
   foreach $rgt (@HRGT) {

      print SPY_OUT "ACTV['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = ", $rgt . "\'"; print SPY_OUT "['%s-%s-reagent'      % (dbse, '$rgt')]\n";
      print SPY_OUT "RXNM['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = ", $rgt . "\'"; print SPY_OUT "dict(zip(ACTV['%s-%s' % (dbse, '$rgt')], [+1]))\n\n";
   }
}
elsif ($route == 2) {
   foreach $rxn (@HRXN) {

      print SPY_OUT "ACTV['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = ", $rxn . "\'"; print SPY_OUT "['%s-%s-reagent'      % (dbse, ''),\n";
      print SPY_OUT "                                                               '%s-%s-reagent'      % (dbse, ''),\n";
      print SPY_OUT "                                                               '%s-%s-reagent'      % (dbse, '') ]\n";
      print SPY_OUT "RXNM['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = ", $rxn . "\'"; print SPY_OUT "dict(zip(ACTV['%s-%s' % (dbse, '$rxn')], []))\n\n";
   }
}
elsif ($route == 3) { 

   print SPY_OUT "ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction\n";
   print SPY_OUT "ACTV_SA = {}  # order of active reagents for non-supermolecular calculations\n";

   print SPY_OUT "for rxn in HRXN:\n\n";

   print SPY_OUT "    RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,\n";
   print SPY_OUT "                                      '%s-%s-monoA-CP'   % (dbse, rxn) : -1,\n";
   print SPY_OUT "                                      '%s-%s-monoB-CP'   % (dbse, rxn) : -1,\n";
   print SPY_OUT "                                      '%s-%s-monoA-unCP' % (dbse, rxn) : -1,\n";
   print SPY_OUT "                                      '%s-%s-monoB-unCP' % (dbse, rxn) : -1 }\n\n";

   print SPY_OUT "    ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]\n\n";

   print SPY_OUT "    ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),\n";
   print SPY_OUT "                                      '%s-%s-monoA-CP'   % (dbse, rxn),\n";
   print SPY_OUT "                                      '%s-%s-monoB-CP'   % (dbse, rxn) ]\n\n";

   print SPY_OUT "    ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),\n";
   print SPY_OUT "                                      '%s-%s-monoA-unCP' % (dbse, rxn),\n";
   print SPY_OUT "                                      '%s-%s-monoB-unCP' % (dbse, rxn) ]\n\n";
}

print SPY_OUT "# <<< Reference Values [kcal/mol] >>>\n";
print SPY_OUT "BIND = {}\n";
print SPY_OUT "nan = float('NaN')\n";
if    ($route == 1) { foreach $rgt (@HRGT) { print SPY_OUT "BIND['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = %8.3f\n", $rgt . "\'", $BINDRGT{$rgt}; } }
elsif ($route == 2) { foreach $rxn (@HRXN) { print SPY_OUT "BIND['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = %8.3f\n", $rxn . "\'", $BINDRXN{$rxn}; } }
elsif ($route == 3) { foreach $rgt (@HRGT) { print SPY_OUT "BIND['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = %8.3f\n", $rgt . "\'", $BINDRGT{$rgt}; } }

print SPY_OUT "\n# <<< Comment Lines >>>\n";
print SPY_OUT "TAGL = {}\n";
if    ($route == 1) { 

   foreach $rgt (@HRGT) { print SPY_OUT "TAGL['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", $TAGLRGT{$rgt};
                          print SPY_OUT "TAGL['%s-%s-reagent'    % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", $TAGLRGT{$rgt};                     }
}
elsif ($route == 2) {

   foreach $rxn (@HRXN) { print SPY_OUT "TAGL['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rxn . "\'", $TAGLRXN{$rxn};                     }
   foreach $rgt (@HRGT) { print SPY_OUT "TAGL['%s-%s-reagent'    % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", $TAGLRGT{$rgt};                     }
}
elsif ($route == 3) {

   foreach $rgt (@HRGT) { print SPY_OUT "TAGL['%s-%s'            % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", $TAGLRGT{$rgt};
                          print SPY_OUT "TAGL['%s-%s-dimer'      % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", "Dimer from " . $TAGLRGT{$rgt};
                          print SPY_OUT "TAGL['%s-%s-monoA-CP'   % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", "Monomer A from " . $TAGLRGT{$rgt};
                          print SPY_OUT "TAGL['%s-%s-monoB-CP'   % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", "Monomer B from " . $TAGLRGT{$rgt};
                          print SPY_OUT "TAGL['%s-%s-monoA-unCP' % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", "Monomer A from " . $TAGLRGT{$rgt};
                          print SPY_OUT "TAGL['%s-%s-monoB-unCP' % (dbse, '";  printf SPY_OUT "%-22s )] = \"\"\"%s \"\"\"\n", $rgt . "\'", "Monomer B from " . $TAGLRGT{$rgt}; }
}

print GPY_OUT "# <<< Geometry Specification Strings >>>\n";
print GPY_OUT "rxnpattern = re.compile(r'^(.+)-(.+)-(.+)\$')\n";
print GPY_OUT "GEOS = {}\n";
print GPY_OUT "for rxn in HRXN:\n";
if ( ($route == 1) || ($route == 2) ) {

   print GPY_OUT   "    for rgt in ACTV['%s-%s' % (dbse, rxn)]:\n\n";

   print GPY_OUT   "        molname = rxnpattern.match(rgt)\n";
   print GPY_OUT   "        GEOS['%s' % (rgt)] = eval('%s_%s' % (dbse, molname.group(2)))\n";
}
elsif ($route == 3) {

   print GPY_OUT "\n    GEOS['%s-%s-dimer'      % (dbse, rxn)] = eval('%s_%s' % (dbse, rxn))\n";
   print GPY_OUT   "    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn))) + monoA_CP\n";
   print GPY_OUT   "    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn))) + monoB_CP\n";
   print GPY_OUT   "    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn))) + monoA_unCP\n";
   print GPY_OUT   "    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = str(eval('%s_%s' % (dbse, rxn))) + monoB_unCP\n";
}


close(LAB_OUT);
close(SPL_OUT);
close(RPL_OUT);
close(GPL_OUT);
close(SPY_OUT);
close(GPY_OUT);

open(PL_OUT, ">resource_$set.pl");
open(PY_OUT, ">$set.py");

system("cat temp1.txt >  resource_$set.pl");
system("cat temp5.txt >> resource_$set.pl");
system("cat temp2.txt >> resource_$set.pl");
system("cat temp3.txt >  $set.py");
system("cat temp4.txt >> $set.py");

close(PL_OUT);
close(PY_OUT);

system("rm -f temp1.txt temp2.txt temp3.txt temp4.txt temp5.txt BFStemp.txt stats.txt");
@groups = split '\s', $(;
$dropresource = 1;
foreach $grp (@groups) { if (getgrgid($grp) eq "sherrill") { $dropresource = 0; } }
if ($dropresource) { system("rm -f resource_$set.pl"); }


print "\n   **  Congratulations, your database file $set.py has been constructed!\n";

print "\n   **  To have a minimally functioning database, do the following:\n\n";

if (($line2 eq "comment") && ($isOS eq "true")) {
    print "       *  If not all neutral singlets, fill in correct charge and multiplicity for all reagents.\n"; }
if (($line2 eq "comment") && ($isOS eq "false")) {
    print "       *  If not all neutral, fill in correct charge for all reagents.\n"; }
if (($route == 3) && ($line2 eq "cgmp")) {
    print "       *  The charge and multiplicity read in from line2 of the xyz files has been assigned to\n";
    print "          fragmentA, leaving fragmentB as a neutral singlet. If this is incorrect for any\n";
    print "          reagents, reapportion the charge and multiplicity correctly between fragments A & B.\n";
}
if (($route == 3) && ($line2 eq "comment")) {
    print "       *  If dimer and both subsystems are not neutral singlets, fill in correct charge and\n";
    print "          multiplicity for each subsystem.\n";
}
if ($route == 2) {
    print "       *  Define the reagents that contribute to reach reaction by filling in the empty single\n";
    print "          quotes in ACTV. Add or delete lines as necessary for each reaction if more or fewer than\n";
    print "          three reagents contribute. See NHTBH.py as an example.\n";
    print "       *  Define the mathematical contribution of reagents to reactions by filling in a number (most\n";
    print "          often +1 or -1) for each reagent to the RXNM of each reaction. See NHTBH.py as an example.\n";
}
print "       *  Move $set.py into \$PSIDATADIR/databases/ so that the PSI4 driver can find it.\n";

print "\n   **  To enhance the functionality/documentation of your database, do the following:\n\n";

print "       *  Rearrange the order of reactions in HRXN, as this will define the order for the database.\n";
print "       *  Fill in the skeleton docstring at top of file, adding sources for geometries and\n";
print "          any reference data.\n";
print "       *  Add a line for your database at the bottom of psi4/doc/userman/sphinx_psithon_doc/source/db.rst\n";
print "          so that it gets slurped into the documentation.\n";
print "       *  Fill in the comment lines of TAGL in plain text. These show up as banners in job output files.\n";
print "       *  Fill in reference values (in kcal/mol) into BIND.\n";
print "       *  If multiple sets of reference values are available, define each in an array BIND_ALTREF so that\n";
print "          they can be called in a psi4 input file as benchmark='ALTREF'. Add the new reference to\n";
print "          the docstring. See S22.py as an example.\n";
print "       *  Fill in the least computationally expensive 2-3 reactions into HRXN_SM and the most expensive\n";
print "          into HRXN_LG so that they can be called in a psi4 input file as subset='small' or subset='large'\n";
print "       *  Define subsets of reactions such as in an array SUBSETARRAY=['reaction', 'reaction'] so that\n";
print "          they can be called in a psi4 input file as subset='SUBSETARRAY'. Add the new subset option to\n";
print "          to the docstring. See NBC10.py for a simple example or CFLOW.py for a complex example.\n";

print "\n";

