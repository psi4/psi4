#!/usr/bin/perl

system ("rm -f merged.gbs");

$numArgs = $#ARGV + 1;
if ($numArgs == 0) { print_usage (); exit; }
if ($numArgs > 4)  { print_usage (); exit; }

$file1 = $ARGV[0];
$file2 = $ARGV[1];
$mode  = $ARGV[2];
$quiet = $ARGV[3];

if (!(-e $file1)) { die "\nERROR: file1 $file1 not present\n\n"; }
if (!(-e $file2)) { die "\nERROR: file2 $file2 not present\n\n"; }

if ( ($mode eq "") || ($mode eq "exclusive") ) { $mode = "exclusive"; }
elsif ($mode eq "inclusive")                   { $mode = "inclusive"; }
elsif ($mode eq "insert_alar")                 { $mode = "alar";      }
elsif ($mode eq "insert_hhe")                  { $mode = "hhe";       }
else                                           { die "\nERROR: improper mode $mode\n\n"; }

if ($quiet eq "quiet") { $quiet = "quiet"; }
else                   { $quiet = "no"; }

if ($quiet eq "no") { print "\nMODE: $mode\n"; }

open(GBS_OUT,">merge-$file1-$file2.gbs");

open IN, $file1;
@text1 = <IN>;
close(IN);

open IN, $file2;
@text2 = <IN>;
close(IN);

# merge comment material
foreach $line (@text1) { if ($line =~ /^!/) { print GBS_OUT $line; } }
print GBS_OUT "!\n";
foreach $line (@text2) { if ($line =~ /^!/) { print GBS_OUT $line; } }
print GBS_OUT "!\n";
my $date = `/bin/date +"%a, %H:%M %b %d, %Y"`;
print GBS_OUT "! merged from $file1 and $file2 at $date";
print GBS_OUT "!\n\n\n\n****\n";

@HELEM = (H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr);
@HLITE = (H,He);
#@HHEAV =      (Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr);
@HALAR =                               (Al,Si,P,S,Cl,Ar);


# merge basis set info
foreach $element (@HELEM) {

   $startline1 = 0;
   $startline2 = 0;
   
   $ii = 0;
   foreach $line (@text1) {

      if ( ($line =~ /^($element)\s+0\s/i) && ($alert eq "yes") ) { $startline1 = $ii; }

      if ($line =~ /\*\*\*\*/) { $alert = "yes"; }
      else                     { $alert = "no";  }

      $ii++;
   }

   $ii = 0;
   foreach $line (@text2) {

      if ( ($line =~ /^($element)\s+0\s/i) && ($alert eq "yes") ) { $startline2 = $ii; }

      if ($line =~ /\*\*\*\*/) { $alert = "yes"; }
      else                     { $alert = "no";  }

      $ii++;
   }

   if ($quiet eq "no") { print "\n$element\t$startline1\t$startline2\t"; }

   if ($mode eq "exclusive") {

      if ( ($startline1 != 0) && ($startline2 != 0) ) {

         if ($quiet eq "no") { print "****"; }
         print GBS_OUT "$element     0\n";

         print_element_file1 ();
         print_element_file2 ();

         print GBS_OUT "****\n";

      }
   }

   elsif ($mode eq "inclusive") {

      if ( ($startline1 != 0) && ($startline2 == 0) ) {

         if ($quiet eq "no") { print "****"; }
         print GBS_OUT "$element     0\n";

         print_element_file1 ();

         print GBS_OUT "****\n";

      }
      elsif ( ($startline1 == 0) && ($startline2 != 0) ) {

         if ($quiet eq "no") { print "****"; }
         print GBS_OUT "$element     0\n";

         print_element_file2 ();

         print GBS_OUT "****\n";

      }
      elsif ( ($startline1 != 0) && ($startline2 != 0) ) {

         print "ERROR: ELEMENT $element in BOTH FILES";
      }
   }

   elsif ($mode eq "hhe") {

      if ( (grep {$_ eq $element} @HLITE) && ($startline2 != 0) ) {

         if ($quiet eq "no") { print "** 2 **"; }
         print GBS_OUT "$element     0\n";

         print_element_file2 ();

         print GBS_OUT "****\n";

      }
      elsif ($startline1 != 0) {

         if ($quiet eq "no") { print "** 1 **"; }
         print GBS_OUT "$element     0\n";

         print_element_file1 ();

         print GBS_OUT "****\n";

      }
   }

   elsif ($mode eq "alar") {

      if ( (grep {$_ eq $element} @HALAR) && ($startline2 != 0) ) {

         if ($quiet eq "no") { print "** 2 **"; }
         print GBS_OUT "$element     0\n";

         print_element_file2 ();

         print GBS_OUT "****\n";

      }
      elsif ($startline1 != 0) {

         if ($quiet eq "no") { print "** 1 **"; }
         print GBS_OUT "$element     0\n";

         print_element_file1 ();

         print GBS_OUT "****\n";

      }
   }

}

if ($quiet eq "no") { print "\n\n"; }
print GBS_OUT "\n";
close(GBS_OUT);

if ($quiet eq "quiet") { system ("mv merge-$file1-$file2.gbs merged.gbs"); }


sub print_element_file1 {

   $ii = 1;
   while ($text1[$startline1+$ii] !~ /\*\*\*\*/) {

      print GBS_OUT "$text1[$startline1+$ii]";
      $ii++;

   }
}


sub print_element_file2 {

   $ii = 1;
   while ($text2[$startline2+$ii] !~ /\*\*\*\*/) {

      print GBS_OUT "$text2[$startline2+$ii]";
      $ii++;

   } 
}


sub print_usage {

   print "\nusage: emsl_manipulate.pl file1 file2 [exclusive,inclusive,insert_hhe,insert_alar] [quiet]\n";
   print   "       where file1 and file2 are downloads of gaussian-format basis sets \n";
   print   "          or diffuse/polarization/tight functions that are to be \n";
   print   "          merged into a single gaussian-format gbs file\n";
   print   "       a third argument of 'exclusive' (default) will only merge \n";
   print   "          elements when the element is available in both files\n";
   print   "       a third argument of 'inclusive' (must be specified) will only \n";
   print   "          merge elements when the element is available in ONE file \n";
   print   "          (useful to combine heavy 3df polar. with H-He 2p polar.)\n";
   print   "       a third argument of 'insert_hhe' uses available elements from the first\n";
   print   "          file, except for h-he which are taken from the second\n";
   print   "       a third argument of 'insert_alar' uses available elements from the first\n";
   print   "          file, except for al-si-p-s-cl-ar which are taken from the second\n";
   print   "       a fourth argument of 'quiet' writes only errors to the screen\n\n";
}

