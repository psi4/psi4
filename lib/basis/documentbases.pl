#!/usr/bin/perl

open(HTML_OUT,">psi4bases.html");
print_html_topfile ();

@HELEM = (H,HE,LI,BE,B,C,N,O,F,NE,NA,MG,AL,SI,P,S,CL,AR,K,CA,SC,TI,V,CR,MN,FE,CO,NI,CU,ZN,GA,GE,AS,SE,BR,KR);

%DELEM = (

   "H"   => "<td bgcolor=#ffdddd>H</td>",
   "HE"  => "<td bgcolor=#ffdddd>He</td>",
   "LI"  => "<td bgcolor=#99ccee>Li</td>",
   "BE"  => "<td bgcolor=#99ccee>Be</td>",
   "B"   => "<td bgcolor=#99ccee>B</td>",
   "C"   => "<td bgcolor=#99ccee>C</td>",
   "N"   => "<td bgcolor=#99ccee>N</td>",
   "O"   => "<td bgcolor=#99ccee>O</td>",
   "F"   => "<td bgcolor=#99ccee>F</td>",
   "NE"  => "<td bgcolor=#99ccee>Ne</td>",
   "NA"  => "<td bgcolor=#ccccee>Na</td>",
   "MG"  => "<td bgcolor=#ccccee>Mg</td>",
   "AL"  => "<td bgcolor=#ccccee>Al</td>",
   "SI"  => "<td bgcolor=#ccccee>Si</td>",
   "P"   => "<td bgcolor=#ccccee>P</td>",
   "S"   => "<td bgcolor=#ccccee>S</td>",
   "CL"  => "<td bgcolor=#ccccee>Cl</td>",
   "AR"  => "<td bgcolor=#ccccee>Ar</td>",
   "K"   => "<td bgcolor=#cceecc>K</td>",
   "CA"  => "<td bgcolor=#cceecc>Ca</td>",
   "SC"  => "<td bgcolor=#cceecc>Sc</td>",
   "TI"  => "<td bgcolor=#cceecc>Ti</td>",
   "V"   => "<td bgcolor=#cceecc>V</td>",
   "CR"  => "<td bgcolor=#cceecc>Cr</td>",
   "MN"  => "<td bgcolor=#cceecc>Mn</td>",
   "FE"  => "<td bgcolor=#cceecc>Fe</td>",
   "CO"  => "<td bgcolor=#cceecc>Co</td>",
   "NI"  => "<td bgcolor=#cceecc>Ni</td>",
   "CU"  => "<td bgcolor=#cceecc>Cu</td>",
   "ZN"  => "<td bgcolor=#cceecc>Zn</td>",
   "GA"  => "<td bgcolor=#cceecc>Ga</td>",
   "GE"  => "<td bgcolor=#cceecc>Ge</td>",
   "AS"  => "<td bgcolor=#cceecc>As</td>",
   "SE"  => "<td bgcolor=#cceecc>Se</td>",
   "BR"  => "<td bgcolor=#cceecc>Br</td>",
   "KR"  => "<td bgcolor=#cceecc>Kr</td>"
);

@HGBS = <*.gbs>;

print HTML_OUT "<tr><td><big><b>Pople</b></big></td></tr>\n\n";
@HPOPLE = ("sto-3g.gbs","3-21g.gbs",<6-31[gp]*.gbs>,<6-311*.gbs>);
foreach $gbs (@HPOPLE) { htmlline($gbs); }

print HTML_OUT "<tr><td><big><b>Dunning D&#950;</b></big></td></tr>\n\n";
@HDDUNNING = (<cc-*dz*.gbs>,<cc-*dpd*.gbs>,<aug-cc-*dz*.gbs>,<aug-cc-*dpd*.gbs>,<d-aug-cc-*dz*.gbs>,<d-aug-cc-*dpd*.gbs>,<heavy-aug-cc-*dz*.gbs>,<heavy-aug-cc-*dpd*.gbs>);
foreach $gbs (@HDDUNNING) { htmlline($gbs); }

print HTML_OUT "<tr><td><big><b>Dunning T&#950;</b></big></td></tr>\n\n";
@HTDUNNING = (<cc-*tz*.gbs>,<cc-*tpd*.gbs>,<aug-cc-*tz*.gbs>,<aug-cc-*tpd*.gbs>,<d-aug-cc-*tz*.gbs>,<d-aug-cc-*tpd*.gbs>,<heavy-aug-cc-*tz*.gbs>,<heavy-aug-cc-*tpd*.gbs>);
foreach $gbs (@HTDUNNING) { htmlline($gbs); }

print HTML_OUT "<tr><td><big><b>Dunning Q&#950;</b></big></td></tr>\n\n";
@HQDUNNING = (<cc-*qz*.gbs>,<cc-*qpd*.gbs>,<aug-cc-*qz*.gbs>,<aug-cc-*qpd*.gbs>,<d-aug-cc-*qz*.gbs>,<d-aug-cc-*qpd*.gbs>,<heavy-aug-cc-*qz*.gbs>,<heavy-aug-cc-*qpd*.gbs>);
foreach $gbs (@HQDUNNING) { htmlline($gbs); }

print HTML_OUT "<tr><td><big><b>Dunning 5&#950;</b></big></td></tr>\n\n";
@H5DUNNING = (<cc-*5z*.gbs>,<cc-*5pd*.gbs>,<aug-cc-*5z*.gbs>,<aug-cc-*5pd*.gbs>,<d-aug-cc-*5z*.gbs>,<d-aug-cc-*5pd*.gbs>,<heavy-aug-cc-*5z*.gbs>,<heavy-aug-cc-*5pd*.gbs>);
foreach $gbs (@H5DUNNING) { htmlline($gbs); }

print HTML_OUT "<tr><td><big><b>Dunning 6&#950;</b></big></td></tr>\n\n";
@H6DUNNING = (<cc-*6z*.gbs>,<cc-*6pd*.gbs>,<aug-cc-*6z*.gbs>,<aug-cc-*6pd*.gbs>,<d-aug-cc-*6z*.gbs>,<d-aug-cc-*6pd*.gbs>,<heavy-aug-cc-*6z*.gbs>,<heavy-aug-cc-*6pd*.gbs>);
foreach $gbs (@H6DUNNING) { htmlline($gbs); }

@HCATBASIS = (@HPOPLE,@HDDUNNING,@HTDUNNING,@HQDUNNING,@H5DUNNING,@H6DUNNING);

@HREST = ();
foreach my $item (@HGBS) {
   if (grep {$_ eq $item} @HCATBASIS) {
   }
   else {
      push(@HREST,$item);
   }
}

print HTML_OUT "<tr><td><big><b>Others</b></big></td></tr>\n\n";
foreach $gbs (@HREST) { htmlline($gbs); }


sub htmlline {
   $basisfile = shift;

   $puream = "";
   foreach $element (@HELEM) { $entryhash{$element} = "<td>&nbsp;</td>"; }

   open IN, $basisfile;
   @text = <IN>;
   close(IN);

   @ltemp = split( /\./, $basisfile);
   $basisfile = @ltemp[0];

   foreach $line (@text) {

      if ($line =~ /cartesian/)    { $puream = "6D/10F"; }
      elsif ($line =~ /spherical/) { $puream = "5D/7F";  }

      if ( ($line =~ /\s0\s/) && ($alert eq "yes") ) {

         @ltemp = split( /\s+/, $line);
         $element = @ltemp[0];
         $element =~ tr/a-z/A-Z/;
         #$entryhash{$element} = "<td>&#10003;</td>";
         #$entryhash{$element} = "<td>$element</td>";
         $entryhash{$element} = "$DELEM{$element}";

      }

      if ($line =~ /\*\*\*\*/) { $alert = "yes"; }
      else                     { $alert = "no";  }

   }

   print HTML_OUT "<tr><td>$basisfile</td><td>$puream</td>";
   foreach $element (@HELEM) { print HTML_OUT "$entryhash{$element}"; }
   print HTML_OUT "</tr>\n\n";

}

close(HTML_OUT);
system("open -a /Applications/Safari.app psi4bases.html");


sub print_html_topfile {

   print HTML_OUT "<table border=\"1\" cellpadding=\"4\" cellspacing=\"3\" width=\"100%\">
   <tr><th>Basis&nbsp;Set&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</th><th>puream</th>
   <th bgcolor=#ffdddd>H</th><th bgcolor=#ffdddd>He</th><th bgcolor=#99ccee>Li</th><th bgcolor=#99ccee>Be</th>
   <th bgcolor=#99ccee>B</th><th bgcolor=#99ccee>C</th><th bgcolor=#99ccee>N</th><th bgcolor=#99ccee>O</th>
   <th bgcolor=#99ccee>F</th><th bgcolor=#99ccee>Ne</th><th bgcolor=#ccccee>Na</th><th bgcolor=#ccccee>Mg</th>
   <th bgcolor=#ccccee>Al</th><th bgcolor=#ccccee>Si</th><th bgcolor=#ccccee>P</th><th bgcolor=#ccccee>S</th>
   <th bgcolor=#ccccee>Cl</th><th bgcolor=#ccccee>Ar</th><th bgcolor=#cceecc>K</th><th bgcolor=#cceecc>Ca</th>
   <th bgcolor=#cceecc>Sc</th><th bgcolor=#cceecc>Ti</th><th bgcolor=#cceecc>V</th><th bgcolor=#cceecc>Cr</th>
   <th bgcolor=#cceecc>Mn</th><th bgcolor=#cceecc>Fe</th><th bgcolor=#cceecc>Co</th><th bgcolor=#cceecc>Ni</th>
   <th bgcolor=#cceecc>Cu</th><th bgcolor=#cceecc>Zn</th><th bgcolor=#cceecc>Ga</th><th bgcolor=#cceecc>Ge</th>
   <th bgcolor=#cceecc>As</th><th bgcolor=#cceecc>Se</th><th bgcolor=#cceecc>Br</th><th bgcolor=#cceecc>Kr</th>
   <tr>\n";
}

