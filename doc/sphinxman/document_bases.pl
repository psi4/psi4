#!/usr/bin/perl

my $DriverPath = "";
if ($#ARGV == 0) { $DriverPath = $ARGV[0] . "/"; }
my $BasisFolder = $DriverPath . "../../lib/basis/";

# Create a HTML table for the Trac page
open(HTML_OUT,">psi4bases.html") or die "I can't write to psi4bases.html\n";
print_html_topfile ();
# Make a LaTeX version for the manual, too
open(LATEX_OUT,">psi4bases.tex") or die "I can't write to psi4bases.tex\n";
print_latex_topfile();
# Make a RST table for the newest incarnation of the manual, too
open(RST_OUT,">source/autodoc_psi4bases.rst") or die "I can't write to psi4basis.rst\n";
print_rst_topfile();

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


%LELEM = (

   "H"   => "\\cellcolor{R0} H  & ",
   "HE"  => "\\cellcolor{R0} He & ",
   "LI"  => "\\cellcolor{R1} Li & ",
   "BE"  => "\\cellcolor{R1} Be & ",
   "B"   => "\\cellcolor{R1} B  & ",
   "C"   => "\\cellcolor{R1} C  & ",
   "N"   => "\\cellcolor{R1} N  & ",
   "O"   => "\\cellcolor{R1} O  & ",
   "F"   => "\\cellcolor{R1} F  & ",
   "NE"  => "\\cellcolor{R1} Ne & ",
   "NA"  => "\\cellcolor{R2} Na & ",
   "MG"  => "\\cellcolor{R2} Mg & ",
   "AL"  => "\\cellcolor{R2} Al & ",
   "SI"  => "\\cellcolor{R2} Si & ",
   "P"   => "\\cellcolor{R2} P  & ",
   "S"   => "\\cellcolor{R2} S  & ",
   "CL"  => "\\cellcolor{R2} Cl & ",
   "AR"  => "\\cellcolor{R2} Ar & ",
   "K"   => "\\cellcolor{R3} K  & ",
   "CA"  => "\\cellcolor{R3} Ca & ",
   "SC"  => "\\cellcolor{R3} Sc & ",
   "TI"  => "\\cellcolor{R3} Ti & ",
   "V"   => "\\cellcolor{R3} V  & ",
   "CR"  => "\\cellcolor{R3} Cr & ",
   "MN"  => "\\cellcolor{R3} Mn & ",
   "FE"  => "\\cellcolor{R3} Fe & ",
   "CO"  => "\\cellcolor{R3} Co & ",
   "NI"  => "\\cellcolor{R3} Ni & ",
   "CU"  => "\\cellcolor{R3} Cu & ",
   "ZN"  => "\\cellcolor{R3} Zn & ",
   "GA"  => "\\cellcolor{R3} Ga & ",
   "GE"  => "\\cellcolor{R3} Ge & ",
   "AS"  => "\\cellcolor{R3} As & ",
   "SE"  => "\\cellcolor{R3} Se & ",
   "BR"  => "\\cellcolor{R3} Br & ",
   "KR"  => "\\cellcolor{R3} Kr & ",
);

%SELEM = (

   "H"   => " |H|  |",
   "HE"  => " |He| |",
   "LI"  => " |Li| |",
   "BE"  => " |Be| |",
   "B"   => " |B|  |",
   "C"   => " |C|  |",
   "N"   => " |N|  |",
   "O"   => " |O|  |",
   "F"   => " |F|  |",
   "NE"  => " |Ne| |",
   "NA"  => " |Na| |",
   "MG"  => " |Mg| |",
   "AL"  => " |Al| |",
   "SI"  => " |Si| |",
   "P"   => " |P|  |",
   "S"   => " |S|  |",
   "CL"  => " |Cl| |",
   "AR"  => " |Ar| |",
   "K"   => " |K|  |",
   "CA"  => " |Ca| |",
   "SC"  => " |Sc| |",
   "TI"  => " |Ti| |",
   "V"   => " |V|  |",
   "CR"  => " |Cr| |",
   "MN"  => " |Mn| |",
   "FE"  => " |Fe| |",
   "CO"  => " |Co| |",
   "NI"  => " |Ni| |",
   "CU"  => " |Cu| |",
   "ZN"  => " |Zn| |",
   "GA"  => " |Ga| |",
   "GE"  => " |Ge| |",
   "AS"  => " |As| |",
   "SE"  => " |Se| |",
   "BR"  => " |Br| |",
   "KR"  => " |Kr| |",
);

opendir(BASISDIR, $BasisFolder) or die "I can't read $BasisFolder\n";
@HGBS = grep { /^.*\.gbs$/ } readdir(BASISDIR);
closedir(BASISDIR);

$divdash = '+' . ('-' x 47) . '+' . ('-' x 8) . ('+------' x ($#HELEM+1)) . '+';
$diveqal = '+' . ('=' x 47) . '+' . ('=' x 8) . ('+======' x ($#HELEM+1)) . '+';
printf RST_OUT "   $divdash\n   | %-45s | %-6s |", "Basis Set", "PUREAM";
foreach $element (@HELEM) { print RST_OUT "$SELEM{$element}"; }
print RST_OUT "\n   $diveqal\n";

print HTML_OUT "<tr><td><big><b>Pople</b></big></td></tr>\n\n";
print LATEX_OUT "\n\n\\\\*\n{\\normalsize\\textbf{Pople}} \\\\*\n";
printf RST_OUT "   $divdash\n   | %-45s | %6s | %249s |\n   $divdash\n", "**Pople**", "", "";
@HPOPLE = ("sto-3g.gbs","3-21g.gbs");
foreach $b (@HGBS) { if ($b =~ /^6-31g.*\.gbs$/) { push(@HPOPLE, $b); } }
foreach $b (@HGBS) { if ($b =~ /^6-31p.*\.gbs$/) { push(@HPOPLE, $b); } }
foreach $b (@HGBS) { if ($b =~ /^6-311.*\.gbs$/) { push(@HPOPLE, $b); } }
foreach $gbs (@HPOPLE) { htmlline($gbs); }
printf "Auto-documenting basis set files %s\n", "Pople";

print HTML_OUT "<tr><td><big><b>Ahlrichs/Karlsruhe</b></big></td></tr>\n\n";
print LATEX_OUT "\n\n\\\\*\n{\\normalsize\\textbf{Ahlrichs/Karlsruhe}} \\\\*\n";
printf RST_OUT "   $divdash\n   | %-45s | %6s | %249s |\n   $divdash\n", "**Ahlrichs/Karlsruhe**", "", "";
@HAHLRICHS = ();
foreach $b (@HGBS) { if ($b =~ /^def2-.*\.gbs$/) { push(@HAHLRICHS, $b); } }
foreach $gbs (@HAHLRICHS) { htmlline($gbs); }
printf "Auto-documenting basis set files %s\n", "Ahlrichs";

print HTML_OUT "<tr><td><big><b>Dunning D&#950;</b></big></td></tr>\n\n";
print LATEX_OUT "\n\n\\\\*\n{\\normalsize\\textbf{Dunning D\$\\bm\\zeta\$}} \\\\*\n";
printf RST_OUT "   $divdash\n   | %-45s | %6s | %249s |\n   $divdash\n", "**Dunning D-zeta**", "", "";
@HDDUNNING = ();
foreach $b (@HGBS) { if ($b =~ /^cc-.*dz.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^cc-.*dpd.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*dz.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*dpd.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*dz.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*dpd.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*dz.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*dpd.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*dz.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*dpd.*\.gbs$/) { push(@HDDUNNING, $b); } }
foreach $gbs (@HDDUNNING) { htmlline($gbs); }
printf "Auto-documenting basis set files %s\n", "Dunning double-zeta";

print HTML_OUT "<tr><td><big><b>Dunning T&#950;</b></big></td></tr>\n\n";
print LATEX_OUT "\n\n\\\\*\n{\\normalsize\\textbf{Dunning T\$\\bm\\zeta\$}} \\\\*\n";
printf RST_OUT "   $divdash\n   | %-45s | %6s | %249s |\n   $divdash\n", "**Dunning T-zeta**", "", "";
@HTDUNNING = ();
foreach $b (@HGBS) { if ($b =~ /^cc-.*tz.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^cc-.*tpd.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*tz.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*tpd.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*tz.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*tpd.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*tz.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*tpd.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*tz.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*tpd.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^may-cc-.*tz.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^may-cc-.*tpd.*\.gbs$/) { push(@HTDUNNING, $b); } }
foreach $gbs (@HTDUNNING) { htmlline($gbs); }
printf "Auto-documenting basis set files %s\n", "Dunning triple-zeta";

print HTML_OUT "<tr><td><big><b>Dunning Q&#950;</b></big></td></tr>\n\n";
print LATEX_OUT "\n\n\\\\*\n{\\normalsize\\textbf{Dunning Q\$\\bm\\zeta\$}} \\\\*\n";
printf RST_OUT "   $divdash\n   | %-45s | %6s | %249s |\n   $divdash\n", "**Dunning Q-zeta**", "", "";
@HQDUNNING = ();
foreach $b (@HGBS) { if ($b =~ /^cc-.*qz.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^cc-.*qpd.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*qz.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*qpd.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*qz.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*qpd.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*qz.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*qpd.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*qz.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*qpd.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^may-cc-.*qz.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^may-cc-.*qpd.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^apr-cc-.*qz.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^apr-cc-.*qpd.*\.gbs$/) { push(@HQDUNNING, $b); } }
foreach $gbs (@HQDUNNING) { htmlline($gbs); }
printf "Auto-documenting basis set files %s\n", "Dunning quadruple-zeta";

print HTML_OUT "<tr><td><big><b>Dunning 5&#950;</b></big></td></tr>\n\n";
print LATEX_OUT "\n\n\\\\*\n{\\normalsize\\textbf{Dunning 5\$\\bm\\zeta\$}} \\\\*\n";
printf RST_OUT "   $divdash\n   | %-45s | %6s | %249s |\n   $divdash\n", "**Dunning 5-zeta**", "", "";
@H5DUNNING = ();
foreach $b (@HGBS) { if ($b =~ /^cc-.*5z.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^cc-.*5pd.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*5z.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*5pd.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*5z.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*5pd.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*5z.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*5pd.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*5z.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*5pd.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^may-cc-.*5z.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^may-cc-.*5pd.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^apr-cc-.*5z.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^apr-cc-.*5pd.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^mar-cc-.*5z.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^mar-cc-.*5pd.*\.gbs$/) { push(@H5DUNNING, $b); } }
foreach $gbs (@H5DUNNING) { htmlline($gbs); }
printf "Auto-documenting basis set files %s\n", "Dunning 5-zeta";

print HTML_OUT "<tr><td><big><b>Dunning 6&#950;</b></big></td></tr>\n\n";
print LATEX_OUT "\n\n\\\\*\n{\\normalsize\\textbf{Dunning 6\$\\bm\\zeta\$}} \\\\*\n";
printf RST_OUT "   $divdash\n   | %-45s | %6s | %249s |\n   $divdash\n", "**Dunning 6-zeta**", "", "";
@H6DUNNING = ();
foreach $b (@HGBS) { if ($b =~ /^cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^aug-cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^d-aug-cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^heavy-aug-cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^jun-cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^may-cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^may-cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^apr-cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^apr-cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^mar-cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^mar-cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^feb-cc-.*6z.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $b (@HGBS) { if ($b =~ /^feb-cc-.*6pd.*\.gbs$/) { push(@H6DUNNING, $b); } }
foreach $gbs (@H6DUNNING) { htmlline($gbs); }
printf "Auto-documenting basis set files %s\n", "Dunning 6-zeta";

@HCATBASIS = (@HPOPLE,@HAHLRICHS,@HDDUNNING,@HTDUNNING,@HQDUNNING,@H5DUNNING,@H6DUNNING);

@HREST = ();
foreach my $item (@HGBS) {
   if (grep {$_ eq $item} @HCATBASIS) {
   }
   else {
      push(@HREST,$item);
   }
}

print HTML_OUT "<tr><td><big><b>Others</b></big></td></tr>\n\n";
print LATEX_OUT "\n\n\\\\*\n{\\normalsize\\textbf{Others}} \\\\*\n";
printf RST_OUT "   $divdash\n   | %-45s | %6s | %249s |\n   $divdash\n", "**Others**", "", "";
foreach $gbs (@HREST) { htmlline($gbs); }
printf "Auto-documenting basis set files %s\n", "others";

sub htmlline {
   $basisfile = shift;

   $puream = "";
   foreach $element (@HELEM) { 
      $entryhash{$element} = "<td>&nbsp;</td>";
      $latexhash{$element} = " & ";
      $rsthash{$element}   = "      |";
   }

   open IN, "$BasisFolder/$basisfile";
   @text = <IN>;
   close(IN);

   @ltemp = split( /\./, $basisfile);
   $basisfile = @ltemp[0];
   $basisfilepl = $basisfile;
   $basisfilepl =~ s/_/\\_/g;

   foreach $line (@text) {

      if ($line =~ /cartesian/)    { $puream = "6D/10F"; }
      elsif ($line =~ /spherical/) { $puream = "5D/7F";  }

      if ( ($line =~ /\s0\s/) && ($alert eq "yes") ) {

         @ltemp = split( /\s+/, $line);
         $element = @ltemp[0];
         $element =~ tr/a-z/A-Z/;
         $entryhash{$element} = "$DELEM{$element}";
         $latexhash{$element} = "$LELEM{$element}";
         $rsthash{$element}   = "$SELEM{$element}";

      }

      if ($line =~ /\*\*\*\*/) { $alert = "yes"; }
      else                     { $alert = "no";  }

   }

   print HTML_OUT "<tr><td>$basisfile</td><td>$puream</td>";
   foreach $element (@HELEM) { print HTML_OUT "$entryhash{$element}"; }
   print HTML_OUT "</tr>\n\n";

   print LATEX_OUT "$basisfilepl & $puream & ";
   foreach $element (@HELEM) { print LATEX_OUT "$latexhash{$element}"; }
   print LATEX_OUT "\\\\\n\n";

   $basislink = ":srcbasis:`" . $basisfile . "`";
   printf RST_OUT "   | %-45s | %-6s |", $basislink, $puream;
   foreach $element (@HELEM) { print RST_OUT "$rsthash{$element}"; }
   print RST_OUT "\n   $divdash\n";
}

print_latex_endfile();
print RST_OUT "\n";
close(RST_OUT);
close(LATEX_OUT);
close(HTML_OUT);
#system("open -a /Applications/Safari.app psi4bases.html");


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


sub print_rst_topfile {

   print RST_OUT ".. role:: pbg\n";
   print RST_OUT ".. raw:: html\n";
   print RST_OUT "\n";
   print RST_OUT "   <script type=\"text/javascript\" src=\"http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js\"></script>\n";
   print RST_OUT "   <script>\n";
   print RST_OUT "     \$(document).ready(function() {\n";
   print RST_OUT "       \$(\".pbg\").parent().addClass(\"pbg-parent\");\n";
   print RST_OUT "     });\n";
   print RST_OUT "   </script>\n";
   print RST_OUT "   <style>\n";
   print RST_OUT "      .pbg-parent {background-color:#ffdddd;}\n";
   print RST_OUT "   </style>\n";
   print RST_OUT "\n";
   print RST_OUT ".. role:: gbg\n";
   print RST_OUT ".. raw:: html\n";
   print RST_OUT "\n";
   print RST_OUT "   <script type=\"text/javascript\" src=\"http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js\"></script>\n";
   print RST_OUT "   <script>\n";
   print RST_OUT "     \$(document).ready(function() {\n";
   print RST_OUT "       \$(\".gbg\").parent().addClass(\"gbg-parent\");\n";
   print RST_OUT "     });\n";
   print RST_OUT "   </script>\n";
   print RST_OUT "   <style>\n";
   print RST_OUT "      .gbg-parent {background-color:#cceecc;}\n";
   print RST_OUT "   </style>\n";
   print RST_OUT "\n";
   print RST_OUT ".. role:: bbg\n";
   print RST_OUT ".. raw:: html\n";
   print RST_OUT "\n";
   print RST_OUT "   <script type=\"text/javascript\" src=\"http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js\"></script>\n";
   print RST_OUT "   <script>\n";
   print RST_OUT "     \$(document).ready(function() {\n";
   print RST_OUT "       \$(\".bbg\").parent().addClass(\"bbg-parent\");\n";
   print RST_OUT "     });\n";
   print RST_OUT "   </script>\n";
   print RST_OUT "   <style>\n";
   print RST_OUT "      .bbg-parent {background-color:#99ccee;}\n";
   print RST_OUT "   </style>\n";
   print RST_OUT "\n";
   print RST_OUT ".. role:: lbg\n";
   print RST_OUT ".. raw:: html\n";
   print RST_OUT "\n";
   print RST_OUT "   <script type=\"text/javascript\" src=\"http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js\"></script>\n";
   print RST_OUT "   <script>\n";
   print RST_OUT "     \$(document).ready(function() {\n";
   print RST_OUT "       \$(\".lbg\").parent().addClass(\"lbg-parent\");\n";
   print RST_OUT "     });\n";
   print RST_OUT "   </script>\n";
   print RST_OUT "   <style>\n";
   print RST_OUT "      .lbg-parent {background-color:#ccccee;}\n";
   print RST_OUT "   </style>\n";
   print RST_OUT "\n";
   print RST_OUT ".. |H|  replace:: :pbg:`H`\n";
   print RST_OUT ".. |He| replace:: :pbg:`He`\n";
   print RST_OUT ".. |Li| replace:: :bbg:`Li`\n";
   print RST_OUT ".. |Be| replace:: :bbg:`Be`\n";
   print RST_OUT ".. |B|  replace:: :bbg:`B`\n";
   print RST_OUT ".. |C|  replace:: :bbg:`C`\n";
   print RST_OUT ".. |N|  replace:: :bbg:`N`\n";
   print RST_OUT ".. |O|  replace:: :bbg:`O`\n";
   print RST_OUT ".. |F|  replace:: :bbg:`F`\n";
   print RST_OUT ".. |Ne| replace:: :bbg:`Ne`\n";
   print RST_OUT ".. |Na| replace:: :lbg:`Na`\n";
   print RST_OUT ".. |Mg| replace:: :lbg:`Mg`\n";
   print RST_OUT ".. |Al| replace:: :lbg:`Al`\n";
   print RST_OUT ".. |Si| replace:: :lbg:`Si`\n";
   print RST_OUT ".. |P|  replace:: :lbg:`P`\n";
   print RST_OUT ".. |S|  replace:: :lbg:`S`\n";
   print RST_OUT ".. |Cl| replace:: :lbg:`Cl`\n";
   print RST_OUT ".. |Ar| replace:: :lbg:`Ar`\n";
   print RST_OUT ".. |K|  replace:: :gbg:`K`\n";
   print RST_OUT ".. |Ca| replace:: :gbg:`Ca`\n";
   print RST_OUT ".. |Sc| replace:: :gbg:`Sc`\n";
   print RST_OUT ".. |Ti| replace:: :gbg:`Ti`\n";
   print RST_OUT ".. |V|  replace:: :gbg:`V`\n";
   print RST_OUT ".. |Cr| replace:: :gbg:`Cr`\n";
   print RST_OUT ".. |Mn| replace:: :gbg:`Mn`\n";
   print RST_OUT ".. |Fe| replace:: :gbg:`Fe`\n";
   print RST_OUT ".. |Co| replace:: :gbg:`Co`\n";
   print RST_OUT ".. |Ni| replace:: :gbg:`Ni`\n";
   print RST_OUT ".. |Cu| replace:: :gbg:`Cu`\n";
   print RST_OUT ".. |Zn| replace:: :gbg:`Zn`\n";
   print RST_OUT ".. |Ga| replace:: :gbg:`Ga`\n";
   print RST_OUT ".. |Ge| replace:: :gbg:`Ge`\n";
   print RST_OUT ".. |As| replace:: :gbg:`As`\n";
   print RST_OUT ".. |Se| replace:: :gbg:`Se`\n";
   print RST_OUT ".. |Br| replace:: :gbg:`Br`\n";
   print RST_OUT ".. |Kr| replace:: :gbg:`Kr`\n";

   print RST_OUT "\n\n.. _`table:basisElement`:\n";
   print RST_OUT "\n.. table:: Basis set availability by element in |PSIfour|\n\n";

}



sub print_latex_topfile {

   print LATEX_OUT "\\section{Basis Set Availability by Element}\\label{basisElement}\n\n";
   print LATEX_OUT "{\n\\tiny\n";
   print LATEX_OUT "\\definecolor{W0}{cmyk}{0.0,0.0,0.0,0.05}\n";
   print LATEX_OUT "\\definecolor{R0}{cmyk}{0.0,0.2,0.2,0.0}\n";
   print LATEX_OUT "\\definecolor{R1}{cmyk}{0.18,0.04,0.0,0.07}\n";
   print LATEX_OUT "\\definecolor{R2}{cmyk}{0.11,0.11,0.0,0.05}\n";
   print LATEX_OUT "\\definecolor{R3}{cmyk}{0.13,0.0,0.16,0.11}\n\n";

   print LATEX_OUT "\\renewcommand{\\arraystretch}{0.5}\n";
   #print LATEX_OUT "\\tiny\n";
   #print LATEX_OUT "\\begin{center}\n";
   print LATEX_OUT "\\begin{landscape}\n";
   print LATEX_OUT "\\begin{longtable}{p{3cm} p{1cm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm}\n";
   print LATEX_OUT "      p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm}\n";
   print LATEX_OUT "      p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm} p{1mm}}\n";
   print LATEX_OUT "\\caption{Element availibility for basis sets built into \\PSIfour.} \\label{table:basisElement} \\\\\n";

   print LATEX_OUT "\\hline\\hline\n";
   print LATEX_OUT "Basis Set & puream & \\cellcolor{R0} H & \\cellcolor{R0} He &\n";
   print LATEX_OUT "\\cellcolor{R1} Li & \\cellcolor{R1} Be & \\cellcolor{R1} B & \\cellcolor{R1} C & \\cellcolor{R1} N & \\cellcolor{R1} O & \n";
   print LATEX_OUT "\\cellcolor{R1} F & \\cellcolor{R1} Ne & \\cellcolor{R2} Na & \\cellcolor{R2} Mg & \\cellcolor{R2} Al & \\cellcolor{R2} Si & \n";
   print LATEX_OUT "\\cellcolor{R2} P & \\cellcolor{R2} S & \\cellcolor{R2} Cl & \\cellcolor{R2} Ar & \\cellcolor{R3} K & \\cellcolor{R3} Ca & \n";
   print LATEX_OUT "\\cellcolor{R3} Sc & \\cellcolor{R3} Ti & \\cellcolor{R3} V & \\cellcolor{R3} Cr & \\cellcolor{R3} Mn & \\cellcolor{R3} Fe & \n";
   print LATEX_OUT "\\cellcolor{R3} Co & \\cellcolor{R3} Ni & \\cellcolor{R3} Cu & \\cellcolor{R3} Zn & \\cellcolor{R3} Ga & \\cellcolor{R3} Ge & \n";
   print LATEX_OUT "\\cellcolor{R3} As & \\cellcolor{R3} Se \\cellcolor{R3} & \\cellcolor{R3} Br & \\cellcolor{R3} Kr & \\\\\n";
   print LATEX_OUT "\\hline\n";
   print LATEX_OUT "\\endfirsthead\n\n";
   print LATEX_OUT "\\hline\n";
   print LATEX_OUT "\\endhead\n\n";
   print LATEX_OUT "\\hline\\hline\n";
   print LATEX_OUT "\\endfoot\n\n";
}


sub print_latex_endfile {

   print LATEX_OUT "\\hline\\hline\n";
   print LATEX_OUT "\\end{longtable}\n";
   print LATEX_OUT "\\end{landscape}\n";
   #print LATEX_OUT "\\end{center}\n";
   #print LATEX_OUT "\\tiny\n";
   print LATEX_OUT "\\renewcommand{\\arraystretch}{1.0}\n";
   print LATEX_OUT "}\n";
   print LATEX_OUT "\\newpage\n\n";

}

