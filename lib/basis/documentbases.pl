#!/usr/bin/perl

open(HTML_OUT,">psi4bases.html");
print_html_topfile ();

@HGBS = <*.gbs>;
#@HELEM = (H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr);
@HELEM = (H,HE,LI,BE,B,C,N,O,F,NE,NA,MG,AL,SI,P,S,CL,AR,K,CA,SC,TI,V,CR,MN,FE,CO,NI,CU,ZN,GA,GE,AS,SE,BR,KR);

foreach $basisfile (@HGBS) {

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
         #print "$alert**$element**\n";
         $entryhash{$element} = "<td>&#10003;</td>";

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



