#!/usr/bin/perl

open(HTML_OUT,">psi4bases.html");
print_html_topfile ();

@HGBS = <*.gbs>;
#@HELEM = (H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr);
@HELEM = (H,HE,LI,BE,B,C,N,O,F,NE,NA,MG,AL,SI,P,S,CL,AR,K,CA,SC,TI,V,CR,MN,FE,CO,NI,CU,ZN,GA,GE,AS,SE,BR,KR);

foreach $basisfile (@HGBS) {

   foreach $element (@HELEM) { $entryhash{$element} = "<td>&nbsp;</td>"; }

   open IN, $basisfile;
   @text = <IN>;
   close(IN);

   @ltemp = split( /\./, $basisfile);
   $basisfile = @ltemp[0];

   foreach $line (@text) {

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

   print HTML_OUT "<tr><td>$basisfile</td>";
   foreach $element (@HELEM) { print HTML_OUT "$entryhash{$element}"; }
   print HTML_OUT "</tr>\n\n";

}

close(HTML_OUT);
system("open -a /Applications/Safari.app psi4bases.html");


sub print_html_topfile {

   print HTML_OUT "<table border=\"1\" cellpadding=\"5\" cellspacing=\"5\" width=\"100%\">
<tr>
<th>Basis&nbsp;Set&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</th>
<th>H</th><th>He</th>
<th>Li</th><th>Be</th><th>B</th><th>C</th><th>N</th><th>O</th><th>F</th><th>Ne</th>
<th>Na</th><th>Mg</th><th>Al</th><th>Si</th><th>P</th><th>S</th><th>Cl</th><th>Ar</th>
<th>K</th><th>Ca</th><th>Sc</th><th>Ti</th><th>V</th><th>Cr</th><th>Mn</th><th>Fe</th><th>Co</th>
<th>Ni</th><th>Cu</th><th>Zn</th><th>Ga</th><th>Ge</th><th>As</th><th>Se</th><th>Br</th><th>Kr</th>
</tr>
<tr>
\n";

}



