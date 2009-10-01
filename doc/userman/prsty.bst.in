%%% ======================================================================
%%%  @BibTeX-style-file{
%%%     filename        = "prsty.bst",
%%%     version         = "3.0",
%%%     date            = "November 10, 1992",
%%%     ISO-date        = "1992.11.10",
%%%     time            = "15:41:54.18 EST",
%%%     author          = "American Physical Society",
%%%     contact         = "Christopher B. Hamlin",
%%%     address         = "APS Publications Liaison Office
%%%                        500 Sunnyside Blvd.
%%%                        Woodbury, NY 11797",
%%%     telephone       = "(516) 576-2390",
%%%     FAX             = "(516) 349-7817",
%%%     email           = "mis@aps.org (Internet)",
%%%     supported       = "yes",
%%%     archived        = "pinet.aip.org/pub/revtex,
%%%                        Niord.SHSU.edu:[FILESERV.REVTEX]",
%%%     keywords        = "REVTeX, version 3.0, Physical Review BibTeX
%%%                        style, American Physical Society",
%%%     codetable       = "ISO/ASCII",
%%%     checksum        = "55300 926 2632 20018",
%%%     docstring       = "This is the Physical Review--style BibTeX style
%%%                        file for the American Physical Society under
%%%                        REVTeX 3.0 (release of November 10, 1992).
%%%
%%%
%%%                        The checksum field above contains a CRC-16
%%%                        checksum as the first value, followed by the
%%%                        equivalent of the standard UNIX wc (word
%%%                        count) utility output of lines, words, and
%%%                        characters.  This is produced by Robert
%%%                        Solovay's checksum utility."
%%% }
%%% ======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           begin prsty.bst                                 %
%                             version 3.0                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%   This file is part of the APS files in the REVTeX 3.0 distribution.      %
%   Version 3.0 of REVTeX, November 10, 1992.                               %
%                                                                           %
%   Copyright (c) 1992 The American Physical Society.                       %
%                                                                           %
%   See the REVTeX 3.0 README file for restrictions and more information.   %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       This is a bibtex style file to produce Physical-Review--style       %
%            references for inclusion in LaTeX/REVTeX documents.            %
%                Main problem: crossref handling is terrible                %
%               Some documentation is available (e-mail below)              %
%               Contact: C. Hamlin, mis@aps.org  (Internet)                 %
%                                   mis@apsedoff (bitnet)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Added "newmisc" entry so that notes at the end of an entry will look
%% right.  Hopefully this won't break other things.
%% -Daniel Crawford, Januray, 1997
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
  }
  {}
  { label }


INTEGERS { vol.used i.temp }

STRINGS { s t }



FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'write$
  if$
}

% Note: nothing left on stack for future processing in bibitem.
FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  newline$
}

FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\em " swap$ * "}" * }
  if$
}

FUNCTION {embolden}
{ duplicate$ empty$
    { pop$ "" }
    { "{\bf " swap$ * "}" * }
  if$
}

FUNCTION {paren}
{ duplicate$ empty$
    { pop$ "" }
    { "(" swap$ * ")" * }
  if$
}

INTEGERS { nameptr namesleft numnames }
INTEGERS { etal }

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames #100 >
  s numnames "{ll}" format.name$ "others" = numnames #1 > and
  or 'etal :=
  etal
    { #1 #1 + 'namesleft := }
    { numnames 'namesleft := }
  if$
    { namesleft #0 > }
    { s nameptr "{f. }{vv~}{ll}{, jj}" format.name$ 't :=
      nameptr #1 >
        { namesleft #1 >
            { ", " * t * }
            { nameptr #2 >
                { "," * }
                'skip$
              if$
              t "others" =
              etal or
                { " {\it et~al.}" * }
                { " and " * t * }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { author format.names }
  if$
}

FUNCTION {format.edited}
{ editor empty$
    { "" }
    { "edited by " editor format.names * }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    'title
  if$
}

FUNCTION {lc.first.letter}
{ 't :=
  ""
  t #1 #1 substring$
  "l" change.case$ *
  t #2 global.max$ substring$
  *
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "--" *
              t #2 global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}


FUNCTION {first.page}
{ 't :=
  ""
    {  t empty$ not t #1 #1 substring$ "-" = not and }
    { t #1 #1 substring$ *
      t #2 global.max$ substring$ 't :=
    }
  while$
}

FUNCTION {format.date}
{ year empty$
    { "" }
    'year
  if$
}

FUNCTION {format.btitle}
{ title emphasize }

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}


FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
        { "pp.\ " pages n.dashify * }
        { "p.\ " pages * }
      if$
    }
  if$
}

FUNCTION {format.pages.a}
{ pages empty$
    { "" }
    { "page" pages first.page tie.or.space.connect }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume field.or.null embolden
  " " swap$ * *
  pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages.a }
        { ", " * pages first.page * }
      if$
    }
  if$
}

FUNCTION {if.comma}
{   { ", " * }
    'skip$
  if$
}

FUNCTION {if.period}
{   { ". " * }
    'skip$
  if$
}

% pushes "" if value passed is empty$
FUNCTION {push.string}
{ 'i.temp :=
  duplicate$ empty$
    { pop$ "" }
    { i.temp if.comma }
  if$
}

% pushes "" if value passed is empty$
FUNCTION {push.string.check}
{ 'i.temp :=
  duplicate$ empty$
    { pop$ "Warning: " swap$ * " missing in entry " cite$ "." * * * warning$
      ""
     }
    { swap$ pop$ i.temp if.comma }
  if$
}

% adds nothing if value passed is empty$
FUNCTION {add.string}
{ 'i.temp :=
  duplicate$ empty$
    { pop$ }
    { * i.temp if.comma }
  if$
}

% adds nothing if value passed is empty$
FUNCTION {new.add.period}
{ 'i.temp :=
  duplicate$ empty$
    { pop$ }
    { * i.temp if.period }
  if$
}

% adds nothing if value passed is empty$
FUNCTION {add.string.check}
{ 'i.temp :=
  duplicate$ empty$
    { pop$ "Warning: " swap$ * " missing in entry " cite$ "." * * * warning$ }
    { swap$ pop$ *
      i.temp if.comma
    }
  if$
}

FUNCTION {format.vol.page}
{ "volume" volume embolden #1 push.string.check
  " " *
  "pages" pages field.or.null first.page #0 add.string.check
  " " *
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  year empty$ note empty$ and and and and
    { "All relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { type empty$
        { " No.~" * }
        { " " * }
      if$
      number *
    }
  if$
}


FUNCTION {format.etc}
{ vol.used
    'skip$
    { number empty$
        { volume empty$
            'skip$
            { "Vol.~" volume * *
              chapter empty$ pages empty$ note empty$ and and
                  'skip$
                  { ", " * }
              if$
            }
          if$
        }
        { "No.~" number * *
          chapter empty$ pages empty$ note empty$ and and
              'skip$
              { ", " * }
          if$
        }
      if$
    }
  if$
  chapter empty$
    'skip$
    { type empty$
        { "Chap.~" chapter * * }
        { type "u" change.case$ "l" change.case$ " " chapter * * * }
      if$
      pages empty$ note empty$ and
        'skip$
        { ", " * }
      if$
    }
  if$
  pages empty$
    'skip$
    { format.pages *
      note empty$
        'skip$
        { ", " * }
      if$
    }
  if$
  note empty$
    'skip$
    { note lc.first.letter * }
  if$
}


FUNCTION {format.pub.add.year.etc}
{ " ("
  publisher empty$
    { "PUBLISHER" "No publisher for entry " cite$ "." * * warning$ }
    { publisher }
  if$
  ", " * *
  address empty$
    { "ADDRESS" "Missing address in entry " cite$ * "." * warning$ }
    { address }
  if$
  ", " * *
  year empty$
   { "YEAR" }
   'year
  if$
  ")" * *
  volume empty$ not number empty$ not or title empty$ series empty$ or and
  chapter empty$ not or
  pages empty$ not or type empty$ not or note empty$ not or
    { ", " *
      format.etc
    }
    'skip$
  if$
}

FUNCTION {format.bedition}
{ edition empty$
    { "" }
    { edition "l" change.case$
      duplicate$ "first" =
        { pop$ "1st " }
        { duplicate$ "second" =
            { pop$ "2nd " }
            { duplicate$ "third" =
                { pop$ "3rd " }
                'skip$
              if$
            }
          if$
        }
      if$
      " ed." *
    }
  if$
}

FUNCTION {format.title.and.series.and.org}
{ #0 'vol.used :=
  ""
  duplicate$ 's :=
  title empty$
    'skip$
    { title 's := }
  if$
  booktitle empty$
    'skip$
    { booktitle 's := }
  if$
  s empty$
    'skip$
    { type$ "inproceedings" = type$ "proceedings" = or type$ "inbook" = or
        { " in " * }
        { author empty$ not editor empty$ not and
            { " in " * }
            'skip$
          if$
        }
      if$
      s emphasize
        series empty$ editor empty$ edition empty$ and and not
        type$ "inproceedings" =   organization empty$ not    and
        type$ "proceedings"   =   organization empty$ not    and or     or
          add.string
    }
  if$
  series empty$
    'skip$
    { s empty$
         'skip$
         { number empty$
             { volume empty$
                 'skip$
                 { "Vol.~" volume " of " * * * #1 'vol.used := }
               if$
             }
             { "No.~" number " in " * * * #1 'vol.used := }
           if$
         }
      if$
      series emphasize
        edition empty$ not
        editor empty$ not
        type$ "inproceedings" =   organization empty$ not    and
        type$ "proceedings" =   organization empty$ not      and
        or   or   or
          add.string
    }
  if$
  type$ "inproceedings" =    type$ "proceedings" =    or
    { organization empty$
        'skip$
        { organization editor empty$ edition empty$ and not add.string }
      if$
    }
    'skip$
  if$
}


FUNCTION {article}
{ output.bibitem
  "author" format.authors #1 push.string.check
  "journal" journal #0 add.string.check " " *
  output
  format.vol.page #0 push.string
  " " *
  "date" format.date paren #0 add.string.check
  note empty$
    'skip$
    { ", " * note lc.first.letter * }
  if$
  fin.entry
  crossref empty$
    'skip$
    { "See Ref.\ \cite{" crossref "}." * * write$ newline$ }
  if$
}

FUNCTION {book}
{ output.bibitem
  format.authors #1 push.string
  "title" format.title.and.series.and.org #0 add.string.check
  format.bedition editor empty$ not add.string
  output
  format.edited #0 push.string
  format.pub.add.year.etc #0 add.string
  fin.entry
  crossref empty$
    'skip$
    { "See Ref.\ \cite{" crossref "}." * * write$ newline$ }
  if$
}

FUNCTION {booklet}
{ output.bibitem
  format.authors #1 push.string
  "title" format.title #0 add.string.check
  howpublished empty$
    'skip$
    { ", " howpublished * * }
  if$
  address empty$
    'skip$
    { ", " address * * }
  if$
  year empty$
    'skip$
    { " (" year ")" * * * }
  if$
  note empty$
    'skip$
    { ", " note  lc.first.letter * * }
  if$
  fin.entry
}

FUNCTION {inbook} { book }

FUNCTION {incollection} { book }

FUNCTION {inproceedings}{ book }

FUNCTION {proceedings} { book }

FUNCTION {conference} { book }



FUNCTION {manual}
{ output.bibitem
  ""
  author empty$
    'skip$
    { format.authors #1 add.string }
  if$
  format.btitle
    edition empty$ address empty$ organization empty$ year empty$ note empty$
                           and and and and not
      add.string
  edition empty$
    'skip$
    { format.bedition
        organization empty$ address empty$ year empty$ note empty$
                           and and and not
          add.string
    }
  if$
  organization empty$
    'skip$
    { organization
        address empty$ year empty$ note empty$ and and not
          add.string
    }
  if$
  address empty$
    'skip$
    { address   year empty$ note empty$ and not   add.string }
  if$
  year empty$
    'skip$
    { year   note empty$ not   add.string }
  if$
  note empty$
    'skip$
    { note lc.first.letter #0 add.string }
  if$
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  "author" format.authors #1 push.string.check
  title empty$
    'skip$
    { "We don't use thesis titles in Phys. Rev, see \cite{"
      cite$ "}" * * warning$
    }
  if$
  "Master's thesis" format.thesis.type   #1   add.string
  "school"   school   #1   add.string.check
  address #1 add.string
  "year" format.date   note empty$ not   add.string.check
  note empty$
    'skip$
    { note lc.first.letter * }
  if$
  fin.entry
}

FUNCTION {phdthesis}
{ output.bibitem
  "author" format.authors #1 push.string.check
  title empty$
    'skip$
    { "We don't use thesis titles in Phys. Rev, see \cite{"
      cite$ "}" * * warning$
    }
  if$
  "Ph.D. thesis" format.thesis.type   #1   add.string
  "school"   school   #1   add.string.check
  address #1 add.string
  "year" format.date   note empty$ not   add.string.check
  note empty$
    'skip$
    { note lc.first.letter * }
  if$
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors
    title empty$ howpublished empty$
      year empty$ note empty$ and and and not push.string
  format.title
    howpublished empty$ year empty$ note empty$
      and and not add.string
  howpublished
    year empty$ note empty$ and not add.string
  format.date note empty$ not add.string
  note empty$
    'skip$
    { note
      author empty$ title empty$ howpublished empty$
      year empty$ and and and
        'skip$
        'lc.first.letter
      if$
      #0 add.string
    }
  if$
  fin.entry
  empty.misc.check
}

FUNCTION {newmisc}
{ output.bibitem
  format.authors
    title empty$ howpublished empty$
      year empty$ note empty$ and and and not push.string
  format.title
    howpublished empty$ year empty$ note empty$
      and and not add.string
  howpublished
    year empty$ note empty$ and not add.string
  format.date note empty$ not new.add.period
  note empty$
    'skip$
    { note
      author empty$ title empty$ howpublished empty$
      year empty$ and and and
        'skip$
        'lc.first.letter
      if$
      #0 add.string
    }
  if$
  fin.entry
  empty.misc.check
}

FUNCTION {techreport}
{ output.bibitem
  "author" format.authors #1 push.string.check
  title empty$
    'skip$
    { "We chuck the title of techreports in Phys. Rev.---\cite{"
       cite$ "}." * * warning$
    }
  if$
  format.tr.number
    institution empty$ address empty$ note empty$ and and not add.string
  "institution"  institution  address empty$ note empty$ and not
    add.string.check
  address #0 add.string
  " (unpublished)" note empty$ not add.string
  note empty$
    'skip$
    { note lc.first.letter * }
  if$
  fin.entry
}



FUNCTION {unpublished}
{ output.bibitem
  "author" format.authors   note empty$ not   push.string.check
  note empty$
    'skip$
    { note lc.first.letter #0 add.string }
  if$
  " (unpublished)" *
  fin.entry
}

FUNCTION {default.type} { misc }




MACRO {jan} {"January"}
MACRO {feb} {"February"}
MACRO {mar} {"March"}
MACRO {apr} {"April"}
MACRO {may} {"May"}
MACRO {jun} {"June"}
MACRO {jul} {"July"}
MACRO {aug} {"August"}
MACRO {sep} {"September"}
MACRO {oct} {"October"}
MACRO {nov} {"November"}
MACRO {dec} {"December"}






MACRO {acmcs} {"ACM Comput. Surv."}
MACRO {acta} {"Acta Inf."}
MACRO {applopt} {"Appl. Opt."}
MACRO {cacm} {"Commun. ACM"}
MACRO {ibmjrd} {"IBM J. Res. Dev."}
MACRO {ibmsj} {"IBM Syst.~J."}
MACRO {ieeese} {"IEEE Trans. Softw. Eng."}
MACRO {ieeetc} {"IEEE Trans. Comput."}
MACRO {ieeetcad}
 {"IEEE Trans. Comput.-Aided Design Integrated Circuits"}
MACRO {ipl} {"Inf. Process. Lett."}
MACRO {jacm} {"J.~ACM"}
MACRO {jcss} {"J.~Comput. Syst. Sci."}
MACRO {scp} {"Sci. Comput. Programming"}
MACRO {sicomp} {"SIAM J. Comput."}
MACRO {tocs} {"ACM Trans. Comput. Syst."}
MACRO {tods} {"ACM Trans. Database Syst."}
MACRO {tog} {"ACM Trans. Gr."}
MACRO {toms} {"ACM Trans. Math. Softw."}
MACRO {toois} {"ACM Trans. Office Inf. Syst."}
MACRO {toplas} {"ACM Trans. Prog. Lang. Syst."}
MACRO {tcs} {"Theor. Comput. Sci."}
MACRO {anp} {"Adv. Nucl. Phys."}
MACRO {ap} {"Ann. Phys."}
MACRO {arnpc} {"Annu. Rev. Nucl. Part. Sci."}
MACRO {baps} {"Bull. Am. Phys. Soc."}
MACRO {cp} {"Computers Phys."}
MACRO {jp} {"J. Phys."}
MACRO {josa} {"J. Opt. Soc. Am."}
MACRO {josaa} {"J. Opt. Soc. Am. A"}
MACRO {josab} {"J. Opt. Soc. Am. B"}
MACRO {nim} {"Nucl. Instrum. Meth."}
MACRO {np} {"Nucl. Phys."}
MACRO {rmp} {"Rev. Mod. Phys."}
MACRO {physica} {"Physica"}
MACRO {pl} {"Phys. Lett."}
MACRO {plb} {"Phys. Lett. B"}
MACRO {pr} {"Phys. Rev."}
MACRO {pra} {"Phys. Rev. A"}
MACRO {prb} {"Phys. Rev. B"}
MACRO {prc} {"Phys. Rev. C"}
MACRO {prd} {"Phys. Rev. D"}
MACRO {pre} {"Phys. Rev. E"}
MACRO {prl} {"Phys. Rev. Lett."}
MACRO {prep} {"Phys. Rep."}
MACRO {ps} {"Phys. Scr."}
MACRO {rpp} {"Rep. Prog. Phys."}
MACRO {sjnp} {"Sov. J. Nucl. Phys."}
MACRO {springmp} {"Springer Tracts in Modern Physics"}
MACRO {yadfiz} {"Yad. Fiz."}
MACRO {zp} {"Z. Phys."}



READ

STRINGS { longest.label }
INTEGERS { number.label longest.label.width }
FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}
FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}
EXECUTE {initialize.longest.label}
ITERATE {longest.label.pass}




FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{"  longest.label  * "}" * write$ newline$
}

EXECUTE {begin.bib}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}
EXECUTE {end.bib}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             end prsty.bst                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
