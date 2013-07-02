%{
#include <cstdlib>
#include <string>
#include "mra-parser.hh"
#include "mra-driver.hh"

#define TOK(X) yylval->s=new std::string(yytext); return token::X;
#define RET(X) return token::X;

#undef yywrap
#define yywrap() 1
#define yyterminate() return token::END
%}

%option noyywrap nounput batch debug

COMMENT  #[^\n]*[\n]
WHITESPACE [ \t]*
EOL [\n]+

PLUS     \+
MINUS    \-
RARROW   \-\>
TIMES    "*"
DIVIDE   \/
POWER    \^
LBRACE   "{"
RBRACE   "}"
LBRACK   \[
RBRACK   \]
LPAREN   "("
RPAREN   ")"
MLINE    \|
MMLINE   \|\|
LANGLE   <
RANGLE   >
SEMICOLON ;
COMMA     ,
PERIOD   "."
ASSIGN    "="
SUBSCRIPT _

F        F
FUNCTION abs|norm|normf|sqrt|exp|sin|cos|tan|cot|sec|arcsin|arccos|arctan|sinh|cosh|tanh|coth|erf|log|ln
PRINT    print
PLOT     plot
READ     read
WRITE    write
WHILE    [wW][hH][iI][lL][eE]
FOR      [fF][oO][rR]
IF       [iI][fF]
ELIF     [eE][lL][iI][fF]
ELSE     [eE][lL][sS][eE]
SUM      [sS][uU][mM]
DEL      [dD][eD][lL]

LET  [lL][eE][tT]
IN   [iI][nN]
END  [eE][nN][dD]

OR       [oO][rR]
AND      [aA][nN][dD]
NOT      [nN][oO][tT]
NE       "!="
EQ       "=="
LE       "<="
GE       ">="

GREEK Delta|Gamma|Lambda|Phi|Pi|Psi|Sigma|Theta|Upsilon|Xi|alpha|beta|chi|delta|digamma|epsilon|eta|gamma|iota|kappa|lambda|mu|nu|omega|phi|pi|psi|rho|sigma|tau|theta|upsilon|varepsilon|varkappa|varphi|varpi|varrho|varsigma|vartheta|xi|zeta

OMEGA Omega

NABLA nabla

NAME     [a-zA-Z][a-zA-Z0-9]*
DIGIT    [0-9]
INTEGER  [0-9]+
FIELD    [RCZ]

FLOAT0   [0-9]*"."[0-9]+
FLOAT1   [0-9]+"."[0-9]*
FLOAT2   [0-9]+[eE][\+\-]?[0-9]+
FLOAT3   [0-9]*"."[0-9]+[eE][\+\-]?[0-9]+
FLOAT    {FLOAT0}|{FLOAT1}|{FLOAT2}|{FLOAT3}

STRING   \"[^\"\n]*\"

%{
# define YY_USER_ACTION  yylloc->columns (yyleng);
%}

%%

%{
  yylloc->step ();
%}

%{
  typedef yy::mra_parser::token token;
%}


{COMMENT}     yylloc->lines(1); yylloc->step(); RET(EOL);
{WHITESPACE}  yylloc->step();
{EOL}         yylloc->lines(yyleng); yylloc->step(); RET(EOL);


{FLOAT} {TOK(FLOAT)}
{INTEGER} {TOK(INTEGER)}

{F}  {RET(F)}

{FIELD} {TOK(FIELD)}

{OR} {RET(OR)}
{AND} {RET(AND)}
{NE} {RET(NE)}
{EQ} {RET(EQ)}
{NOT} {RET(NOT)}
{LE}  {RET(LE)}
{GE}  {RET(GE)}

{READ} {RET(READ)}
{WRITE} {RET(WRITE)}
{PRINT} {RET(PRINT)}
{PLOT} {RET(PLOT)}
{WHILE} {RET(WHILE)}
{FOR} {RET(FOR)}
{IF}    {RET(IF)}
{ELIF}  {RET(ELIF)}
{ELSE}  {RET(ELSE)}
{SUM}   {RET(SUM)}
{DEL}   {RET(DEL)}
{FUNCTION} {TOK(FUNCTION)}
{GREEK} {TOK(GREEK)}
{OMEGA} {RET(OMEGA)}

{LET} {RET(LET)}
{IN}  {RET(IN)}
{END} {RET(ENDPROG)}

{NAME} {TOK(NAME)}

{SEMICOLON} {RET(SEMICOLON)}
{SUBSCRIPT} {RET(SUBSCRIPT)}
{COMMA} {RET(COMMA)}
{PERIOD} {RET(PERIOD)}

{STRING} {TOK(STRING)}

{PLUS} {RET(PLUS)}
{MINUS} {RET(MINUS)}
{TIMES} {RET(TIMES)}
{DIVIDE} {RET(DIVIDE)}
{POWER} {RET(POWER)}

{LBRACE} {RET(LBRACE)}
{RBRACE} {RET(RBRACE)}
{LBRACK} {RET(LBRACK)}
{RBRACK} {RET(RBRACK)}
{LPAREN} {RET(LPAREN)}
{RPAREN} {RET(RPAREN)}
{LANGLE} {RET(LANGLE)}
{RANGLE} {RET(RANGLE)}

{ASSIGN} {RET(ASSIGN)}
{RARROW} {RET(RARROW)}
{MMLINE} {RET(MMLINE)}
{MLINE}  {RET(MLINE)}

%%

void mra_driver::scan_begin ()
{
  yy_flex_debug = trace_scanning;
  if (file == "-")
    yyin = stdin;
  else if (!(yyin = fopen (file.c_str (), "r")))
    {
      error (std::string ("cannot open ") + file);
      exit (1);
    }
}

void mra_driver::scan_end ()
{
  fclose (yyin);
}
