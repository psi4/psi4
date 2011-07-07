%skeleton "lalr1.cc"
%defines
%define parser_class_name "mra_parser"


%code requires {
# include <iostream>
# include <string>
# include <cmath>
# include <cstdio>

class mra_driver;
struct Exp;
}

%union
{
  std::string *s;
  Exp *e;
};

%parse-param { mra_driver& driver }
%lex-param   { mra_driver& driver }
%locations
%initial-action
{
  @$.begin.filename = @$.end.filename = &driver.file;
};

%debug
%error-verbose

%code {
# include "mra-driver.hh"
}

%type  <e> exp
%type  <e> statement
%type  <e> letstatement
%type  <e> assignment
%type  <e> paramlist
%type  <e> arglist
%type  <e> varlist
%type  <e> variable
%type  <e> summation
%type  <e> derivative
%type  <e> laplacian
%type  <e> bsh
%type  <e> braket
%type  <e> mline
%type  <e> norm2
%type  <e> normf

%token <s> NAME
%token <s> GREEK
%token <s> FLOAT
%token <s> INTEGER
%token <s> STRING
%token <s> FUNCTION
%token <s> FIELD

%token LPAREN
%token RPAREN
%token LBRACE
%token RBRACE
%token LBRACK
%token RBRACK
%token LANGLE
%token RANGLE
%token MLINE
%token MMLINE

%token F
%token OMEGA
%token PRINT
%token PLOT
%token READ
%token WRITE
%token WHILE
%token FOR
%token IF
%token ELIF
%token ELSE
%token SUM
%token DEL

%token SEMICOLON
%token COMMA
%token PERIOD

%token ASSIGN
%token RARROW

%token EOL
%token END 0 

%token LET
%token IN
%token ENDPROG

%left  OR
%left  AND
%left  LT GT LE GE EQ NE
%left  LANGLE RANGLE MLINE
%left  PLUS MINUS
%left  SUMMATION
%left  TIMES DIVIDE
%left  DERIVATIVE 
%left  NOT NEGATE
%nonassoc BRAKET
%right SUBSCRIPT POWER

%printer    { debug_stream () << $$.s; } <>

%%

program: LET EOL                                   { driver.add(new Exp("Let")); } 
            letstatements 
         IN EOL                                    { driver.add(new Exp("In")); } 
            statementlist 

letstatements: /*empty*/
        | letstatements letstatement

letstatement:
          EOL                                      {}
        | varlist  IN FIELD EOL                    { driver.add(new Exp("decl",$1,new Exp(*$3))); driver.addsyms($1,*$3); }
        | variable ASSIGN exp EOL                  { driver.add(new Exp("=",$1, $3, $3->type()));
                                                     if ($1->type() == Exp::UNKNOWN) {
                                                        $1->set_type($3->type());
                                                        driver.insert_sym("inferring", $1);
                                                     }
                                                   }   
        | variable ASSIGN paramlist RARROW exp EOL { driver.add(new Exp("def",$1,$3,$5,$5->type())); 
                                                     if ($1->type() == Exp::UNKNOWN) {
                                                        $1->set_type($5->type());
                                                        driver.insert_sym("inferring", $1);
                                                     }
                                                   }
        | OMEGA ASSIGN LBRACK exp COMMA exp RBRACK POWER INTEGER EOL
                                                   { driver.add(new Exp("domain",$4,$6,new Exp(*$9,Exp::INTEGER))); }

statementlist: 
          /*empty*/ 
        | statementlist statement

statement:
          EOL                                      {}
        | assignment EOL                           { driver.add($1); }
        | exp EOL                                  { driver.add($1); }
        | PRINT arglist EOL                        { driver.add(new Exp("print",$2)); }
        | PLOT arglist EOL                         { driver.add(new Exp("plot",$2)); }
        | WHILE exp EOL                            { driver.add(new Exp("while",$2)); driver.scopedepth++; } 
        | IF    exp EOL                            { driver.add(new Exp("if",$2)); driver.scopedepth++; } 
        | ELIF  exp EOL                            { driver.add(new Exp("elif",$2)); } 
        | ELSE  EOL                                { driver.add(new Exp("else")); } 
        | FOR variable                             { $2->set_type(Exp::INTEGER); driver.insert_sym("inferring",$2);}
          IN LBRACK exp COMMA exp RBRACK EOL       { driver.add(new Exp("for",$2,$6,$8)); driver.scopedepth++; }
        | ENDPROG EOL                              { if(driver.scopedepth==0) {
                                                       driver.add(new Exp("End")); // End of program
                                                     }
                                                     else {
                                                       driver.scopedepth--;
                                                       driver.add(new Exp("end")); // End of block
                                                     }
                                                   }
        
assignment:
          variable ASSIGN exp                      { $$ = new Exp("=",$1, $3, $3->type());
                                                     if ($1->type() == Exp::UNKNOWN) {
                                                        $1->set_type($3->type());
                                                        driver.insert_sym("inferring", $1);
                                                     }
                                                   }  
summation:
          SUM SUBSCRIPT variable ASSIGN exp POWER exp { 
                                                        $3->set_type(Exp::INTEGER); 
                                                        driver.insert_sym("inferring",$3); 
                                                        $$ = new Exp("sum", $3, $5, $7);                                                        
                                                      } 

derivative:
          DEL SUBSCRIPT exp                        { $$ = new Exp("del",$3); }    

laplacian: 
          DEL POWER INTEGER                        { $$ = new Exp("laplacian"); if (*$3 != "2") driver.error("bad power for laplacian?"); }

bsh:
          LPAREN exp MINUS DEL POWER INTEGER RPAREN POWER MINUS INTEGER  { $$ = new Exp("bsh",$2); }
        | DEL POWER MINUS INTEGER                  { $$ = new Exp("bsh",new Exp("0.0",Exp::REAL)); if (*$4 != "2") driver.error("bad power for inv laplacian?"); }



braket:
          LANGLE mline RANGLE                      { $$ = $2; }
        | LANGLE exp RANGLE                        { $$ = new Exp("trace",$2,Exp::brakettype($2->type(),$2->type())); }


normf:
          MMLINE exp MMLINE SUBSCRIPT F            { $$ = new Exp("mathfun", new Exp("normf"), $2, Exp::REAL); }

norm2:
          MMLINE exp MMLINE                        { $$ = new Exp("mathfun", new Exp("norm"), $2, Exp::REAL); }

mline:
          exp MLINE exp                            { $$ = new Exp("braket", $1, $3, Exp::brakettype($1->type(),$3->type())); }

exp:
          INTEGER                                  { $$ = new Exp(*$1,Exp::INTEGER); }
        | STRING                                   { $$ = new Exp(*$1,Exp::STRING); }
        | FLOAT                                    { $$ = new Exp(*$1,Exp::REAL); }
        | variable                                 { $$ = $1; } 
        | braket %prec BRAKET                      { $$ = $1; }
        | exp PLUS exp                             { $$ = new Exp("+",$1,$3); }
        | exp MINUS exp                            { $$ = new Exp("-",$1,$3); }
        | exp TIMES exp                            { $$ = new Exp("*",$1,$3); }
        | exp DIVIDE exp                           { $$ = new Exp("/",$1,$3); }
        | exp POWER exp SUBSCRIPT exp              { $$ = new Exp("^_",$1,$3,$5,Exp::subscrtype($1->type())); }
        | exp SUBSCRIPT exp POWER exp              { $$ = new Exp("^_",$1,$5,$3,Exp::subscrtype($1->type())); }
        | exp POWER exp                            { $$ = new Exp("^",$1,$3,$1->type()); }
        | exp SUBSCRIPT exp                        { $$ = new Exp("_",$1,$3,Exp::subscrtype($1->type())); }
        | exp OR exp                               { $$ = new Exp("or",$1,$3); }
        | exp AND exp                              { $$ = new Exp("and",$1,$3); }
        | exp LANGLE exp %prec LT                  { $$ = new Exp("lt",$1,$3); }
        | exp RANGLE exp %prec GT                  { $$ = new Exp("gt",$1,$3); }
        | exp LE exp                               { $$ = new Exp("le",$1,$3); }
        | exp GE exp                               { $$ = new Exp("ge",$1,$3); }
        | exp EQ exp                               { $$ = new Exp("eq",$1,$3); }
        | exp NE exp                               { $$ = new Exp("ne",$1,$3); }
        | NOT exp                                  { $$ = new Exp("not",$2); }
        | MINUS exp                                { $$ = new Exp("-",$2,$2->type()); }
        | LPAREN exp RPAREN                        { $$ = new Exp("()",$2,$2->type()); }
        | F            variable                    { $$ = new Exp("F", $2, Exp::funtype($2->type())); }
        | F     LPAREN variable RPAREN             { $$ = new Exp("F", $3, Exp::funtype($3->type())); }
        |       MLINE  exp MLINE                   { $$ = new Exp("mathfun", new Exp("abs"), $2, Exp::REAL); }
        | normf                                    { $$ = $1; }
        | norm2                                    { $$ = $1; }
        | FUNCTION LPAREN exp RPAREN               { $$ = new Exp("mathfun", new Exp(*$1), $3); }
        | variable LPAREN arglist RPAREN           { $$ = new Exp("call", $1, $3, $1->type()); }
        | derivative exp %prec DERIVATIVE          { $$ = $1; $$->set_type($2->type()); $$->add($2); }
        | summation  exp %prec SUMMATION           { $$ = $1; $$->set_type($2->type()); $$->add($2); }
        | laplacian  exp %prec SUMMATION           { $$ = $1; $$->set_type($2->type()); $$->add($2); }
        | bsh        exp %prec SUMMATION           { $$ = $1; $$->set_type($2->type()); $$->add($2); }

variable:
          NAME                          { $$ = new Exp(*$1, driver.lookup_type(*$1)); } 
        | F                             { $$ = new Exp("F", driver.lookup_type("F")); } 
        | FIELD                         { $$ = new Exp(*$1, driver.lookup_type(*$1)); } 
        | GREEK                         { $$ = new Exp("greek", new Exp(*$1), driver.lookup_type(*$1)); }

/* comma separated list of variables used in type declaration */
varlist:
          variable                      { $$ = new Exp("varlist", $1); }
        | varlist COMMA variable        { $$ = $1; $1->add($3); }

/* comma separated list of variables used to define function parameters */
paramlist:
          variable                      { $$ = new Exp("paramlist", $1); }
        | paramlist COMMA variable      { $$ = $1; $1->add($3); }

/* comma separated list of expressions used to pass arguments to functions */
arglist:
        exp                             { $$ = new Exp("arglist", $1); }
        | arglist COMMA exp             { $$ = $1; $1->add($3); }

%%

void yy::mra_parser::error (const yy::mra_parser::location_type& l,
                            const std::string& m)
{
  driver.error (l, m);
}
