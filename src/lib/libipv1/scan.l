%{
#include <string.h>
#include <tmpl.h>
#include "ip_types.h"
#include "ip_global.h"
#include "scan.gbl"
#include "scan.lcl"
#include "y.tab.h"
int lineno;
%}
string  [A-Za-z0-9$_\.+-]*
qstring \"[^"\n]+\"
%%
{string}		{ yylval.str = (char *)malloc(strlen(yytext)+1);
				  if (!yylval.str) {
				    perror("{string} rule");
				    ip_error("{string} rule: malloc failed");
				    }
				  strcpy(yylval.str,yytext);
                  if (ip_uppercase) cvs_toupper(yylval.str);
				  return(T_STRING);
				  }
{qstring}		{ yylval.str = (char *)malloc(strlen(yytext));
				  if (!yylval.str) {
				    perror("{qstring} rule");
				    ip_error("{qstring} rule: malloc failed");
				    }
				  strcpy(yylval.str,&yytext[1]);
				  yylval.str[strlen(yylval.str)-1] = '\0';
				  return(T_STRING);
				  }
[ \t]+		;
"\n"            lineno++;
%.*$			;
#               { return(0); }
[();=:]			{ return((int) yytext[0]); }
.				{ ip_error("Illegal character"); }
%%

/* Convert a string to uppercase. */
LOCAL_FUNCTION void cvs_toupper(char *s)
{
  for (; *s!='\0'; s++) {
    if (*s>='a' && *s <='z') *s = *s + 'A' - 'a';
    }
  }

/* Show position. */
void showpos(void)
{
  printf("error occurred at line number %d\n",lineno);
  }

