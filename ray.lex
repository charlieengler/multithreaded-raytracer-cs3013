%{

#include <stdio.h>
#include "ray_ast.h"
#include "ray.yacc.generated_h"

%}

%option reentrant
%option bison-bridge
%option bison-locations
%option yylineno

%option header-file="ray.lex.generated_h"

%%

[ \t\n]+			{ ; }
[-+]?([0-9]+)(\.[0-9]+)?(e-?[0-9]+)?	{ sscanf(yytext, "%lf", &yylval->dblval); return FLOAT; }
sphere				{ return SPHERE; }
plane				{ return PLANE; }
\{				{ return LBRACE; }
\}				{ return RBRACE; }
pos				{ return POS; }
radius				{ return RADIUS; }
normal				{ return NORMAL; }
rgba				{ return RGBA; }
reflectance			{ return REFLECTANCE; }
color				{ return COLOR; }
light				{ return LIGHT; }
velocity			{ return VELOCITY; }
<<EOF>>				{ return YYEOF; }
.		{ fprintf(stderr, "bad input character '%s' at line %d\n", yytext, yylineno); return YYEOF; }

%%
