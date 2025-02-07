%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ray_ast.h"
#include "ray.yacc.generated_h"
#include "ray.lex.generated_h"


%}

%define api.pure full
%define parse.error detailed
%locations
%parse-param { struct context *context }
%parse-param { yyscan_t yyscanner }
%lex-param { yyscan_t yyscanner }

%code provides {

void yyerror (YYLTYPE *y, struct context *context, yyscan_t yyscanner, char const *s);

}

%start scene_graph


%token FLOAT SPHERE PLANE COLOR LBRACE RBRACE RGBA REFLECTANCE POS NORMAL RADIUS LIGHT VELOCITY

%union {
	struct context *context;
	sphere sphere;
	spherespec *spherespec;
	spherespecs *spherespecs;
	plane plane;
	planespec *planespec;
	planespecs *planespecs;
	light light;
	lightspec *lightspec;
	lightspecs *lightspecs;	
	color *color;
	colorspec *colorspec;
	colorspecs *colorspecs;
	pt4 *pt4;
	pt3 *pt3;
	double dblval;
}

%type <context> scene_graph
%type <sphere> sphere
%type <spherespec> spherespec
%type <spherespecs> spherespecs
%type <plane> plane
%type <planespec> planespec
%type <planespecs> planespecs
%type <light> light
%type <lightspec> lightspec
%type <lightspecs> lightspecs
%type <color> color
%type <colorspec> colorspec
%type <colorspecs> colorspecs
%type <pt3> pt3
%type <pt4> pt4
%type <dblval> dblval FLOAT;


%%                   /* beginning of rules section */

scene_graph:					{ }
	|	scene_graph sphere		{ /* added in sphere code */ }
	|	scene_graph plane		{ context_add_plane(context, $2); }
	|	scene_graph light		{ context_add_light(context, $2); }
	;

plane:		PLANE LBRACE planespecs RBRACE	{ memset(&$$, 0, sizeof($$)); for (planespec *ps = $3->first; ps; ps = ps->next) { apply_planespec(&$$, ps); } free_planespecs($3); }
	;

sphere:		SPHERE LBRACE spherespecs RBRACE {
							sphere archetype = {0};
							for (spherespec *ss = $3->first; ss; ss = ss->next) {
								if (ss->position == NULL) {
									apply_spherespec(&archetype, ss);
								}
							}
							for (spherespec *ss = $3->first; ss; ss = ss->next) {
								if (ss->position != NULL) {
									sphere s = archetype;
									apply_spherespec(&s, ss);
									context_add_sphere(context, s);
								}
							}
							free_spherespecs($3);
						}
	;

light:		LIGHT LBRACE lightspecs RBRACE	{ memset(&$$, 0, sizeof($$)); for (lightspec *ls = $3->first; ls; ls = ls->next) { apply_lightspec(&$$, ls); } free_lightspecs($3); }
	;

color:		COLOR LBRACE colorspecs RBRACE	{ $$ = new_color(); for(colorspec *cs = $3->first; cs; cs = cs->next) { apply_colorspec($$, cs); } free_colorspecs($3); }
	;

planespecs:					{ $$ = new_planespecs(); }
	|	planespecs planespec		{ $$ = $1; append_ll($$, $2); }
	;

spherespecs:					{ $$ = new_spherespecs(); }
	|	spherespecs spherespec		{ $$ = $1; append_ll($$, $2); }
	;

lightspecs:					{ $$ = new_lightspecs(); }
	|	lightspecs lightspec		{ $$ = $1; append_ll($$, $2); }
	;

colorspecs:					{ $$ = new_colorspecs(); }
	|	colorspecs colorspec		{ $$ = $1; append_ll($$, $2); }
	;

planespec:	POS pt3				{ $$ = new_planespec(); $$->position = $2; }
	|	NORMAL pt3			{ $$ = new_planespec(); $$->normal = $2; }
	|	color				{ $$ = new_planespec(); $$->color = $1; }
	;

spherespec:	POS pt3				{ $$ = new_spherespec(); $$->position = $2; }
	|	RADIUS dblval			{ $$ = new_spherespec(); $$->radius = $2; }
	|	VELOCITY pt3			{ $$ = new_spherespec(); $$->velocity = $2; }
	|	color				{ $$ = new_spherespec(); $$->color = $1; }
	;

lightspec:	POS pt3				{ $$ = new_lightspec(); $$->position = $2; }
	|	color				{ $$ = new_lightspec(); $$->color = $1; }
	;

colorspec:	RGBA pt4			{ $$ = new_colorspec(); $$->rgba = $2; }
	|	REFLECTANCE dblval		{ $$ = new_colorspec(); $$->reflectance = $2; }
	;

pt4:		LBRACE dblval dblval dblval dblval RBRACE	{ $$ = new_pt4(); $$->v[0] = $2; $$->v[1] = $3; $$->v[2] = $4; $$->v[3] = $5; }
	;

pt3:		LBRACE dblval dblval dblval RBRACE	{ $$ = new_pt3(); $$->v[0] = $2; $$->v[1] = $3; $$->v[2] = $4; }
	;

dblval:		FLOAT				{ $$ = $1; }
	;

%%

void yyerror (YYLTYPE *y, struct context *context, yyscan_t yyscanner, char const *s) {
	fprintf(stderr, "%s at line %d\n", s, yyget_lineno(yyscanner)); 
}

int yywrap(yyscan_t yyscanner)
{
	return 1;
}
