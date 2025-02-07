
CFLAGS += -Wpedantic -pedantic-errors
CFLAGS += -Werror
CFLAGS += -Wall
CFLAGS += -Wextra
#CFLAGS += -Waggregate-return
CFLAGS += -Wbad-function-cast
CFLAGS += -Wcast-align
CFLAGS += -Wno-cast-qual	# free() should accept const pointers
CFLAGS += -Wno-declaration-after-statement
CFLAGS += -Wfloat-equal
CFLAGS += -Wformat=2
CFLAGS += -Wlogical-op
CFLAGS += -Wmissing-include-dirs
CFLAGS += -Wno-missing-declarations
CFLAGS += -Wno-missing-prototypes
CFLAGS += -Wnested-externs
CFLAGS += -Wpointer-arith
CFLAGS += -Wredundant-decls
CFLAGS += -Wsequence-point
CFLAGS += -Wshadow
CFLAGS += -Wno-strict-prototypes
CFLAGS += -Wswitch
CFLAGS += -Wundef
CFLAGS += -Wunreachable-code
CFLAGS += -Wunused-but-set-parameter
CFLAGS += -Wno-unused-parameter
CFLAGS += -Wno-maybe-uninitialized
CFLAGS += -Wwrite-strings

BINARIES += ray

all: $(BINARIES) $(EXPECTED_SH)

OPT = -O3

ray: ray.yacc.generated.o ray.lex.generated.o ray.o ray_console.o ray_ast.o ray_math.o ray_render.o ray_bmp.o ray_physics.o
	gcc -g $(OPT) $^ -lpthread -lm -o $@

%.generated.o: %.generated_c
	gcc -g $(OPT) -x c $< -DYYDEBUG=1 -c -o $@ -MD -MF $(@:.o=.d)

%.o: %.c Makefile
	gcc -g $(OPT) -x c $(CFLAGS) $< -DYYDEBUG=1 -c -o $@ -MD -MF $(@:.o=.d)

%.yacc.generated_c: %.yacc Makefile
	bison -Wconflicts-sr -Wcounterexamples --locations --language=c --header=$$(echo $@ | sed 's/c$$/h/') -o $@ $<

%.lex.generated_c: %.lex Makefile
	flex -o $@ $<

# Force lex/yacc (flex/bison) runs before regular source compilation since we depend on generated headers.
ray.o: ray.lex.generated_c ray.yacc.generated_c
ray.yacc.generated_c: ray.lex.generated_c

project2.zip: FORCE
	rm -rf $@ project2/ && mkdir project2/
	cp *.c *.h *.lex *.yacc Makefile setup.sh scene*.txt project2/
	zip -r $@ project2/

clean:
	rm -f *.o *.d *.generated[_.][chdo] project2.zip project2_starter.zip $(BINARIES) expected_section?.txt

.PHONY: all clean submission_zip expected FORCE

-include *.d

