#acc2.c
#f4_control.c
#f4_dynamics.c
#fuzz_ex.c
#ga_example.c
#py_example.py
#simple_example2.c
#simple_example.c
#simple_fuzzy_example2.c
#simple_fuzzy_example.c
#spring_mass_ode_control.c
#spring_mass_ode_solver.c
CC = gcc
CFLAGS = -I../src/include

libs = -L../lib -lfuzzy -lm

bindir = ../bin

#.PHONY: all
#all: simple_example simple_example2 ga_example

$(bindir)/simple_fuzzy_example: simple_fuzzy_example.c ../lib/libfuzzy.a
	mkdir -p $(bindir)
	$(CC) $(CFLAGS) $< $(libs) -o $@

clean:
	$(rm) $(bindir)/simple_fuzzy_example
