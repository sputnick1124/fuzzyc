AM_CFLAGS = -O4 -Wall -fopenmp
AM_LDFLAGS = -lm

#lib_LIBRARIES = libfuzzy.a libga.a

#libfuzzy_a_SOURCES = fuzzy.c

#libga_a_SOURCES = fuzzy.c ga.c



ga.o: ga.c fuzzy.o
	gcc $CFLAGS -o ga.o fuzzy.o -c ga.c $LDFLAGS

fuzzy.o: fuzzy.c
	gcc $CFLAGS -o fuzzy.o -c fuzzy.c $LDFLAGS

