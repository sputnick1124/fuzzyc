AM_CFLAGS = -O4 -Wall -fopenmp
AM_LDFLAGS = -lm

ga.o: ga.c fuzzy.o
	gcc $CFLAGS -o ga.o fuzzy.o -c ga.c $LDFLAGS

fuzzy.o: fuzzy.c
	gcc $CFLAGS -o fuzzy.o -c fuzzy.c $LDFLAGS

