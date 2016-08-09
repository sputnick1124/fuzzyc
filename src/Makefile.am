SUBDIRS = example test
CFLAGS = -O4 -Wall

fuzzy.o: fuzzy.c
	gcc $(CFLAGS) -o $@ -c $+

ga.o: ga.c fuzzy.o
	gcc $(CFLAGS) -o $@ -c $+
