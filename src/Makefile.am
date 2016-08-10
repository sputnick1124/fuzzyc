SUBDIRS = example test
CFLAGS = -O4 -Wall -I../include -DDEBUG

#fuzzy.o: fuzzy.c
#	gcc $(CFLAGS) -o $@ -c $+

#ga.o: ga.c fuzzy.o
#	gcc $(CFLAGS) -o $@ -c $+

fuzzy_SOURCES = fuzzy.c
ga_SOURCES =  fuzzy.c ga.c
gnuplot_i_sources = gnuplot_i.c ../include/gnuplot_i.h
