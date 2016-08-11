CFLAGS = -g

bin_PROGRAMS = simple_example2 ga_example spring_mass_ode_control
simple_example2_SOURCES = simple_example2.c ../fuzzy.c
ga_example_SOURCES = ga_example.c ../fuzzy.c ../ga.c ../gnuplot_i.c
spring_mass_ode_control_SOURCES = spring_mass_ode_control.c ../gnuplot_i.c ../fuzzy.c ../ga.c

#spring_mass_ode_control.o: spring_mass_ode_control.c ../fuzzy.c ../ga.c
#	gcc -O4 -Wall -o spring_mass_ode_control.o -c spring_mass_ode_control.c 

#spring_mass_ode_control: spring_mass_ode_control.o ../gnuplot_i.o
#	gcc -o spring_mass_ode_control spring_mass_ode_control.o ../gnuplot_i.o ../fuzzy.o ../ga.o -lgsl -lcblas -latlas -lm
spring_mass_ode_control_LDADD = -lgsl -lgslcblas -latlas -lm

