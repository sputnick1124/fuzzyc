AM_CFLAGS = -fopenmp -g -Wall

bin_PROGRAMS = simple_example2 ga_example spring_mass_ode_control f4_dynamics f4_control
simple_example2_SOURCES = simple_example2.c ../fuzzy.c
ga_example_SOURCES = ga_example.c ../fuzzy.c ../ga.c ../gnuplot_i.c
spring_mass_ode_control_SOURCES = spring_mass_ode_control.c ../gnuplot_i.c ../fuzzy.c ../ga.c
f4_dynamics_SOURCES = f4_dynamics.c ../gnuplot_i.c ../fuzzy.c ../ga.c
f4_control_SOURCES = f4_control.c ../gnuplot_i.c ../fuzzy.c ../ga.c

spring_mass_ode_control_LDADD = -lgsl -lgslcblas -latlas
f4_dynamics_LDADD = -lgsl -lgslcblas -latlas
f4_control_LDADD = -lgsl -lgslcblas -latlas

