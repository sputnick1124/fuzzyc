bin_PROGRAMS = simple_example2 ga_example #spring_mass_ode_solver
simple_example2_SOURCES = simple_example2.c ../fuzzy.c
ga_example_SOURCES = ga_example.c ../fuzzy.c ../ga.c ../gnuplot_i.c
spring_mass_ode_solver_SOURCES = spring_mass_ode_solver.c ../gnuplot_i.c

#spring_mass_ode_solver.o: spring_mass_ode_solver.c
#	gcc -O4 -Wall -o spring_mass_ode_solver.o -c spring_mass_ode_solver.c 

#spring_mass_ode_solver: spring_mass_ode_solver.o ../gnuplot_i.o
#	gcc -o spring_mass_ode_solver spring_mass_ode_solver.o ../gnuplot_i.o -lgsl -lcblas -latlas -lm

CFLAGS = -g
