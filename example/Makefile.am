AM_CFLAGS = -O4 -fopenmp -I../src

noinst_PROGRAMS = simple_example2 ga_example spring_mass_ode_control f4_dynamics f4_control


simple_example2_SOURCES = simple_example2.c ../src/fuzzy.c

ga_example_SOURCES = ga_example.c ../src/fuzzy.c ../src/ga.c ../src/gnuplot_i.c

spring_mass_ode_control_SOURCES = spring_mass_ode_control.c ../src/gnuplot_i.c ../src/fuzzy.c ../src/ga.c

f4_dynamics_SOURCES = f4_dynamics.c ../src/gnuplot_i.c ../src/fuzzy.c ../src/ga.c

f4_control_SOURCES = f4_control.c ../src/gnuplot_i.c ../src/fuzzy.c ../src/ga.c

spring_mass_ode_control_LDADD =  -lgsl -lgslcblas -latlas
f4_dynamics_LDADD =  -lgsl -lgslcblas -latlas
f4_control_LDADD =  -lgsl -lgslcblas -latlas

