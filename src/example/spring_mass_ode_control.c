#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "../../include/fuzzy.h"
#include "../../include/ga.h"
#include "../../include/gnuplot_i.h"

const double om = 0.15;


int spring_mass_dyn(double t, const double y[], double f[],
					void *params) {
	(void)(t);
	double * u = *(double **)params;
	f[0] = y[1];
	f[1] = -2*u[0]*om*y[1] - om*om*y[0];
	return GSL_SUCCESS;
}

double
spring_fitness(struct Fis * fis)
{
	double * u = malloc(sizeof(double));

	gsl_odeiv2_system sp_mass = {spring_mass_dyn, NULL, 2, &u};

	gsl_odeiv2_driver * d =
	  gsl_odeiv2_driver_alloc_y_new (&sp_mass, gsl_odeiv2_step_rk8pd,
									1e-6, 1e-6, 0.0);

	int i;
	double t = 0.0, t1 = 100.0;
	double y[2] = {0.2, 0.0};
//	double tt[1000], y0[1000], y1[1000];
	double cost = 0;

	for (i = 1; i < 1000; i++) {
		evalfis(u,y,fis);
		double ti = i * t1 / 1000.0;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

		if (status != GSL_SUCCESS) {
			printf("error, return value=%d\n",status);
			break;
		}
//		tt[i] = t;
//		y0[i] = y[0];
//		y1[i] = y[1];
		cost += fabs(y[0]) + fabs(y[1]);
	}
	gsl_odeiv2_driver_free (d);
	return cost;
}

void
plot_fis_control(struct Fis * fis)
{
	gnuplot_ctrl * h1;
	h1 = gnuplot_init();
	gnuplot_setstyle(h1, "lines");

	double * u = malloc(sizeof(double));

	gsl_odeiv2_system sp_mass = {spring_mass_dyn, NULL, 2, &u};

	gsl_odeiv2_driver * d =
	  gsl_odeiv2_driver_alloc_y_new (&sp_mass, gsl_odeiv2_step_rk8pd,
									1e-6, 1e-6, 0.0);

	int i, j;
	double t, t1;
	double y[2];
	double tt[1000], y0[1000], y1[1000];

	for (j = 0; j < 2; j++) {
		y[0] = 0.2; y[1] = 0.0;
		t = 0.0; t1 = 100.0;
		for (i = 1; i < 1000; i++) {
			if (j < 1) {evalfis(u,y,fis);}
			else {u[0] = 2;}
			double ti = i * t1 / 1000.0;
			int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

			if (status != GSL_SUCCESS) {
				printf("error, return value=%d\n",status);
				break;
			}
			tt[i] = t;
			y0[i] = y[0];
			y1[i] = y[1];
		}
		gsl_odeiv2_driver_reset (d);
		gnuplot_plot_xy(h1, tt, y0, 1000, "Mass position");

		gnuplot_plot_xy(h1, tt, y1, 1000, "Mass velocity");

	}
	gsl_odeiv2_driver_free (d);

	printf("Press any key to exit program\n");
	getchar();

	gnuplot_close(h1);

}

int main(int argc, char *argv[]) {
	srand((long int)time(NULL));
	srand48(rand());

	int num_in = 2;
	int in_mfs[2] = {5,5};
	int num_out = 1;
	int out_mfs[1] = {5};
	struct Specs * spcs = specs_set(num_in, in_mfs, num_out, out_mfs);
	struct HyperParams * hp = malloc(sizeof(struct HyperParams));

	hp->pop_size = 100;
	hp->elite = 0.05;
	hp->crossover = 0.5;
	hp->mutate = 0.25;
	hp->max_gen = 200;

	struct Fis * bestfis = run_ga(spcs, hp, spring_fitness);

	plot_fis_control(bestfis);

	return 0;
}
