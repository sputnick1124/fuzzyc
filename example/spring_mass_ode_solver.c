#include <stdio.h>
#include <unistd.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "gnuplot_i.h"

const double om = 0.15;

int spring_mass_dyn(double t, const double y[], double f[],
					void *params) {
	(void)(t);
	f[0] = y[1];
	f[1] = -2*om*y[1] - om*om*y[0];
	return GSL_SUCCESS;
}

int main(int argc, char *argv[]) {
	gnuplot_ctrl * h1;
	h1 = gnuplot_init();
	gnuplot_setstyle(h1, "points");

	gnuplot_ctrl * h2;
	h2 = gnuplot_init();
	gnuplot_setstyle(h2, "points");

	gsl_odeiv2_system sp_mass = {spring_mass_dyn, NULL, 2, NULL};

	gsl_odeiv2_driver * d =
	  gsl_odeiv2_driver_alloc_y_new (&sp_mass, gsl_odeiv2_step_rk8pd,
									1e-6, 1e-6, 0.0);

	int i;
	double t = 0.0, t1 = 100.0;
	double y[2] = {0.2, 0.0};
	double tt[100000], y0[100000], y1[100000];

	for (i = 1; i < 100000; i++) {
		double ti = i * t1 / 100000.0;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

		if (status != GSL_SUCCESS) {
			printf("error, return value=%d\n",status);
			break;
		}
		tt[i] = t;
		y0[i] = y[0];
		y1[i] = y[1];
	}

	gnuplot_plot_xy(h1, tt, y0, 100000, "Mass position");

	gnuplot_plot_xy(h2, tt, y1, 100000, "Mass velocity");

	printf("Press any key to exit program\n");
	getchar();

	gnuplot_close(h1);
	gnuplot_close(h2);

	gsl_odeiv2_driver_free (d);
	return 0;
}
