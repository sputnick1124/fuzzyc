#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "../../include/fuzzy.h"
#include "../../include/ga.h"
#include "../../include/gnuplot_i.h"

static const int nsteps = 10000;
static const double t1 = 20.0;

double stepinfo(int nsteps, double t[nsteps], double x[nsteps]) {
	double max = 0;
	double ts = 0;
	int i;
	for (i = 0; i < nsteps; i++) {
		max = (x[i] > max ? x[i] : max);
		ts = (fabs(x[i] - 1) > 0.02 ? t[i] : ts);
	}
	double os = (max - 1) * 100;
	if (max > 1000) {return 101.0;}
	return ts;
}

int f4_dyn(double t, const double y[], double f[],
					void *params) {
	(void)(t);
	double * u = *(double **)params;
	f[0] = -10.875975715524719 * y[0] -
			 9.371205550737207 * y[1] -
			 6.097137901127494 * y[2] -
			 0.273373807458803 * y[3] -
			 0.138811795316565 * y[4] + u[0];
	f[1] = y[0];
	f[2] = y[1];
	f[3] = y[2];
	f[4] = y[3];
	return GSL_SUCCESS;
}

double f4_fitness(struct Fis * fis) {
	double * u = malloc(sizeof(double));;

	gsl_odeiv2_system f4 = {f4_dyn, NULL, 5, &u};


	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new(&f4, gsl_odeiv2_step_rk8pd,
									1e-3, 1e-3, 0.0);

	int i;
	double t = 0.0;
	double y[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
	double tt[nsteps], y0[nsteps];
	double kp, ki, kd;
	double ep[nsteps], ei[nsteps], ed[nsteps];
	double e[3], out[3];
	ep[0] = 0; ei[0] = 0; ed[0] = 0;

	tt[0] = 0; y0[0] = 0;
	u[0] = 0;
	for (i = 1; i < nsteps; i++) {
		double ti = i * t1 / nsteps;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

		if (status != GSL_SUCCESS) {
			printf("error, return value=%d\n",status);
			break;
		}

		tt[i] = t;
		y0[i] = 14.575021682567217 * y[2] +
				 5.884648742411102 * y[3] +
				 0.443191673894189 * y[4];
		ep[i] = (1 - y0[i]);
		ei[i] = ei[i - 1] + (t1 / nsteps) * ep[i];
		ed[i] = (ep[i] - ep[i - 1]) / (t1 / nsteps);
		e[0] = (ep[i] + 10) / 10.0;
		e[1] = (ei[i] + 10) / 10.0;
		e[2] = (ed[i] + 10) / 10.0;
		evalfis(out,e,fis);
//		u[0] = (out[0] - 0.5) * 10; //If single output FIS directly controls force
		kp = out[0] * 20; ki = out[1] * 5; kd = out[2] * 10; //If fuzzy-PID
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	return stepinfo(nsteps, tt, y0);
}

double PID_fitness(double kp, double ki, double kd, gnuplot_ctrl * h) {
	double * u = malloc(sizeof(double));;

	gsl_odeiv2_system f4 = {f4_dyn, NULL, 5, &u};


	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new(&f4, gsl_odeiv2_step_rk8pd,
									1e-3, 1e-3, 0.0);

	int i;
	double t = 0.0;
	double y[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
	double tt[nsteps], y0[nsteps];
	double ep[nsteps], ei[nsteps], ed[nsteps];
	ep[0] = 1; ei[0] = 0; ed[0] = 0;
	u[0] = 0;
	tt[0] = 0; y0[0] = 0;
	u[0] = 0;

	for (i = 1; i < nsteps; i++) {
		double ti = i * t1 / nsteps;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

		if (status != GSL_SUCCESS) {
			printf("error, return value=%d\n",status);
			break;
		}

		tt[i] = t;
		y0[i] = 14.575021682567217 * y[2] +
				 5.884648742411102 * y[3] +
				 0.443191673894189 * y[4];
		ep[i] = (1 - y0[i]);
		ei[i] = ei[i - 1] + (t1 / nsteps) * ep[i];
		ed[i] = (ep[i] - ep[i - 1]) / (t1 / nsteps);
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	gnuplot_plot_xy(h,tt,y0,nsteps,"PID Response");
	double cost = stepinfo(nsteps, tt, y0);
	printf("%f\n",cost);
	return cost;
}

void f4_fis_plot(struct Fis * fis) {
	double * u = malloc(sizeof(double));;

	gsl_odeiv2_system f4 = {f4_dyn, NULL, 5, &u};


	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new(&f4, gsl_odeiv2_step_rk8pd,
									1e-3, 1e-3, 0.0);

	int i;
	double t = 0.0;
	double y[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
	double tt[nsteps], y0[nsteps];
	double kp, ki, kd;
	double ep[nsteps], ei[nsteps], ed[nsteps];
	double e[3], out[3];
	ep[0] = 1; ei[0] = 0; ed[0] = 0;

	tt[0] = 0; y0[0] = 0;
	u[0] = 0;

	for (i = 1; i < nsteps; i++) {
		double ti = i * t1 / nsteps;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

		if (status != GSL_SUCCESS) {
			printf("error, return value=%d\n",status);
			break;
		}

		tt[i] = t;
		y0[i] = 14.575021682567217 * y[2] +
				 5.884648742411102 * y[3] +
				 0.443191673894189 * y[4];
		ep[i] = (1 - y0[i]);
		ei[i] = ei[i - 1] + (t1 / nsteps) * ep[i];
		ed[i] = (ep[i] - ep[i - 1]) / (t1 / nsteps);
		e[0] = (ep[i] + 10) / 10.0;
		e[1] = (ei[i] + 10) / 10.0;
		e[2] = (ed[i] + 10) / 10.0;
		evalfis(out,e,fis);
//		u[0] = (out[0] - 0.5) * 10; //If single output FIS directly controls force
		kp = out[0] * 20; ki = out[1] * 5; kd = out[2] * 10; //If fuzzy-PID
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	gnuplot_ctrl * h1;
	h1 = gnuplot_init();
	gnuplot_setstyle(h1,"lines");
	gnuplot_plot_xy(h1, tt, y0, nsteps, "FIS-controlled response");
	PID_fitness(5.8, 0.47, 2.85,h1);
	getchar();
	gnuplot_close(h1);
	printf("min(ep) = %f\nmax(ep) = %f\n",minimum(nsteps,ep),maximum(nsteps,ep));
	printf("min(ei) = %f\nmax(ei) = %f\n",minimum(nsteps,ei),maximum(nsteps,ei));
	printf("min(ed) = %f\nmax(ed) = %f\n",minimum(nsteps,ed),maximum(nsteps,ed));
}


int main(int argc, char *argv[]) {
	srand(time(NULL));
	srand48(rand());
	int num_in = 3;
	int in_mfs[3] = {5,5,5};
//	int num_out = 1;
//	int out_mfs[1] = {5};
	int num_out = 3;
	int out_mfs[3] = {5,5,5};

	struct Specs * spcs = specs_set(num_in, in_mfs, num_out, out_mfs);
	struct HyperParams * hp = malloc(sizeof(struct HyperParams));

	hp->pop_size = 100;
	hp->elite = 0.05;
	hp->crossover = 0.5;
	hp->mutate = 0.25;
	hp->max_gen = 200;

	struct Fis * bestfis = run_ga(spcs, hp, f4_fitness);

	f4_fis_plot(bestfis);

	fis_destroy(bestfis);

//	printf("%f\n",PID_fitness(5.8, 0.47, 2.85));

	specs_clear(spcs);
	free(hp);

	return 0;
}
