#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "fuzzy.h"
#include "ga.h"
#include "gnuplot_i.h"

//#define FUZZY

static const int nsteps = 2000;
static const double t1 = 20.0;
static const double scalar = 20.0;
static const double abs_err = 1e-10;
static const double rel_err = 1e-10;
static const int nout = 3;

double stepinfo(int nsteps, double t[nsteps], double x[nsteps], int flag) {
	double max = 0;
	double ts = 0;
	double osc = 0.0;
	int i;
	int der = 1; //check for change in derivative (the idea being to penalize oscillations)
	for (i = 0; i < nsteps; i++) {
		max = (x[i] > max ? x[i] : max);
		ts = (fabs(x[i] - 1) > 0.02 ? t[i] : ts);
		if (flag && !flag) {
			printf("t=%0.2f\tx(t) = %0.4f\tfabs(x(t)-1)=%0.4f\tts=%0.2f\n",t[i], x[i], fabs(x[i] - 1.0), ts);
		} else {}
		if (i > 0) {
			if ( (der && (x[i] < x[i - 1])) | (!der && (x[i] > x[i-1]))) {
				der ^= 1;
				osc += 1.0;
//				printf("x[%d] = %f, x[%d] = %f, der = %d, osc = %0.0f\n", i,x[i],i-1,x[i-1],der,osc);
			}
		}
	}
	double os = (max > 1 ? max - 1 : 0);
	if (max > 1000) {return 101.0;}
	if (flag) {printf("ts,os,osc=[%f,%f,%f]\n",ts,os,osc);}
	return (ts + os);
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

int f4_deg_dyn(double t, const double y[], double f[],
					void *params) {
	(void)(t);
	double * u = *(double **)params;
	f[0] = -10.719861231569817 * y[0] -
			 7.528187337380746 * y[1] -
			 3.229835212489159 * y[2] -
			 0.160104076322637 * y[3] -
			 0.069384215091067 * y[4] + u[0];
	f[1] = y[0];
	f[2] = y[1];
	f[3] = y[2];
	f[4] = y[3];
	return GSL_SUCCESS;
}

int f4_sub_dyn(double t, const double y[], double f[],
					void *params) {
	(void)(t);
	double * u = *(double **)params;
	f[0] = -11.262534184138559 * y[0] -
			20.738377392889699 * y[1] -
			81.061987237921599 * y[2] +
			  0.069724247948952 * y[3] +
			  0.130013673655424 * y[4] + u[0];
	f[1] = y[0];
	f[2] = y[1];
	f[3] = y[2];
	f[4] = y[3];
	return GSL_SUCCESS;
}

int f4_sup_dyn(double t, const double y[], double f[],
					void *params) {
	(void)(t);
	double * u = *(double **)params;
	f[0] = -10.5051664753157 * y[0] -
			20.3960964408726 * y[1] -
			153.6739380022962 * y[2] -
			0.7158438576349 * y[3] -
			0.0706084959816 * y[4] + u[0];
	f[1] = y[0];
	f[2] = y[1];
	f[3] = y[2];
	f[4] = y[3];
	return GSL_SUCCESS;
}

double
generic_fitness_comp(struct Fis * fis,
		int (* fun)(double t, const double y[], double dydt[], void * params),
		double C[3],
		gnuplot_ctrl * h)
{
	double * u = malloc(sizeof(double));;
	int nsteps = 100000;
	gsl_odeiv2_system f4 = {fun, NULL, 5, &u};


	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new(&f4, gsl_odeiv2_step_rk8pd,
									abs_err, rel_err, 0.0);

	int i;
	double t = 0.0;
	double y[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double tt[nsteps], y0[nsteps], u_t[nsteps];
	double kp, ki, kd;
	double ep[nsteps], ei[nsteps], ed[nsteps];
	double e[3], out[nout];
	ep[0] = 1; ei[0] = 0; ed[0] = 0;

	tt[0] = 0; y0[0] = 0;
	u[0] = 0;
	u_t[0] = 0;
	for (i = 1; i < nsteps; i++) {
		double ti = i * t1 / nsteps;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

		if (status != GSL_SUCCESS) {
			printf("error, return value=%d\n",status);
			break;
		}

		tt[i] = t;
		y0[i] = C[0] * y[2] +
				 C[1] * y[3] +
				 C[2] * y[4];
		ep[i] = (1 - y0[i]);
		ei[i] = ei[i - 1] + (t1 / nsteps) * ep[i];
		ed[i] = (ep[i] - ep[i - 1]) / (t1 / nsteps);
		e[0] = (ep[i] + scalar / 2.0) / scalar;
		e[1] = (ei[i] + scalar / 2.0) / scalar;
		e[2] = (ed[i] + scalar / 2.0) / scalar;
		evalfis(out,e,fis);
		#ifdef FUZZY
		u[0] = out[0] * 10 - 5; //If single output FIS directly controls force
		#endif
		#ifndef FUZZY
		kp = out[0] * 20; ki = out[1] * 5; kd = out[2] * 10; //If fuzzy-PID
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
		#endif
		u_t[i] = u[0];
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	gnuplot_plot_xy(h, tt, y0, nsteps, "Fuzzy PID");
/*	gnuplot_ctrl * force_plot;
	force_plot = gnuplot_init();
	gnuplot_setstyle(force_plot, "lines");
	gnuplot_plot_xy(force_plot, tt, u_t, nsteps, "Control Force");
	getchar();
	gnuplot_close(force_plot);*/
	return stepinfo(nsteps, tt, y0, 1);
}

double f4_fitness(struct Fis * fis) {
	double * u = malloc(sizeof(double));;

	gsl_odeiv2_system f4 = {f4_dyn, NULL, 5, &u};


	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new(&f4, gsl_odeiv2_step_rk8pd,
									abs_err, rel_err, 0.0);

	int i;
	double t = 0.0;
	double y[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double tt[nsteps], y0[nsteps];
	double kp, ki, kd;
	double ep[nsteps], ei[nsteps], ed[nsteps];
	double e[3], out[nout];
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
	double c3[3] = {113.8331814038286, 58.1700091157703, 0.7103463992707};
		y0[i] = 14.575021682567217 * y[2] +
				 5.884648742411102 * y[3] +
				 0.443191673894189 * y[4];
/*		y0[i] = 113.8331814038286 * y[2] +
				 58.1700091157703 * y[3] +
				 0.7103463992707 * y[4];*/
		ep[i] = (1 - y0[i]);
		ei[i] = ei[i - 1] + (t1 / nsteps) * ep[i];
		ed[i] = (ep[i] - ep[i - 1]) / (t1 / nsteps);
		e[0] = (ep[i] + scalar / 2.0) / scalar;
		e[1] = (ei[i] + scalar / 2.0) / scalar;
		e[2] = (ed[i] + scalar / 2.0) / scalar;
		evalfis(out,e,fis);
		#ifdef FUZZY
		u[0] = out[0] * 10 - 5; //If single output FIS directly controls force
		#endif
		#ifndef FUZZY
		kp = out[0] * 20; ki = out[1] * 5; kd = out[2] * 10; //If fuzzy-PID
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
		#endif
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	return stepinfo(nsteps, tt, y0, 0);
}

double PID_fitness(double kp, double ki, double kd,
		int (* fun)(double t, const double y[], double dydt[], void * params),
		double C[3],
		gnuplot_ctrl * h) {
	double * u = malloc(sizeof(double));;
	int nsteps = 100000;
	gsl_odeiv2_system f4 = {fun, NULL, 5, &u};


	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new(&f4, gsl_odeiv2_step_rk8pd,
									abs_err, rel_err, 0.0);

	int i;
	double t = 0.0;
	double y[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
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
		y0[i] = C[0] * y[2] +
				 C[1] * y[3] +
				 C[2] * y[4];
		ep[i] = (1 - y0[i]);
		ei[i] = ei[i - 1] + (t1 / nsteps) * ep[i];
		ed[i] = (ep[i] - ep[i - 1]) / (t1 / nsteps);
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	gnuplot_plot_xy(h,tt,y0,nsteps,"PID Response");
	double cost = stepinfo(nsteps, tt, y0, 1);
//	printf("%f\n",cost);
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
	double y[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double tt[nsteps], y0[nsteps];
	double kp, ki, kd;
	double ep[nsteps], ei[nsteps], ed[nsteps];
	double e[3], out[nout];
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
		e[0] = (ep[i] + scalar / 2.0) / scalar;
		e[1] = (ei[i] + scalar / 2.0) / scalar;
		e[2] = (ed[i] + scalar / 2.0) / scalar;
		evalfis(out,e,fis);
		#ifdef FUZZY
		u[0] = out[0] * 10 - 5; //If single output FIS directly controls force
		#endif
		#ifndef FUZZY
		kp = out[0] * 20; ki = out[1] * 5; kd = out[2] * 10; //If fuzzy-PID
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
		#endif
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	gnuplot_ctrl * h1;
	h1 = gnuplot_init();
	gnuplot_setstyle(h1,"lines");
	gnuplot_plot_xy(h1, tt, y0, nsteps, "FIS-controlled response");
	printf("FIS data:\n");
	stepinfo(nsteps, tt, y0, 0);
//	printf("PID data:\n");
//	PID_fitness(5.8, 0.47, 2.85,h1);
	getchar();
	gnuplot_close(h1);
//	printf("min(ep) = %f\nmax(ep) = %f\n",minimum(nsteps,ep),maximum(nsteps,ep));
//	printf("min(ei) = %f\nmax(ei) = %f\n",minimum(nsteps,ei),maximum(nsteps,ei));
//	printf("min(ed) = %f\nmax(ed) = %f\n",minimum(nsteps,ed),maximum(nsteps,ed));
}


int main(int argc, char *argv[]) {
	srand((long long int)time(NULL));
	srand48(rand());
	int num_in = 3;
	int in_mfs[3] = {3,3,3};
	#ifndef FUZZY
	int num_out = 3;
	int out_mfs[3] = {3,3,3};
	#endif
	#ifdef FUZZY
	int num_out = 1;
	int out_mfs[1] = {3};
	#endif

	struct Specs * spcs = specs_set(num_in, in_mfs, num_out, out_mfs);
	struct HyperParams * hp = malloc(sizeof(struct HyperParams));

	hp->pop_size = 100;
	hp->elite = 0.05;
	hp->crossover = 0.5;
	hp->mutate = 0.25;
	hp->max_gen = 200;

	FILE * fd = fopen("output.fis", "w");
	struct Fis * bestfis = run_ga(spcs, hp, f4_fitness, fd);
	fclose(fd);

/** Plot the comparison graphs between PID- and Fuzzy-controlled systems**/
	double c1[3] = {14.575021682567217, 5.884648742411102, 0.443191673894189};
	double c2[3] = {14.575021682567217, 5.949696444058977, 0.455333911535126};
	double c3[3] = {113.8331814038286, 58.1700091157703, 0.7103463992707};
	double c4[3] = {74.626865671641795, 12.531572904707232, -0.017962112514351};


//	f4_fis_plot(bestfis);
//	fis_print(bestfis,NULL);

	gnuplot_ctrl * h1;
	gnuplot_ctrl * h2;
	gnuplot_ctrl * h3;
	gnuplot_ctrl * h4;
	h1 = gnuplot_init();
	h2 = gnuplot_init();
	h3 = gnuplot_init();
	h4 = gnuplot_init();
	gnuplot_setstyle(h1,"lines");
	gnuplot_setstyle(h2,"lines");
	gnuplot_setstyle(h3,"lines");
	gnuplot_setstyle(h4,"lines");
	gnuplot_set_xlabel(h1, "Time, s");
	gnuplot_set_xlabel(h2, "Time, s");
	gnuplot_set_xlabel(h3, "Time, s");
	gnuplot_set_xlabel(h4, "Time, s");
	gnuplot_set_ylabel(h1, "Theta, rad");
	gnuplot_set_ylabel(h2, "Theta, rad");
	gnuplot_set_ylabel(h3, "Theta, rad");
	gnuplot_set_ylabel(h4, "Theta, rad");
	printf("F4 dyn\n");
	generic_fitness_comp(bestfis, f4_dyn, c1, h1);
	PID_fitness(5.8, 0.47, 2.85,f4_dyn, c1,h1);
	printf("F4 deg50 dyn\n");
	generic_fitness_comp(bestfis, f4_deg_dyn, c2, h2);
	PID_fitness(5.8, 0.47, 2.85,f4_deg_dyn, c2,h2);
	printf("F4 sub dyn\n");
	generic_fitness_comp(bestfis, f4_sub_dyn, c3, h3);
	PID_fitness(5.8, 0.47, 2.85,f4_sub_dyn, c3,h3);
	printf("F4 sup dyn\n");
	generic_fitness_comp(bestfis, f4_sup_dyn, c4, h4);
	PID_fitness(5.8, 0.47, 2.85,f4_sup_dyn, c4,h4);
	getchar();
	gnuplot_cmd(h1, "set terminal postscript color eps");
	gnuplot_cmd(h2, "set terminal postscript colour eps");
	gnuplot_cmd(h3, "set terminal postscript colour eps");
	gnuplot_cmd(h4, "set terminal postscript colour eps");
	gnuplot_cmd(h1, "set output \"f4_approach.eps\"");
	gnuplot_cmd(h2, "set output \"f4_deg.eps\"");
	gnuplot_cmd(h3, "set output \"f4_sub.eps\"");
	gnuplot_cmd(h4, "set output \"f4_sup.eps\"");
	gnuplot_cmd(h1, "replot");
	gnuplot_cmd(h2, "replot");
	gnuplot_cmd(h3, "replot");
	gnuplot_cmd(h4, "replot");
	gnuplot_close(h1);
	gnuplot_close(h2);
	gnuplot_close(h3);
	gnuplot_close(h4);

	/*Generate lookup tables for MATLAB*/
	#ifndef FUZZY
	double i, j, k;
	int in, out;
	char * filename = "lookup_table";
	char new_fn[100];
	char outinnum[13];
	double x[3], y[3];
	double sc[3] = {20,5,10};

	for (out = 0; out < num_out; out++) {
		for (k = 0.00; k < 1.01; k += 0.01) {
			x[2] = k;
			strcpy(new_fn,filename);
			sprintf(outinnum,"%d-%0.2f.csv",out, k);
			strcat(new_fn, outinnum);
			fd = fopen(new_fn, "w");
//			fd = stdout;
			for (j = 0; j < 1.01; j += 0.01) {
				x[1] = j;
				for (i = 0; i < 1.01; i += 0.01) {
					x[0] = i;
					evalfis(y,x,bestfis);
					fprintf(fd,"%f",sc[out] * y[out]);
					if (i < 1) {fprintf(fd,", ");}
				}
				fprintf(fd,"\n");
			}
			fclose(fd);
		}
	}
	#endif

	fis_destroy(bestfis);
	specs_clear(spcs);
	free(hp);

	return 0;
}

