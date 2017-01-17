#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "fuzzy.h"
#include "ga.h"
#include "gnuplot_i.h"

const double om = 0.5;
const int nsteps = 2500;
const double t1 = 30.0;

int spring_mass_dyn(double t, const double y[], double f[],
					void *params) {
	(void)(t);
	double * p = *(double **)params;
    double u = p[0];
    double om = p[1];
	f[0] = y[1];
	f[1] = -2*u*om*y[1] - om*om*y[0];
	return GSL_SUCCESS;
}

double
spring_fitness(struct Fis * fis)
{
	double * p = malloc(sizeof(double) * 2);

	gsl_odeiv2_system sp_mass = {spring_mass_dyn, NULL, 2, &p};

	gsl_odeiv2_driver * d =
	  gsl_odeiv2_driver_alloc_y_new (&sp_mass, gsl_odeiv2_step_rk8pd,
									1e-6, 1e-6, 0.0);

	int i;
	double t = 0.0;
	double y[2] = {0.2, 0.0};
    double ts = 0.0;
    double u[1];
//	double tt[1000], y0[1000], y1[1000];
	double cost = 0;
    p[1] = 0.5;

	for (i = 1; i < nsteps; i++) {
		evalfis(u,y,fis);
        p[0] = u[0];
		double ti = i * t1 / nsteps;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

		if (status != GSL_SUCCESS) {
			printf("error, return value=%d\n",status);
			break;
		}
//		tt[i] = t;
//		y0[i] = y[0];
//		y1[i] = y[1];
//		cost += fabs(y[0]) + fabs(y[1]);
        ts = (fabs(y[0]) > 0.003 || fabs(y[1]) > 0.003) ? t : ts;
	}
	gsl_odeiv2_driver_free (d);
	free(p);
	return ts;
}

void
plot_fis_control(struct Fis * fis, double ic[2], double om,
                gnuplot_ctrl * h1, int sra, double *tss, FILE * logfile)
{
	double * p = malloc(sizeof(double) * 2);

	gsl_odeiv2_system sp_mass = {spring_mass_dyn, NULL, 2, &p};

	gsl_odeiv2_driver * d =
	  gsl_odeiv2_driver_alloc_y_new (&sp_mass, gsl_odeiv2_step_rk8pd,
									1e-6, 1e-6, 0.0);

	int i, j;
	double t;
	double y[2];
    double u[1];
	double tt[nsteps*10], y0[nsteps*10], y1[nsteps*10];
    double ut[nsteps*10];
    double ts = 0;
    p[1] = om;

	tt[0] = 0; y0[0] = 0; y1[0] = 0;
	for (j = 0; j < 2; j++) {
		y[0] = ic[0]; y[1] = ic[1];
		t = 0.0;
        int stop = 0;
//		for (i = 1; stop < 1; i++) {
		for (i = 1; i < nsteps * 10; i++) {
			if (j < 1) {evalfis(u,y,fis); p[0] = u[0];}
			else if (sra > 0) {p[0] = 2*om;}
            else {p[0] = 2;}
			double ti = i * t1 / nsteps;
			int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

			if (status != GSL_SUCCESS) {
				printf("error, return value=%d\n",status);
				break;
			}
			tt[i] = t;
			y0[i] = y[0];
			y1[i] = y[1];
            ut[i] = u[0];
            ts = ((fabs(y[0]) > 0.003 || fabs(y[1]) > 0.003)) ? t : ts;
            if (sra > 0 && j < 1 && i < nsteps) {fprintf(logfile,"%f\t%f\t%f\n",t,y[0],y[1]);}
		}
        fprintf(logfile,"\n");
        tss[j] = ts;
        ts = 0.0;
		gsl_odeiv2_driver_reset (d);
        if (j < 1){
            if (sra > 0 && i < nsteps) {}
            else {
                if (sra > 0) {}
                else {
            		gnuplot_plot_xy(h1, tt, y0, nsteps, "Fuzzy Position");

            		gnuplot_plot_xy(h1, tt, y1, nsteps, "Fuzzy Velocity");
                }
            }
        } else {
            if (sra > 0) {}
            else {
        		gnuplot_plot_xy(h1, tt, y0, nsteps, "Critical position");

        		gnuplot_plot_xy(h1, tt, y1, nsteps, "Critical velocity");
            }
        }
	}
	gsl_odeiv2_driver_free (d);

	free(p);
//    printf("om:%0.3f, ts1: %0.2f, ts2: %0.2f\n",om,tss[0],tss[1]);
}

void
log_fis_control(struct Fis * fis, double ic[2], double om,
               FILE * logfile, int sra, double *tss)
{
	double * p = malloc(sizeof(double) * 2);

	gsl_odeiv2_system sp_mass = {spring_mass_dyn, NULL, 2, &p};

	gsl_odeiv2_driver * d =
	  gsl_odeiv2_driver_alloc_y_new (&sp_mass, gsl_odeiv2_step_rk8pd,
									1e-6, 1e-6, 0.0);

	int i;
	double t;
	double y[2];
    double u[1];
	double tt[nsteps*10], y0[nsteps*10], y1[nsteps*10];
    double ut[nsteps*10];
    double ts = 0.0;
    p[1] = om;

	tt[0] = 0; y0[0] = 0; y1[0] = 0;
	y[0] = ic[0]; y[1] = ic[1];
	t = 0.0;
    int stop = 0;
//		for (i = 1; stop < 1; i++) {
	for (i = 1; i < nsteps * 10; i++) {
		evalfis(u,y,fis);
        p[0] = u[0];
		double ti = i * t1 / nsteps;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

		if (status != GSL_SUCCESS) {
			printf("error, return value=%d\n",status);
		}
		tt[i] = t;
		y0[i] = y[0];
	    y1[i] = y[1];
        ut[i] = u[0];
        ts = ((fabs(y[0]) > 0.003) || (fabs(y[1]) > 0.003)) ? t : ts;
        fprintf(logfile,"%f\t%f\t%f\n",t,y[0],y[1]);
	}
    fprintf(logfile,"\n");
    ts = 0.0;
    gsl_odeiv2_driver_reset (d);
	gsl_odeiv2_driver_free (d);


	free(p);
//    printf("om:%0.3f, ts1: %0.2f, ts2: %0.2f\n",om,tss[0],tss[1]);
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

	hp->pop_size = 50;
	hp->elite = 0.05;
	hp->crossover = 0.55;
	hp->mutate = 0.35;
	hp->max_gen = 500;

    FILE * fisfile = fopen("sra.fis","w");
	struct Fis * bestfis = run_ga(spcs, hp, spring_fitness, fisfile);
    fclose(fisfile);

	gnuplot_ctrl * h1;
	h1 = gnuplot_init();
	gnuplot_setstyle(h1, "lines");
    gnuplot_cmd(h1,"set xlabel font 'Verdana,20'");
    gnuplot_cmd(h1,"set xtics font 'Verdana,20'");
    gnuplot_cmd(h1,"set ylabel font 'Verdana,20'");
    gnuplot_cmd(h1,"set ytics font 'Verdana,20'");
    gnuplot_cmd(h1,"set key font 'Verdana,20'");
    gnuplot_set_xlabel(h1,"Time (s)");

	gnuplot_ctrl * h2;
	h2 = gnuplot_init();
	gnuplot_setstyle(h2, "lines");
    gnuplot_cmd(h2,"set xlabel font 'Verdana,20'");
    gnuplot_cmd(h2,"set xtics font 'Verdana,20'");
    gnuplot_cmd(h2,"set ylabel font 'Verdana,20'");
    gnuplot_cmd(h2,"set ytics font 'Verdana,20'");
    gnuplot_cmd(h2,"set key font 'Verdana,20'");
    gnuplot_set_xlabel(h2,"Time (s)");

	gnuplot_ctrl * h3;
	h3 = gnuplot_init();
	gnuplot_setstyle(h3, "lines");
    gnuplot_cmd(h3,"set xlabel font 'Verdana,20'");
    gnuplot_cmd(h3,"set xtics font 'Verdana,20'");
    gnuplot_cmd(h3,"set ylabel font 'Verdana,20'");
    gnuplot_cmd(h3,"set ytics font 'Verdana,20'");
    gnuplot_cmd(h3,"set key font 'Verdana,20'");
    gnuplot_set_xlabel(h3,"Time (s)");

	gnuplot_ctrl * h4;
	h4 = gnuplot_init();
	gnuplot_setstyle(h4, "lines");
    gnuplot_cmd(h4,"set xlabel font 'Verdana,20'");
    gnuplot_cmd(h4,"set xtics font 'Verdana,20'");
    gnuplot_cmd(h4,"set ylabel font 'Verdana,20'");
    gnuplot_cmd(h4,"set ytics font 'Verdana,20'");
    gnuplot_cmd(h4,"set key font 'Verdana,20'");
    gnuplot_set_xlabel(h4,"$\\\\omega$ \\\\si{\\\\radian\\\\per\\\\s}");
    gnuplot_set_ylabel(h4,"Settling Time (s)");


    double ic1[2] = {0.2,0.0};
    double ic2[2] = {0.0,0.2};
    double tss[2];
	plot_fis_control(bestfis,ic1,0.5,h1,0,tss,stdout);
    plot_fis_control(bestfis,ic2,0.5,h2,0,tss,stdout);

    /*Stochastic robustness analysis*/
    FILE * logfile = fopen("sra.txt","w");
    FILE * datafile = fopen("sra_ts.txt","w");
    int i;
    for (i = 0; i < 100; i++) {
        double om = drand48()*0.75 + 0.25;
        plot_fis_control(bestfis,ic1,om,h1,1,tss,logfile);
        fprintf(datafile,"%f\t%f\t%f\n",om,tss[0],tss[1]);
    }
    fclose(datafile);

    gnuplot_cmd(h3,"plot \"sra.txt\" using 1:2 with lines title \"Position\", \"sra.txt\" using 1:3 with lines title \"Velocity\"");
    gnuplot_cmd(h4,"plot \"sra_ts.txt\" using 1:2 title \"Fuzzy Damping\", \"sra_ts.txt\" using 1:3 title \"Critical Damping\"");

	fis_destroy(bestfis);
	specs_clear(spcs);
	free(hp);
    gnuplot_cmd(h1,"set terminal epslatex color");
    gnuplot_cmd(h1,"set output \"ic1.tex\"");
    gnuplot_cmd(h1,"replot");
    gnuplot_cmd(h2,"set terminal epslatex color");
    gnuplot_cmd(h2,"set output \"ic2.tex\"");
    gnuplot_cmd(h2,"replot");
    gnuplot_cmd(h3,"set terminal epslatex color");
    gnuplot_cmd(h3,"set output \"sra.tex\"");
    gnuplot_cmd(h3,"replot");
    gnuplot_cmd(h4,"set terminal epslatex color");
    gnuplot_cmd(h4,"set output \"sra_ts.tex\"");
    gnuplot_cmd(h4,"replot");
	printf("Press any key to exit program\n");
    getchar();
	gnuplot_close(h1);
	gnuplot_close(h2);
	gnuplot_close(h3);
	gnuplot_close(h4);

	return 0;
}
