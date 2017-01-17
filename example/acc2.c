#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "fuzzy.h"
#include "ga.h"
#include "gnuplot_i.h"

const int nsteps = 2000;

int acc_dyn(double t, const double y[], double f[],
            void *params) {
    (void)(t);
    double * u = *(double **)params;
    const double om = 0.5;
    f[0] = y[1];
    f[1] = -2*om*u[0]*y[1] - om*om*y[0];
    return GSL_SUCCESS;
}

double stepinfo(int nsteps, double t[nsteps], double x[nsteps]) {
    double ts = 0;
    int i;
    for (i = 0; i < nsteps; i++) {
        ts = fabs(x[i]) >= 0.003 ? t[i] : ts;
    }
    return ts;
}

void run_dynamics(struct Fis * fis, double tt[nsteps], double y0[nsteps],
                    double y1[nsteps], double ic[2]) {
    double * u = malloc(sizeof(double));

    gsl_odeiv2_system acc = {acc_dyn, NULL, 2, &u};

    gsl_odeiv2_driver * d =
        gsl_odeiv2_driver_alloc_y_new(&acc, gsl_odeiv2_step_rk8pd,
                                    1e-7, 1e-7, 0.0);

    int i;
    double t = 0.0;
//    double tt[nsteps], y0[nsteps];
    double y[2];
    y[0] = ic[0];
    y[1] = ic[1];
//    double y[2] = {0.0, 0.2];
    y0[0] = y[0];
    y1[0] = y[1];

    double x[2];
    tt[0] = 0.0;
    int pi;
    for (pi = 0; pi < 2; pi++) {
        printf("y[%d]: %0.2f\t",pi, y[pi]);
    }
    printf("\n");
    for(i = 1; 1 < nsteps; i++) {
        double ti = i * 100.0 / nsteps;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n",status);
            break;
        }

    for (pi = 0; pi < 2; pi++) {
        printf("y[%d]: %0.4f\t",pi, y[pi]);
    }
    printf("\n");
        tt[i] = t;
        y0[i] = y[0];
        y1[i] = y[1];
        x[0] = (y[0] + 1) / 2.0;
        x[1] = (y[1] + 1) / 2.0;
        evalfis(u,x,fis);
        u[0] = u[0]*2;
    }
    free(u);
    gsl_odeiv2_driver_free(d);

}

double acc_fitness(struct Fis * fis) {
    double tt[nsteps], y0[nsteps], y1[nsteps];
    double ic[2] = {0.2, 0.0};
//    double ic[2] = {0.0, 0.2};
    run_dynamics(fis,tt,y0,y1,ic);
    return stepinfo(nsteps,tt,y0);
}


int main(int argc, char *argv[]) {
    srand((long long int)time(NULL));
    srand48(rand());
    int num_in = 2;
    int in_mfs[2] = {3,3};
    int num_out = 1;
    int out_mfs[1] = {5};

    struct Specs * spcs = specs_set(num_in, in_mfs, num_out, out_mfs);
    struct HyperParams * hp = malloc(sizeof(struct HyperParams));

    hp->pop_size = 100;
    hp->elite = 0.05;
    hp->crossover = 0.5;
    hp->mutate = 0.25;
    hp->max_gen = 200;

    FILE * fd = fopen("acc.fis","w");
    struct Fis * bestfis = run_ga(spcs, hp, acc_fitness, fd);
    fclose(fd);

    gnuplot_ctrl * h1;
    h1 = gnuplot_init();
    gnuplot_setstyle(h1,"lines");
    double tt[nsteps], y0[nsteps], y1[nsteps];
    double ic[2] = {0.2,0.0};

    run_dynamics(bestfis, tt, y0, y1, ic);

    gnuplot_plot_xy(h1, tt, y0, nsteps, "Displacement");
    gnuplot_plot_xy(h1, tt, y1, nsteps, "Velocity");
    getchar();
    gnuplot_close(h1);
}
