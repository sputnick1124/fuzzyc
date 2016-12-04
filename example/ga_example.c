#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ga.h"
#include "fuzzy.h"
#include "gnuplot_i.h"

const int pi = 3.14159;
const int scale = 2;
const int shift = -1;

double
r_squared(
	int n,
	double y[n],
	double y_a[n])
{
	int i;
	double SS_t = 0, SS_r = 0, y_m;
	y_m = sum_d(y,n) / (double)n;
	for (i = 0; i < n; i++) {
		SS_t += pow(y[i] - y_m, 2);
		SS_r += pow(y[i] - y_a[i], 2);
	}
	return SS_r / SS_t;
}

double
fit_line(struct Fis * fis)
{
    int i;
    double dx = 0.01;
    int max = (int) (1.0 / dx);
	double y[max], y_a[max];
    double x[1];
    double out[1];
    for (i = 0; i < max; i++) {
        x[0] = (double)i * dx;
//        x[1] = (double)i * dx;
        evalfis(out,x,fis);
		y[i] = (double)i * dx;
		y_a[i] = out[0];
    }
    return r_squared(max,y,y_a);
}

double
fit_sine(struct Fis * fis)
{
    int i;
    double dx = 0.01;
    int max = (int) (1.0 / dx);
	double y[max], y_a[max];
    double x[1];
    double out[1];
    for (i = 0; i < max; i++) {
        x[0] = (double)i * dx;
//        x[1] = (double)i * dx;
        evalfis(out,x,fis);
		y[i] = sin(x[0] * 4 * pi);
		y_a[i] = scale * out[0] + shift;
    }
    return r_squared(max,y,y_a);
}

void
plot_line(struct Fis * fis)
{
    int i;
    double dx = 0.01;
    int max = (int) (1.0 / dx);
    double x[1];
    double y_a[max];
    double out[1];
	double y[max], x_i[max];

	gnuplot_ctrl * h1;
	h1 = gnuplot_init();
	gnuplot_setstyle(h1,"lines");
    for (i = 0; i < max; i++) {
        x[0] = (double)i * dx;
//        x[1] = (double)i * dx;
		x_i[i] = x[0];
        evalfis(out,x,fis);
		y[i] = scale * out[0] + shift;
		y_a[i] = sin(x_i[i] * 4 * pi);
//		y_a[i] = x_i[i];
    }
	gnuplot_plot_xy(h1, x_i, y, max, "Fuzzy Output");
	gnuplot_plot_xy(h1, x_i, y_a, max, "Expected Output");
//    gnuplot_cmd(h1, "set terminal postscript color eps");
//    gnuplot_cmd(h1, "set output \"line.eps\"");
//    gnuplot_cmd(h1, "replot");
	printf("Press any key to close window\n");
	getchar();
	gnuplot_close(h1);
}


int
main(void)
{
	srand((long int)time(NULL));
	srand48(rand());
	int num_in = 1;
	int num_out = 1;
	int in_mfs[1] = {48};
	int out_mfs[1] = {48};
	struct Specs * spcs = specs_set(num_in, in_mfs, num_out, out_mfs);
	struct HyperParams * hp = malloc(sizeof(struct HyperParams));

	hp->pop_size = 100;
	hp->elite = 0.05;
	hp->crossover = 0.5;
	hp->mutate = 0.25;
	hp->max_gen = 10000;

	struct Fis * bestfis = run_ga(spcs, hp, fit_sine, NULL);

	int r;
	for (r = 0; r < spcs->num_rule; r++) {
		printf("rule output mf: %f, %f, %f\n",
				bestfis->rule_list[r]->output[0][0],
				bestfis->rule_list[r]->output[0][1],
				bestfis->rule_list[r]->output[0][2]);
	}
	plot_line(bestfis);
	free(hp);


	specs_clear(spcs);

//	double x[1], out[1];
/*	while (1) {
		scanf("%f\n",&x[0]);
		evalfis(out,x,bestfis);
		printf("result = %f\n",out[0]);
	}*/

	fis_destroy(bestfis);
	return 0;
}
