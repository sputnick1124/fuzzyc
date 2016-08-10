#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "../../include/gnuplot_i.h"
#include "../../include/fuzzy.h"
#include "../../include/ga.h"

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
		printf("%f\t%f\t%f\t%f\n", y[i]-y_m, y[i]-y_a[i], y[i], y_a[i]);
    }
    return 1 - SS_r / SS_t;
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
    double cost = 0;
    for (i = 0; i < max; i++) {
        x[0] = (double)i * dx;
//        x[1] = (double)i * dx;
        evalfis(out,x,fis);
        y[i] = (double)i * dx;
        y_a[i] = out[0];
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
        y[i] = out[0];
    }
    gnuplot_plot_xy(h1, x_i, y, max, "Fuzzy Output");
    gnuplot_plot_xy(h1, x_i, x_i, max, "Expected Output");
    printf("Press any key to close window\n");
    getchar();
    gnuplot_close(h1);
}

int main(void)
{
	srand(clock());
	srand48(clock());
	int num_in = 1;
	int in_mfs[1] = {3};
	int num_out = 1;
	int out_mfs[1] = {3};

	int num_rule = prod_i(in_mfs, num_in);
	int rules[num_rule * (num_in + num_out)];
	int consequents[num_rule * num_out];

	int num_params = 3 * (sum_i(in_mfs,num_in) + sum_i(out_mfs,num_out));
	double params[num_params];

	init_params(params, num_in, in_mfs, num_out, out_mfs);
	consequents[0] = 0;
	consequents[1] = 1;
	consequents[2] = 2;
//	rand_params(params, num_in, in_mfs, num_out, out_mfs);
//	rand_consequents(num_out, num_rule, consequents, out_mfs);

	struct Specs * spcs = specs_set(num_in, in_mfs, num_out, out_mfs);
	struct Individual * ind = individual_create(num_params, params, num_rule, consequents);

	struct Fis * fis = individual_to_fis(ind, spcs);

//	double x[2] = {0.2, 0.25};
	double x[1] = {0.3};
	double out[1];

	evalfis(out,x,fis);
	printf("%f\n",out[0]);
	printf("fitness = %f\n",fit_line(fis));
	plot_line(fis);
	while (1) {
//		scanf("%lf %lf", &x[0], &x[1]);
//		printf("x = [%f, %f]\n",x[0],x[1]);
		scanf("%lf", &x[0]);
		printf("x = [%f]\n", x[0]);
		evalfis(out,x,fis);
		printf("%f\n",out[0]);
	}


	fis_destroy(fis);
	specs_clear(spcs);
	individual_destroy(ind);

	return 0;
}
