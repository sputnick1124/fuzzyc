#ifndef _FUZZLIB_H_
#define _FUZZLIB_H_

double minimum(int n_vals, double values[n_vals]);

double maximum(int n_vals, double values[n_vals]);

double inTriMF(double params[3], double x);

void outTriMF(double x[][4], double * params[3], int num_out, double y);

void defuzzWMeanOfCent(double * out,
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firingStrengths[num_rule]);

void defuzzMeanOfCent(double * out,
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firingStrengths[num_rule]);

struct Rule {
	double ** input;
	int num_in;
	double ** output;
	int num_out;
};

struct Rule *create_rule(double * input[3], int num_in, double * output[3], int num_out);

void print_rule(struct Rule *rule);

void destroy_rule(struct Rule *rule);

void evalfis(double * out, double * x, struct Rule ** rules, int num_rule);

#endif
