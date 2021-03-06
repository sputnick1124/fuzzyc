#ifndef _FUZZLIB_H_
#define _FUZZLIB_H_

double minimum(int n_vals, double values[n_vals]);

double maximum(int n_vals, double values[n_vals]);

double inTriMF(double params[3], double x);

void outTriMF(double x[][4], double * params[3], int num_out, double y);

void defuzzWeightedMV(double * out,
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firingStrengths[num_rule]);

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

void get_fis(struct Rule ** rule_list,
	double params[],
	int num_in,
	int num_out,
	int num_rule,
	int rules[num_rule][num_in + num_out],
	int in_mfs[],
	int out_mfs[]);

struct Rule *create_rule(double input[][3], int num_in, double output[][3], int num_out);

struct Fis {
	int num_rule;
	struct Rule ** rule_list;
};

struct Fis * fis_create(
	double params[],
	int num_in,
	int num_out,
	int num_rule,
	int rules[num_rule][num_in + num_out],
	int in_mfs[],
	int out_mfs[]);

void
fis_destroy(struct Fis * fis);

void
fis_print(struct Fis * fis,
	FILE * fd);

void print_rule(struct Rule * rule, FILE * fd);

void destroy_rule(struct Rule * rule);

void destroy_rules(struct Rule ** rules, int num_rule);

void evalrules(double * out, double * x, struct Rule ** rules, int num_rule);

void
evalfis(
	double * out,
	double * x,
	struct Fis * fis);

#endif
