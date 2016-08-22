#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "engine.h"
#include "mex.h"
#include "matrix.h"

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


struct Individual {
	double * params;
	int * consequents;
};

struct Individual *individual_create(
	int num_params,
	double params[],
	int num_rule,
	int num_out,
	int consequents[]);

void individual_copy(
	struct Individual * ind1,
	struct Individual * ind2,
	int num_params,
	int num_rule,
	int num_out);

void individual_destroy(struct Individual * ind);

void individuals_destroy(struct Individual ** ind, int num_ind);

struct Specs {
	int num_in;
	int * in_mfs;
	int num_out;
	int * out_mfs;
	int num_rule;
	int * rules;
	int num_params;
	double * ranges;
};

struct Specs *specs_set(
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[]);

struct Specs *specs_copy(struct Specs *spcs);

void specs_clear(struct Specs *spcs);

void specs_print(struct Specs * spcs, FILE * fd);

void individual_print(struct Individual * ind, struct Specs * spcs, FILE * fd);

struct HyperParams
{
	int pop_size;
	float elite;
	float crossover;
	float mutate;
	int max_gen;
};

struct HyperParams *
set_default_hp();

struct Fis *
individual_to_fis(
	struct Individual * ind,
	struct Specs * spcs);

typedef double (*fitness_fcn)(struct Fis * fis, Engine *engp);

int sum_i(int adds[], size_t len);

double sum_d(double adds[], size_t len);

int prod_i(int mults[], size_t len);

int rand_i(unsigned int max);

void init_partition(double params[], int num_sets);

void rand_partition(double params[], int num_sets);

void init_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]);

void rand_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]);

void init_antecedents(int num_in, int num_out, int rules[], int in_mfs[]);

void rand_consequents(int num_out, int num_rule, int consequents[num_rule * num_out], int out_mfs[]);

void
add_consequents(
	int num_in,
	int num_out,
	int num_rule,
	int rules[num_rule * (num_in + num_out)],
	int consequents[num_rule * num_out]);

void
sp_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[]);

void
sp_c_crossover(
	int num_params,
	int child1[],
	int child2[],
	int parent1[],
	int parent2[]);

void
tp_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[]);

void
tp_c_crossover(
	int num_params,
	int child1[],
	int child2[],
	int parent1[],
	int parent2[]);

void
blx_a_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[]);

void
individual_crossover(
    int num_params,
    int num_rules,
    int num_out,
    struct Individual *p_ind1,
    struct Individual *p_ind2,
    struct Individual *c_ind1,
    struct Individual *c_ind2);

void
param_range(
	int num_params,
	double ranges[num_params * 2],
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[]);

void
r_mutation(
	int num_params,
	double ranges[num_params * 2],
	double chromosome[num_params],
	int num_genes);

void
rb_mutation(
	int num_params,
	double ranges[num_params * 2],
	double chromosome[num_params],
	int cur_gen,
	int max_gen,
	double b,
	int num_genes);

void
consequent_mutate(
    int num_rules,
    int num_out,
    int consequent[],
    int out_mfs[],
    int num_genes);

void
individual_mutate(
    struct Individual *ind,
    int num_params,
    double ranges[num_params * 2],
    int cur_gen,
    int max_gen,
    double b,
    int num_rules,
    int num_out,
    int out_mfs[],
    int num_genes);

void
population_init(
    int pop_size,
    struct Individual ** population,
    int num_in,
    int in_mfs[],
    int num_out,
    int out_mfs[],
	int num_params,
    int num_rules,
    int rules[num_rules * (num_in + num_out)]);

void
population_iter(
    struct Individual ** pop_now,
    struct Individual ** pop_next,
    int rank[],
    int cur_gen,
	struct HyperParams * hp,
    struct Specs * spcs);

void
population_switch(
	struct Individual *** pop1,
	struct Individual *** pop2);

void
population_rank(
	int pop_size,
	int rank[pop_size],
	struct Individual ** population,
	struct Specs * spcs,
	double (*fitness_fcn)(struct Fis * fis, Engine *engp),
	double *fit_min,
	Engine *engp);

struct Fis *
run_ga(
    struct Specs * spcs,
    struct HyperParams * hp,
    fitness_fcn fit_fcn,
	FILE * fis_log,
	Engine * engp);


double minimum(int n_val, double values[n_val]) {
	/*Find minimum value in an array of doubles*/
	double min = HUGE_VAL;
	int i;
	for (i = 0; i < n_val; i++) {
		min = (values[i] < min ? values[i] : min);
	}
	return min;
}

double maximum(int n_val, double values[n_val]) {
	/*Find maximum value in an array of doubles*/
	double max = 0;
	int i;
	for (i = 0; i < n_val; i++) {
		max = (values[i] > max ? values[i] : max);
	}
	return max;
}

double inTriMF(double params[3], double x) {
	/**/
	double a = params[0];
	double x_star = params[1];
	double b = params[2];
	double m1, b1, m2, b2;
	#ifdef DEBUG
	printf("inTriMF: x = %0.2f\n",x);
	#endif
	if (a != x_star) {
		m1 = 1 / (x_star - a);
		b1 = -m1 * a;
	} else {
		m1 = 0;
		b1 = 0;
	}

	if (b != x_star) {
		m2 = -1 / (b - x_star);
		b2 = 1 - (m2 * x_star);
	} else {
		m2 = 0;
		b2 = 0;
	}

	if (x > b) {
		x = b;
	} else if (x < a) {
		x = a;
	}

	if (x < x_star) {
		return fmax((m1 * x) + b1, 0);
	} else if ( x > x_star) {
		return fmax((m2 * x) + b2, 0);
	} else {
		return 1;
	}
}

void outTriMF(double x[][4], double * params[], int num_out, double y) {
	int o;
	for (o = 0; o < num_out; o++) {
		double a = params[o][0];
		double x_star = params[o][1];
		double b = params[o][2];
		double m1, b1, m2, b2;
		if (a != x_star) {
			m1 = 1 / (x_star - a);
			b1 = -m1 * a;
		} else {
			m1 = 0;
			b1 = 0;
		}

		if (b != x_star) {
			m2 = -1 / (b - x_star);
			b2 = 1 - (m2 * x_star);
		} else {
			m2 = 0;
			b2 = 0;
		}

		x[o][0] = a;
		x[o][1] = (m1 ? ((y - b1) / m1) : a);
		x[o][2] = (m2 ? ((y - b2) / m2) : b);
		x[o][3] = b;
	}
}

void defuzzWeightedMV(double out[],
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firing_strengths[num_rule]) {

	int r, o;
	double hmv, ht;
	ht = 0;
	hmv = 0;
	for (o = 0; o < num_out; o++) {
		for (r = 0; r < num_rule; r++) {
			if (!firing_strengths[r]) {continue;}
			ht += firing_strengths[r];
			hmv += firing_strengths[r] * (traps[r][o][2] + traps[r][o][1]) / 2;
			#ifdef DEBUG
			printf("num_rule = %d\n",num_rule);
			printf("%f, %f, %f\n", firing_strengths[r], traps[r][o][1], traps[r][o][2]);
			printf("defuzzWeightedMV: rule %d: out %d\nht = %f\nhmv = %f\n",
				r, o, ht, hmv);
			#endif
		}
		out[o] = (ht ? hmv / ht : 0);
	}
}


void defuzzWMeanOfCent(double out[],
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firing_strengths[num_rule]) {

	int r, o;
	double mom, ht;
	double a, b, c, d;
	double c_temp, d_temp;
	double s_temp;
	mom = 0;
	ht = 0;
	for (o = 0; o < num_out; o++) {
		for (r = 0; r < num_rule; r++) {
			a = traps[r][o][2] - traps[r][o][1];
			b = traps[r][o][3] - traps[r][o][0];
			c_temp = (traps[r][o][1] - traps[r][o][0]);
			d_temp = (traps[r][o][3] - traps[r][o][2]);
			s_temp = firing_strengths[r]*firing_strengths[r];
			c = sqrt(c_temp*c_temp + s_temp);
			d = sqrt(d_temp*d_temp + s_temp);
			mom += firing_strengths[r] * (b/2 + (2*a + b)*(c*c - d*d) / (6*(b*b - a*a)));
			ht += firing_strengths[r];
			#ifdef DEBUG
			printf("trap[%d] = [%0.2f, %0.2f, %0.2f, %0.2f]\n",
					r,
					traps[r][o][0],
					traps[r][o][1],
					traps[r][o][2],
					traps[r][o][3]);
			printf("[a,b,c,d] = [%0.2f, %0.2f, %0.2f, %0.2f]\n",
					a,b,c,d);
			printf("[c,d,s_temp] = [%0.4f, %0.4f, %0.4f]\n",
					c_temp, d_temp, s_temp);
			printf("mom = %0.4f\nht = %0.3f\n",
					mom, ht);
			#endif
		}
		out[o] = mom / ht;
	}
}

void defuzzMeanOfCent(double out[],
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firing_strengths[num_rule]) {

	int r, o;
	double mom;
	double a, b, c, d;
	double c_temp, d_temp;
	double h;
	mom = 0;
	for (o = 0; o < num_out; o++) {
		for (r = 0; r < num_rule; r++) {
			a = traps[r][o][2] - traps[r][o][1];
			b = traps[r][o][3] - traps[r][o][0];
			c_temp = (traps[r][o][1] - traps[r][o][0]);
			d_temp = (traps[r][o][3] - traps[r][o][2]);
			h = firing_strengths[r];
			c = (c_temp > 0 ? sqrt(c_temp*c_temp + h*h) : h);
			d = (d_temp > 0 ? sqrt(d_temp*d_temp + h*h) : h);
			mom += (b/2 + (2*a + b)*(c*c - d*d) / (6*(b*b - a*a)));
			#ifdef DEBUG
			printf("trap[%d] = [%0.2f, %0.2f, %0.2f, %0.2f]\n",
					r,
					traps[r][o][0],
					traps[r][o][1],
					traps[r][o][2],
					traps[r][o][3]);
			printf("[a,b,c,d] = [%0.2f, %0.2f, %0.2f, %0.2f]\n",
					a,b,c,d);
			printf("[c_temp,d_temp,h] = [%0.4f, %0.4f, %0.4f]\n",
					c_temp, d_temp, h);
			printf("mom = %0.4f\n",
					mom);
			#endif
		}
		out[o] = mom / (double) num_rule;
	}
}

struct Rule *create_rule(double input[][3], int num_in, double output[][3], int num_out) {
    int i, p;
    struct Rule *rule = malloc(sizeof(struct Rule));
    assert(rule != NULL);

	rule->input = malloc(num_in * sizeof(double *));

    rule->num_in = num_in;
    for (i = 0; i < num_in; i ++) {
		rule->input[i] = malloc(3 * sizeof(double));
		if (rule->input[i] == NULL) {
			for ( ; i >= 0; i--) {
				free(rule->input[i]);
			}
			free(rule);
			return NULL;
		}
		for (p = 0; p < 3; p++) {
        	rule->input[i][p] = input[i][p];
		}
    }

	rule->output = malloc(num_out * sizeof(double *));

    rule->num_out = num_out;
    for (i = 0; i < num_out; i ++) {
		rule->output[i] = malloc(3 * sizeof(double));
		if (rule->output[i] == NULL) {
			for ( ; i >= 0; i--) {
				free(rule->output[i]);
			}
			free(rule);
			return NULL;
		}
		for (p = 0; p < 3; p++) {
        	rule->output[i][p] = output[i][p];
		}
    }

    return rule;
}

void print_rule(struct Rule *rule, FILE * fd) {
    int i, p;

	if (fd == NULL) {
		fd = stdout;
	}

    fprintf(fd, "%d inputs:\n",rule->num_in);
    for (i = 0; i < rule->num_in; i++) {
        for (p = 0; p < 3; p++) {
            fprintf(fd, "%0.2f ", rule->input[i][p]);
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "%d outputs:\n",rule->num_out);
    for (i = 0; i < rule->num_out; i++) {
        for (p = 0; p < 3; p++) {
            fprintf(fd, "%0.2f ", rule->output[i][p]);
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "\n\n");
}

void destroy_rules(struct Rule ** rule_list, int num_rule) {
	int i;
	for (i = 0; i < num_rule; i++) {
		destroy_rule(rule_list[i]);
	}
}

void destroy_rule(struct Rule *rule) {
	int i;
	assert(rule != NULL);
	for (i = 0; i < rule->num_in; i++) {
		free(rule->input[i]);
	}
	free(rule->input);
	for (i = 0; i < rule->num_out; i++) {
		free(rule->output[i]);
	}
	free(rule->output);
	free(rule);
}

struct Fis * fis_create(
    double params[],
    int num_in,
    int num_out,
    int num_rule,
    int rules[num_rule][num_in + num_out],
    int in_mfs[],
    int out_mfs[])
{
	struct Fis *fis = malloc(sizeof(struct Fis));
	assert(fis != NULL);

	fis->num_rule = num_rule;
	fis->rule_list = malloc(num_rule * sizeof(struct Rule *));
	if (fis->rule_list == NULL) {
		free(fis);
		return NULL;
	}

	get_fis(fis->rule_list, params, num_in, num_out, num_rule, rules, in_mfs, out_mfs);
	return fis;
}

void
fis_destroy(struct Fis *fis)
{
	destroy_rules(fis->rule_list,fis->num_rule);
	free(fis->rule_list);
	free(fis);
}

void
evalrules(
	double * out,
	double * x,
	struct Rule ** rules,
	int num_rule)
{
	int r, in;
	double firing_strengths[num_rule];
	double s_temp[rules[0]->num_in];
	double out_vals[num_rule][rules[0]->num_out][4];
	#ifdef DEBUG
	printf("num_rule = %d\n",num_rule);
	#endif
	for (r = 0; r < num_rule; r++) {
		#ifdef DEBUG
		print_rule(rules[r]);
		#endif
		for (in = 0; in < rules[r]->num_in; in++) {
			s_temp[in] = inTriMF(rules[r]->input[in], x[in]);
			#ifdef DEBUG
			printf("evalfis: x[%d] = %0.2f\n",in,x[in]);
			printf("inTriMF([%0.2f, %0.2f, %0.2f], %0.2f) = %0.2f\n",
					rules[r]->input[in][0],
					rules[r]->input[in][1],
					rules[r]->input[in][2],
					x[in],
					s_temp[in]);
			#endif
		}
		firing_strengths[r] = minimum(rules[r]->num_in,s_temp);
		outTriMF(out_vals[r],
				 rules[r]->output,
				 rules[0]->num_out,
				 firing_strengths[r]);
		#ifdef DEBUG
		int d, o;
		printf("firing_strength = min(inTriMF(x)) = %0.2f\n",firing_strengths[r]);
		for (o = 0; o < rules[0]->num_out; o++) {
			printf("outTriMF([%0.2f, %0.2f, %0.2f], %0.2f) = [",
					rules[r]->output[o][0],
					rules[r]->output[o][1],
					rules[r]->output[o][2],
					firing_strengths[r]);
			for (d = 0; d < 4; d++) {
				printf("%0.5f ",out_vals[r][o][d]);
			}
		}
		printf("]\n\n\n");
		#endif
	}
	defuzzWeightedMV(out,
		num_rule,
		rules[0]->num_out,
		out_vals,
		firing_strengths);
}

void
evalfis(
	double * out,
	double * x,
	struct Fis * fis)
{
	evalrules(out, x, fis->rule_list, fis->num_rule);
}


void get_fis(struct Rule ** rule_list,
	double params[],
	int num_in,
	int num_out,
	int num_rule,
	int rules[num_rule][num_in + num_out],
	int in_mfs[],
	int out_mfs[]) {

	int in, out, rule;
	int pin, pout;

	int ppoint, out_start;
	int _b;

	for (rule = 0; rule < num_rule; rule++) {
		#ifdef DEBUG
			printf("RULE%d\n",rule);
		#endif
		double input_list[num_in][3];
		double output_list[num_out][3];
		for (in = 0; in < num_in; in++) {
			ppoint = 0;
			if (in > 0) {
				for (_b = in-1; _b >= 0; _b--) {
					ppoint += 3 * in_mfs[_b];
				}
			} else { };
			ppoint += 3 * rules[rule][in];
			for (pin = 0; pin < 3; pin++) {
				#ifdef DEBUG
					printf("get_fis: INPUT:");
					printf("in = %d, pin = %d, ppoint = %d\n",in,pin,ppoint);
				#endif
				input_list[in][pin] = params[ppoint];
				ppoint += 1;
			}
		}

		out_start = 0;
		for (in = 0; in < num_in; in++) {
			out_start += 3 * in_mfs[in];
		}

		for (out = 0; out < num_out; out++) {
			ppoint = out_start;
			if (out > 0) {
				for (_b = out-1; _b >= 0; _b--) {
					ppoint += 3 * out_mfs[_b];
				}
			} else { };
			ppoint += 3 * rules[rule][num_in + out];
			for (pout = 0; pout < 3; pout++) {
				#ifdef DEBUG
					printf("get_fis: OUTPUT:");
					printf("out = %d, pout = %d, ppoint = %d\n",out, pout, ppoint);
				#endif
				output_list[out][pout] = params[ppoint];
				ppoint += 1;
			}
		}

		rule_list[rule] = create_rule(input_list, num_in, output_list, num_out);
	}
}

void
fis_print(struct Fis * fis, FILE * fd)
{
	int rule, in, out;

	if (fd == NULL) {
		fd = stdout;
	}

	for (rule = 0; rule < fis->num_rule; rule++) {
		fprintf(fd, "INPUT %d:\n\t[", rule);
		for (in = 0; in < fis->rule_list[rule]->num_in; in++) {
			print_rule(fis->rule_list[rule], fd);
		}
		fprintf(fd, "]\nOUTPUT %d:\n\t[", rule);
		for (out = 0; out < fis->rule_list[rule]->num_out; out++) {
			print_rule(fis->rule_list[rule], fd);
		}
		fprintf(fd, "]\n");
	}
}


/** Some utility functions to take care of the mechanics **/
int
prod_i(int mults[], size_t len)
{
	int retval = 1;
	int i;
	for (i = 0; i < len; i++) {
		retval *= mults[i];
	}
	return retval;
}

int
sum_i(int adds[], size_t len)
{
	int retval = 0;
	int i;
	for (i = 0; i < len; i++) {
		retval += adds[i];
	}
	return retval;
}

double
sum_d(double adds[], size_t len)
{
	double retval = 0;
	int i;
	for (i = 0; i < len; i++) {
		retval += adds[i];
	}
	return retval;
}

int
rand_i(unsigned int max)
{
	/* return random integer in the range [0,max)*/
	int r;
	const unsigned int buckets = RAND_MAX / max;
	const unsigned int limit = buckets * max;

	do {
		r = rand();
	} while (r >= limit);

	/* cast assumes that max is much less than RAND_MAX*/
	return (int) (r / buckets);
}

int
rand_tri_i(int max)
{
	/* return integer in the range [0, max) from a triangular distribution
	with a,c = 0, and b = max
	see https://en.wikipedia.org/wiki/Triangular_distribution*/
    double u = drand48();

	/* implicit cast to int is intentional*/
    return max - sqrt((1 - u) * max * max);
}

static int
cmpdouble_p(const void *pp1, const void *pp2)
{
	double ** p1 = (double **)pp1;
	double ** p2 = (double **)pp2;
	if (**p1 > **p2) {
	return 1;
	} else if (**p1 < **p2) {
		return -1;
	} else {
		return 0;
	}
}

struct Specs *
specs_set(
	int num_in,
    int in_mfs[],
    int num_out,
    int out_mfs[])
{
	int in, out;
	struct Specs * spcs = malloc(sizeof(struct Specs));
	assert(spcs != NULL);

	spcs->num_in = num_in;
	spcs->in_mfs = malloc(num_in * sizeof(int));
	if (spcs->in_mfs == NULL) {
		free(spcs);
		return NULL;
	}

	spcs->num_out = num_out;
	spcs->out_mfs = malloc(num_out * sizeof(int));
	if (spcs->out_mfs == NULL) {
		free(spcs->in_mfs);
		free(spcs);
		return NULL;
	}

	spcs->num_rule = prod_i(in_mfs, num_in);
	spcs->rules = malloc(spcs->num_rule * (num_in + num_out) * sizeof(int));
	if (spcs->rules == NULL) {
		free(spcs->out_mfs);
		free(spcs->in_mfs);
		free(spcs);
		return NULL;
	}

	spcs->num_params = 3 * (sum_i(in_mfs, num_in) + sum_i(out_mfs, num_out));
	spcs->ranges = malloc(spcs->num_params * 2 * sizeof(double));
	if (spcs->ranges == NULL) {
		free(spcs->rules);
		free(spcs->out_mfs);
		free(spcs->in_mfs);
		free(spcs);
		return NULL;
	}

	for (in = 0; in < num_in; in++) {
		spcs->in_mfs[in] = in_mfs[in];
	}
	for (out = 0; out < num_out; out++) {
		spcs->out_mfs[out] = out_mfs[out];
	}

	init_antecedents(
		spcs->num_in,
		spcs->num_out,
		spcs->rules,
		spcs->in_mfs);
	param_range(
		spcs->num_params,
		spcs->ranges,
		spcs->num_in,
		spcs->in_mfs,
		spcs->num_out,
		spcs->out_mfs);

	return spcs;
}

struct Specs *
specs_copy(struct Specs * oldspcs)
{
	struct Specs * spcs = malloc(sizeof(struct Specs));
	assert(spcs != NULL);

	spcs->num_in = oldspcs->num_in;
	spcs->in_mfs = malloc(spcs->num_in * sizeof(int));
	if (spcs->in_mfs == NULL) {
		free(spcs);
		return NULL;
	}

	spcs->num_out = oldspcs->num_out;
	spcs->out_mfs = malloc(spcs->num_out * sizeof(int));
	if (spcs->in_mfs == NULL) {
		free(spcs->in_mfs);
		free(spcs);
		return NULL;
	}

	spcs->num_rule = oldspcs->num_rule;
	spcs->rules = malloc(spcs->num_rule * (spcs->num_in + spcs->num_out) * sizeof(int));
	if (spcs->rules == NULL) {
		free(spcs->out_mfs);
		free(spcs->in_mfs);
		free(spcs);
		return NULL;
	}

	spcs->num_params = oldspcs->num_params;
	spcs->ranges = malloc(spcs->num_params * 2 * sizeof(double));
	if (spcs->ranges == NULL) {
		free(spcs->rules);
		free(spcs->out_mfs);
		free(spcs->in_mfs);
		free(spcs);
		return NULL;
	}

	memcpy(spcs->in_mfs, oldspcs->in_mfs, spcs->num_in * sizeof(int));

	memcpy(spcs->out_mfs, oldspcs->out_mfs, spcs->num_out * sizeof(int));

	memcpy(spcs->rules, oldspcs->rules, (spcs->num_rule * (spcs->num_in + spcs->num_out) * sizeof(int)));

	memcpy(spcs->ranges, oldspcs->ranges, spcs->num_params * 2 * sizeof(double));

	return spcs;
}

void
specs_clear(struct Specs *spcs)
{
	free(spcs->ranges);
	free(spcs->rules);
	free(spcs->out_mfs);
	free(spcs->in_mfs);
	free(spcs);
}

/** Initialization functions for starting off the GA population**/
struct Individual *
individual_create(
	int num_params,
	double params[],
	int num_rule,
	int num_out,
	int consequents[])
{
	int p, r;
	struct Individual *ind = malloc(sizeof(struct Individual));
	assert(ind != NULL);

	ind->params = malloc(num_params * sizeof(double));
	if (ind->params == NULL) {
		free(ind);
		return NULL;
	}
	for (p = 0; p < num_params; p++) {
		ind->params[p] = params[p];
	}

	ind->consequents = malloc(num_rule * num_out * sizeof(int));
	if (ind->consequents == NULL) {
		free(ind->params);
		free(ind);
		return NULL;
	}
	for (r = 0; r < num_rule * num_out; r++) {
		ind->consequents[r] = consequents[r];
	}

	return ind;
}

void
individual_copy(
	struct Individual * ind1,
	struct Individual * ind2,
	int num_params,
	int num_rule,
	int num_out)
{
	memcpy(ind2->params, ind1->params, num_params * sizeof(double));

	memcpy(ind2->consequents, ind1->consequents, num_rule * num_out * sizeof(int));
}

void
individual_destroy(struct Individual * ind)
{
	free(ind->params);
	free(ind->consequents);
	free(ind);
}

void
individuals_destroy(struct Individual ** inds, int num_ind)
{
	int i;
	for (i = 0; i < num_ind; i++) {
		individual_destroy(inds[i]);
	}
}

void
init_partition(double params[], int num_sets) {
	int p;
	int num_params = num_sets * 3;
	params[0] = 0; params[1] = 0;
	params[num_params-1] = 1;
	params[num_params-2] = 1;
	double I = 1 / ((double) num_sets - 1);
	for (p = 2; p < num_params - 2; p++) {
		if (p % 3 == 0) {
			params[p] = ((p / 3) - 1) * I;
		} else if (p % 3 == 1) {
			params[p] = ((p - 1) / 3) * I;
		} else {
			params[p] = (((p - 2) / 3) + 1) * I;
		}
	}
}

void
rand_partition(double params[], int num_sets)
{
	int p;
	double r;
	int num_params = num_sets * 3;
	params[0] = 0; params[1] = 0;
	params[num_params-1] = 1;
	params[num_params-2] = 1;
	double I = 1 / ((double) num_sets - 1);

	for (p = 2; p < num_params - 2; p++) {
		r = drand48() * I;
		if (p == 3) {
			params[p] = r / 2;
		} else if ( p == num_params - 4) {
			params[p] = (((p - 2) / 3) + 1) * I - r / 2;
		} else if (p % 3 == 0) {
			params[p] = ((p / 3) - 1) * I + r - (I / 2);
		} else if (p % 3 == 1) {
			params[p] = ((p - 1) / 3) * I + r - (I / 2);
		} else {
			params[p] = (((p - 2) / 3) + 1) * I + r - (I / 2);
		}
	}
}

void
init_params(
	double params[],
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[])
{
	int in, out, p;

	p = in_mfs[0];
	init_partition(params, p);
	for (in = 1; in < num_in; in++) {
		init_partition(&params[p * 3], in_mfs[in]);
		p += in_mfs[in];
	}

	for (out = 0; out < num_out; out++) {
		init_partition(&params[p * 3], out_mfs[out]);
		p += out_mfs[out];
	}
}

void
rand_params(
	double params[],
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[])
{
	int in, out, p;

	p = in_mfs[0];
	rand_partition(params, p);
	for (in = 1; in < num_in; in++) {
		rand_partition(&params[p * 3], in_mfs[in]);
		p += in_mfs[in];
	}

	for (out = 0; out < num_out; out++) {
		rand_partition(&params[p * 3], out_mfs[out]);
		p += out_mfs[out];
	}
}

void
init_antecedents(
    int num_in,
    int num_out,
    int *rules,
    int in_mfs[])
{
    int in, rule;
    unsigned int flag = 0;
    int num_rule = 1;
    int width = num_in + num_out;
    for (in = 0; in < num_in; in++) {
        rules[0 + in] = 0;
        num_rule *= in_mfs[in];
    }
    for (rule = 1; rule < num_rule; rule++) {
        if (rules[(rule-1) * width + 0] == (in_mfs[0] - 1)) {
            rules[rule * width + 0] = 0;
            flag ^= (1 << 1);
        } else {
            rules[rule * width + 0] = rules[(rule-1) * width + 0] + 1;
        }
        for (in = 1; in < num_in; in++) {
            if (flag & (1 << in)) {
                if (rules[(rule-1) * width + in] == (in_mfs[in]-1)) { 
                    rules[rule * width + in] = 0; 
                    flag ^= (1 << (in + 1)); 
                } else {
                    rules[rule * width + in] = rules[(rule-1) * width + in] + 1;
                }
                flag ^= (1 << in);
            } else {
                rules[rule * width + in] = rules[(rule-1) * width + in];
            }
        }
    }
}

void
rand_consequents(
	int num_out,
	int num_rule,
	int consequents[num_rule * num_out],
	int out_mfs[])
{
	int out, rule;
	for (rule = 0; rule < num_rule; rule++) {
		for (out = 0; out < num_out; out++) {
			consequents[rule * num_out + out] = rand_i(out_mfs[out]);
		}
	}
}

void
add_consequents(
	int num_in,
	int num_out,
	int num_rule,
	int *rules,
	int consequents[num_rule * num_out])
{
	int out, rule;
	int width = num_in + num_out;
	for (rule = 0; rule < num_rule; rule++) {
		for (out = 0; out < num_out; out++) {
			rules[rule * width + num_in + out] = consequents[rule * num_out + out];
		}
	}
}

/** Crossover functions **/
void
sp_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[])
{
/* Single-point crossover*/
	/*We definitely don't care about the first or last two params*/
	int i = 2 + rand_i(num_params - 4);
	int p;
	for (p = 0; p < num_params; p++) {
		child1[p] = (p <= i ? parent1[p] : parent2[p]);
		child2[p] = (p <= i ? parent2[p] : parent1[p]);
	}
}

void
sp_c_crossover(
	int num_params,
	int child1[],
	int child2[],
	int parent1[],
	int parent2[])
{
/* Single-point crossover*/
	int i = rand_i(num_params);
	int p;
	printf("%d\n",i);
	for (p = 0; p < num_params; p++) {
		child1[p] = (p <= i ? parent1[p] : parent2[p]);
		child2[p] = (p <= i ? parent2[p] : parent1[p]);
	}
}

void
tp_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[])
{
/* Two-point crossover*/
	int i = 2 + rand_i(num_params - 2);
	int j = (i < num_params - 2 ? i + rand_i(num_params - 2 - i) : i);
	int p;
	for (p = 0; p <= j; p++) {
		child1[p] = (p <= i ? parent1[p] : parent2[p]);
		child2[p] = (p <= i ? parent2[p] : parent1[p]);
	}
	for (p = j; p < num_params; p++) {
		child1[p] = parent1[p];
		child2[p] = parent2[p];
	}
}

void
tp_c_crossover(
	int num_params,
	int child1[],
	int child2[],
	int parent1[],
	int parent2[])
{
/* Two-point crossover*/
	int i = rand_i(num_params);
	int j = (i < num_params ? i + rand_i(num_params - i) : i);
	int p;
	for (p = 0; p <= j; p++) {
		child1[p] = (p <= i ? parent1[p] : parent2[p]);
		child2[p] = (p <= i ? parent2[p] : parent1[p]);
	}
	for (p = j; p < num_params; p++) {
		child1[p] = parent1[p];
		child2[p] = parent2[p];
	}
}

void
blx_a_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[])
{
/* Blended crossover with alpha exploration factor*/
	double I;
	double xmin, xmax;
	double new_xmin, new_xmax;
	double alpha = 0.1;
	double lambda;
	int p;
	for (p = 0; p < num_params - 2; p++) {
		xmin = fmin(parent1[p], parent2[p]);
		xmax = fmax(parent1[p], parent2[p]);
		I = xmax - xmin;
		lambda = drand48();
		new_xmin = xmin - I * alpha;
		new_xmin = (new_xmin > 0 ? new_xmin : 0);
		new_xmax = xmax + I * alpha;
		new_xmax = ( new_xmax < 1 ? new_xmax : 1);
		child1[p] = lambda * new_xmin + (1 - lambda) * new_xmax;
		child2[p] = (1 - lambda) * new_xmin + lambda * new_xmax;
		if ((child1[p] < 0) | (child1[p] > 1)) {
			printf("p1,p2,I,c1,L = [%0.4f,%0.4f,%0.4f,%0.4f,%0.4f] \n",parent1[p],parent2[p],I,child1[p],lambda);
			exit(EXIT_FAILURE);
		}
		if ((child2[p] < 0) | (child2[p] > 1)) {
			printf("p1,p2,I,c2,L = [%0.4f,%0.4f,%0.4f,%0.4f,%0.4f] \n",parent1[p],parent2[p],I,child1[p],lambda);
			exit(EXIT_FAILURE);
		}
	}
	child1[num_params-2] = parent1[num_params-2];
	child1[num_params-1] = parent1[num_params-1];
	child2[num_params-2] = parent2[num_params-2];
	child2[num_params-1] = parent2[num_params-1];
}

void
individual_crossover(
	int num_params,
	int num_rules,
	int num_out,
	struct Individual *p_ind1,
	struct Individual *p_ind2,
	struct Individual *c_ind1,
	struct Individual *c_ind2)
{
	blx_a_crossover(
		num_params,
		c_ind1->params,
		c_ind2->params,
		p_ind1->params,
		p_ind2->params);
	tp_c_crossover(
		num_rules * num_out,
		c_ind1->consequents,
		c_ind2->consequents,
		p_ind1->consequents,
		p_ind2->consequents);
}

/** Mutation functions **/
void
param_range(
	int num_params,
	double ranges[num_params * 2],
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[])
{
	int p = 0;
	int in, out, mf, mfp;
	int mfs[num_in + num_out];
	double I;

	for (in = 0; in < num_in; in++) {
		mfs[in] = in_mfs[in];
	}
	for (out = 0; out < num_out; out++) {
		mfs[num_in + out] = out_mfs[out];
	}

	for (mf = 0; mf < num_in + num_out; mf++) {
		I = 1 / ((double) mfs[mf] - 1);
		for (mfp = 0 ; mfp < mfs[mf] * 3; mfp++, p++) {
			if ((mfp == 0) | (mfp == 1)
				| (mfp == mfs[mf] * 3 - 1)
				| (mfp == mfs[mf] * 3 - 2)) {
				ranges[p * 2 + 0] = -1;
				ranges[p * 2 + 1] = -1;
			} else if (mfp == 3) {
				ranges[p * 2 + 0] = 0;
				ranges[p * 2 + 1] = I / 2;
			} else if (mfp == (mfs[mf] * 3 - 4)) {
				ranges[p * 2 + 0] = 1 - I / 2;
				ranges[p * 2 + 1] = 1;
			} else if ((mfp % 3) == 0) { 
				ranges[p * 2 + 0] = ((mfp / 3) - 1) * I - I / 2;
				ranges[p * 2 + 1] = ((mfp / 3) - 1) * I + I / 2;
			} else if ((mfp % 3) == 1) { 
				ranges[p * 2 + 0] = (((mfp - 1) / 3)) * I - I / 2;
				ranges[p * 2 + 1] = (((mfp - 1) / 3)) * I + I / 2;
			} else { 
				ranges[p * 2 + 0] = (((mfp - 2) / 3) + 1) * I - I / 2;
				ranges[p * 2 + 1] = (((mfp - 2) / 3) + 1) * I + I / 2;
			}
		}
	}
}

void
r_mutation(
	int num_params,
	double ranges[num_params * 2],
	double chromosome[num_params],
	int num_genes)
{
	int p, i;
	for (i = 0; i < num_genes; i++) {
		do {
			p = rand_i(num_params);
		} while (ranges[p * 2] == -1);
		chromosome[p] = ranges[p * 2] + rand() * (ranges[p * 2 + 1] - ranges[p * 2]);
	}
}

void
rb_mutation(
	int num_params,
	double ranges[num_params * 2],
	double chromosome[num_params],
	int cur_gen,
	int max_gen,
	double b,
	int num_genes)
{
	int tau;
	double r, del;
	int p, i;
	for (i = 0; i < num_genes; i++) {
		tau = rand_i(2);
		r = drand48();
		do {
			p = rand_i(num_params);
		} while (ranges[p * 2] == -1);
		if (tau) {
			del = (chromosome[p] - ranges[p * 2]) * (1 - r * (1 - pow((cur_gen / max_gen),b)));
			chromosome[p] -= del;
		} else {
			del = (ranges[p * 2 + 1] - chromosome[p]) * (1 - r * (1 - pow((cur_gen / max_gen),b)));
			chromosome[p] += del;
		}
	}
}

void
consequent_mutate(
	int num_rules,
	int num_out,
	int consequents[],
	int out_mfs[],
	int num_genes)
{
	int gen, mf, i;
	int mut[] = {-1, 1};
	int tau = rand_i(2);
	for (i = 0; i < num_genes; i++) {
		gen = rand_i(num_rules * num_out);
		mf = gen / (num_rules * num_out);
		if (consequents[gen] < out_mfs[mf] - 1) {
			if (consequents[gen] != 0) {
				consequents[gen] += mut[tau];
			} else {
				consequents[gen] += 1;
			}
		} else {
			consequents[gen] -= 1;
		}
	}
}

void
individual_mutate(
	struct Individual *ind,
	int num_params,
	double ranges[num_params * 2],
	int cur_gen,
	int max_gen,
	double b,
	int num_rules,
	int num_out,
	int out_mfs[],
	int num_genes)
{
	rb_mutation(
		num_params,
		ranges,
		ind->params,
		cur_gen,
		max_gen,
		b,
		num_genes);
	consequent_mutate(
		num_rules,
		num_out,
		ind->consequents,
		out_mfs,
		num_genes);

}

/**GA execution**/
void
population_init(
	int pop_size,
	struct Individual ** population,
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[],
	int num_params,
	int num_rules,
	int rules[num_rules * (num_in + num_out)])
{
	int i;
	double tmp_params[num_params];
	int tmp_consequents[num_rules * num_out];
	for (i = 0; i < pop_size; i++) {
		rand_params(tmp_params, num_in, in_mfs, num_out, out_mfs);
		rand_consequents(num_out, num_rules, tmp_consequents, out_mfs);
		population[i] = individual_create(num_params, tmp_params, num_rules, num_out, tmp_consequents);
	}
}

void
population_iter(
	struct Individual ** pop_now,
	struct Individual ** pop_next,
	int rank[],
	int cur_gen,
	struct HyperParams * hp,
	struct Specs * spcs)
{
	int ind, select1, select2;
	int pop_size = hp->pop_size;
	int num_elite = pop_size * hp->elite;
	int num_cross = pop_size * hp->crossover;
	int num_mut = pop_size * hp->mutate;
	for (ind = 0; ind < pop_size; ind++) {
		if (ind < num_elite) {
			individual_copy(
				pop_now[rank[ind]],
				pop_next[ind],
				spcs->num_params,
				spcs->num_rule,
				spcs->num_out);
		} else if (ind < num_elite + num_cross + num_mut) {
			select1 = rand_tri_i(pop_size);
			select2 = rand_tri_i(pop_size);
			individual_crossover(
				spcs->num_params,
				spcs->num_rule,
				spcs->num_out,
				pop_now[rank[select1]],
				pop_now[rank[select2]],
				pop_next[ind],
				pop_next[ind + 1]);
				ind++;
			if (ind >= num_elite + num_cross) {
				individual_mutate(
					pop_next[ind],
					spcs->num_params,
					spcs->ranges,
					cur_gen,
					hp->max_gen,
					1.5,
					spcs->num_rule,
					spcs->num_out,
					spcs->out_mfs,
					4);
				individual_mutate(
					pop_next[ind + 1],
					spcs->num_params,
					spcs->ranges,
					cur_gen,
					hp->max_gen,
					1.5,
					spcs->num_rule,
					spcs->num_out,
					spcs->out_mfs,
					4);
			}
		} else {
			rand_params(
				pop_next[ind]->params,
				spcs->num_in,
				spcs->in_mfs,
				spcs->num_out,
				spcs->out_mfs);
			rand_consequents(
				spcs->num_out,
				spcs->num_rule,
				pop_next[ind]->consequents,
				spcs->out_mfs);
		}
	}
}

struct Fis *
individual_to_fis(
	struct Individual * ind,
	struct Specs * spcs)
{
	int rule, arg;
	int width = spcs->num_in + spcs->num_out;
	int rules[spcs->num_rule][width];
	add_consequents(
		spcs->num_in,
		spcs->num_out,
		spcs->num_rule,
		spcs->rules,
		ind->consequents);
	for (rule = 0; rule < spcs->num_rule; rule++) {
		for (arg = 0; arg < width; arg++) {
			rules[rule][arg] = spcs->rules[rule * width + arg];
		}
	}
	struct Fis * fis = fis_create(
		ind->params,
		spcs->num_in,
		spcs->num_out,
		spcs->num_rule,
		rules,
		spcs->in_mfs,
		spcs->out_mfs);

	return fis;
}

void
population_rank(
	int pop_size,
	int rank[pop_size],
	struct Individual ** population,
	struct Specs * spcs,
	fitness_fcn fit_fcn,
	double * fit_min,
	Engine * engp)
{
	int ind;
	double fitness[pop_size];
	double * ranked_fitness[pop_size];
	double tmp_fit;
	struct Fis * tmp_fis;
	struct Specs * tmp_spcs;
	{
		for (ind = 0; ind < pop_size; ind++) {
			tmp_spcs = specs_copy(spcs);
			tmp_fis = individual_to_fis(population[ind],tmp_spcs);
			tmp_fit = fit_fcn(tmp_fis, engp);
			fitness[ind] = tmp_fit;
			ranked_fitness[ind] = &fitness[ind];
			fis_destroy(tmp_fis);
			specs_clear(tmp_spcs);
		}
	} 
	*fit_min = minimum(pop_size, fitness);
	qsort(ranked_fitness, pop_size, sizeof(double *), cmpdouble_p);
	for (ind = 0; ind < pop_size; ind++) {
		rank[ind] = (int)(ranked_fitness[ind] - &fitness[0]);
	}
}

void
population_switch(
	struct Individual *** pop1,
	struct Individual *** pop2)
{
	struct Individual ** tmp = *pop1;
	*pop1 = *pop2;
	*pop2 = tmp;
}

struct Fis *
run_ga(
	struct Specs * spcs,
	struct HyperParams * hp,
	fitness_fcn fit_fcn,
	FILE * fis_log,
	Engine *engp)
{
	struct Individual ** pop1 = malloc(hp->pop_size * sizeof(struct Individual *));
	struct Individual ** pop2 = malloc(hp->pop_size * sizeof(struct Individual *));
	double fitness_hist[hp->max_gen];
	int rank[hp->pop_size];
	int gen;
	population_init(
		hp->pop_size,
		pop1,
		spcs->num_in,
		spcs->in_mfs,
		spcs->num_out,
		spcs->out_mfs,
		spcs->num_params,
		spcs->num_rule,
		spcs->rules);
	population_init(
		hp->pop_size,
		pop2,
		spcs->num_in,
		spcs->in_mfs,
		spcs->num_out,
		spcs->out_mfs,
		spcs->num_params,
		spcs->num_rule,
		spcs->rules);

	for (gen = 0; gen < hp->max_gen; gen++) {
		population_rank(hp->pop_size, rank, pop1, spcs, fit_fcn, &fitness_hist[gen], engp);
		mexPrintf("Ind[%d] Fitness: %f\n",rank[0],fitness_hist[gen]);
		if ( (gen > 10) && (fabs(sum_d(&fitness_hist[gen - 10], 10)/10.0 - fitness_hist[gen]) < 1e-17) ) {
			struct Fis * ret_fis =  individual_to_fis(pop1[rank[0]],spcs);
			printf("Best fitness:%f\n",fit_fcn(ret_fis, engp));
			individual_print(pop1[rank[0]], spcs, NULL);
			individuals_destroy(pop1, hp->pop_size);
			individuals_destroy(pop2, hp->pop_size);
			free(pop1);
			free(pop2);
			printf("%d Generations\n",gen);
			return ret_fis;
		} else if (gen == hp->max_gen - 1) {
			break;
		}
		population_iter(pop1, pop2, rank, gen, hp, spcs);
		population_switch(&pop1, &pop2);
	}
	individual_print(pop1[rank[0]], spcs, fis_log);
	struct Fis * ret_fis =  individual_to_fis(pop1[rank[0]],spcs);
	individuals_destroy(pop1, hp->pop_size);
	individuals_destroy(pop2, hp->pop_size);
	free(pop1);
	free(pop2);
	return ret_fis;
}

void
specs_print(struct Specs * spcs, FILE * fd)
{
	if (fd == NULL) {
		return;
	}

	int in, out, rule, width;
	fprintf(fd, "Specs:\n");
	fprintf(fd, "num_in: %d\tin_mfs: {",spcs->num_in);
	for (in = 0; in < spcs->num_in; in++) {
		fprintf(fd, "%d ", spcs->in_mfs[in]);
	}
	fprintf(fd, "}\nnum_out: %d\tout_mfs: {",spcs->num_out);
	for (out = 0; out < spcs->num_out; out++) {
		fprintf(fd, "%d ", spcs->out_mfs[out]);
	}

	fprintf(fd,"}\nnum_rule: %d\n", spcs->num_rule);
	width = spcs->num_in + spcs->num_out;
	for (rule = 0; rule < spcs->num_rule; rule++) {
		fprintf(fd, "\t{");
		for (in = 0; in < spcs->num_in; in++) {
			fprintf(fd, "%d ", spcs->rules[rule * width + in]);
		}
		fprintf(fd, "}\n");
	}
}

void
individual_print(struct Individual * ind, struct Specs * spcs, FILE * fd)
{
	if (fd == NULL) {
		return;
	}

	int p, in, out;
	int rule, arg;
	int width = spcs->num_in + spcs->num_out;
	add_consequents(
		spcs->num_in,
		spcs->num_out,
		spcs->num_rule,
		spcs->rules,
		ind->consequents);

	fprintf(fd, "num_in: %d\tin_mfs: {",spcs->num_in);
	for (in = 0; in < spcs->num_in; in++) {
		fprintf(fd, "%d ", spcs->in_mfs[in]);
	}

	fprintf(fd, "}\nnum_out: %d\tout_mfs: {",spcs->num_out);
	for (out = 0; out < spcs->num_out; out++) {
		fprintf(fd, "%d ", spcs->out_mfs[out]);
	}
	fprintf(fd, "}\n");

	int mf = 0;
	fprintf(fd, "MF Parameters:\n");
	for (in = 0; in < spcs->num_in; in++) {
		fprintf(fd, "\tInput %d:\n", in);
		for (p = 0; p < spcs->in_mfs[in]; p++) {
			fprintf(fd, "\t{%f %f %f}\n",
				ind->params[3 * mf + 3 * p],
				ind->params[3 * mf + 3 * p + 1],
				ind->params[3 * mf + 3 * p + 2]);
		}
		mf += spcs->in_mfs[in];
	}
	for (out = 0; out < spcs->num_out; out++) {
		fprintf(fd, "\tOutput %d:\n", out);
		for (p = 0; p < spcs->out_mfs[out]; p++) {
			fprintf(fd, "\t{%f %f %f}\n",
				ind->params[3 * mf + 3 * p],
				ind->params[3 * mf + 3 * p + 1],
				ind->params[3 * mf + 3 * p + 2]);
		}
		mf += spcs->out_mfs[out];
	}


	fprintf(fd, "Rule Matrix:\n");

	for (rule = 0; rule < spcs->num_rule; rule++) {
		fprintf(fd, "{");
		for (arg = 0; arg < width; arg++) {
			fprintf(fd, "%d ", spcs->rules[rule * width + arg]);
		}
		fprintf(fd, "}\n");
	}
}


static const int nsteps = 2000;
static const double t1 = 20.0;
static const double scalar = 30.0;
static const double abs_err = 1e-10;
static const double rel_err = 1e-10;
static const int nout = 3;

static char * filenames[4] = {"approach.tex", "deg.tex", "sub.tex", "sup.tex"};


double stepinfo(int nsteps, double t[nsteps], double x[nsteps], int flag, Engine *engp) {
	double xf = x[nsteps - 1];
	double *ts;
	double *os;
	double *a;
	mxArray *T = NULL, *X = NULL, *Ts = NULL, *Os = NULL, *A = NULL;
	T = mxCreateDoubleMatrix(1,nsteps,mxREAL);
	X = mxCreateDoubleMatrix(1,nsteps,mxREAL);
	memcpy((void *)mxGetPr(T), (void *)t, nsteps*sizeof(t));
	memcpy((void *)mxGetPr(X), (void *)x, nsteps*sizeof(x));
	engPutVariable(engp,"T",T);
	engPutVariable(engp,"X",X);
	engEvalString(engp,"si = stepinfo(X,T);");
	engEvalString(engp,"ts = si.SettlingTime;");
	engEvalString(engp,"os = si.Overshoot;");
	engEvalString(engp,"A = 5.0;");
	Ts = engGetVariable(engp,"ts");
	Os = engGetVariable(engp,"os");
	A = engGetVariable(engp,"A");
	ts = mxGetPr(Ts);
	os = mxGetPr(Os);
	a = mxGetPr(A);
	/*mexPrintf("%f\t%f\t%f\n",*ts,*os,*a);*/
	mxDestroyArray(T);
	mxDestroyArray(X);
	mxDestroyArray(Ts);
	mxDestroyArray(Os);
	return (*ts + *os / 20.0);
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
		int flag,
		Engine *engp)
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
			mexPrintf("error, return value=%d\n",status);
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
		u[0] = out[0] * 10 - 5;
		#endif
		#ifndef FUZZY
		kp = out[0] * 20; ki = out[1] * 5; kd = out[2] * 10;
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
		#endif
		u_t[i] = u[0];
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	return stepinfo(nsteps, tt, y0, flag, engp);
}

double f4_fitness(struct Fis * fis, Engine *engp) {
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
			mexPrintf("error, return value=%d\n",status);
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
		u[0] = out[0] * 10 - 5;
		#endif
		#ifndef FUZZY
		kp = out[0] * 20; ki = out[1] * 5; kd = out[2] * 10;
		u[0] = kp*ep[i] + ki*ei[i] + kd*ed[i];
		#endif
	}
	free(u);
	gsl_odeiv2_driver_free (d);
	return stepinfo(nsteps, tt, y0, -1, engp);
}

double PID_fitness(double kp, double ki, double kd,
		int (* fun)(double t, const double y[], double dydt[], void * params),
		double C[3],
		Engine *engp)
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
			mexPrintf("error, return value=%d\n",status);
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
	double cost = stepinfo(nsteps, tt, y0, 4,engp);
	return cost;
}


/*void mexFunction(int nlhs, mxArray *plhs[],
			int nrhs, const mxArray *prhs[]) {*/
int main()
{
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
	mexPrintf("I've Started!\n");
	hp->pop_size = 100;
	hp->elite = 0.05;
	hp->crossover = 0.5;
	hp->mutate = 0.25;
	hp->max_gen = 200;

	Engine *engp;
	if(!(engp = engOpen(""))) {
		mexPrintf("\nCan't start MATLAB engine!!\n");
		return EXIT_FAILURE;
	}

	FILE * fd = fopen("output.fis", "w");
	mexPrintf("Hopefully this isn't the last output you see from me...\n");
	struct Fis * bestfis = run_ga(spcs, hp, f4_fitness, NULL, engp);
	fclose(fd);

/** Plot the comparison graphs between PID- and Fuzzy-controlled systems**/
	double c1[3] = {14.575021682567217, 5.884648742411102, 0.443191673894189};
	double c2[3] = {14.575021682567217, 5.949696444058977, 0.455333911535126};
	double c3[3] = {113.8331814038286, 58.1700091157703, 0.7103463992707};
	double c4[3] = {74.626865671641795, 12.531572904707232, -0.017962112514351};


	char outdir[100];
	time_t T;
	time(&T);
	struct tm * timeinfo = localtime(&T);
	strftime(outdir,sizeof(outdir),"%H-%M-%S/",timeinfo);

	mexPrintf("F4 dyn\n");
	generic_fitness_comp(bestfis, f4_dyn, c1, 0,engp);
	PID_fitness(5.8, 0.47, 2.85,f4_dyn, c1,engp);
	mexPrintf("F4 deg50 dyn\n");
	generic_fitness_comp(bestfis, f4_deg_dyn, c2, 1,engp);
	PID_fitness(5.8, 0.47, 2.85,f4_deg_dyn, c2,engp);
	mexPrintf("F4 sub dyn\n");
	generic_fitness_comp(bestfis, f4_sub_dyn, c3, 2,engp);
	PID_fitness(5.8, 0.47, 2.85,f4_sub_dyn, c3,engp);
	mexPrintf("F4 sup dyn\n");
	generic_fitness_comp(bestfis, f4_sup_dyn, c4, 3,engp);
	PID_fitness(5.8, 0.47, 2.85,f4_sup_dyn, c4,engp);
	mexPrintf("Save plots and lookup tables? (Y/n):\n");

	/*Generate lookup tables for MATLAB*/
	#ifndef FUZZY
	double i, j, k;
	int in, out;
	char * filename = "lookup_table";
	char new_fn[100];
	char outinnum[13];
	double x[3], y[3];
	double sc[3] = {20,5,10};

/*	for (out = 0; out < num_out; out++) {
		for (k = 0.00; k < 1.01; k += 0.01) {
			x[2] = k;
			strcpy(new_fn,filename);
			sprintf(outinnum,"%d-%0.2f.csv",out, k);
			strcat(new_fn, outinnum);
			fd = fopen(new_fn, "w");
			for (j = 0; j < 1.01; j += 0.01) {
				x[1] = j;
				for (i = 0; i < 1.01; i += 0.01) {
					x[0] = i;
					evalfis(y,x,bestfis);
					mexPrintf("%f",sc[out] * y[out]);
					if (i < 1) {mexPrintf(", ");}
				}
				mexPrintf("\n");
			}
			fclose(fd);
		}
	}*/
	#endif

	engClose(engp);
	fis_destroy(bestfis);
	specs_clear(spcs);
	free(hp);
	return;
}

