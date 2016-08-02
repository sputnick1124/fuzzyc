#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../include/ga.h"
#include "../include/fuzzy.h"

/** Some utility functions to take care of the mechanics **/
int prod_i(int mults[], size_t len) {
	int retval = 1;
	int i;
	for (i = 0; i < len; i++) {
		retval *= mults[i];
	}
	return retval;
}

int sum_i(int adds[], size_t len) {
	int retval = 0;
	int i;
	for (i = 0; i < len; i++) {
		retval += adds[i];
	}
	return retval;
}

int rand_i(unsigned int max) {
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

/** Initialization functions for starting off the GA population**/
void init_partition(double params[], int num_sets) {
	int p;
	int num_params = num_sets * 3;
	params[0] = 0; params[1] = 0; //Set lower bounds
	params[num_params-1] = 1;	//Set upper bounds
	params[num_params-2] = 1;
	double I = 1 / ((double) num_sets - 1); //Interval to divide up the fuzzy partition
	for (p = 2; p < num_params - 2; p++) {
		if (p % 3 == 0) { //a_i
			params[p] = ((p / 3) - 1) * I;
		} else if (p % 3 == 1) { //b_i
			params[p] = ((p - 1) / 3) * I;
		} else { //c_i
			params[p] = (((p - 2) / 3) + 1) * I;
		}
	}
}

void rand_partition(double params[], int num_sets) {
	int p;
	double r;
	double tmp;
	int num_params = num_sets * 3;
	params[0] = 0; params[1] = 0; //Set lower bounds
	params[num_params-1] = 1;	//Set upper bounds
	params[num_params-2] = 1;
	double I = 1 / ((double) num_sets - 1); //Interval to divide up the fuzzy partition

	for (p = 2; p < num_params - 2; p++) {
		r = drand48() * I; //random shift from partition division
		if (p == 3) {
			params[p] = r / 2; //Don't want the first element slipping below 0...
		} else if ( p == num_params - 4) {
			params[p] = (((p - 2) / 3) + 1) * I - r / 2; //...or the last past 1
		} else if (p % 3 == 0) { //a_i
			params[p] = ((p / 3) - 1) * I + r - (I / 2);
		} else if (p % 3 == 1) { //b_i
			params[p] = ((p - 1) / 3) * I + r - (I / 2);
		} else { //c_i
			params[p] = (((p - 2) / 3) + 1) * I + r - (I / 2);
		}
	}
}

void init_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]) {
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

void rand_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]) {
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

void init_antecedents(int num_in,
	int num_out,
	int rules[][num_in + num_out],
	int in_mfs[]) {

	int in, rule;
	unsigned int flag = 0;
	int num_rule = 1;
	//Initialize first rule to zeros and compute num_rule
	for (in = 0; in < num_in; in++) {
		rules[0][in] = 0;
		num_rule *= in_mfs[in];
	}
	for (rule = 1; rule < num_rule; rule++) {
		if (rules[rule-1][0] == (in_mfs[0] - 1)) {
			rules[rule][0] = 0;
			flag ^= (1 << 1);//Tell next element to increment
		} else {
			rules[rule][0] = rules[rule-1][0] + 1;
		}
		for (in = 1; in < num_in; in++) {
			if (flag & (1 << in)) { //If this element's bit is set
				if (rules[rule-1][in] == (in_mfs[in]-1)) { //And previous element is max
					rules[rule][in] = 0; //Start back at zero
					flag ^= (1 << (in + 1)); //Set next element's bit
				} else {
					rules[rule][in] = rules[rule-1][in] + 1;
				}
				flag ^= (1 << in);//Unset this element's bit regardless
			} else {
				rules[rule][in] = rules[rule-1][in];
			}
		}
	}
}

void rand_consequents(
	int num_in,
	int num_out,
	int num_rule,
	int rules[][num_in + num_out],
	int out_mfs[]) {

	int out, rule;
	for (rule = 0; rule < num_rule; rule++) {
		for (out = 0; out < num_out; out++) {
			rules[rule][num_in + out] = rand_i(out_mfs[out]);
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
blx_a_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[])
{
/* Blended crossover with alpha exploration factor*/
	double I;
	double x1, x2;
	double xmin, xmax;
	double alpha = 0.5;
	int p;
	for (p = 0; p < num_params - 2; p++) {
		xmin = fmin(parent1[p], parent2[p]);
		xmax = fmax(parent1[p], parent2[p]);
		I = xmax - xmin;
		child1[p] = xmin - I * alpha + drand48() * I * 2 * alpha;
		child2[p] = xmin - I * alpha + drand48() * I * 2 * alpha;
	}
	child1[num_params-2] = parent1[num_params-2];
	child1[num_params-1] = parent1[num_params-1];
	child2[num_params-2] = parent2[num_params-2];
	child2[num_params-1] = parent2[num_params-1];
}

/** Mutation functions **/
void
param_range(
	int num_params,
	double ranges[num_params][2],
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
//		printf("%d\t%d\n",in_mfs[in],mfs[in]);
	}
	for (out = 0; out < num_out; out++) {
		mfs[num_in-1 + out] = out_mfs[out];
//		printf("%d\t%d\n",out_mfs[out],mfs[num_in-1+out]);
	}

	for (mf = 0; mf < num_in + num_out; mf++) {
		for (mfp = 0 ; mfp < mfs[mf] * 3 - 2; p++, mfp++) {
			I = 1 / ((double) mfs[mf] - 1);
			ranges[p][0] = 0;
			ranges[p][0] = I / 2;
			if (p % 3 == 0) { //a_i
				ranges[p][0] = ((p / 3) - 1) * I - I / 2;
				ranges[p][1] = ((p / 3) - 1) * I + I / 2;
			} else if (p % 3 == 1) { //b_i
				ranges[p][0] = (((p - 1) / 3)) * I - I / 2;
				ranges[p][1] = (((p - 1) / 3)) * I + I / 2;
			} else { //c_i
				ranges[p][0] = (((p - 2) / 3) + 1) * I - I / 2;
				ranges[p][1] = (((p - 2) / 3) + 1) * I + I / 2;
			}
		}
	}
}

