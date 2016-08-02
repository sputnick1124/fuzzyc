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
			params[p] = ((p-1) / 3) * I;
		} else { //c_i
			params[p] = (((p-2) / 3) + 1) * I;
		}
	}
}

void rand_partition(double params[], int num_sets) {
	int p;
	double r;
	double tmp;
	int num_params = num_sets * 3;
	srand48((long int)clock());
	params[0] = 0; params[1] = 0; //Set lower bounds
	params[num_params-1] = 1;	//Set upper bounds
	params[num_params-2] = 1;
	double I = 1 / ((double) num_sets - 1); //Interval to divide up the fuzzy partition

	for (p = 2; p < num_params - 2; p++) {
		r = drand48() * I; //random shift from partition division
		printf("I = %0.2f\t r = %0.2f\n",I, r);
		if (p == 3) {
			params[p] = r; //Don't want the first element slipping below 0...
			printf("p = r = %f\n",r);
		} else if ( p == num_params - 4) {
			params[p] = (((p - 2) / 3) + 1) * I - r; //...or the last past 1
		} else if (p % 3 == 0) { //a_i
			params[p] = ((p / 3) - 1) * I + r - (I / 2);
		} else if (p % 3 == 1) { //b_i
			params[p] = ((p-1) / 3) * I + r - (I / 2);
		} else { //c_i
			params[p] = (((p-2) / 3) + 1) * I + r - (I / 2);
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

/** Crossover and mutation functions **/
void
sp_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[])
{
/* Single-point crossover
