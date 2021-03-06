#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#include "ga.h"
#include "fuzzy.h"
#include "debug.h"

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
//	printf("specs_set1: ");
//	ranges_print_test(spcs->num_params,ranges);
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
//	printf("specs_set2: ");
//	ranges_print_test(spcs->num_params,spcs->ranges);

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
    /*Everything assumes a normalized partition on [0,1]. This evenly divides
    th partition into num_sets divisions and places mfs in it with ramps at the ends.*/
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

void
rand_partition(double params[], int num_sets)
{
    /*Divides the partition the same as init_partition, but considers the even
    divisions as buckets within which to place the mf params randomly. This is
    coverage is guaranteed because there will always be overlap (at least at creation).*/
	int p;
	double r;
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
    /*The whole idea here is to mimic a rotary counting mechanism, but we can't
    make the assumption that the counter is in decimal. In fact, each digit may be
    in its own number base, so we have to keep track of each digit place and
    increment appropriately*/
    int in, rule;
    unsigned int flag = 0;
    int num_rule = 1;
    int width = num_in + num_out;
    //Initialize first rule to zeros and compute num_rule
    for (in = 0; in < num_in; in++) {
        rules[0 + in] = 0;
        num_rule *= in_mfs[in];
    }
    for (rule = 1; rule < num_rule; rule++) {
        if (rules[(rule-1) * width + 0] == (in_mfs[0] - 1)) {
            rules[rule * width + 0] = 0;
            flag ^= (1 << 1);//Tell next element to increment
        } else {
            rules[rule * width + 0] = rules[(rule-1) * width + 0] + 1;
        }
        for (in = 1; in < num_in; in++) {
            if (flag & (1 << in)) { //If this element's bit is set
                if (rules[(rule-1) * width + in] == (in_mfs[in]-1)) { //And previous element is max
                    rules[rule * width + in] = 0; //Start back at zero
                    flag ^= (1 << (in + 1)); //Set next element's bit
                } else {
                    rules[rule * width + in] = rules[(rule-1) * width + in] + 1;
                }
                flag ^= (1 << in);//Unset this element's bit regardless
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
//			printf("%d ",consequents[c]);
		}
//		printf("\n");
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
	int i = 2 + rand_i(num_params - 4);
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
//	printf("%d %d \n",i,j);
	for (p = 0; p <= j; p++) {
		child1[p] = (p <= i ? parent1[p] : parent2[p]);
		child2[p] = (p <= i ? parent2[p] : parent1[p]);
//		printf("|%d|%d",child1[p],child2[p]);
	}
	for (p = j; p < num_params; p++) {
		child1[p] = parent1[p];
		child2[p] = parent2[p];
//		printf("|%d|%d",child1[p],child2[p]);
	}
//	printf("\n");
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
	double alpha = 0;
	double lambda;
	int p;
	for (p = 2; p < num_params - 2; p++) {
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
	child1[0] = parent1[0];
	child1[1] = parent1[1];
	child2[0] = parent2[0];
	child2[1] = parent2[1];
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
//		printf("%d\t%d\n",in_mfs[in],mfs[in]);
	}
	for (out = 0; out < num_out; out++) {
		mfs[num_in + out] = out_mfs[out];
//		printf("%d\t%d\n",out_mfs[out],mfs[num_in+out]);
	}

	for (mf = 0; mf < num_in + num_out; mf++) {
		I = 1 / ((double) mfs[mf] - 1);
		for (mfp = 0 ; mfp < mfs[mf] * 3; mfp++, p++) {
//			printf("%d\t%d\t%d\t%d\t%d\t%d\n",p,mfp, mfp%3, mfp%3==1, mf, mfs[mf]*3);
			if ((mfp == 0) | (mfp == 1)
				| (mfp == mfs[mf] * 3 - 1)
				| (mfp == mfs[mf] * 3 - 2)) {
//				printf("mfp in [0,1,%d,%d]\n",mfs[mf]*3-2, mfs[mf]*3-1);
				ranges[p * 2 + 0] = -1;
				ranges[p * 2 + 1] = -1;
			} else if (mfp == 3) {
				ranges[p * 2 + 0] = 0;
				ranges[p * 2 + 1] = I / 2;
			} else if (mfp == (mfs[mf] * 3 - 4)) {
				ranges[p * 2 + 0] = 1 - I / 2;
				ranges[p * 2 + 1] = 1;
			} else if ((mfp % 3) == 0) { //a_i
//				printf("mfpMOD 3 is 0\n");
				ranges[p * 2 + 0] = ((mfp / 3) - 1) * I - I / 2;
				ranges[p * 2 + 1] = ((mfp / 3) - 1) * I + I / 2;
			} else if ((mfp % 3) == 1) { //b_i
//				printf("mfp MOD 3 is 1\n");
				ranges[p * 2 + 0] = (((mfp - 1) / 3)) * I - I / 2;
				ranges[p * 2 + 1] = (((mfp - 1) / 3)) * I + I / 2;
			} else { //c_i
//				printf("mfp MOD 3 is 2\n");
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
//		printf("r_mutation: range = [%f, %f]\n",ranges[p][0],ranges[p][1]);
	}
}

void
rb_mutation(
	int num_params,
//	double ** ranges,
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
//				printf("UP\n");
			}
		} else {
			consequents[gen] -= 1;
//			printf("DOWN\n");
		}
	}
}

void
individual_mutate(
	struct Individual *ind,
	int num_params,
//	double ** ranges,
	double ranges[num_params * 2],
	int cur_gen,
	int max_gen,
	double b,
	int num_rules,
	int num_out,
	int out_mfs[],
	int num_genes)
{
//	ranges_print_test(num_params,ranges);
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
//	init_params(tmp_params, num_in, in_mfs, num_out, out_mfs);
//	rand_consequents(num_out, num_rules, tmp_consequents, out_mfs);
//	population[0] = individual_create(num_params, tmp_params, num_rules, tmp_consequents);
	for (i = 0; i < pop_size; i++) {
		rand_params(tmp_params, num_in, in_mfs, num_out, out_mfs);
		rand_consequents(num_out, num_rules, tmp_consequents, out_mfs);
		population[i] = individual_create(num_params, tmp_params, num_rules, num_out, tmp_consequents);
	}
}


void
population_iter_cascade(
        int num_pop,
        struct Individual ** pop_now[num_pop],
        struct Individual ** pop_next[num_pop],
        int rank[],
        int cur_gen,
        struct HyperParams * hp,
        struct Specs ** spcs)
{
	int pop, ind;
	for (pop = 0; pop < num_pop; pop++) {
		population_iter(
				pop_now[pop],
				pop_next[pop],
				rank,
				cur_gen,
				hp,
				spcs[pop]);
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
	double * fit_min)
{
	int ind;
	double fitness[pop_size];
	double * ranked_fitness[pop_size];
	double tmp_fit;
	struct Fis * tmp_fis;
	struct Specs * tmp_spcs;
	// Parallelize with omp directives
	#pragma omp parallel private(tmp_fis, tmp_spcs)
	{
		#pragma omp for
		for (ind = 0; ind < pop_size; ind++) {
			tmp_spcs = specs_copy(spcs);
			tmp_fis = individual_to_fis(population[ind],tmp_spcs);
			tmp_fit = fit_fcn(tmp_fis);
			fitness[ind] = tmp_fit;
			ranked_fitness[ind] = &fitness[ind];
			fis_destroy(tmp_fis);
			specs_clear(tmp_spcs);
		}
	} //end omp parallelization
	*fit_min = minimum(pop_size, fitness);
	qsort(ranked_fitness, pop_size, sizeof(double *), cmpdouble_p);
	for (ind = 0; ind < pop_size; ind++) {
		rank[ind] = (int)(ranked_fitness[ind] - &fitness[0]);
	}
}

void
population_rank_cascade(
	int num_pop,
    int pop_size,
	int rank[pop_size],
	struct Individual ** populations[num_pop],
	struct Specs * spcs[num_pop],
	fitness_fcn_cascade fit_fcn,
	double *fit_min)
{
    int ind, pop;
    double fitness[pop_size];
    double * ranked_fitness[pop_size];
    double tmp_fit;
    struct Fis * tmp_fis[num_pop];
    struct Specs * tmp_spcs[num_pop];
    // Parallelize with omp directives
    #pragma omp parallel private(tmp_fis, tmp_spcs, pop)
    {
        #pragma omp for
        for (ind = 0; ind < pop_size; ind++) {
            for (pop = 0; pop < num_pop; pop++) {
                tmp_spcs[pop] = specs_copy(spcs[pop]);
                tmp_fis[pop] = individual_to_fis(populations[pop][ind],tmp_spcs[pop]);
            }
            tmp_fit = fit_fcn(num_pop, tmp_fis);
            fitness[ind] = tmp_fit;
            ranked_fitness[ind] = &fitness[ind];
            for (pop = 0; pop < num_pop; pop++) {
                fis_destroy(tmp_fis[pop]);
                specs_clear(tmp_spcs[pop]);
            }
        }
    } //end omp parallelization
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

void
population_switch_cascade(
        struct Individual **** pop1s,
        struct Individual **** pop2s)
{
	struct Individual *** tmp = *pop1s;
	*pop1s = *pop2s;
	*pop2s = tmp;
}

struct Fis *
run_ga(
	struct Specs * spcs,
	struct HyperParams * hp,
	fitness_fcn fit_fcn,
	FILE * fis_log)
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
		population_rank(hp->pop_size, rank, pop1, spcs, fit_fcn, &fitness_hist[gen]);
		printf("Ind[%d] Fitness: %f\n",rank[0],fitness_hist[gen]);
        int stag =50;
		if ( (gen > stag) && (fabs(sum_d(&fitness_hist[gen - stag], stag)/(double)stag - fitness_hist[gen]) < 1e-7)) {
			struct Fis * ret_fis =  individual_to_fis(pop1[rank[0]],spcs);
			printf("Best fitness:%f\n",fit_fcn(ret_fis));
			individual_print(pop1[rank[0]], spcs, fis_log);
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
    printf("Best fitness:%f\n",fit_fcn(ret_fis));
	printf("%d Generations\n",gen);
	individuals_destroy(pop1, hp->pop_size);
	individuals_destroy(pop2, hp->pop_size);
	free(pop1);
	free(pop2);
	return ret_fis;
}

void
run_cascade_ga(
	int num_pop,
	struct Fis * fis_list[num_pop],
	struct Specs * spcs[num_pop],
	struct HyperParams * hp,
	fitness_fcn_cascade fit_fcn,
	FILE * fis_log)
{
	int pop, poptmp;
	struct Individual *** pop1s = malloc(num_pop * sizeof(struct Individual **));
	struct Individual *** pop2s = malloc(num_pop * sizeof(struct Individual **));
	double fitness_hist[hp->max_gen];
	int rank[hp->pop_size];
	int gen;
	for (pop = 0; pop < num_pop; pop++) {
		pop1s[pop] = malloc(hp->pop_size * sizeof(struct Individual *));
		pop2s[pop] = malloc(hp->pop_size * sizeof(struct Individual *));
		population_init(
			hp->pop_size,
			pop1s[pop],
			spcs[pop]->num_in,
			spcs[pop]->in_mfs,
			spcs[pop]->num_out,
			spcs[pop]->out_mfs,
			spcs[pop]->num_params,
			spcs[pop]->num_rule,
			spcs[pop]->rules);
		population_init(
			hp->pop_size,
			pop2s[pop],
			spcs[pop]->num_in,
			spcs[pop]->in_mfs,
			spcs[pop]->num_out,
			spcs[pop]->out_mfs,
			spcs[pop]->num_params,
			spcs[pop]->num_rule,
			spcs[pop]->rules);
	}

	for (gen = 0; gen < hp->max_gen; gen++) {
		population_rank_cascade(num_pop, hp->pop_size, rank, pop1s, spcs, fit_fcn, &fitness_hist[gen]);
		printf("Ind[%d] Fitness: %f\n",rank[0],fitness_hist[gen]);
//		individual_print(pop1[rank[0]], ga_log);
		if ( (gen > 10) && (fabs(sum_d(&fitness_hist[gen - 10], 10)/10.0 - fitness_hist[gen]) < 1e-17) ) {
//			printf("Best fitness:%f\n",fit_fcn(ret_fis));
			for (poptmp = 0; poptmp < num_pop; poptmp++) {
				fis_list[poptmp] = individual_to_fis(pop1s[poptmp][rank[0]],spcs[poptmp]);
				individuals_destroy(pop1s[poptmp], hp->pop_size);
				individuals_destroy(pop2s[poptmp], hp->pop_size);
				free(pop1s[poptmp]);
				free(pop2s[poptmp]);
			}
			free(pop1s);
			free(pop2s);
			printf("%d Generations\n",gen + 1);
			return;
		} else if (gen == hp->max_gen - 1) {
			printf("%d Generations\n",gen + 1);
			break;
		}
		population_iter_cascade(num_pop, pop1s, pop2s, rank, gen, hp, spcs);
		population_switch_cascade(&pop1s, &pop2s);
	}
	for (poptmp = 0; poptmp < num_pop; poptmp++) {
		fis_list[poptmp] = individual_to_fis(pop1s[poptmp][rank[0]],spcs[poptmp]);
		individuals_destroy(pop1s[poptmp], hp->pop_size);
		individuals_destroy(pop2s[poptmp], hp->pop_size);
		free(pop1s[poptmp]);
		free(pop2s[poptmp]);
	}
	free(pop1s);
	free(pop2s);
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

	fprintf(fd, "num_in: {%d}\nin_mfs: {",spcs->num_in);
	for (in = 0; in < spcs->num_in; in++) {
		fprintf(fd, "%d ", spcs->in_mfs[in]);
	}

	fprintf(fd, "}\nnum_out: {%d}\nout_mfs: {",spcs->num_out);
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
