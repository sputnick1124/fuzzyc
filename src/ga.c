#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "../include/ga.h"
#include "../include/fuzzy.h"


void ranges_print_test(int num_params, double * ranges) {
    int p, limit;
	for (p = 0; p < num_params; p++) {
		printf("%d: ",p);
		for (limit = 0; limit < 2; limit++) {
	        printf("%0.3f\t", ranges[p * 2 + limit]);
		}
		printf("\n");
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
    int retval;
    double u = drand48();

	/* implicit cast to int is intentional*/
    return max - sqrt((1 - u) * max * max);
}

struct Specs *
specs_set(
	int num_in,
    int in_mfs[],
    int num_out,
    int out_mfs[],
	int rules[],
	double ranges[])
{
	int in, out, rule, ant, range, limit;
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
	int in, out, rule, ant, range, limit;
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

	for (in = 0; in < spcs->num_in; in++) {
		spcs->in_mfs[in] = oldspcs->in_mfs[in];
	}
	for (out = 0; out < spcs->num_out; out++) {
		spcs->out_mfs[out] = oldspcs->out_mfs[out];
	}

	for (rule = 0; rule < spcs->num_rule; rule++) {
		for (ant = 0; ant < spcs->num_in; ant++) {
			spcs->rules[rule * spcs->num_in + ant] = oldspcs->rules[rule * spcs->num_in + ant];
		}
	}

	for (range = 0; range < spcs->num_params; range++) {
		for (limit = 0; limit < 2; limit++) {
			spcs->ranges[range * 2 + limit] = oldspcs->ranges[range * 2 + limit];
		}
	}

	return spcs;
}

void
specs_clear(struct Specs *spcs)
{
	int range, rule;
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

	ind->consequents = malloc(num_rule * sizeof(int));
	if (ind->consequents == NULL) {
		free(ind->params);
		free(ind);
		return NULL;
	}
	for (r = 0; r < num_rule; r++) {
		ind->consequents[r] = consequents[r];
	}

	return ind;
}

void
individual_copy(
	struct Individual * ind1,
	struct Individual * ind2,
	int num_params,
	int num_rule)
{
	int rule, p;
	for (p = 0; p < num_params; p++) {
		ind2->params[p] = ind1->params[p];
	}
	for (rule = 0; rule < num_rule; rule++) {
		ind2->consequents[rule] = ind1->consequents[rule];
	}
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
	int out, rule, c;
	c = 0;
	for (out = 0; out < num_out; out++) {
		for (rule = 0; rule < num_rule; c++, rule++) {
			consequents[c] = rand_i(out_mfs[out]);
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
	int c = 0;
	int width = num_in + num_out;
	for (rule = 0; rule < num_rule; c++, rule++) {
		for (out = 0; out < num_out; out++) {
			rules[rule * width + num_in + out] = consequents[c];
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
consequent_crossover(
	int num_params,
	int child1[],
	int child2[],
	int parent1[],
	int parent2[])
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
	consequent_crossover(
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
			if (mfp == 0 | mfp == 1
				| mfp == mfs[mf] * 3 - 1
				| mfp == mfs[mf] * 3 - 2) {
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
	int con, out, gen, mf, i;
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

	for (i = 0; i < pop_size; i++) {
		rand_params(tmp_params, num_in, in_mfs, num_out, out_mfs);
		rand_consequents(num_out, num_rules, tmp_consequents, out_mfs);
		population[i] = individual_create(num_params, tmp_params, num_rules, tmp_consequents);
	}
}

void
population_iter(
	struct Individual ** pop_now,
	struct Individual ** pop_next,
	int pop_size,
	int rank[],
	float elite,
	float crossover,
	float mutate,
	int cur_gen,
	int max_gen,
	struct Specs * spcs)
{
	int ind, select1, select2;
	int num_elite = pop_size * elite;
	int num_cross = pop_size * crossover;
	int num_mut = pop_size * mutate;
	int num_rand = pop_size - (num_elite + num_cross + num_mut);
	for (ind = 0; ind < pop_size; ind++) {
		if (ind < num_elite) {
			individual_copy(
				pop_now[rank[ind]],
				pop_next[ind],
				spcs->num_params,
				spcs->num_rule);
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
					max_gen,
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
					max_gen,
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
