#ifndef _GALIB_H_
#define _GALIB_H_

struct Individual {
	double * params;
	int * consequents;
};

struct Individual *individual_create(
	int num_params,
	double params[],
	int num_rule,
	int consequents[]);

void individual_copy(
	struct Individual * ind1,
	struct Individual * ind2,
	int num_params,
	int num_rule);

void individual_destroy(struct Individual * ind);

void individuals_destroy(struct Individual ** ind, int num_ind);

struct Specs {
	int num_in;
	int * in_mfs;
	int num_out;
	int * out_mfs;
	int num_rule;
	int ** rules;
	int num_params;
	double ** ranges;
};

struct Specs *specs_set(
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[],
	int rules[][num_in + num_out],
	double ranges[][2]);

struct Specs *specs_copy(struct Specs *spcs);

void specs_clear(struct Specs *spcs);

int sum_i(int adds[], size_t len);

int prod_i(int mults[], size_t len);

int rand_i(unsigned int max);

void init_partition(double params[], int num_sets);

void rand_partition(double params[], int num_sets);

void init_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]);

void rand_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]);

void init_antecedents(int num_in, int num_out, int rules[][num_in], int in_mfs[]);

void rand_consequents(int num_out, int num_rule, int consequents[num_rule * num_out], int out_mfs[]);

void
add_consequents(
	int num_in,
	int num_out,
	int num_rule,
	int rules[num_rule][num_in + num_out],
	int consequents[num_rule * num_out]);

void
sp_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[]);

void
tp_crossover(
	int num_params,
	double child1[],
	double child2[],
	double parent1[],
	double parent2[]);

void
consequent_crossover(
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
	double ranges[num_params][2],
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[]);

void
r_mutation(
	int num_params,
	double ranges[num_params][2],
	double chromosome[num_params],
	int num_genes);

void
rb_mutation(
	int num_params,
//	double **ranges,
	double ranges[num_params][2],
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
//    double **ranges,
    double ranges[][2],
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
    int rules[][num_in + num_out]);

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
    struct Specs * spcs);

#endif
