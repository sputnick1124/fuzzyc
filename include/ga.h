#ifndef _GALIB_H_
#define _GALIB_H_

int sum_i(int adds[], size_t len);

int prod_i(int mults[], size_t len);

int rand_i(unsigned int max);

void init_partition(double params[], int num_sets);

void rand_partition(double params[], int num_sets);

void init_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]);

void rand_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]);

void init_antecedents(int num_in, int num_out, int rules[][num_in], int in_mfs[]);

void rand_consequents(int num_out, int num_rule, int consequents[num_rule * num_out], int out_mfs[]);

void add_consequents(
	int num_in,
	int num_out,
	int num_rule,
	int rules[num_rule][num_in + num_out],
	int consequents[num_rule * num_out]);

void sp_crossover(int num_params, double child1[], double child2[], double parent1[], double parent2[]);

void tp_crossover(int num_params, double child1[], double child2[], double parent1[], double parent2[]);

void blx_a_crossover(int num_params, double child1[], double child2[], double parent1[], double parent2[]);

void param_range(
	int num_params,
	double ranges[num_params][2],
	int num_in,
	int in_mfs[],
	int num_out,
	int out_mfs[]);

void r_mutation(
	int num_params,
	double ranges[num_params][2],
	double chromosome[num_params],
	int num_genes);

void rb_mutation(
	int num_params,
	double ranges[num_params][2],
	double chromosome[num_params],
	int cur_gen,
	int max_gen,
	double b,
	int num_genes);

struct Individual {
	double * params;
	double * consequents;
};

struct Individual *individual_create(
	int num_params,
	double params[],
	int num_rule,
	int consequents[]);

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

void individual_destroy(struct Individual *ind);


#endif
