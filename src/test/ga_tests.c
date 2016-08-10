#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../../include/ga.h"
#include "../../include/fuzzy.h"

int test_chromo(int num_params, double chromo[]) {
	int p;
	for (p = 2; p < num_params; p += 3) {
		if (chromo[p-2] > chromo[p-1]) {
			return 0;
		} else if ( chromo[p-1] > chromo[p]) {
			return 0;
		} else {
			return 1;
		}
	}
}

int test_consequents(int num_out, int num_rule, int consequents[], int out_mfs[]) {
	int out, rule;
	int c = 0;
	for (out = 0; out < num_out; out++) {
		for (rule = 0; rule < num_rule; c++, rule++) {
			if (consequents[c] >= out_mfs[out]) {
				return 0;
			}
		}
	}
	return 1;
}

int individual_test(struct Individual *ind, int num_params, int num_out, int num_rule, int out_mfs[]) {
	int r1, r2;
	r1 = test_chromo(num_params, ind->params);
	r2 = test_consequents(num_out, num_rule, ind->consequents, out_mfs);
	return (r1 + r2) == 2;
}

void chromosome_print(int num_params, double chromo[]) {
	int p;
	printf("[");
	for (p = 0; p < num_params; p++) {
		if (p % 3 == 0) {
			printf("|");
		}
		printf("%f ",chromo[p]);
	}
	printf("]\n");
}


void cons_print(int num_rules, int consequents[]) {
	int r;
	for (r = 0; r < num_rules; r++) {
		printf("| %d |", consequents[r]);
	}
}

double
fit_line(struct Fis * fis)
{
	int i;
	double dx = 0.01;
	int max = (int) (1.0 / dx);
	double x[1];
	double out[1];
	double cost = 0;
	for (i = 0; i < max; i++) {
		x[0] = (double)i * dx;
		evalfis(out,x,fis);
		cost += dx * (out[0] - (double)i) * (out[0] - (double)i);
	}
	return cost;
}

int main(int argc, char * argv[]) {
	int num_in = 3;
	int in_mfs[] = {3, 4, 2};
	int num_out = 1;
	int out_mfs[] = {3};
	int num_params = 3 * (sum_i(in_mfs, num_in) + sum_i(out_mfs,num_out));
	int num_rules = prod_i(in_mfs,num_in);
	int rules[num_rules * (num_in + num_out)];
	int consequents[num_rules * num_out];
	double params[num_params];
	double testp1[num_params];
	double testp2[num_params];
	double testc1[num_params];
	double testc2[num_params];
	double testc3[num_params];
	double testc4[num_params];
	double ranges[num_params * 2];
	int ind;
	int in, test;
	int fails = 0;
	int test_num = 1000000;

	srand48((long int) time(NULL));
	srand((long int) time(NULL));
//	init_partition(params, in_mfs[0]);
//	init_partition(&params[in_mfs[0]*3], in_mfs[1]);
//	init_partition(&params[sum_i(in_mfs,2)*3], in_mfs[2]);
	rand_params(params, num_in, in_mfs, num_out, out_mfs);

	init_antecedents(num_in, num_out, rules, in_mfs);
	rand_consequents(num_out, num_rules, consequents, out_mfs);
	add_consequents(num_in, num_out, num_rules, rules, consequents);

	param_range(num_params, ranges, num_in, in_mfs, num_out, out_mfs);

	/*Test that rand_params returns a valid list of params*/
	printf("Testing rand_params\n");
	for (test = 0; test < test_num; test++) {
		rand_params(params, num_in, in_mfs, num_out, out_mfs);
		if (!test_chromo(num_params, params)) {
			printf("Test %d failed:\n",test);
			chromosome_print(num_params, params);
			fails++;
			//return 1;
		}
	}
	printf("%d tests failed\n",fails);

	/*Test that sp_crossover returns a valid list of params*/
	printf("Testing Single-point crossover\n");
	fails = 0;
	for (test = 0; test < test_num; test++) {
		rand_params(testp1, num_in, in_mfs, num_out, out_mfs);
		rand_params(testp2, num_in, in_mfs, num_out, out_mfs);
		sp_crossover(num_params, testc1, testc2, testp1, testp2);
		if (!(test_chromo(num_params, testc1) & test_chromo(num_params, testc2))) {
			printf("Test %d failed:\n",test);
			chromosome_print(num_params, testp1);
			chromosome_print(num_params, testp2);
			chromosome_print(num_params, testc1);
			chromosome_print(num_params, testc2);
			return 1;
			fails++;
		}
	}
	printf("%d tests failed\n",fails);

	/*Test that tp_crossover returns a valid list of params*/
	printf("Testing Two-point crossover\n");
	fails = 0;
	for (test = 0; test < test_num; test++) {
		rand_params(testp1, num_in, in_mfs, num_out, out_mfs);
		rand_params(testp2, num_in, in_mfs, num_out, out_mfs);
		tp_crossover(num_params, testc1, testc2, testp1, testp2);
		if (!(test_chromo(num_params, testc1) & test_chromo(num_params, testc2))) {
			/*chromosome_print(num_params, testc1);
			chromosome_print(num_params, testc2);*/
			//return 1;
			fails++;
		}
	}
	printf("%d tests failed\n",fails);

	/*Test that blx_a_crossover returns a valid list of params*/
	printf("Testing BLX-alpha crossover\n");
	fails = 0;
	for (test = 0; test < test_num; test++) {
		rand_params(testp1, num_in, in_mfs, num_out, out_mfs);
		rand_params(testp2, num_in, in_mfs, num_out, out_mfs);
		blx_a_crossover(num_params, testc1, testc2, testp1, testp2);
		if (!(test_chromo(num_params, testc1) & test_chromo(num_params, testc2))) {
			/*printf("Test %d failed:\n",test);
			chromosome_print(num_params, testc1);
			chromosome_print(num_params, testc2);*/
			//return 1;
			fails++;
		}
	}
	printf("%d tests failed\n",fails);

	/*Test that repeated blx_a_crossover returns a valid list of params*/
	printf("Testing recursive BLX-alpha crossover\n");
	fails = 0;
	rand_params(testp1, num_in, in_mfs, num_out, out_mfs);
	rand_params(testp2, num_in, in_mfs, num_out, out_mfs);
	blx_a_crossover(num_params, testc1, testc2, testp1, testp2);
	for (test = 0; test < test_num; test++) {
		blx_a_crossover(num_params, testc3, testc4, testc1, testc2);
		if (!(test_chromo(num_params, testc1) & test_chromo(num_params, testc2))) {
			/*printf("Test %d failed:\n",test);
			chromosome_print(num_params, testc1);
			chromosome_print(num_params, testc2);*/
			//return 1;
			fails++;
		}
		*testc1 = *testc3;
		*testc2 = *testc4;
	}
	printf("%d tests failed\n",fails);


	/*Test that mutations result in valid chromosomes*/
	printf("Testing random mutation\n");
	fails = 0;
	for (test = 0; test < test_num; test++) {
		rand_params(testp1, num_in, in_mfs, num_out, out_mfs);
		r_mutation(num_params, ranges, testp1, 5);
		if (!test_chromo(num_params, testp1)) {
			/*printf("Test %d failed:\n",test);
			chromosome_print(num_params, testc1);
			chromosome_print(num_params, testc2);*/
			//return 1;
			fails++;
		}
		*testc1 = *testc3;
		*testc2 = *testc4;
	}
	printf("%d tests failed\n",fails);

	printf("Testing non-uniform random mutation\n");
	fails = 0;
	for (test = 0; test < test_num; test++) {
		rand_params(testp1, num_in, in_mfs, num_out, out_mfs);
		rb_mutation(num_params, ranges, testp1, test, test_num, 1.5, 5);
		if (!test_chromo(num_params, testp1)) {
			/*printf("Test %d failed:\n",test);
			chromosome_print(num_params, testc1);
			chromosome_print(num_params, testc2);*/
			//return 1;
			fails++;
		}
		*testc1 = *testc3;
		*testc2 = *testc4;
	}
	printf("%d tests failed\n",fails);

	/*Test the creation of a population and the validity of its individuals*/
	int pop_size = 100;
	struct Individual * population[pop_size];
	struct Individual * tmp_pop[pop_size];
	population_init(pop_size, population, num_in, in_mfs, num_out, out_mfs, num_params, num_rules, rules);
	population_init(pop_size, tmp_pop, num_in, in_mfs, num_out, out_mfs, num_params, num_rules, rules);

	printf("Testing the population initialization\n");
	fails = 0;
	for (test = 0; test < pop_size; test++) {
		if (!individual_test(population[test], num_params, num_out, num_rules, out_mfs)) {
			fails++;
		}
	}
	printf("%d individuals are invalid\n",fails);


	printf("Testing individual crossover\n");
	fails = 0;
	int ind1, ind2;
	for (test = 0; test < test_num; test++) {
		ind1 = rand_i(pop_size);
		ind2 = rand_i(pop_size);
		individual_crossover(num_params,
			num_rules,
			num_out,
			population[ind1],
			population[ind2],
			tmp_pop[ind1],
			tmp_pop[ind2]);
		if (!individual_test(tmp_pop[ind1], num_params, num_out, num_rules, out_mfs)
			| !individual_test(tmp_pop[ind1], num_params, num_out, num_rules, out_mfs)) {
			fails++;
		}
	}
	printf("%d tests failed\n",fails);

	struct Specs *spcs =  specs_set(num_in, in_mfs, num_out, out_mfs);

	printf("Testing individual mutation\n");
//	printf("ga_tests.c: num_params = %d\n",num_params);
	fails = 0;
	for (test = 0; test < test_num; test++) {
		ind1 = rand_i(pop_size);
		individual_mutate(
			population[ind1],
			num_params,
			ranges,
			test,
			1000000,
			1.5,
			num_rules,
			num_out,
			out_mfs,
			5);
		if (!individual_test(population[ind1], num_params, num_out, num_rules, out_mfs)) {
			fails++;
		}
	}
	printf("%d tests failed\n",fails);
//	ranges_print(num_params,ranges);
//	ranges_print(spcs->num_params,spcs->ranges);

//	struct Specs *spcs =  specs_set(num_in, in_mfs, num_out, out_mfs, rules,ranges);
	int rank[pop_size];
	for (ind = 0; ind < pop_size; ind++) {
		rank[ind] = pop_size - ind - 1;
//		printf("rank[%d] = %d\n",ind,rank[ind]);
	}

	struct HyperParams * hp = malloc(sizeof(struct HyperParams));
	hp->pop_size = pop_size;
	hp->elite = 0.05;
	hp->crossover = 0.5;
	hp->mutate = 0.25;
	hp->max_gen = 100;
	population_iter(
		population,
		tmp_pop,
		rank,
		100,
		hp,
		spcs);

//	printf("Ranges:\n");
//	ranges_print(num_params,(double (*)[2])spcs->ranges);
//	ranges_print(num_params,ranges);

	printf("Testing generation iteration\n");
	fails = 0;
	for (ind = 0; ind < pop_size; ind++) {
		if (!individual_test(population[ind],num_params, num_out, num_rules, out_mfs)) {
			fails++;
			printf("ind: pop: %d\n",ind);
			chromosome_print(num_params, population[ind]->params);
			printf("\n\n");
		}
		if (!individual_test(tmp_pop[ind],num_params, num_out, num_rules, out_mfs)) {
			fails++;
			printf("ind: pop_new: %d\n",ind);
			chromosome_print(num_params, tmp_pop[ind]->params);
			printf("\n\n");
		}
	}
	printf("%d tests failed\n",fails);


	printf("Testing conversion to fis\n");
	fails = 0;
	individuals_destroy(population, pop_size);
	int o;
	struct Fis * tmp_fis;
	double out[num_out];
	double x[num_in];
	for (in = 0; in < num_in; in++) {x[in] = drand48();}
	for (test = 0; test < test_num / pop_size; test++) {
		population_init(pop_size, population, num_in, in_mfs, num_out, out_mfs, num_params, num_rules, rules);
		for (ind = 0; ind < pop_size; ind++) {
			tmp_fis = individual_to_fis(population[ind],spcs);
			evalfis(out, x, tmp_fis);
			for (o = 0; o < num_out; o++) {
				if (!isfinite(out[o]) | out[o] < 0 | out[o] > 1) {
					fails++;
					printf("test: %d, ind: %d, out[%d] = %f\n",test,ind,o,out[o]);
				}
			}
			fis_destroy(tmp_fis);
		}
		individuals_destroy(population, pop_size);
	}
	printf("%d tests failed\n",fails);

	specs_clear(spcs);
//	individuals_destroy(population, pop_size);
	individuals_destroy(tmp_pop, pop_size);

	int imf[1] = {3};
	int omf[1] = {3};

	spcs = specs_set(1,imf,1,omf);
//	hp->max_gen = 5;
	struct Fis * bestfis = run_ga(spcs,hp,fit_line);

	specs_clear(spcs);
	free(hp);

	double xx;
	double x_vals[1], out_vals[1];
	double y[100];
	int c = 0;
	for (xx = 0.0; xx < 1.01; c++, xx += 0.01) {
		x_vals[0] = xx;
		evalfis(out_vals,x_vals,bestfis);
//		y[c] = out_vals[0];
		printf("%f\n",out_vals[0]);
	}

//	chromosome_print(18,bestfis->params);

	fis_destroy(bestfis);

//	printf("Ranges:\n");
//	ranges_print(num_params, ranges);


	return 0;
}
