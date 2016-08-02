#include <stdio.h>
#include <time.h>
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

int main(int argc, char * argv[]) {
	int num_in = 3;
	int in_mfs[] = {3, 4, 2};
	int num_out = 1;
	int out_mfs[] = {3};
	int num_params = 3 * (sum_i(in_mfs, num_in) + sum_i(out_mfs,num_out));
	int num_rules = prod_i(in_mfs,num_in);
	int rules[num_rules][num_in + num_out];
	double params[num_params];
	double testp1[num_params];
	double testp2[num_params];
	double testc1[num_params];
	double testc2[num_params];
	double testc3[num_params];
	double testc4[num_params];
	const double ranges[num_params][2];
	int in, test;
	int fails = 0;

	srand48((long int) time(NULL));
	srand((long int) time(NULL));
//	init_partition(params, in_mfs[0]);
//	init_partition(&params[in_mfs[0]*3], in_mfs[1]);
//	init_partition(&params[sum_i(in_mfs,2)*3], in_mfs[2]);
	rand_params(params, num_in, in_mfs, num_out, out_mfs);

	init_antecedents(num_in, num_out, rules, in_mfs);
	rand_consequents(num_in, num_out, num_rules, rules, out_mfs);

	param_range(num_params, ranges, num_in, in_mfs, num_out, out_mfs);

	/*Test that rand_params returns a valid list of params*/
	printf("Testing rand_params\n");
	for (test = 0; test < 1000000; test++) {
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
	for (test = 0; test < 1000000; test++) {
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
	for (test = 0; test < 1000000; test++) {
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
	for (test = 0; test < 1000000; test++) {
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
	for (test = 0; test < 1000000; test++) {
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

	return 0;
}
