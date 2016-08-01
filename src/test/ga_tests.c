#include <stdio.h>
#include "../../include/ga.h"
#include "../../include/fuzzy.h"

int main(int argc, char * argv[]) {
	int num_in = 3;
	int in_mfs[] = {3, 4, 2};
	int num_out = 1;
	int out_mfs[] = {3};
	int num_params = 3 * (sum_i(in_mfs, num_in) + sum_i(out_mfs,num_out));
	int num_rules = prod_i(in_mfs,num_in) * prod_i(out_mfs,num_out);
	int rules[num_rules][num_in + num_out];
	double params[num_params];
	int rule, in, p;

//	init_partition(params, in_mfs[0]);
//	init_partition(&params[in_mfs[0]*3], in_mfs[1]);
//	init_partition(&params[sum_i(in_mfs,2)*3], in_mfs[2]);
	rand_params(params, num_in, in_mfs, num_out, out_mfs);

	init_antecedents(num_in, num_out, rules, in_mfs);

	printf("[");
	for (p = 0; p < num_params; p++) {
		if (p % 3 == 0) {
			printf("|");
		}
		printf("%0.2f ",params[p]);
	}
	printf("]\n");

	for (rule = 0; rule < prod_i(in_mfs, num_in); rule++) {
		printf("[");
		for (in = 0; in < num_in + num_out; in++) {
			printf("%d ",rules[rule][in]);
		}
		printf("]\n");
	}
	return 0;
}
