#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fuzzy.h"


int main(int argc, char * argv[]) {
	int i;
	double x[argc-1];
	double retval[1];
	struct Rule * rule1;
	struct Rule * rule2;
	struct Rule * rule3;
	struct Rule * rule4;
	struct Rule * rule_list[4];

	double size_small[3] = {0, 0, 1};
	double size_large[3] = {0, 1, 1};

	double weight_small[3] = {0, 0, 1};
	double weight_large[3] = {0, 1, 1};

	double qual_bad[3] = {0, 0, 0.5};
	double qual_med[3] = {0, 0.5, 1};
	double qual_good[3] = {0.5, 1, 1};

	double *input_list1[] = {size_small, weight_small};
	double *input_list2[] = {size_small, weight_large};
	double *input_list3[] = {size_large, weight_small};
	double *input_list4[] = {size_large, weight_large};

	double *output_list1[] = {qual_bad};
	double *output_list2[] = {qual_med};
	double *output_list3[] = {qual_good};

	rule1 = create_rule(input_list1, 2, output_list1, 1);
	rule2 = create_rule(input_list2, 2, output_list2, 1);
	rule3 = create_rule(input_list3, 2, output_list2, 1);
	rule4 = create_rule(input_list4, 2, output_list3, 1);

	rule_list[0] = rule1;
	rule_list[1] = rule2;
	rule_list[2] = rule3;
	rule_list[3] = rule4;

	for (i = 1; i < argc; i++) {
		x[i - 1] = atof(argv[i]);
	}
	for (i = 1; i < argc; i++) {
		printf("%0.2f\n",x[i - 1]);
	}


	evalfis(retval, x, rule_list, 4);

	printf("%0.5f\n", retval[0]);
	long long unsigned li;
	for (li = 0; li < 100000; li++) {
		x[0] += (1.0 + li) / li;
		evalfis(retval, x, rule_list, 4);
	}
	for (li = 0; li < 100000; li++) {
		x[0] += (1.0 + li) / li;
		evalfis(retval, x, rule_list, 4);
	}
	for (li = 0; li < 100000; li++) {
		x[0] += (1.0 + li) / li;
		evalfis(retval, x, rule_list, 4);
	}
//	evalfis(retval, x, rule_list, 4);

	int r;
	for (r = 0; r < 4; r++) {
		destroy_rule(rule_list[r]);
	}


	return 0;
}

