#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fuzzy.h"


int main(int argc, char * argv[]) {
	struct Rule * rule_list[4];
	int i;
	double x[2];
	double retval[1];

	double params[] = {
	0.0, 0.0, 1.0,
	0.0, 1.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 1.0, 1.0,
	0.0, 0.0, 0.5,
	0.0, 0.5, 1.0,
	0.5, 1.0, 1.0};

	int rules[][3] = {
	{0, 0, 0},
	{0, 1, 1},
	{1, 0, 1},
	{1, 1, 2}};

	int inmfs[] = {2, 2};
	int outmfs[] = {3};

	for (i = 1; i < argc; i++) {
		x[i - 1] = atof(argv[i]);
	}
	for (i = 1; i < argc; i++) {
		printf("%0.2f\n",x[i - 1]);
	}

	get_fis(rule_list, params, 2, 1, 4, rules, inmfs, outmfs);

	evalrules(retval, x, rule_list, 4);

	printf("%0.5f\n", retval[0]);
	long long unsigned li;
	int c = 0;
	for (li = 0; li < 100001; li++) {
		x[0] += (1.0 + li) / li;
		evalrules(retval, x, rule_list, 4);
		if (li == 100000 && c <= 8) {
			c += 1;
			li = 0;
		}
	}

	int r;
	for (r = 0; r < 4; r++) {
		destroy_rule(rule_list[r]);
	}


	return 0;
}

