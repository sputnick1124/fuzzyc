#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fuzzy.h"


int main(int argc, char *argv[]) {
	struct Rule * rule_list[4];
	double retval[1];
	double x = 0, y = 0;
	double xin[2];
	double out[1000 * 1000];
	int i, j;

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

	get_fis(rule_list, params, 2, 1, 4, rules, inmfs, outmfs);


	for (i = 0; i < 1000; i++) {
		for (j = 0; j < 1000; j++) {
			x = (1.0 / 1000.0) * (double)i;
			y = (1.0 / 1000.0) * (double)j;
			xin[0] = x;
			xin[1] = y;
			evalrules(retval, xin, rule_list, 4);
			out[i*j] = retval[0];
		}
	}

//	for (i = 0; i < 1000*1000; i++){
//		printf("%f\n",out[i]);
//	}
	int r;
	for (r = 0; r < 4; r++) {
		destroy_rule(rule_list[r]);
	}
	return 0;
}


