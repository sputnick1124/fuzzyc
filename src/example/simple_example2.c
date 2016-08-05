#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../include/fuzzy.h"


double genfuzex(double x, double y) {
//	struct Rule * rule_list[4];
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

	int in_mfs[] = {2, 2};
	int out_mfs[] = {3};

//	get_fis(rule_list, params, 2, 1, 4, rules, in_mfs, out_mfs);
	struct Fis * fis = fis_create(params, 2,1,4, rules,in_mfs, out_mfs);

	double xin[2];
	xin[0] = x;
	xin[1] = y;
//	evalrules(retval, xin, fis->rule_list, fis->num_rule);
	evalfis(retval, xin, fis);

//	int r;
//	for (r = 0; r < 4; r++) {
//		destroy_rule(rule_list[r]);
//	}

	fis_destroy(fis);
	return retval[0];
}


int main(int argc, char *argv[]) {
	double x, y, z;

	x = atof(argv[1]);
	y = atof(argv[2]);

	z = genfuzex(x,y);
	printf("%f\n",z);
	return 0;
}

