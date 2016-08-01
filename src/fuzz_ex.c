#include <gnuplot_plugin.h>
#include <math.h>
#include "fuzzy.h"


double genfuzex(double x, double y) {
	struct Rule * rule_list[4];
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

	get_fis(rule_list, params, 2, 1, 4, rules, inmfs, outmfs);

	double xin[2];
	xin[0] = x;
	xin[1] = y;
	evalrules(retval, xin, rule_list, 4);

	int r;
	for (r = 0; r < 4; r++) {
		destroy_rule(rule_list[r]);
	}
	return retval[0];
}
