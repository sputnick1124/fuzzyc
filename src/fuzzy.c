#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "fuzzy.h"

double minimum(int n_val, double values[n_val]) {
	/*Find minimum value in an array of doubles*/
	double min = HUGE_VAL;
	int i;
	for (i = 0; i < n_val; i++) {
		min = (values[i] < min ? values[i] : min);
	}
	return min;
}

double maximum(int n_val, double values[n_val]) {
	/*Find maximum value in an array of doubles*/
	double max = 0;
	int i;
	for (i = 0; i < n_val; i++) {
		max = (values[i] > max ? values[i] : max);
	}
	return max;
}

double inTriMF(double params[3], double x) {
	/**/
	double a = params[0];
	double x_star = params[1];
	double b = params[2];
	double m1, b1, m2, b2;
	#ifdef DEBUG
	printf("inTriMF: x = %0.2f\n",x);
	#endif
	if (a != x_star) {
		m1 = 1 / (x_star - a);
		b1 = -m1 * a;
	} else {
		m1 = 0;
		b1 = 0;
	}

	if (b != x_star) {
		m2 = -1 / (b - x_star);
		b2 = 1 - (m2 * x_star);
	} else {
		m2 = 0;
		b2 = 0;
	}

	if (x > b) {
		x = b;
	} else if (x < a) {
		x = a;
	}

	if (x < x_star) {
		return fmax((m1 * x) + b1, 0);
	} else if ( x > x_star) {
		return fmax((m2 * x) + b2, 0);
	} else {
		return 1;
	}
}

void outTriMF(double x[][4], double * params[], int num_out, double y) {
	int o;
//	for (i = 0; i < num_out; i++) {
//		x[i] = malloc(sizeof(double));
//	}
	for (o = 0; o < num_out; o++) {
		double a = params[o][0];
		double x_star = params[o][1];
		double b = params[o][2];
		double m1, b1, m2, b2;
		if (a != x_star) {
			m1 = 1 / (x_star - a);
			b1 = -m1 * a;
		} else {
			m1 = 0;
			b1 = 0;
		}

		if (b != x_star) {
			m2 = -1 / (b - x_star);
			b2 = 1 - (m2 * x_star);
		} else {
			m2 = 0;
			b2 = 0;
		}

		x[o][0] = a;
		x[o][1] = (m1 ? ((y - b1) / m1) : a);
		x[o][2] = (m2 ? ((y - b2) / m2) : b);
		x[o][3] = b;
	}
}

void defuzzWeightedMV(double out[],
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firing_strengths[num_rule]) {

	int r, o;
	double hmv, ht;
	ht = 0;
	hmv = 0;
	for (o = 0; o < num_out; o++) {
		for (r = 0; r < num_rule; r++) {
			if (!firing_strengths[r]) {continue;}
			ht += firing_strengths[r];
			hmv += firing_strengths[r] * (traps[r][o][2] + traps[r][o][1]) / 2;
			#ifdef DEBUG
			printf("num_rule = %d\n",num_rule);
			printf("%f, %f, %f\n", firing_strengths[r], traps[r][o][1], traps[r][o][2]);
			printf("defuzzWeightedMV: rule %d: out %d\nht = %f\nhmv = %f\n",
				r, o, ht, hmv);
			#endif
		}
		out[o] = (ht ? hmv / ht : 0);
	}
//	return out;
}


void defuzzWMeanOfCent(double out[],
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firing_strengths[num_rule]) {

	int r, o;
	double mom, ht;
	double a, b, c, d;
	double c_temp, d_temp;
	double s_temp;
	mom = 0;
	ht = 0;
	for (o = 0; o < num_out; o++) {
		for (r = 0; r < num_rule; r++) {
			a = traps[r][o][2] - traps[r][o][1];
			b = traps[r][o][3] - traps[r][o][0];
			c_temp = (traps[r][o][1] - traps[r][o][0]);
			d_temp = (traps[r][o][3] - traps[r][o][2]);
			s_temp = firing_strengths[r]*firing_strengths[r];
			c = sqrt(c_temp*c_temp + s_temp);
			d = sqrt(d_temp*d_temp + s_temp);
			mom += firing_strengths[r] * (b/2 + (2*a + b)*(c*c - d*d) / (6*(b*b - a*a)));
			ht += firing_strengths[r];
			#ifdef DEBUG
			printf("trap[%d] = [%0.2f, %0.2f, %0.2f, %0.2f]\n",
					r,
					traps[r][o][0],
					traps[r][o][1],
					traps[r][o][2],
					traps[r][o][3]);
			printf("[a,b,c,d] = [%0.2f, %0.2f, %0.2f, %0.2f]\n",
					a,b,c,d);
			printf("[c,d,s_temp] = [%0.4f, %0.4f, %0.4f]\n",
					c_temp, d_temp, s_temp);
			printf("mom = %0.4f\nht = %0.3f\n",
					mom, ht);
			#endif
		}
		out[o] = mom / ht;
	}
//	return out;
}

void defuzzMeanOfCent(double out[],
	int num_rule,
	int num_out,
	double traps[num_rule][num_out][4],
	double firing_strengths[num_rule]) {

	int r, o;
	double mom;
	double a, b, c, d;
	double c_temp, d_temp;
	double h;
	mom = 0;
	for (o = 0; o < num_out; o++) {
		for (r = 0; r < num_rule; r++) {
			a = traps[r][o][2] - traps[r][o][1];
			b = traps[r][o][3] - traps[r][o][0];
			c_temp = (traps[r][o][1] - traps[r][o][0]);
			d_temp = (traps[r][o][3] - traps[r][o][2]);
			h = firing_strengths[r];
			c = (c_temp > 0 ? sqrt(c_temp*c_temp + h*h) : h);
			d = (d_temp > 0 ? sqrt(d_temp*d_temp + h*h) : h);
			mom += (b/2 + (2*a + b)*(c*c - d*d) / (6*(b*b - a*a)));
			#ifdef DEBUG
			printf("trap[%d] = [%0.2f, %0.2f, %0.2f, %0.2f]\n",
					r,
					traps[r][o][0],
					traps[r][o][1],
					traps[r][o][2],
					traps[r][o][3]);
			printf("[a,b,c,d] = [%0.2f, %0.2f, %0.2f, %0.2f]\n",
					a,b,c,d);
			printf("[c_temp,d_temp,h] = [%0.4f, %0.4f, %0.4f]\n",
					c_temp, d_temp, h);
			printf("mom = %0.4f\n",
					mom);
			#endif
		}
		out[o] = mom / (double) num_rule;
	}
//	return out;
}

struct Rule *create_rule(double input[][3], int num_in, double output[][3], int num_out) {
    int i, p;
    struct Rule *rule = malloc(sizeof(struct Rule));
    assert(rule != NULL);

	rule->input = malloc(num_in * sizeof(double *));

    rule->num_in = num_in;
    for (i = 0; i < num_in; i ++) {
		rule->input[i] = malloc(3 * sizeof(double));
		if (rule->input[i] == NULL) {
			for ( ; i >= 0; i--) {
				free(rule->input[i]);
			}
			free(rule);
			return NULL;
		}
		for (p = 0; p < 3; p++) {
        	rule->input[i][p] = input[i][p];
		}
    }

	rule->output = malloc(num_out * sizeof(double *));

    rule->num_out = num_out;
    for (i = 0; i < num_out; i ++) {
		rule->output[i] = malloc(3 * sizeof(double));
		if (rule->output[i] == NULL) {
			for ( ; i >= 0; i--) {
				free(rule->output[i]);
			}
			free(rule);
			return NULL;
		}
		for (p = 0; p < 3; p++) {
        	rule->output[i][p] = output[i][p];
		}
    }

    return rule;
}

void print_rule(struct Rule *rule, FILE * fd) {
    int i, p;

	if (fd == NULL) {
		fd = stdout;
	}
//	fprintf(fd, "rule is located at %p\n",rule);

    fprintf(fd, "%d inputs:\n",rule->num_in);
    for (i = 0; i < rule->num_in; i++) {
        for (p = 0; p < 3; p++) {
            fprintf(fd, "%0.2f ", rule->input[i][p]);
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "%d outputs:\n",rule->num_out);
    for (i = 0; i < rule->num_out; i++) {
        for (p = 0; p < 3; p++) {
            fprintf(fd, "%0.2f ", rule->output[i][p]);
        }
        fprintf(fd, "\n");
    }
    fprintf(fd, "\n\n");
}

void destroy_rules(struct Rule ** rule_list, int num_rule) {
	int i;
	for (i = 0; i < num_rule; i++) {
		destroy_rule(rule_list[i]);
	}
}

void destroy_rule(struct Rule *rule) {
	int i;
	assert(rule != NULL);
	for (i = 0; i < rule->num_in; i++) {
		free(rule->input[i]);
	}
	free(rule->input);
	for (i = 0; i < rule->num_out; i++) {
		free(rule->output[i]);
	}
	free(rule->output);
	free(rule);
}

struct Fis * fis_create(
    double params[],
    int num_in,
    int num_out,
    int num_rule,
    int rules[num_rule][num_in + num_out],
    int in_mfs[],
    int out_mfs[])
{
	struct Fis *fis = malloc(sizeof(struct Fis));
	assert(fis != NULL);

	fis->num_rule = num_rule;
	fis->rule_list = malloc(num_rule * sizeof(struct Rule *));
	if (fis->rule_list == NULL) {
		free(fis);
		return NULL;
	}

	get_fis(fis->rule_list, params, num_in, num_out, num_rule, rules, in_mfs, out_mfs);
	return fis;
}

void
fis_destroy(struct Fis *fis)
{
	destroy_rules(fis->rule_list,fis->num_rule);
	free(fis->rule_list);
	free(fis);
}

void
evalrules(
	double * out,
	double * x,
	struct Rule ** rules,
	int num_rule)
{
	int r, in;
	double firing_strengths[num_rule];
	double s_temp[rules[0]->num_in];
	double out_vals[num_rule][rules[0]->num_out][4];
	#ifdef DEBUG
	printf("num_rule = %d\n",num_rule);
	#endif
	for (r = 0; r < num_rule; r++) {
		#ifdef DEBUG
		print_rule(rules[r]);
		#endif
		for (in = 0; in < rules[r]->num_in; in++) {
			s_temp[in] = inTriMF(rules[r]->input[in], x[in]);
			#ifdef DEBUG
			printf("evalfis: x[%d] = %0.2f\n",in,x[in]);
			printf("inTriMF([%0.2f, %0.2f, %0.2f], %0.2f) = %0.2f\n",
					rules[r]->input[in][0],
					rules[r]->input[in][1],
					rules[r]->input[in][2],
					x[in],
					s_temp[in]);
			#endif
		}
		firing_strengths[r] = minimum(rules[r]->num_in,s_temp);
		outTriMF(out_vals[r],
				 rules[r]->output,
				 rules[0]->num_out,
				 firing_strengths[r]);
		#ifdef DEBUG
		int d, o;
		printf("firing_strength = min(inTriMF(x)) = %0.2f\n",firing_strengths[r]);
		for (o = 0; o < rules[0]->num_out; o++) {
			printf("outTriMF([%0.2f, %0.2f, %0.2f], %0.2f) = [",
					rules[r]->output[o][0],
					rules[r]->output[o][1],
					rules[r]->output[o][2],
					firing_strengths[r]);
			for (d = 0; d < 4; d++) {
				printf("%0.5f ",out_vals[r][o][d]);
			}
		}
		printf("]\n\n\n");
		#endif
	}
	defuzzWeightedMV(out,
		num_rule,
		rules[0]->num_out,
		out_vals,
		firing_strengths);
}

void
evalfis(
	double * out,
	double * x,
	struct Fis * fis)
{
	evalrules(out, x, fis->rule_list, fis->num_rule);
}


void get_fis(struct Rule ** rule_list,
	double params[],
	int num_in,
	int num_out,
	int num_rule,
	int rules[num_rule][num_in + num_out],
	int in_mfs[],
	int out_mfs[]) {

	int in, out, rule;
	int pin, pout;

	int ppoint, out_start;
	int _b;

	for (rule = 0; rule < num_rule; rule++) {
		#ifdef DEBUG
			printf("RULE%d\n",rule);
		#endif
		double input_list[num_in][3];
		double output_list[num_out][3];
		for (in = 0; in < num_in; in++) {
			ppoint = 0;
			if (in > 0) {
				for (_b = in-1; _b >= 0; _b--) {
					ppoint += 3 * in_mfs[_b];
				}
			} else { };
			ppoint += 3 * rules[rule][in];
			for (pin = 0; pin < 3; pin++) {
				#ifdef DEBUG
					printf("get_fis: INPUT:");
					printf("in = %d, pin = %d, ppoint = %d\n",in,pin,ppoint);
				#endif
				input_list[in][pin] = params[ppoint];
				ppoint += 1;
			}
		}

		out_start = 0;
		for (in = 0; in < num_in; in++) {
			out_start += 3 * in_mfs[in];
		}

		for (out = 0; out < num_out; out++) {
			ppoint = out_start;
			if (out > 0) {
				for (_b = out-1; _b >= 0; _b--) {
					ppoint += 3 * out_mfs[_b];
				}
			} else { };
			ppoint += 3 * rules[rule][num_in + out];
			for (pout = 0; pout < 3; pout++) {
				#ifdef DEBUG
					printf("get_fis: OUTPUT:");
					printf("out = %d, pout = %d, ppoint = %d\n",out, pout, ppoint);
				#endif
				output_list[out][pout] = params[ppoint];
				ppoint += 1;
			}
		}

		rule_list[rule] = create_rule(input_list, num_in, output_list, num_out);
	}
}

void
fis_print(struct Fis * fis, FILE * fd)
{
	int rule, in, out;

	if (fd == NULL) {
		fd = stdout;
	}

	for (rule = 0; rule < fis->num_rule; rule++) {
		fprintf(fd, "INPUT %d:\n\t[", rule);
		for (in = 0; in < fis->rule_list[rule]->num_in; in++) {
			print_rule(fis->rule_list[rule], fd);
		}
		fprintf(fd, "]\nOUTPUT %d:\n\t[", rule);
		for (out = 0; out < fis->rule_list[rule]->num_out; out++) {
			print_rule(fis->rule_list[rule], fd);
		}
		fprintf(fd, "]\n");
	}
}
