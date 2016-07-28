#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "fuzzy.h"

double minimum(int n_val, double values[n_val]) {
	double min = HUGE_VAL;
	int i;
	for (i = 0; i < n_val; i++) {
		min = (values[i] < min ? values[i] : min);
	}
	return min;
}

double maximum(int n_val, double values[n_val]) {
	double max = 0;
	int i;
	for (i = 0; i < n_val; i++) {
		max = (values[i] > max ? values[i] : max);
	}
	return max;
}

double inTriMF(double params[3], double x) {
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
	out = 0;
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

struct Rule *create_rule(double * input[3], int num_in, double * output[3], int num_out) {
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

void print_rule(struct Rule *rule) {
    int i, p;

	printf("rule is located at %p\n",rule);

    printf("%d inputs:\n",rule->num_in);
    for (i = 0; i < rule->num_in; i++) {
        for (p = 0; p < 3; p++) {
            printf("%0.2f ", rule->input[i][p]);
        }
        printf("\n");
    }
    printf("%d outputs:\n",rule->num_out);
    for (i = 0; i < rule->num_out; i++) {
        for (p = 0; p < 3; p++) {
            printf("%0.2f ", rule->output[i][p]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void destroy_rule(struct Rule *rule) {
	int i;
	assert(rule != NULL);
	for (i = 0; i < rule->num_in; i++) {
		free(rule->input[i]);
	}
	free(rule->input);
	for (i = 0; i < rule->num_in; i++) {
		free(rule->output[i]);
	}
	free(rule->output);
	free(rule);
}

void evalfis(double * out, double * x, struct Rule ** rules, int num_rule) {
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
	defuzzMeanOfCent(out,
		num_rule,
		rules[0]->num_out,
		out_vals,
		firing_strengths);
}

