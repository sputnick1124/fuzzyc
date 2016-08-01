#include <stdio.h>

int prod_i(int mults[], size_t len) {
	int retval = 1;
	int i;
	for (i = 0; i < len; i++) {
		retval *= mults[i];
	}
	return retval;
}

int sum_i(int adds[], size_t len) {
	int retval = 0;
	int i;
	for (i = 0; i < len; i++) {
		retval += adds[i];
	}
	return retval;
}

void init_partition(double params[], int num_sets) {
	int p;
	int num_params = num_sets * 3;
	params[0] = 0; params[1] = 0; //Set lower bounds
	params[num_params-1] = 1;	//Set upper bounds
	params[num_params-2] = 1;
	double I = 1 / ((double) num_sets - 1); //Interval to divide up the fuzzy partition
	for (p = 2; p < num_params - 2; p++) {
		if (p % 3 == 0) { //a_i
			params[p] = ((p / 3) - 1) * I;
		} else if (p % 3 == 1) { //b_i
			params[p] = ((p-1) / 3) * I;
		} else { //c_i
			params[p] = (((p-2) / 3) + 1) * I;
		}
	}
}

void init_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]) {
	int in, out, p;

	p = in_mfs[0];
	init_partition(params, p);
	for (in = 1; in < num_in; in++) {
		init_partition(&params[p * 3], in_mfs[in]);
		p += in_mfs[in];
	}

	for (out = 0; out < num_out; out++) {
		init_partition(&params[p * 3], out_mfs[out]);
		p += out_mfs[out];
	}
}

void init_rules(int num_in, int rules[][num_in], int in_mfs[]) {
	int in, rule;
	unsigned int flag = 0;
	int num_rule = 1;
	//Initialize first rule to zeros and compute num_rule
	for (in = 0; in < num_in; in++) {
		rules[0][in] = 0;
		num_rule *= in_mfs[in];
	}
	for (rule = 1; rule < num_rule; rule++) {
		if (rules[rule-1][0] == (in_mfs[0] - 1)) {
			rules[rule][0] = 0;
			flag ^= (1 << 1);//Tell next element to increment
		} else {
			rules[rule][0] = rules[rule-1][0] + 1;
		}
		for (in = 1; in < num_in; in++) {
			if (flag & (1 << in)) { //If this element's bit is set
				if (rules[rule-1][in] == (in_mfs[in]-1)) { //And previous element is max
					rules[rule][in] = 0; //Start back at zero
					flag ^= (1 << (in + 1)); //Set next element's bit
				} else {
					rules[rule][in] = rules[rule-1][in] + 1;
				}
				flag ^= (1 << in);//Unset this element's bit regardless
			} else {
				rules[rule][in] = rules[rule-1][in];
			}
		}
	}
}

int main(int argc, char * argv[]) {
	int in_mfs[3] = {3, 4, 2};
	int num_in = 3;
	int out_mfs[1] = {3};
	int num_out = 1;
	int num_params = 3 * (sum_i(in_mfs, num_in) + sum_i(out_mfs,num_out));
	int num_rules = prod_i(in_mfs,num_in) * prod_i(out_mfs,num_out);
	int rules[num_rules][3];
	double params[num_params];
	int rule, in, p;

//	init_partition(params, in_mfs[0]);
//	init_partition(&params[in_mfs[0]*3], in_mfs[1]);
//	init_partition(&params[sum_i(in_mfs,2)*3], in_mfs[2]);
	init_params(params, num_in, in_mfs, num_out, out_mfs);

	init_rules(num_in, rules, in_mfs);

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
		for (in = 0; in < num_in; in++) {
			printf("%d ",rules[rule][in]);
		}
		printf("]\n");
	}
	return 0;
}
