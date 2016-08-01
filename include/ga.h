#ifndef _GALIB_H_
#define _GALIB_H_

int sum_i(int adds[], size_t len);

int prod_i(int mults[], size_t len);

void init_partition(double params[], int num_sets);

void rand_partition(double params[], int num_sets);

void init_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]);

void rand_params(double params[], int num_in, int in_mfs[], int num_out, int out_mfs[]);

void init_antecedents(int num_in, int num_out, int rules[][num_in], int in_mfs[]);

#endif
