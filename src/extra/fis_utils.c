#include <stdio.h>
#include <string.h>

#include "fuzzy.h"

struct Fis * fis_parse(FILE * fisfile) {
    char * tmp_chunk = NULL;
    char * chunk = NULL;
    char *token, *save_ptr;
    size_t len = 0;
    int i, j;

    int num_in, num_out;
    // Find first chunk of text between {} delimiters. This is num_in.
    getdelim(&chunk,&len,'{',fisfile);
    tmp_chunk = chunk;
    getdelim(&chunk,&len,'}',fisfile);
    free(tmp_chunk);
    tmp_chunk = chunk;
    sscanf(chunk,"%d",&num_in);

    // Next chunk of text in delimiters contains num_in integers.
    // Need to tokenize and iterate through to extract them.
    int in_mfs[num_in];
    getdelim(&chunk,&len,'{',fisfile);
    free(tmp_chunk);
    tmp_chunk = chunk;
    getdelim(&chunk,&len,'}',fisfile);
    free(tmp_chunk);
    tmp_chunk = chunk;
    for (i = 0; i < num_in; i++, chunk=NULL) {
        token = strtok_r(chunk," ", &save_ptr);
        sscanf(token,"%d",&in_mfs[i]);
    }

    // Find next chunk of text between {} delimiters. This is num_out.
    getdelim(&chunk,&len,'{',fisfile);
    free(tmp_chunk);
    tmp_chunk = chunk;
    getdelim(&chunk,&len,'}',fisfile);
    free(tmp_chunk);
    tmp_chunk = chunk;
    sscanf(chunk,"%d",&num_out);

    // Next chunk of text in delimiters contains num_out integers.
    // Need to tokenize and iterate through to extract them.
    int out_mfs[num_out];
    getdelim(&chunk,&len,'{',fisfile);
    free(tmp_chunk);
    tmp_chunk = chunk;
    getdelim(&chunk,&len,'}',fisfile);
    free(tmp_chunk);
    tmp_chunk = chunk;
    for (i = 0; i < num_out; i++, chunk=NULL) {
        token = strtok_r(chunk," ",&save_ptr);
        sscanf(token,"%d",&out_mfs[i]);
    }

    // Here we can assume that the doubles come in triplets (since we're using triangles)
    int num_mfs = sum_i(in_mfs) + sum_i(out_mfs);
    int p = 0;
    double params[num_mfs * 3];
    for (i = 0; i < num_in; i++) {
        for (j = 0; j < in_mfs[i]; j++) {
            getdelim(&chunk,&len,'{',fisfile);
            free(tmp_chunk);
            tmp_chunk = chunk;
            getdelim(&chunk,&len,'}',fisfile);
            free(tmp_chunk);
            tmp_chunk = chunk;
            sscanf(chunk,"%lf %lf %lf",&params[p],&params[p+1],&params[p+2]);
            p += 3;
        }
    }
    for (i = 0; i < num_out; i++) {
        for (j = 0; j < out_mfs[i]; j++) {
            getdelim(&chunk,&len,'{',fisfile);
            free(tmp_chunk);
            tmp_chunk = chunk;
            getdelim(&chunk,&len,'}',fisfile);
            free(tmp_chunk);
            tmp_chunk = chunk;
            sscanf(chunk,"%lf %lf %lf",&params[p],&params[p+1],&params[p+2]);
            p += 3;
        }
    }

    int num_rule = prod_i(in_mfs);
    int num_args = num_in + num_out;
    int rules[num_rules][num_args];
    for (i = 0; i < num_rule; i++) {
        getdelim(&chunk,&len,'{',fisfile);
        free(tmp_chunk);
        tmp_chunk = chunk;
        getdelim(&chunk,&len,'}',fisfile);
        free(tmp_chunk);
        tmp_chunk = chunk;
        for (j = 0; j < num_args; j++, chunk = NULL) {
            token = strtok_r(chunk," ",&save_ptr);
            sscanf(token,"%d",&rules[i][j]);
        }
    }

    struct Fis * fis = fis_create(
        params,
        num_in,
        num_out,
        num_rule,
        rules,
        in_mfs,
        out_mfs);
    return fis;
}


