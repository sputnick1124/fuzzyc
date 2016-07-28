#include <stdio.h>
#include "fuzzy.h"

int main(int argc, char * argv[]) {
	double input[1][3] = {0,0,1.0};
	double output[3] = {0,0.5,1};
	struct Rule *rule = create_rule(output, 1, input);
	printf("%d\n",(int)sizeof(rule));
	return 0;
}
