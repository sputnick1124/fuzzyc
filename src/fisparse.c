#include <stdio.h>
#include <yaml.h>

#include "fuzzy.h"

int main(int argc, char *argv[]) {
	char * file = argv[1];
	FILE *fh = fopen(file, "r");
	yaml_parser_t parser;
	yaml_event_t  event;   /* New variable */

	/* Initialize parser */
	if(!yaml_parser_initialize(&parser)) {
		fputs("Failed to initialize parser\n", stderr);
        return 1;
	}
	if(fh == NULL) {
		fputs("Failed to open file\n",stderr);
        return 1;
	}

	/* Set input file */
	yaml_parser_set_input_file(&parser, fh);
	int done = 0;
	for (;!done;) {
		if (!yaml_parser_parse(&parser, &event)) {
			char *err_msg;
			sprintf(err_msg,"Parser error %d\n", parser.error);
			fputs(err_msg,stderr);
			exit(EXIT_FAILURE);
		}

		switch(event.type) {
			case YAML_NO_EVENT: printf("No event!"); break;

			/* Stream start/end */
			case YAML_STREAM_START_EVENT: puts("STREAM START"); break;
			case YAML_STREAM_END_EVENT:   puts("STREAM END");   done = 1; break;
			/* Block delimeters */
			case YAML_DOCUMENT_START_EVENT: puts("<b>Start Document</b>"); break;
			case YAML_DOCUMENT_END_EVENT:   puts("<b>End Document</b>");   break;
			case YAML_SEQUENCE_START_EVENT: puts("<b>Start Sequence</b>"); break;
			case YAML_SEQUENCE_END_EVENT:   puts("<b>End Sequence</b>");   break;
			case YAML_MAPPING_START_EVENT:  puts("<b>Start Mapping</b>");  break;
			case YAML_MAPPING_END_EVENT:    puts("<b>End Mapping</b>");    break;
			/* Data */
			case YAML_ALIAS_EVENT:  printf("Got alias (anchor %s)\n", event.data.alias.anchor); break;
			case YAML_SCALAR_EVENT: printf("Got scalar (value %s)\n", event.data.scalar.value); break;
		}
		if(event.type != YAML_STREAM_END_EVENT) {
			yaml_event_delete(&event);
		}
	}
	yaml_event_delete(&event);
	/* END new code */

	/* Cleanup */
	yaml_parser_delete(&parser);
	fclose(fh);
	return 0;
}

