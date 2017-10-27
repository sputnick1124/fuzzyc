#include <Python.h>
#include <stdio.h>
#include "fuzzy.h"

#ifdef PYTHON2
#pragma message("Using PYTHON2")
#endif
#ifdef PYTHON3
#pragma message("Using PYTHON3")
#endif

static PyObject *
fuzzy_evalrules(PyObject *self, PyObject *args) {
	PyObject *x;
	PyObject *params;
	int num_in;
	int num_out;
	int num_rule;
	PyObject *rules;
	PyObject *inmfs;
	PyObject *outmfs;
	if (!PyArg_ParseTuple(args,
		"OOiiiOOO",
		&x,
		&params,
		&num_in,
		&num_out,
		&num_rule,
		&rules,
		&inmfs,
		&outmfs)) {
		//error
		printf("Something is terribly wrong");
	}

	int i, j;
	int num_params = 3 * (num_in + num_out + num_rule);
	double c_params[num_params];
	int c_rules[num_rule][(num_in + num_out)];
	int c_inmfs[num_in];
	int c_outmfs[num_out];
	double c_x[num_in];

	PyObject *x_input = PyObject_GetIter(x);
//	printf("x = [");
	for (i = 0; i < num_in; i++) {
		PyObject *next = PyIter_Next(x_input);
		c_x[i] = PyFloat_AsDouble(next);
//		printf("%0.2f",c_x[i]);
	}
//	printf("]\n\n");

	/*Get the param list into a double[]*/
	PyObject *param_list = PyObject_GetIter(params);
	if (!param_list) {
		//error not iterator
	}

//	printf("params = [");
	for (i = 0; i < num_params; i++) {
		PyObject *next = PyIter_Next(param_list);
		if (!next) {
			//nothing left in the iterator
			break;
		}

		if (!PyFloat_Check(next)) {
			//error, we were expecting a floating point value
		}

		c_params[i] = PyFloat_AsDouble(next);
//		printf("%0.2f ",c_params[i]);
	}
//	printf("]\n\n");

//	printf("rules lives at %p\n",rules);
	PyObject *rule_list = PyObject_GetIter(rules);
	if (!rule_list) {
		//error not iterator
//		printf("something is wrong with the rules\n");
	}

//	printf("rules = [");
	for (i = 0; i < num_rule; i++) {
		for (j = 0; j < num_in + num_out; j++) {
			PyObject *next = PyIter_Next(rule_list);
			if (!next) {
				//nothing left in rule
				break;
			}
			#ifdef PYTHON2
			c_rules[i][j] = (int) PyInt_AsLong(next);
			#endif
			#ifdef PYTHON3
			c_rules[i][j] = (int) PyLong_AsLong(next);
			#endif
//			printf("%d",c_rules[i][j]);
		}
	}
//	printf("]\n");


	PyObject *inmf_list = PyObject_GetIter(inmfs);
	if (!inmf_list) {
		//error not iterator
	}

//	printf("inmfs = [");
	for (i = 0; i < num_in; i++) {
		PyObject *next = PyIter_Next(inmf_list);
		if (!next) {
			//nothing left
			break;
		}

		#ifdef PYTHON2
		c_inmfs[i] = (int) PyInt_AsLong(next);
		#endif
		#ifdef PYTHON3
		c_inmfs[i] = (int) PyLong_AsLong(next);
		#endif

//		printf("%d ",c_inmfs[i]);
	}
//	printf("]\n\n");

	PyObject *outmf_list = PyObject_GetIter(outmfs);
	if (!outmf_list) {
		//error not iterator
	}

//	printf("outmfs = [");
	for (i = 0; i < num_out; i++) {
		PyObject *next = PyIter_Next(outmf_list);
		if (!next) {
			//nothing left
			break;
		}

		#ifdef PYTHON2
		c_outmfs[i] = (int) PyInt_AsLong(next);
		#endif
		#ifdef PYTHON3
		c_outmfs[i] = (int) PyLong_AsLong(next);
		#endif

//		printf("%d ",c_outmfs[i]);
	}
//	printf("]\n\n");

	struct Rule * c_rule_list[num_rule];
	double retval[num_out];

	get_fis(c_rule_list, c_params, num_in, num_out, num_rule, c_rules, c_inmfs, c_outmfs);
	evalrules(retval, c_x, c_rule_list, num_rule);

	PyObject *RetList = PyList_New(num_out);
	for (i = 0; i < num_out; i++) {
		PyList_SetItem(RetList, i, PyFloat_FromDouble(retval[i]));
	}

	for (i = 0; i < 4; i++) {
		destroy_rule(c_rule_list[i]);
	}

	return RetList;
}

static PyMethodDef FuzzyMethods[] = {
    {"evalrules",  fuzzy_evalrules, METH_VARARGS,
     "Evaluate FIS represented by rules and MF params."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

#ifdef PYTHON2
PyMODINIT_FUNC
initfuzzy(void)
{
    (void) Py_InitModule("fuzzy", FuzzyMethods);
}
#endif

#ifdef PYTHON3
static struct PyModuleDef fuzzymodule = {
	PyModuleDef_HEAD_INIT,
	"fuzzy",   /* name of module */
	NULL, /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
	FuzzyMethods
};

PyMODINIT_FUNC
PyInit_fuzzy(void)
{
	return PyModule_Create(&fuzzymodule);
}
#endif


