import fuzzy

size_small = [0, 0, 1]
size_large = [0, 1, 1]

weight_small = [0, 0, 1]
weight_large = [0, 1, 1]

qual_bad = [0.0, 0.0, 0.5]
qual_med = [0.0, 0.5, 1.0]
qual_good = [0.5, 1.0, 1.0]

params = []
params.extend(size_small)
params.extend(size_large)
params.extend(weight_small)
params.extend(weight_large)
params.extend(qual_bad)
params.extend(qual_med)
params.extend(qual_good)

num_in = 2
inmfs = [2, 2]

num_out = 1
outmfs = [3]

num_rule = 4
rules = [
	0, 0, 0,
	0, 1, 1,
	1, 0, 1,
	1, 1, 2]

x = [0.2, 0.25]
retval = fuzzy.evalrules(x, params, num_in, num_out, num_rule, rules, inmfs, outmfs)

print(retval)

