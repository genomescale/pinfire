#!/usr/bin/python2.7

import os
import sys

fix_folder = sys.argv[1]
original_filenames = sorted(os.listdir(fix_folder))

for original_fn in original_filenames:
	original_parameters = original_fn.split("_")
	fixed_parameters = []
	for op in original_parameters:
		op_value = op[:-1]
		op_key = op[-1]
		if op_key == "l" and op_value.isdigit():
			fixed_op = "%03dl" % int(op_value)
			fixed_parameters.append(fixed_op)
		else:
			fixed_parameters.append(op)

	fixed_fn = "_".join(fixed_parameters)
	
	if original_fn != fixed_fn:
		original_path = os.path.join(fix_folder, original_fn)
		fixed_path = os.path.join(fix_folder, fixed_fn)
		os.rename(original_path, fixed_path)