#!/usr/bin/env python
import pylab


file = 'distparams'
models = "model_names.dat"

file_root = "file_root = chains/"
params = 'parameter_names = chains/'



with open(models) as inf:
    for line in inf:
        parts = line.split()
	modelnm = parts[8].replace(".paramnames", "")

	file1 = open(file + '.ini', 'r')
	lines = file1.readlines()
	lines[0] = file_root + modelnm + '\n'
	lines[1] = params + modelnm + '.paramnames' + '\n'
	file1.close()

	output = file + '_' + modelnm + '.ini'
	file2 = open(output, 'w')
	for line in lines:
		file2.write(line)
	file2.close()
	print './getdiest' + output
