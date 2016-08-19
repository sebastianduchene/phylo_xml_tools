import os, re, sys
import dendropy as dp
import numpy as np

input_tree_file, output_file = sys.argv[1:3]

phylogram = dp.Tree.get_from_path(input_tree_file, 'newick', preserve_underscores=True)
out_file = open(output_file, 'w')

dates = list()
for i in phylogram.taxon_namespace:
    split_name = re.split('_', i.label)
    taxon_name = '_'.join(split_name[0:-1])
    taxon_date = split_name[-1]
    dates.append([i.label,'\t',taxon_date,'\n'])

dates_array = np.array(dates)

try:
    sys.argv[3]
    dates_array[:, 2] = dates_array[np.random.uniform(size = dates_array.shape[0]).argsort(),2]
except:
    pass

out_file.write(str(len(dates))+'\n')

for d in dates_array:
    out_file.write(''.join(list(d)))

print 'Date file is:'

print dates_array

out_file.close()

