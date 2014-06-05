import re, sys, os, random, copy

file_name, n_rand, out_name = sys.argv[1:]

print file_name, n_rand, out_name

#file_name = raw_input("Please input the path to your file, or drag it here: ")
#n_rand = int(raw_input("How many data sets with randomised dates should I produce? "))
#out_name = raw_input("What should be the name of the output file(s)? ")

data_lines = open(re.sub(r' ', '', file_name)  , 'r').readlines()

# find the location of the dates. This is the original position

dates_location = []
for i in range(len(data_lines)):
    if re.search(r'<date', data_lines[i]) == None:
        continue
    else:
        dates_location.append(i)

# extract the dates

dates_blocks = [data_lines[i] for i in dates_location]

# shuffle the order of the dates n times

shuffle_index = []

for i in range(n_rand):
    shuffle_index.append(random.sample(range(len(dates_location)), len(dates_location)))

# insert in the original position

def change_out_names(data_file, new_name):
    new_lines = [re.sub(r'fileName=.+log' , 'fileName="'+out_name+str(i+1)+'.log', m) for m in data_file]
    new_lines = [re.sub(r'fileName=.+trees', 'fileName="'+out_name+str(i+1)+'.trees', m) for m in new_lines]
    new_lines = [re.sub(r'operatorAnalysis=.+ops', 'operatorAnalysis="'+out_name+str(i+1)+'.ops', m) for m in new_lines]
    return new_lines


for i in range(len(shuffle_index)):
    shuff_temp = shuffle_index[i]
    lines_temp = copy.copy(data_lines)

    for k in range(len(dates_blocks)):
        print "The true date is %s " %lines_temp[dates_location[k]]  
        print "The new assigned date is %s" %dates_blocks[shuff_temp[k]]
        lines_temp[dates_location[k]] = dates_blocks[shuff_temp[k]]
    lines_temp = change_out_names(lines_temp, out_name)

    print "I am writting replicate %s in %s, out of %s replicates" %(i + 1, out_name+str(i+1)+".xml", n_rand)
    open(out_name+str(i+1)+".xml", 'w').writelines(lines_temp)

print "I have saved the files in %s" %os.getcwd()
