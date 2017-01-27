import re, os, sys
from random import shuffle
import xml.etree.ElementTree as ET

print """
To run, use:
python randomise_dates_beast2.py Number_of_randomisations xml_file_name
"""

n_rand = int(sys.argv[1])
input_name = sys.argv[2]

tree = ET.parse(input_name)
root = tree.getroot()

# Go down the tree to find the trait:
run = root.findall('run')[0]
state = run.findall('state')[0]
tr = state.findall('tree')[0]
trait = tr.findall('trait')[0]
dates_text = trait.text
dates_split = re.split(',', dates_text)

# Repeat from here down for each randomisation:
for rand_rep in range(n_rand):

    randomisation_name = re.sub('[.]xml','_randomised_'+str(rand_rep), input_name)
    print 'Saving files in '+randomisation_name
    dates = list()
    taxa_names = list()

    for d in dates_split:
        temp = re.split('=', d)
        taxa_names.append(temp[0])
        dates.append(temp[1])

    shuffle(dates)
    new_dates = list()

    for i in range(len(dates)):
        new_dates.append(taxa_names[i]+'='+dates[i])
    
    dates_text = ','.join(new_dates)
    trait.text = dates_text
    logs = run.findall('logger')

    for log in logs:
        log_keys = log.keys()
        if 'fileName' in log_keys:
            log.set('fileName', re.sub('.+[.]', randomisation_name+'.', log.get('fileName')))

    tree.write(randomisation_name+'.xml')
