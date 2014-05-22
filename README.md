Date randomisation for BEAST xml files
======================================

Sebastian Duchene

May 22 2014


In this example I show how to randomise the dates for BEAST XML files with the python script randomise_dates.py. The script has been tested in Unix machines, but it should work on windows computers with a python installation. 

Downloaad randomise_dates.py and the example_file or this complete repository.

The first requirement is an xml file with heterochronous data, such as for samples from viruses or ancient DNA.

Follow these instructions for Mac or Unix machines. For windows you can open python and run the script as a module, then follow these steps from (4).:

1- Open a terminal window (you can find it using spotlight).

2- type in "cd", followed by a space and drag in the folder where you want the xml files with the randomised dates:

```js
git_testing sebastianduchene$ cd /Users/sebastianduchene/Desktop/deprate_sims2/git_testing/example_files/
```

3- At the terminal window type in "python", followed by a space. Now drag the script, which is called randomise_dates.py. The terminal window should look somethin like this:

```js
git_testing sebastianduchene$ python /Users/sebastianduchene/Desktop/deprate_sims2/git_testing/example_files/randomise_dates.py
```

4- Type enter. You will be prompted to drag your input file. This should be an xml file generated in BEAST with sampling times for the tips:

```js
Please input the path to your file, or drag it here:
```

5- Now drag your xml file. In this case I will use the example file in this repository:

```js
Please input the path to your file, or drag it here: /Users/sebastianduchene/Desktop/deprate_sims2/git_testing/example_files/hcv31_ns5b.xml
```

6- Type enter. You will be prompted for a number of randomised data sets to generate. Type in a number. Typically 5 is a good number:

```js
How many data sets with randomised dates should I produce? 5
```

7- Type enter. You will be prompted to name the output files. I will name them "random_test" to distinguish them from my original file:

```js
What should be the name of the output file(s)? random_test
``` 

8- Type enter. The program will output some lines, where you can confirm the the correct lines have been changed in your xml file. The last line tells you where it saved the files.

```js
The true date is <date value="2002" direction="forwards" units="years"/>

The new assigned date is <date value="2006" direction="forwards" units="years"/>

The true date is <date value="2003" direction="forwards" units="years"/>

The new assigned date is <date value="2007" direction="forwards" units="years"/>

I am writting replicate 5 in out_test_25.xml, out of 5 replicates

I have saved the files in /Users/sebastianduchene/Desktop/deprate_sims2/git_testing/example_files/
```

9- Run these files in BEAST. The output for the log, trees, and ops files are the same as those of the randomised files. In thsi example it would be random_test1.log, random_test1.trees and random_test1.ops,