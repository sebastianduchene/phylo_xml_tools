# Summarise trees

## Map bootstrap replicates on a treee

- Make sure R is installed. If not download from [here](http://www.r-project.org)

- Open a terminal window.

- Type 'Rscript' and drag 'map_bootstrap.R', the targer tree, the bootstrap trees. At the end, type the name of the output file where you want to save the annotated tree. The prompt should look something like this like this:

```
your_computer_name$ Rscript map_bootstrap.R my_target.tree my_bootstrap_reps.trees my_annotated.tree
```

- Hit enter.

- Open the annotated tree (my_annotated.tree in the example above) in figtree and find the bootstrap values using the nodelabels menu.

- Note that the example trees here are in Newick format. But this should also work for trees in Nexus format.
