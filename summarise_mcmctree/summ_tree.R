library(ape)

dir()
tree_files <- dir(pattern = '[.]tre')

all_trees_lines <- vector()
for(i in 1:length(tree_files)){
    tree_lines <- readLines(tree_files[i])
    all_trees_lines[i] <- gsub('[[](&|[0-9]|%|=| |[{]|[}]|[.]|,|[]])+', '', tree_lines[4])
}

all_trees <- read.tree(text = all_trees_lines)

write.nexus(all_trees, file = 'all_trees.trees')
