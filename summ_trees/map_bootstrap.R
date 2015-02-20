args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3) stop('Arguments missing. Please include the target tree, the bootstrap trees, and the output file name')

if(!(require('ape', character.only = TRUE))){
  install.packages('ape')
}

tr <- tryCatch(read.tree(args[1]), error = function(x) read.nexus(args[1]))
boot_trees <- tryCatch(read.tree(args[2]), error = function(x) read.nexus(args[2]))
out_name = args[3]

tr_labeled <- makeNodeLabel(tr)

tr$node.label <- prop.clades(phy = tr, part = prop.part(boot_trees))
write.tree(tr, file = out_name)

cat('\nI saved the tree with bootstrap values in: ', getwd(), '/', out_name, '\n\n', sep = '')

