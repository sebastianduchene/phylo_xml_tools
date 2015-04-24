
library(phangorn)

## pval, tree length, tree topology

get_gc_summary <- function(gc_output){
    p_vals <- sapply(1:length(gc_output), function(x) gc_output[[x]][[1]][[3]])

    get_TL_diff <- function(x){
        true_tree_length <- sum(gc_output[[x]][[2]]$edge.length)
        est_tree_length <- sum(gc_output[[x]][[1]][[4]]$tree$edge.length)
        return( (est_tree_length - true_tree_length) / true_tree_length)
    }

    get_TP_diff <- function(x) dist.topo(gc_output[[x]][[2]], gc_output[[x]][[1]][[4]]$tree)

    tl_diffs <- sapply(1:length(gc_output), function(x) get_TL_diff(x))
    tp_diffs <- sapply(1:length(gc_output), function(x) get_TP_diff(x))

    return(list(p_vals, tl_diffs, tp_diffs))
}


######################
######################
load('jc_jc.Rdata')
jc_jc_summary <- get_gc_summary(jc_jc)

load('jc_gtr.Rdata')
jc_gtr_summary <- get_gc_summary(jc_gtr)

load('jc_gtr_g.Rdata')
jc_gtrg_summary <- get_gc_summary(jc_gtr_g)

load('gtr_jc.Rdata')
gtr_jc_summary <- get_gc_summary(gtr_jc)

load('gtr_gtr.Rdata')
gtr_gtr_summary <- get_gc_summary(gtr_gtr)

load('gtr_gtr_g.Rdata')
gtr_gtrg_summary <- get_gc_summary(gtr_gtr_g)

load('long_tree_gtrg.Rdata')
long_tree_gtrt_summary <- get_gc_summary(gtrg_long_tree)
