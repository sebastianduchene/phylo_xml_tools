
simulate_hetero <- function(ntax = 50, slen = 1000){


    # 10 taxa have different base frequencies
    tr1 <- rtree(ntax - 11)
    tr1$edge.length <- rlnorm(length(tr1$edge.length), meanlog = -4.5, sd = 0.3)
    s1 <- simSeq(tr1, l = slen, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3))
    
    tr2 <- rtree(10)
    tr2$edge.length <- rlnorm(length(tr2$edge.length), meanlog = -4.5, sd = 0.3)
    s2 <- simSeq(tr2, l = slen, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.3, 0.4, 0.2, 0.1))

    tr3 <- bind.tree(tr1, tr2, where = 1)
    tr3$tip.label <- paste0('t', 1:length(tr3$tip.label))

    names(s1) <- paste0('t', 1:length(s1))
    names(s1) <- paste0('t', length(s1):(length(s1)-1 + length(s2)))

    s3 <- rbind(as.DNAbin(s1)[-length(s1), ], as.DNAbin(s2))
 # Print multinomial likelihood to show that there are differences in the number of site patterns.   
#    print(c(multlik(s3), multlik(as.DNAbin(simSeq(tr3, bf = c(0.1, 0.2, 0.4, 0.3))))))
    
    return(list(tr3, s3))
}


gc_test <- function(dna_data, parallel = F, nsims = 10, rm_gaps = TRUE, model = NULL){
    require(phangorn)

    print(dna_data)

    concat_list <- function(c_list){
        if(length(c_list) == 2 ){
            return(c(c_list[[1]], c_list[[2]]))
        }else if(length(c_list) > 2){
            return(c(c_list[[1]], concat_list(c_list[-1])))
        }else{
            return(c_list)
        }
    }

    multlik <- function(al){
        if(class(al) != 'DNAbin') al <- as.DNAbin(al)
        if(!is.matrix(al)) stop('Please supply the sequences as a matrix')
        nsites <- ncol(al)
        al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
        return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
    }

    rem_gaps <- function(dna_data){
        if(!is.matrix(dna_data)){
            dna_data <- as.DNAbin(dna_data)
        }
        has_gap <- function(x){
            return(!(any(c('a', 'c', 'g', 't') %in% x) & !any(c('-', 'n') %in% x) & length(unique(x)) <= 4))
        }
        gap_sites <- sapply(1:ncol(dna_data), function(d) has_gap(as.character(dna_data[, d])))
        return(dna_data[, !gap_sites])
    }

    start_tree <- nj(dist.ml(phyDat(dna_data)))

    if(rm_gaps){
        dna_data <- phyDat(rem_gaps(dna_data))
    }

    if(is.null(model)) stop('please specify a substitution model. It sholud be one of JC, GTR, GTR+G, GTR+G+I')

    if(model == 'JC'){
        opt_data <- function(dna_data) optim.pml(pml(tree = start_tree, data = dna_data, k = 1, model = 'JC'), optNni = T, optBf = F, optQ = F, model = 'JC')
    }else if(model == 'GTR'){
        opt_data <- function(dna_data) optim.pml(pml(tree = start_tree, data = dna_data, k = 1, model = 'GTR'), optNni = T, optBf = T, optQ = T, model = 'GTR')
    }else if(model == 'GTR+G'){
        opt_data <- function(dna_data) optim.pml(pml(tree = start_tree, data = dna_data, k = 4, model = 'GTR'), optNni = T, optBf = T, optQ = T, optGamma = T, model = 'GTR')
    }else  if(model == 'GTR+G+I'){
        opt_data <- function(dna_data) optim.pml(pml(tree = start_tree, data = dna_data, k = 4, model = 'GTR'), optNni = T, optBf = T, optQ = T, optGamma = T, model = 'GTR', optInv = T)
    }else{
        stop('please specify a substitution model. It sholud be one of JC, GTR, GTR+G, GTR+G+I')
    }

    s_len <- length(as.DNAbin(dna_data)[1, ])

    get_sim_rep <- function(mle){
        if(model == 'JC'){
            sim_dat <- simSeq(mle$tree, l =  s_len)
            sim_opt <- opt_data(sim_dat)
        }else if(model == 'GTR'){
            sim_dat <- simSeq(mle$tree, l = s_len, bf = mle$bf, mle$Q)
            sim_opt <- opt_data(sim_dat)
        }else if(model == 'GTR+G'){
                    rates = phangorn:::discrete.gamma(mle$shape, k = 4)
                    sim_dat_all<- lapply(rates, function(r) simSeq(mle$tree, l = round(s_len/4, 0), Q = mle$Q, bf = mle$bf, rate = r))
                    sim_dat <- concat_list(sim_dat_all)
                    sim_opt <- opt_data(sim_dat)
                }else if(model == 'GTR+G+I'){
                    rates = phangorn:::discrete.gamma(mle$shape, k = 4)
                    ninv <- mle$inv * s_len
                    nvar <- s_len - ninv
#
                    print(c(ninv, nvar))
#
                    sim_dat_var <- lapply(rates, function(r) simSeq(mle$tree, l = round(nvar/6, 0), Q = mle$Q, bf = mle$bf, rate = r))
                    sim_dat_inv <- simSeq(mle$tree, l = round(ninv, 0), Q = mle$Q, bf = mle$bf, rate = 0.0005)
                    sim_dat_all <- c(sim_dat_var, list(sim_dat_inv))
                    sim_dat <- concat_list(sim_dat_all)


                    sim_opt <- opt_data(sim_dat)
                }else{
                    stop('please specify a substitution model. It sholud be one of JC, GTR, GTR+G, GTR+G+I')
                }
 #
        print(sim_dat)
        print(sim_opt)
#
        return(multlik(sim_dat) - sim_opt$logLik)
    }

    mle_dna_data <- opt_data(dna_data)
    unc_lik_dna_data <- multlik(dna_data)
    con_lik_dna_data <- mle_dna_data$logLik
    delta_stat <- unc_lik_dna_data - con_lik_dna_data

#
    print(mle_dna_data)
#

    if(parallel){
        require(foreach)
        require(doParallel)
        cl <- makeCluster(6)
        registerDoParallel(cl)
        delta_sims <- foreach(x = 1:nsims, .packages = c('phangorn', 'ape'), .combine = c) %dopar% get_sim_rep(mle_dna_data)
        stopCluster(cl)
    }else{
        delta_sims<- vector()
        for(i in 1:nsims){
            delta_sims[i] <- get_sim_rep(mle_dna_data)
        }
    }

    return(list(delta_stat, delta_sims, sum(delta_stat > delta_sims) / nsims, mle_dna_data))
}





# SImulate a toy data set under GTR+G


    concat_list <- function(c_list){
        if(length(c_list) == 2 ){
            return(c(c_list[[1]], c_list[[2]]))
        }else if(length(c_list) > 2){
            return(c(c_list[[1]], concat_list(c_list[-1])))
        }else{
            return(c_list)
        }
    }


    multlik <- function(al){
        if(class(al) != 'DNAbin') al <- as.DNAbin(al)
        if(!is.matrix(al)) stop('Please supply the sequences as a matrix')
        nsites <- ncol(al)
        al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
        return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
    }


library(phangorn)
library(methods)

if(F){
#Simulate JC
jc_jc <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    sim_data <- simSeq(fix_tree, l = 1000)
    jc_jc[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'JC'), fix_tree)
    cat('Completed simulation replicate', i, 'for JC_JC\n')
}
save(jc_jc, file = 'jc_jc.Rdata')

jc_gtr <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    sim_data <- simSeq(fix_tree, l = 1000)
    jc_gtr[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'GTR'), fix_tree)
    cat('Completed simulation replicate', i, 'for JC_GTR\n')
}
save(jc_gtr, file = 'jc_gtr.Rdata')

jc_gtr_g <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    sim_data <- simSeq(fix_tree, l = 1000)
    jc_gtr_g[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'GTR+G'), fix_tree)
    cat('Completed simulation replicate', i, 'for JC_GTR+G\n')
}
save(jc_gtr_g, file = 'jc_gtr_g.Rdata')

# Simulate GTR
gtr_jc <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    sim_data <- simSeq(fix_tree, l = 1000, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3))
    gtr_jc[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'JC'), fix_tree)
    cat('Completed simulation replicate', i, 'for GTR_JC\n')

}
save(gtr_jc, file = 'gtr_jc.Rdata')

gtr_gtr <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    sim_data <- simSeq(fix_tree, l = 1000, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3))
    gtr_gtr[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'GTR'), fix_tree)
    cat('Completed simulation replicate', i, 'for GTR_GTR\n')
}
save(gtr_gtr, file = 'gtr_gtr.Rdata')


gtr_gtr_g <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    sim_data <- simSeq(fix_tree, l = 1000, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3))
    gtr_gtr_g[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'GTR+G'), fix_tree)
    cat('Completed simulation replicate', i, 'for GTR_GTR+G\n')
}
save(gtr_gtr_g, file = 'gtr_gtr_g.Rdata')


# Simulate GTR+G
gtrg_jc <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    gamma_cats <- phangorn:::discrete.gamma(1, 4)
    sim_data <- concat_list(lapply(gamma_cats, function(r) simSeq(fix_tree, l = 250, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3), rate = r)))
    gtrg_jc[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'JC'), fix_tree)
    cat('Completed simulation replicate', i, 'for GTR+G_JC\n')
}
save(gtrg_jc, file = 'gtrg_jc.Rdata')

gtrg_gtr <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    gamma_cats <- phangorn:::discrete.gamma(1, 4)
    sim_data <- concat_list(lapply(gamma_cats, function(r) simSeq(fix_tree, l = 250, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3), rate = r)))
    gtrg_gtr[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'GTR'), fix_tree)
    cat('Completed simulation replicate', i, 'for GTR+G_GTR\n')
}
save(gtrg_gtr, file = 'gtrg_gtr.Rdata')

gtrg_gtrg <- list()
for(i in 1:10){
    fix_tree <- rtree(50)
    fix_tree$edge.length <- rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    gamma_cats <- phangorn:::discrete.gamma(1, 4)
    sim_data <- concat_list(lapply(gamma_cats, function(r) simSeq(fix_tree, l = 250, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3), rate = r)))
    gtrg_gtrg[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 100, model = 'GTR+G'), fix_tree)
    cat('Completed simulation replicate', i, 'for GTR+G_GTR+G\n')
}
save(gtrg_gtrg, file = 'gtrg_gtrg.Rdata')


}

long_tree_gtrg <- list()
for(i in 1){
    fix_tree <- rtree(50)
    fix_tree$edge.length <-  rlnorm(n = length(fix_tree$edge.length), meanlog = -4.5, sd = 0.3)
    gamma_cats <- phangorn:::discrete.gamma(0.05, 4)
    sim_data <- concat_list(lapply(gamma_cats, function(r) simSeq(fix_tree, l = 250, Q = c(1.34, 4.81, 0.93, 1.24, 5.56, 1), bf = c(0.1, 0.2, 0.4, 0.3), rate = r)))

    long_tree_gtrg[[i]] <- list(gc_test(sim_data, parallel = F, nsims = 10, model = 'GTR+G'), fix_tree)
    cat('Completed simulation replicate', i, 'for GTR+G+long tree_GTR+G\n')
}
save(long_tree_gtrg, file = 'gtrg_gtrg.Rdata')

stop('run long tree')

hetero_gtrg <- list()
for(i in 1:3){
    sim_het <- simulate_hetero(ntax = 50, slen = 1000)
    fix_tree <- sim_het[[1]]
    sim_data <- sim_het[[2]]
    hetero_gtrg[[i]] <- list(gc_test(sim_data, parallel = T, nsims = 10, model = 'GTR+G'), fix_tree)
    cat('Completed simulation replicate', i, 'for hetero GTR+G\n')
} 
