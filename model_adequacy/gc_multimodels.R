
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
        opt_data <- function(dna_data) optim.pml(pml(tree = start_tree, data = dna_data, k = 4, model = 'GTR'), optNni = T, optBf = T, optQ = T, model = 'GTR')
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
                    pinv <- mle$inv
                    pvar <- 1 - pinv
                    sim_dat_var <- lapply(rates, function(r) simSeq(mle$tree, l = round((pvar*s_len)/4, 0), Q = mle$Q, bf = mle$bf, rate = r))
                    sim_dat_inv <- simSeq(mle$tree, l = round(s_len*pinv, 0), Q = mle$Q, bf = mle$bf, rate = 0.0005)
                    sim_dat_all <- c(sim_dat_var, list(sim_dat_inv))
                    sim_dat <- concat_list(sim_dat_all)
                    sim_opt <- opt_data(sim_dat)
                }else{
                    stop('please specify a substitution model. It sholud be one of JC, GTR, GTR+G, GTR+G+I')
                }

        return(multlik(sim_dat) - sim_opt$logLik)
    }

    mle_dna_data <- opt_data(dna_data)
    unc_lik_dna_data <- multlik(dna_data)
    con_lik_dna_data <- mle_dna_data$logLik
    delta_stat <- unc_lik_dna_data - con_lik_dna_data

    if(parallel){
        require(foreach)
        require(doParallel)
        cl <- makeCluster(4)
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

fix_tree <- rtree(10)
fix_tree$edge.length <- 5 * fix_tree$edge.length / sum(fix_tree$edge.length)
gamma_cats <- phangorn:::discrete.gamma(alpha = 1, k = 4)

slen <- 1000
s1_list <- lapply(gamma_cats, function(x) simSeq(fix_tree, l = round(slen / length(gamma_cats), 0), rate = x, model = 'GTR', bf = c(0.2, 0.1, 0.3, 0.4)))
s2_list <- c(s1_list, list(simSeq(fix_tree, l = 100, rate = 0.0005)))
dna_data <- concat_list(s2_list)

t_fun <- gc_test(dna_data, parallel= T, nsims = 100, model = 'GTR+G+I')
#hist(t_fun[[2]])

#####
