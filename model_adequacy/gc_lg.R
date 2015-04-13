
gc_test_lg_g <- function(aa_data, parallel = F, nsims = 10, rm_gaps = TRUE){
    require(phangorn)

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
        if(!is.matrix(al)) al <- as.matrix(al)
        al <- as.character(al)
        nsites <- ncol(al)
        al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
        return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
    }

    rem_gaps <- function(aa_data){
        if(!is.matrix(aa_data)){
            aa_data <- as.character(aa_data)
        }
        has_gap <- function(x){
            return(any(c('-', 'n') %in% x))
        }
        gap_sites <- sapply(1:ncol(aa_data), function(d) has_gap(as.character(aa_data[, d])))
        return(aa_data[, !gap_sites])
    }

    start_tree <- nj(dist.ml(phyDat(aa_data, type = 'AA')))

    if(rm_gaps){
        aa_data <- phyDat(rem_gaps(aa_data), type = 'AA')
    }

    opt_data <- function(aa_data) optim.pml(pml(tree = start_tree, data = aa_data, k = 4, model = 'LG'), optNni = T, optBf = T, optQ = F, optGamma = T, model = 'LG')

    get_sim_rep <- function(mle, data_phydat){
        rates = phangorn:::discrete.gamma(mle$shape, k = 4)
        sim_dat_all<- lapply(rates, function(r) simSeq(mle$tree, l = round(length(data_phydat[[1]])/4, 0), bf = mle$bf, rate = r, model = 'LG', type = 'AA'))
        sim_dat <- concat_list(sim_dat_all)
        sim_opt <- opt_data(sim_dat)
        return(multlik(sim_dat) - sim_opt$logLik)
    }

    print('started on empirical data')
    mle_aa_data <- opt_data(aa_data)
    unc_lik_aa_data <- multlik(aa_data)
    con_lik_aa_data <- mle_aa_data$logLik
    delta_stat <- unc_lik_aa_data - con_lik_aa_data

    print('started simulations')
    if(parallel){
        require(foreach)
        require(doParallel)
        cl <- makeCluster(4)
        registerDoParallel(cl)
        delta_sims <- foreach(x = 1:nsims, .packages = c('phangorn', 'ape'), .combine = c) %dopar% tryCatch(get_sim_rep(mle_aa_data, aa_data), error = function(x) return(0))
        stopCluster(cl)
    }else{
        delta_sims<- vector()
        for(i in 1:nsims){
            delta_sims[i] <- tryCatch(get_sim_rep(mle_aa_data, aa_data), error = function(x) return(0))
            print(paste('completed sim', i))
        }
    }

    return(list(delta_stat, delta_sims, sum(delta_stat > delta_sims) / nsims))
}





# SImulate a toy data set under LG+G



    concat_list <- function(c_list){
        if(length(c_list) == 2 ){
            return(c(c_list[[1]], c_list[[2]]))
        }else if(length(c_list) > 2){
            return(c(c_list[[1]], concat_list(c_list[-1])))
        }else{
            return(c_list)
        }
    }



library(phangorn)

fix_tree <- rtree(5)
fix_tree$edge.length <- fix_tree$edge.length / 10
gamma_cats <- phangorn:::discrete.gamma(alpha = 100, k = 4) / 2

slen <- 40

s1_list <- lapply(gamma_cats, function(x) simSeq(fix_tree, l = round(slen / length(gamma_cats), 0), rate = x, model = 'LG', type = 'AA' ))
aa_data <- concat_list(s1_list)

t1 <- gc_test_lg_g(aa_data, rm_gaps = T, parallel = T, nsims = 50)

#opt_1 <- optim.pml(pml(fix_tree, data = aa_data, k = 4, model = 'LG'), optGamma = T, optBf = T, optNni = T)

#t_fun <- gc_test_gtr_g(aa_data, parallel= T, nsims = 100)
#hist(t_fun[[2]])
#####
