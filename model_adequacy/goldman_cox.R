# Using SRD06 model. This is two partitions for 1+2 and 3 codon postions. Both use HKY+G with 4 categories
sim_goldman_cox <- function(dna_data, parallel = F, nsims = 10, rm_gaps = TRUE){
    require(phangorn)

    multlik <- function(al){
        if(class(al) != 'DNAbin') al <- as.DNAbin(al)
        if(!is.matrix(al)) stop('Please supply the sequences as a matrix')
        nsites <- ncol(al)
        al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
        return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
    }

    rem_gaps <- function(dna_data){
    	     has_gap <- function(x){
	        return(!(any(c('a', 'c', 'g', 't') %in% x) & !any(c('-', 'n') %in% x) & length(unique(x)) <= 4))
	     }

	     gap_sites <- sapply(1:ncol(dna_data), function(d) has_gap(as.character(dna_data[, d])))
	     return(dna_data[, !gap_sites])
    }

    start_tree <- nj(dist.ml(phyDat(dna_data)))

# Excluding gaps, as suggested by Ripplinger and Sullivan (). This should be done after selecting codon positions to avoid shifts in the readong frame.
    if(rm_gaps){
	cp12 <- phyDat(rem_gaps(dna_data[, -(seq(from = 3, to = ncol(dna_data), by = 3))]) )
    	cp3 <- phyDat(rem_gaps(dna_data[, seq(from = 3, to = ncol(dna_data), by = 3)]) )
    }else{
	cp12 <- phyDat(dna_data[, -(seq(from = 3, to = ncol(dna_data), by = 3))])
	cp3 <- phyDat(dna_data[, seq(from = 3, to = ncol(dna_data), by = 3)])
    }

    opt_data <- function(dna_data) optim.pml(pml(tree = start_tree, data = dna_data, k = 4, model = 'HKY'), optNni = T, optBf = T, optQ = T, optGamma = T, model = 'HKY')

    get_sim_rep <- function(mle, data_phydat){
        rates = phangorn:::discrete.gamma(mle$shape, k = 4)
        sim_dat_all<- lapply(rates, function(r) simSeq(mle$tree, l = round(length(data_phydat[[1]])/4, 0), Q = mle$Q, bf = mle$bf, rate = r))
        sim_dat <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
        sim_opt <- opt_data(sim_dat)
        return(multlik(sim_dat) - sim_opt$logLik)
    }

# Test models for cp12
    ml_cp12 <- opt_data(cp12)
    unc_lik_cp12 <- multlik(cp12)
    con_lik_cp12 <- ml_cp12$logLik

    if(parallel){
        require(foreach)
        require(doParallel)
        cl <- makeCluster(10)
        registerDoParallel(cl)
        sims_cp12 <- foreach(x = 1:nsims, .packages = c('phangorn', 'ape'), .combine = c) %dopar% get_sim_rep(ml_cp12, cp12)
        stopCluster(cl)
    }else{
        sims_cp12 <- vector()
        for(i in 1:nsims){
            sims_cp12[i] <- get_sim_rep(ml_cp12, cp12)
        }
    }
# Then do the same for cp3
    ml_cp3 <- opt_data(cp3)
    unc_lik_cp3 <- multlik(cp3)
    con_lik_cp3 <- ml_cp3$logLik

    if(parallel){
        cl <- makeCluster(10)
        require(doParallel)
        registerDoParallel(cl)
        sims_cp3 <- foreach(x = 1:nsims, .packages  = c('phangorn', 'ape'), .combine = c) %dopar% get_sim_rep(ml_cp3, cp3)
	stopCluster(cl)
    }else{
        sims_cp3 <- vector()
        for(i in 1:nsims){
            sims_cp3[i] <- get_sim_rep(ml_cp3, cp3)
        }
    }
    return(list(unc_lik_cp12 - con_lik_cp12, sims_cp12, unc_lik_cp3 - con_lik_cp3, sims_cp3))
}



#dna_data<- read.dna('Heinze_dna.fasta', format = 'fasta')
#tr <- rtree(10)
#sim_1 <- as.DNAbin(simSeq(tr, 500, rate = 1))
#sim_2 <- as.DNAbin(simSeq(tr, 500, rate = 1))

#pb_test1<- sim_goldman_cox(cbind(sim_1, sim_2), parallel = T)

