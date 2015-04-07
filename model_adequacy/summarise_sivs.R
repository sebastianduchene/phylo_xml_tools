library(phangorn)
library(ape)
require(methods)

get_mle <- function(mle_file_name, beast_path){

    #Define command to esitmate marginal likelihoods using BEAST
    mle_cmd <- '<?xml version=\"1.0\" standalone=\"yes\"?>
    <beast>
    <pathSamplingAnalysis fileName=\"temp_mle.log\">
    <likelihoodColumn name=\"pathLikelihood.delta\"/>
    <thetaColumn name=\"pathLikelihood.theta\"/>
    </pathSamplingAnalysis>
    <steppingStoneSamplingAnalysis fileName=\"temp_mle.log\">
    <likelihoodColumn name=\"pathLikelihood.delta\"/>
    <thetaColumn name=\"pathLikelihood.theta\"/>
    </steppingStoneSamplingAnalysis>
    </beast>
    '

    #Call BEAST to estimate marginal likelihood using a xml template defined above
    system(paste('cp', mle_file_name, 'temp_mle.log'))
    cat(mle_cmd, file = 'temp_mle.xml')
    mle_est <- system(paste(beast_path, 'temp_mle.xml'), intern = TRUE)

    #Pull out marginal likelihoods from output in terminal
    path_sampling <- grep('^log marginal.+path sampling', mle_est, value = TRUE)
    path_sampling <- gsub('[A-Z]|[a-z]|[(]|[)]| |[a-z][.][a-z]|=', '', path_sampling)

    step_stone_sampling <- grep('^log marginal.+stepping stone sampling', mle_est, value = TRUE)
    step_stone_sampling <- gsub('[A-Z]|[a-z]|[(]|[)]| |[a-z][.][a-z]|=', '', step_stone_sampling)

    #Clean up
    system('rm temp_mle.xml temp_mle.log')

    #Make vector with marginal likelihods estimated using the two methods
    mles <- as.numeric(c(path_sampling, step_stone_sampling))
    names(mles) <- c('path_sampling', 'stepping_stone')
    return(mles)
}

######GET RATOGRAMS
getRatogB <- function(trees, out.file = "ratogs.tree", bug_beast1 = T){

    if(bug_beast1){
        trees_raw <- readLines(trees)
        trees_mod <- gsub('[]]([0-9]|[.])+', '', trees_raw)
        trees_mod <- gsub('[[]&rate=', '', trees_mod)
        writeLines(trees_mod, con = out.file)

    }else{
	  trs <- readLines(trees)
	  trs <- gsub("[[]&([a-z])+[=]", ":", trs)
	  trs <- gsub("[]]:([0-9]|[.])+", "", trs)
	  writeLines(trs, out.file)
      }
}

#SIMULATE DATA

simPhylogsBL <- function(treesf, logdat, N = 100, subsmod = "JC", l = 1000, sample = T, ratogindir = F, ...){
	trees <- read.nexus(treesf)
	logdat <- read.table(logdat, header = T, comment = "#" )#, sep = ",")
	samp <- sample(round(length(trees)*0.5, 0):length(trees), N)
	trees <- trees[samp]
	logdat <- logdat[samp,]
	sim <- list()
	if(any(c("ucldMean", "meanClockRate", "ucld.mean") %in% colnames(logdat))){
	if(ratogindir == T){
            ratogs <- read.nexus("ratogs.tree")[samp]
            print('I read the ratogram from the directory')
	} else {
	      getRatogB(treesf)
	      ratogs <- read.nexus("ratogs.tree")[samp]
	}
	}
	for(i in 1:nrow(logdat)){
            print(paste('simulating data set ', i))
            tr <- trees[[i]]

            if("clock.rate" %in% colnames(logdat)){
                print('in first if')
	      	     tr$edge.length <- tr$edge.length * logdat[i, "clock.rate"]
		     sim[[i]] <- list(phylogram = tr)
            } else if("ucld.mean" %in% colnames(logdat) | "rateChangeCount" %in% colnames(logdat)){
                print('in second if')
                trr <- ratogs[[i]]
                print('These are the branch rates')
                print(trr$edge.length)
                print('These are the branch times')
                print(tr$edge.length)
                trr$edge.length <- trr$edge.length * tr$edge.length
                print('These are the phylogram branch lengths')
                print(trr$edge.length)
                sim[[i]] <- list(phylogram = trr)
                print(sim[[i]])
                print('processed phylogram correctly')
	      }
	      if(subsmod == "JC"){
	      	     sim[[i]][[3]] <- simSeq(sim[[i]][[1]], l = l)
	      } else if(subsmod == "GTR+G"){
                  Q <- c(logdat$ac[i], logdat$ag[i], logdat$at[i], logdat$cg[i], 1, logdat$gt[i])
                  bf <- c(logdat$frequencies1[i], logdat$frequencies2[i], logdat$frequencies3[i], logdat$frequencies4[i])
                  gamma_cats <- phangorn:::discrete.gamma(alpha = logdat$alpha[i], k = 4)
                  print('Using GTR+G')
                  print('loaded model parameters as follows:')
                  print(c(Q, bf, gamma_cats))
                  s_len <- round(l / 4, 0)
                  aln_1 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = s_len, rate = gamma_cats[1])
                  aln_2 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = s_len, rate = gamma_cats[2])
                  aln_3 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = s_len, rate = gamma_cats[3])
                  aln_4 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = s_len, rate = gamma_cats[4])
	      	  sim[[i]][[3]] <- c(aln_1, aln_2, aln_3, aln_4)

	      }else if(subsmod == 'HKY+I'){
                  Q <- c(1, 2*logdat$kappa[i], 1, 1, 2*logdat$kappa[i], 1)
                  bf <- c(logdat$frequencies1[i], logdat$frequencies2[i], logdat$frequencies3[i], logdat$frequencies4[i])
                  pinv <- logdat$pInv[i]
                  print('Using HKY+I')
                  print('loaded parameters as follows:')
                  print(c(Q, bf, pinv))
                  inv_len <- round(l * pinv, 0)
                  var_len<- l - inv_len
                  print('invariable sites')
                  print(inv_len)
                  print('variable sites')
                  print(var_len)
#                  return(list(sim[[i]][[1]], Q, bf, inv_len))
#                  break
                  aln_1 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = inv_len, rate = 0.05)
                  print('simulated invariable sites')
                  aln_2 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = var_len, rate = 1)
                  sim[[i]][[3]] <- c(aln_1, aln_2)
              }
	}
	return(sim)
}



#GET MULTINOMIAL LIKELIHOODS

multlik <- function(al){
	if(class(al) != "DNAbin"){ al <- as.list(as.DNAbin(al)) } else { al <- as.list(al) }
	mat <- as.character(as.list(as.matrix(al))[[1]])
	for(i in 2:length(as.list(al))){
	      mat <- rbind(mat, as.character(as.list(as.matrix(al))[[i]]))
	}
	al <- mat
	nsites <- ncol(al)
	usites <- unique(al, MARGIN = 2)
	liks <- 0
	for(i in 1:ncol(usites)){
	      insts <- sum(apply(al, 2, identical, usites[, i]))
	      liks <- liks + log((insts / nsites)^insts)
	      #print(liks)
	}
	return(liks)
}



# This function takes a chronogram and a phylogram, re-roots the phylogram appropriately, and returns both ladderized.
rootnlad <- function(chron, phyl){

	 rootnode <- length(chron$tip.label) + 1
	 outgr <- chron$tip.label[Descendants(chron, Descendants(chron, rootnode, "children")[1], "tips")[[1]]]
	 chron <- ladderize(root(chron, outgroup = outgr))
	 phyl <- ladderize(root(phyl, outgroup = outgr))

	 return(list(chron, phyl))
}

# This function processes simulated datasets created with the function simPhylogsBL.R by taking posterior phylograms and simulated alignmetns, and creating posterior predictive simulated phylograms for each.
prepsimsBL <- function(sims, empdat, chron, writedat = F){
	 empdat <- as.phyDat(read.dna(empdat, format = 'fasta'))
	 topo <- read.tree(chron)
	 topo$edge.length <- rep(1, length(topo$edge.length))
	 emphy <- optim.pml(pml(topo, empdat, model = 'GTR', k = 4), optQ = T, optGamma = T, optBf = T, control = pml.control(maxit = 50))

         print(paste(emphy$Q, emphy$shape, sep = '_'))
	 empmlik <- multlik(empdat)
	 empbl <- emphy$tree$edge.length
	 simbl <- emphy$tree$edge.length
	 simultlik <- vector()
	 for(i in 1:length(sims)){
             cat('running simluation', i, '\n')
                                        #print(class(sims[[i]][[3]]))
             #NOTE THAT SOME SIMULATIONS DON'T WORK. SKIP ITERATIONS IF THIS HAPPENS
             ppsphy <- tryCatch(optim.pml(pml(topo, sims[[i]][[3]] , model = 'GTR'), optQ = T, optBf = T,  control = pml.control(maxit = 50)), error = function(e) NULL)
             if(is.null(ppsphy)){
                 print('skipping')
                 next
             }

             print(paste(emphy$Q, sep = '_'))
	       #trs <- rootnlad(sims[[i]][[1]], ppsphy$tree)
	       #sims[[i]][[1]] <- trs[[1]]
	       #sims[[i]][[2]] <- trs[[2]]
	       if(writedat == T){
                   write.tree(sims[[i]][[1]], paste0("postphylog", i, ".tre"))
                   write.csv(sims[[i]][[2]], paste0("ppsphylog", i, ".tre"))
                   write.nexus.data(as.list(as.DNAbin(sims[[i]][[3]])), paste0("ppsdat", i, ".nex"))
                   dat <- readLines(paste0("ppsdat", i, ".nex"))
                   dat <- gsub("DATATYPED", "DATATYPE=D", dat)
                   writeLines(dat, paste0("ppsdat", i, ".nex"))
	       }
	       simbl <- rbind(simbl, ppsphy$tree$edge.length)
	       simultlik[i] <- multlik(sims[[i]][[3]])
	       print(paste("finished sim", i))

	 }
	 simbl <- simbl[2:nrow(simbl),]
	 return(list(empbl, simbl, empmlik, simultlik))

}


# GETT SUMMARY STATISTICS
get_summ_stats <- function(pps){
    pvals <- sapply(1:length(pps[[1]]), function(x) sum(pps[[2]][, x] < pps[[1]][x]) / 100)
    pvals_out<- 1 - (sum(pvals > 0.95 | pvals < 0.05) / length(pvals))
    precision <- sapply(1:ncol(pps[[2]]), function(x) diff(quantile(pps[[2]][, x], c(0.05, 0.95))) / mean(pps[[2]][, x]))
    precision_out <- mean(precision)
    multLik <- sum(pps[[3]] > pps[[4]] | pps[[3]] < pps[[4]]) / length(pps[[4]])
    return(c(pvals = pvals_out, precision = precision_out, multLik = multLik))
}




##############
# End defining functions
##############



## SETTING UP...



#if(F){
## SIV PPS ANALYSES
    setwd("/Users/sebastianduchene/Desktop/modad/HIVs/runs")
    mles_siv <- matrix(NA, 5, 3)
    colnames(mles_siv) <- c('strict', 'rlc', 'ucld')
    rownames(mles_siv) <- c('ps', 'ss', 'accuracy', 'precision', 'mult_lik')

    mles_siv[1:2, 1] <- get_mle('siv_2_strict.mle.log', '~/Desktop/phylo_programs/beast181/bin/beast')
    mles_siv[1:2, 2] <- get_mle('siv_2_rlc.mle.log', '~/Desktop/phylo_programs/beast181/bin/beast')
    mles_siv[1:2, 3] <- get_mle('siv_2_ucl.mle.log', '~/Desktop/phylo_programs/beast181/bin/beast')

    sim_strict <- simPhylogsBL(treesf = 'siv_2_strict.trees', logdat = 'siv_2_strict.log', subsmod = 'GTR+G', )
   pps_strict <- prepsimsBL(sims = sim_strict, empdat = 'siv_2.fasta', chron = 'out_tre.tree')
    pvals_strict <- sapply(1:length(pps_strict[[1]]), function(x) sum(pps_strict[[2]][, x] < pps_strict[[1]][x]) / 100)

    mles_siv[3:5, 1] <- get_summ_stats(pps_strict)



sim_rlc <- simPhylogsBL(treesf = 'siv_2_rlc.trees', logdat = 'siv_2_rlc.log', subsmod = 'GTR+G')
   pps_rlc <-  prepsimsBL(sims = sim_rlc, empdat = 'siv_2.fasta', chron = 'out_tre.tree')
    pvals_rlc <- sapply(1:length(pps_rlc[[1]]), function(x) sum(pps_rlc[[2]][, x] < pps_rlc[[1]][x]) / 100)

mles_siv[3:5, 2] <- get_summ_stats(pps_rlc)

    sim_ucl <- simPhylogsBL(treesf = 'siv_2_ucl.trees', logdat = 'siv_2_ucl.log', ratogindir = F, subsmod = 'GTR+G')
    pps_ucl <-  prepsimsBL(sims = sim_ucl, empdat = 'siv_2.fasta', chron = 'out_tre.tree')
    pvals_ucl <- sapply(1:length(pps_ucl[[1]]), function(x) sum(pps_ucl[[2]][, x] < pps_ucl[[1]][x]) / 100)

    mles_siv[3:5, 3] <- get_summ_stats(pps_ucl)

    write.table(mles_siv, file = 'summ_stats_siv.txt', row.names = F)
write.table(rbind(pvals_strict, pvals_rlc, pvals_ucl), file = 'siv_pvals.txt', row.names = F)

#}
