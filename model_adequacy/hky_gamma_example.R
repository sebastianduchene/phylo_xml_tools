
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

# EXAMPLE OF USING GTR+G
  # THE Q MATRIX CAN BE EXTRACTED FROM THE TRANSITION RATES.
                  Q <- c(logdat$ac[i], logdat$ag[i], logdat$at[i], logdat$cg[i], 1, logdat$gt[i])
  # BASE FREQUENCIES ARE ALSO EXTRACTED
                  bf <- c(logdat$frequencies1[i], logdat$frequencies2[i], logdat$frequencies3[i], logdat$frequencies4[i])
  # TO OBTAIN THE GAMMA CATETORIES USE THE phangorn:::discrete.gamma function. THIS RETURS THE 
  # NUMBER OF VALUES SPECIFIED IN k.
                  gamma_cats <- phangorn:::discrete.gamma(alpha = logdat$alpha[i], k = 4)
                  print('Using GTR+G')
                  print('loaded model parameters as follows:')
                  print(c(Q, bf, gamma_cats))
                  s_len <- round(l / 4, 0)
  # IF THE NUMBER OF GAMMA CATEGORIES IS 4, THEN SIMULATE 4 ALIGNMENTS WITH 1/4 THE LENGTH OF THE ORIGINAL
  # ALIGNMENT. NOTE THE IF ANY OF THESE ARE 0, THE simSeq FUNCTION WILL CRASH.
                  aln_1 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = s_len, rate = gamma_cats[1])
                  aln_2 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = s_len, rate = gamma_cats[2])
                  aln_3 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = s_len, rate = gamma_cats[3])
                  aln_4 <- simSeq(sim[[i]][[1]], Q = Q, bf = bf, l = s_len, rate = gamma_cats[4])
  # CONCATENATE THE AIGNMENTS INTO AS phyDat OBJECT.
	      	  sim[[i]][[3]] <- c(aln_1, aln_2, aln_3, aln_4)

	      }else if(subsmod == 'HKY+I'){
# EXAMPLE WITH HKY
  #THE FOLLOWING CODE SHOWS HOW TO SPECIFY Q FROM KAPPA.
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
