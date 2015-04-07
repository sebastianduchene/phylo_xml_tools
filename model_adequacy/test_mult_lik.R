
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
#Bollback
            liks <- (log(insts) * insts) + liks
	}
#Bollback
        liks <- liks - (nsites*log(nsites))
	return(liks)
    }

#Faster unconstrained likelihood calculation. Need to add error messages
multlik_alt <- function(al){
    nsites <- ncol(al)
    al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
    return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
}


library(ape)
library(phangorn)
al <- as.DNAbin(simSeq(rtree(50), 1000, rate = 0.1))

system.time(t1 <- multlik(al))
print(t1)

system.time(t2 <- multlik_alt(al))
print(t2)



al_patterns <- table(sapply(1:ncol(al), function(x) paste(al[, x], collapse = '')))

sum(sapply(al_patterns, function(x) (log(x) * x))) - (50*log(50))

source('goldman_cox.R')

aln1 <- read.dna('~/Desktop/hiv/sc1.POL.fasta', format = 'fasta')
sim_1 <- sim_goldman_cox(aln1)
