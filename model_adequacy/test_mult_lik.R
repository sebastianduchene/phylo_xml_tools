    multlik <- function(al){
	if(class(al) != "DNAbin"){ al <- as.list(as.DNAbin(al)) } else { al <- as.list(al) }
	mat <- as.character(as.list(as.matrix(al))[[1]])
	for(i in 2:length(as.list(al))){
            mat <- rbind(mat, as.character(as.list(as.matrix(al))[[i]]))
	}

        no_gaps <- which(sapply(1:ncol(mat), function(x) if('-' %in% as.character(mat[,x])){ FALSE}else{TRUE}))
          if(length(no_gaps) >= ncol(mat)) stop('all sites have gaps')
        mat <- mat[, no_gaps]

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

library(ape)
source('goldman_cox.R')

aln1 <- read.dna('~/Desktop/hiv/sc1.POL.fasta', format = 'fasta')
sim_goldman_cox(aln1)
