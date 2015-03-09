
#These data sholud fail the test according to CR2 but not CR1
#true_data <- rnorm(100, 0.05, 0.5)
#reps_1 <- lapply(1:5, function(x) rnorm(100, 0.01, 0.01))
#reps_2 <- lapply(1:5, function(x) rnorm(100, 0.1, 0.001))
#reps_3 <- lapply(1:5, function(x) rnorm(100, 0.2, 0.001))
#list_rands<- c(reps_1, reps_2, reps_3)

calc_drt <- function(true_data, list_rands){

    quant_reps<- sapply(list_rands, function(x) quantile(x, c(0.025, 0.975))) # Should return a 2 X length(list_rands) matrix
    cr1 <- sapply(1:ncol(quant_reps), function(x) mean(true_data) > quant_reps[2, x] | mean(true_data) < quant_reps[1, x] )
    cr2 <- sapply(1:ncol(quant_reps), function(x) quantile(true_data, 0.025) > quant_reps[2, x] | quantile(true_data, 0.975) < quant_reps[1, x])

    return( list( drt_results = c(CR1 = all(cr1), CR2 = all(cr2)), true_data = c(mean(true_data), quantile(true_data, c(0.025, 0.975))), rands = quant_reps))
}

#calc_drt(true_data = true_data, list_rands = rep_list)

#test_list <- list(rnorm(10), rnorm(10, 1), rnorm(10, 3), rnorm(10, 55))
c_list <- function(mat_list){
    if(length(mat_list) == 2){
        return(cbind(mat_list[[1]], mat_list[[2]]))
    }else{
        return(cbind(mat_list[[1]],  c_list(mat_list[-1])))
    }
}
#c_list(test_list)

# To run: 
	# Make a list with the columns of interest from the randomised analyses. 
	# Pull out a vector of the column of interest from the true data set.
	
	