
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

plot_drt <- function(true_data, rand_log_files, out_plot_name){
    require(ggplot2)

    rand_log_tails <- list()
    for(i in 1:length(rand_log_files)){
        rand_data_temp <- rand_log_files[[i]]
        rand_log_tails[[i]]<- rand_data_temp[(length(rand_data_temp) - 1000):length(rand_data_temp)]
    }

    complete_frame <- data.frame(cbind(true_dat[(length(true_dat)-1000):length(true_dat)], c_list(rand_log_tails)))
    colnames(complete_frame) <- c('true_dat', paste0('randomisation_', 1:length(rand_log_tails)))

    plot_data <- rbind(mean_vals = colMeans(complete_frame), sapply(1:ncol(complete_frame), function(x) quantile(complete_frame[, x], c(0.025, 0.975))),  loc_vals = c(1, seq(from = 2.1, to = 3, by = 1/ length(rand_log_tails))) )
    plot_data <- t(plot_data)
    plot_data <- data.frame(plot_data, cols = c('black', rep('grey', nrow(plot_data) -1 )))
    colnames(plot_data) <- c('mean_vals', 'lower', 'upper', 'loc_vals', 'cols')

    plot_dat_1 <- ggplot(plot_data, aes(x = loc_vals, y = mean_vals, colour = cols)) + geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper, width = .01)) + theme_bw() + scale_color_manual(values = c('black', 'grey')) + guides(colour = FALSE) + ylab(expression(paste('Substitution rate (', log[10], ' subs/site/year)')))
    print(plot_dat_1)

    if(!missing(out_plot_name)){
        pdf(out_plot_name, useDingbats = F)
        plot_dat_1
        dev.off()
    }
}


#rand_log_files <- lapply(dir('.', pattern = 'random.+log'), function(x) log10(read.table(x, head = T)$ucld.mean))
#true_dat <- log10(read.table('Heinze_TBF.log', head = T)$ucld.mean)

#plot_drt(true_dat, rand_log_files)
