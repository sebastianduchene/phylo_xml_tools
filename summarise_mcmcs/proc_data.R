

# Funciton diagnose_mcmc to get the ess for all parameters quantiles, means and others
diagnose_mcmc <- function(mcmc_run){
    require(coda)
    if(is.mcmc(mcmc_run)) mcmc_run <- as.data.frame(mcmc_run)
    res_matrix <- matrix(NA, ncol(mcmc_run) - 1, 4)

    for(i in 1:(ncol(mcmc_run) - 1)){
        res_matrix[i, ] <- c(mean(mcmc_run[, i + 1]), quantile(mcmc_run[, i + 1], c(0.025, 0.975)), effectiveSize(mcmc_run[, i + 1]))
    }
    rownames(res_matrix) <- colnames(mcmc_run)[2:ncol(mcmc_run)]
    colnames(res_matrix) <- c('Mean', 'Lower95%', 'Upper95%', 'ESS')
    return(res_matrix)
}

#Function to get margninal likelihoods from mle log files from beast
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

concat_list <- function(mat_list){
  if(length(mat_list) == 1){
    return(mat_list[[1]])
  }else if(length(mat_list) == 2){
    return(rbind(mat_list[[1]], mat_list[[2]]))
  }else{
    return(rbind(mat_list[[1]], concat_list(mat_list[-1])))
  }  
}


