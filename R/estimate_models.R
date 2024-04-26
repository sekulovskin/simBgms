#' Estimate Bayesian Graphical Models
#'
#' This function estimates Bayesian Graphical Models (BGM) using either the `bgms` or `BDGraph` R packages, depending on the level of measurement specified. It supports parallelization for efficient computation on multiple CPU cores.
#'
#' @param level Character vector specifying the type of data. It can be either "Gaussian" or "Discrete".
#' @param variable_type What kind of variables should be simulated? Currently this can be a single character string specifying the variable type of all p variables at once. Currently, bgm supports “ordinal” and “blume-capel”. Binary variables are automatically treated as “ordinal’’. Defaults to variable_type = "ordinal". Only relevant if level = "Discrete".
#' @param reference_category An integer vector of length no_variables specifying the reference_category category that is used for the Blume-Capel model (details below). Can be any integer value between 0 and no_categories (or \code{no_categories[i]}).
#' @param repetitions  Numeric. The number of repetitions in the simulation.
#' @param data List. The simulated data sets.
#' @param no_observations Numeric. The sample size of the simulated datasets.
#' @param no_variables Numeric. The number of variables in the simulated data.
#' @param no_categories Numeric. The number of ordinal categories (only relevant if level = "Discrete").
#' @param interaction_scale Numeric. The interaction_scale parameter for the interaction prior (only relevant if level = "Discrete").
#' @param iter Numeric. The number of iterations for the estimation algorithm.
#' @param induce_sparsity Logical. If TRUE then it induces sparsity in the network structure based on the values provided in the \code{density} argument.
#' @param density Numeric. The density/sparsity of the adjacency matrix used to generate the data sets (only relevant if \code{induce_sparsity = TRUE}).
#' @param no_cores The number of CPU cores to use for parallel computation.
#' @param edge_prior The inclusion or exclusion of individual edges in the network is modeled with binary indicator variables that capture the structure of the network. The argument edge_prior is used to set a prior distribution for the edge indicator variables, i.e., the structure of the network. Currently, two options are implemented: The Bernoulli model edge_prior = "Bernoulli" assumes that the probability that an edge between two variables is included is equal to inclusion_probability and independent of other edges or variables. When inclusion_probability = 0.5, this means that each possible network structure is given the same prior weight. The Beta-Bernoulli model edge_prior = "Beta-Bernoulli" assumes a beta prior for the unknown inclusion probability with shape parameters beta_bernoulli_alpha and beta_bernoulli_beta. If beta_bernoulli_alpha = 1 and beta_bernoulli_beta = 1, this means that networks with the same complexity (number of edges) get the same prior weight. The default is edge_prior = "Bernoulli"
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime prior density for the threshold parameters. Must be positive values. If the two values are equal, the prior density is symmetric about zero. If threshold_beta is greater than threshold_alpha, the distribution is skewed to the left, and if threshold_beta is less than threshold_alpha, it is skewed to the right. Smaller values tend to lead to more diffuse prior distributions.
#' @param edge_selection Should the function perform Bayesian edge selection on the edges of the MRF in addition to estimating its parameters (edge_selection = TRUE), or should it just estimate the parameters (edge_selection = FALSE)? The default is edge_selection = TRUE.
#' @param inclusion_probability	 The prior edge inclusion probability for the Bernoulli model. Can be a single probability, or a matrix of p rows and p columns specifying an inclusion probability for each edge pair. The default is inclusion_probability = 0.5. Corresponds to the g.prior argument in the BDgraph package.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta The two shape parameters of the Beta prior density for the Bernoulli inclusion probability. Must be positive numbers. Defaults to beta_bernoulli_alpha = 1 and beta_bernoulli_beta = 1.
#' @param edge_prior Character. The prior distribution for the edge parameters (only relevant if level = "Discrete").
#' @param save Logical. If TRUE, the samples fro the MCMC algorithm are saved from all the iterations \code{iter}. Only possible when level = "Discrete".
#' @param df.prior Numeric. The degrees of freedom for the prior distribution of the precision matrix (only relevant if level = "Gaussian").
#' @param algorithm Character. The algorithm to use for estimating the Gaussian Graphical Model (only relevant if level = "Gaussian").
#' @param burnin Numeric. The number of burn-in iterations for the MCMC algorithm
#' @return A list containing the estimated Bayesian Graphical Models and the associated parameter gird.
#'
#' @details This function estimates BGMs using either the `bgms` or `BDGraph` R packages, depending on the specified level. It parallelizes the estimation process for efficient computation.
#'
#' @export
#'
#' @importFrom foreach %dopar%
#'
#'


estimate_models <- function(level = c("Gaussian", "Discrete"),
                            variable_type = "ordinal",
                            reference_category = 1,
                            repetitions,
                            data,
                            no_observations,
                            no_variables,
                            no_categories,
                            interaction_scale = 2.5,
                            iter = 1e4,
                            induce_sparsity = TRUE,
                            density = 0.5,
                            no_cores = NULL,
                            save = FALSE,
                            edge_prior = "Bernoulli",
                            threshold_alpha = 0.5,
                            threshold_beta = 0.5, edge_selection = TRUE,
                            inclusion_probability = 0.5,
                            beta_bernoulli_alpha = 1,
                            beta_bernoulli_beta = 1,
                            df.prior = 3,
                            algorithm = "bdmcmc",
                            burnin = 1e3) {



  # Set default number of cores if not specified by user
  if(is.null(no_cores)) {
    no_cores <- parallel::detectCores()
  }

  if(induce_sparsity == TRUE){
    density <- density
  } else {
    density <- 1
  }

  cl <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)

  if(level == "Discrete"){


    param_grid <- expand.grid(repetitions  = 1:repetitions,
                              s = 1:length(interaction_scale),
                              l = 1:length(no_categories),
                              k = 1:length(density),
                              i = 1:length(no_observations),
                              j = 1:length(no_variables)
    )


    est <- foreach::foreach(idx = seq_len(nrow(param_grid))) %dopar% {
      row_params <- param_grid[idx,]
      h <- row_params[["repetitions"]]
      l <- row_params[["l"]]
      k <- row_params[["k"]]
      i <- row_params[["i"]]
      j <- row_params[["j"]]
      s <- row_params[["s"]]

      # Estimate
      tryCatch({
        bgms::bgm(data[[h]][[i]][[j]][[k]][[l]],
                  variable_type = variable_type,
                  reference_category = reference_category,
                  interaction_scale = interaction_scale[s],
                  edge_prior = edge_prior,
                  iter = iter,
                  save = save,
                  threshold_alpha = threshold_alpha,
                  threshold_beta = threshold_beta,
                  edge_selection = edge_selection,
                  inclusion_probability = inclusion_probability,
                  beta_bernoulli_alpha = beta_bernoulli_alpha,
                  beta_bernoulli_beta = beta_bernoulli_beta,
                  burnin = burnin)
      }, error = function(e) {
        # If an error occurs, store the error message in the corresponding position of the "est" list
        paste0("Error occurred while estimating model ", idx, ": ", conditionMessage(e))
      })
    }
  } else{


    param_grid <- expand.grid(repetitions  = 1:repetitions ,
                              k = 1:length(density),
                              i = 1:length(no_observations),
                              j = 1:length(no_variables)
    )

    est <- foreach::foreach(idx = seq_len(nrow(param_grid))) %dopar% {
      row_params <- param_grid[idx,]
      h <- row_params[["repetitions"]]
      k <- row_params[["k"]]
      i <- row_params[["i"]]
      j <- row_params[["j"]]

      # Estimate
      tryCatch({
        BDgraph::bdgraph(data[[h]][[i]][[j]][[k]],
                         g.prior = inclusion_probability,
                         iter = iter,
                         algorithm = algorithm,
                         burnin = burnin,
                         df.prior = df.prior,
                         save = save)
      }, error = function(e) {
        # If an error occurs, store the error message in the corresponding position of the "est" list
        paste0("Error occurred while estimating model ", idx, ": ", conditionMessage(e))
      })
    }
  }
  # Stop parallel backend
  parallel::stopCluster(cl)



  # Return estimated models
  return(list(est, param_grid))
}
