#' Run Simulation Studies for Bayesian Graphical Models
#'
#' This is the main function of the package. The function simulates data from either a Gaussian or Discrete (ordinal or binary) Markov Random Field (MRF) and estimates the model parameters using Bayesian model averaging methods by utilizing the \code{R} packages \code{BDgraph} or \code{bgms} (depending on the variable type). The function provides flexibility in specifying the data characteristics as well as the edge weight and network structure priors. The function is able to parallelize the estimation of the models using multiple cores, providing computational efficiency when running large simulation studies. The output is summarized in a ready-to-plot data frame.
#'
#' @param level The type of data to be simulated. It can be either "Gaussian" or "Discrete". The choice of which model (and package) for model estimation depends on this argument. If \code{level = Gaussian} then the package \code{BDgraph} is used, and if \code{level = Discrete} then the package \code{bgms} is used.
#' @param variable_type The type of ordinal variables to be simulated if \code{level = "Discrete"}? This can be a single character string specifying the variable type of all p variables at once. Currently, the \code{R} package \code{bgms} supports \code{“ordinal”} and \code{“blume-capel”} variable types. Binary variables are automatically treated as \code{“ordinal”}. Defaults to \code{variable_type = "ordinal"}. For more details see the documentation of the \code{R} package \code{bgms}.
#' @param reference_category An integer vector of length \code{no_variables} specifying the \code{reference_category} that is used in case \code{variable_type = 'blume-capel'}. Can be any integer value between 0 and no_categories (or \code{no_categories[i]}).
#' @param repetitions  The number of repetitions in the simulation. It is recommended to always run at least 10 repetitions in order to ensure stability of the results.
#' @param no_observations The sample size of the simulated datasets.
#' @param no_variables The number of variables in the simulated data.
#' @param no_categories The number of ordinal categories (only relevant if \code{level = "Discrete"}).
#' @param ordinal_sim_type The mechanism for obtaining the parameters for simulating ordinal data. It can be \code{irt}, so that maximum pseudolikelihood estimates are obtained from an IRT data set (this is the recommended option if you would like to control for different sparsity levels of the graph), \code{dataset}, so that the user provides a data set to obtain the maximum pseudolikelihood estimates, or \code{parameters}, so that the user provides the thresholds and interactions to simulate the data (specified through the \code{thrseholds} and \code{interactions} arguments. Defaults to \code{irt}.
#' @param induce_sparsity  If \code{TRUE} the sparsity is induced on the network structure based on the values provided in the \code{density} argument. Only relevant if \code{ordinal_sim_type = 'irt'}.
#' @param density The density/sparsity of the adjacency matrix used to generate the data sets (only relevant if \code{induce_sparsity = TRUE}).
#' @param interaction_scale  The interaction_scale parameter for the interaction and edge priors.
#' @param no_cores  The number of CPU cores to use for parallel computation. The parallel computation is only used for the Bayesian analysis, and not for simulating the data. Defaults to the maximum number of cores available on the user's machine.
#' @param edge_prior The inclusion or exclusion of individual edges in the network is modeled with binary indicator variables that capture the structure of the network. The argument \code{edge_prior} is used to set a prior distribution for the edge indicator variables, i.e., the structure of the network. Currently, two options are implemented: The Bernoulli model \code{edge_prior = 'Bernoulli'} assumes that the probability that an edge between two variables is included is equal to inclusion_probability and independent of other edges or variables. When inclusion_probability = 0.5, this means that each possible network structure is given the same prior weight. The Beta-Bernoulli model \code{edge_prior = 'Beta-Bernoulli'} assumes a beta prior for the unknown inclusion probability with shape parameters \code{beta_bernoulli_alpha} and \code{beta_bernoulli_beta}. If \code{beta_bernoulli_alpha = 1} and \code{beta_bernoulli_beta = 1}, this means that networks with the same complexity (number of edges) get the same prior weight. The default is \code{edge_prior = 'Bernoulli'}. For more details see the documentation of the \code{R} package \code{bgms}.
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime prior density for the threshold parameters. Must be positive values. If the two values are equal, the prior density is symmetric about zero. If threshold_beta is greater than \code{threshold_alpha}, the distribution is skewed to the left, and if threshold_beta is less than \code{threshold_alpha}, it is skewed to the right. Smaller values tend to lead to more diffuse prior distributions.
#' @param edge_selection Should the function perform Bayesian edge selection on the edges of the MRF in addition to estimating its parameters (\code{edge_selection = TRUE}), or should it just estimate the parameters (\code{edge_selection = FALSE})? The default is \code{edge_selection = TRUE}.
#' @param inclusion_probability	 The prior edge inclusion probability for the Bernoulli model. Can be a single probability, or a matrix of p rows and p columns specifying an inclusion probability for each edge pair. The default is \code{inclusion_probability = 0.5}.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta The two shape parameters of the Beta prior density for the Bernoulli inclusion probability. Must be positive numbers. Defaults to \code{beta_bernoulli_alpha = 1} and \code{beta_bernoulli_beta = 1}.
#' @param dataset Data frame or matrix. Only relevant if \code{level = 'Discrete'} and \code{irt = FALSE}, the user can provide a dataset which will be used to generate the maximum pseudolikelihood parameters. Note that the number of variables in the provided data set needs to be equal to or larger than the highest number provided in the argument \code{no_variables}. Also, the number of ordinal categories in the provided data set must be larger than or equal to the number of ordinal categories provided in \code{no_categories}.
#' @param save If \code{TRUE}, the samples fro the MCMC algorithm are saved from all the iterations \code{iter}. Only avialable when \code{level = 'Discrete'}.
#' @param df.prior The degrees of freedom for the prior distribution of the precision matrix (only relevant if \code{level = 'Gaussian'}).
#' @param thresholds In case level = "Discrete" and \code{ordinal_sim_type = 'parameters'}, the user can provide a matrix containing the thresholds for the ordinal variables. The matrix should have the same number of rows as the maximum number of variables and the same number of columns as the maximum number of categories minus one.
#' @param interactions In case level = "Discrete" and \code{ordinal_sim_type = 'parameters'}, the user can provide a matrix containing the interaction parameters for the ordinal variables. The matrix should be symmetric with zeros on the diagonal.
#' @param burnin The number of burn-in iterations for the MCMC algorithm.
#' @param algorithm The algorithm to use for estimating the Gaussian Graphical Model (only relevant if \code{level = 'Gaussian'}).
#' @param bc_alpha,bc_beta  \code{bc_alpha} the linear contribution of the Blume-Capel model and \code{bc_beta} the quadratic contribution. Defaults to 0.5 and 0.5, respectively.
#' @return A list containing: (i) a data frame with the summarized results (only if \code{save = FALSE}), averaged across \code{repetitions}; (ii) a list with the estimated models; (iii) the simulated data sets; (iv) the maximum pseudolikelihood parameters or precision matrices; (v) the graph adjacency matrices; and (vi) a data frame with the parameter grid for an overview of the design of the simulation.
#' @param iter The number of iterations for the MCMC algorithm.
#' @details This function integrates the simulation and estimation process for Bayesian Graphical Models. It provides options for parallel computation using multiple CPU cores.
#'
#' @export
#'
#' @examples
#' # An example with an Ordinal MRF model
#' sim_bgm(level = "Discrete",
#'         variable_type = "ordinal",
#'         repetitions = 2,
#'         no_observations = 100,
#'         no_variables = c(5, 6),
#'         no_categories = 3,
#'         ordinal_sim_type = "irt",
#'         induce_sparsity = TRUE,
#'         density = 0.5,
#'         interaction_scale = 1,
#'         no_cores = 2,
#'         iter = 1e4)

sim_bgm <- function(level = c("Gaussian", "Discrete"),
                    variable_type = "ordinal",
                    reference_category = 1,
                    repetitions,
                    no_observations,
                    no_variables,
                    no_categories,
                    ordinal_sim_type = c("irt", "dataset", "parameters"),
                    induce_sparsity = TRUE,
                    density = 0.5,
                    dataset = NULL,
                    thresholds = NULL,
                    interactions = NULL,
                    interaction_scale = 2.5,
                    no_cores = NULL,
                    edge_prior = "Bernoulli",
                    threshold_alpha = 0.5,
                    threshold_beta = 0.5,
                    edge_selection = TRUE,
                    inclusion_probability = 0.5,
                    beta_bernoulli_alpha = 1,
                    beta_bernoulli_beta = 1,
                    save = FALSE,
                    df.prior = 3,
                    algorithm = "bdmcmc",
                    bc_alpha = 0.5,
                    bc_beta = -0.5,
                    iter = 1e4,
                    burnin = 1e3,
                    average = TRUE) {

  # Set default number of cores if not specified by user
  if(is.null(no_cores)) {
    no_cores <- parallel::detectCores()
  }


  if(induce_sparsity == TRUE){
    density <- density
  } else {
    density <- 1
  }

  if (!(level %in% c("Gaussian", "Discrete"))) {
    stop("The level of measurement must be either 'Discrete' or 'Gaussian'")
  } else {
    if (level == "Discrete") {

      data <- generate_data(level = "Discrete",
                            variable_type = variable_type,
                            reference_category = reference_category,
                            repetitions  = repetitions,
                            no_observations = no_observations,
                            no_variables = no_variables,
                            no_categories = no_categories,
                            ordinal_sim_type = ordinal_sim_type,
                            induce_sparsity = induce_sparsity,
                            density = density,
                            dataset = dataset,
                            thresholds = thresholds,
                            interactions = interactions,
                            bc_alpha = bc_alpha,
                            bc_beta = bc_beta)

      est <- estimate_models(level = "Discrete",
                             repetitions  = repetitions ,
                             data = data[[1]],
                             no_observations = no_observations,
                             no_variables = no_variables,
                             no_categories = no_categories,
                             interaction_scale = interaction_scale,
                             iter = iter,
                             induce_sparsity = induce_sparsity,
                             density = density,
                             no_cores = no_cores,
                             edge_prior = edge_prior,
                             save = save,
                             threshold_alpha = threshold_alpha,
                             threshold_beta = threshold_beta,
                             edge_selection = edge_selection,
                             inclusion_probability = inclusion_probability,
                             beta_bernoulli_alpha = beta_bernoulli_alpha,
                             beta_bernoulli_beta = beta_bernoulli_beta,
                             burnin = burnin)


      if(save == TRUE){
        warning("The samples from the MCMC algorithm are saved from all the iterations,
                 therefore, the results will not be summarized. If you wish to obtain
                 summarized results, set the argument save = FALSE.")

        output <- list("estimated models" = est[[1]],
                       "simulated data" = data[[1]],
                       "mple parameters" = data[[2]],
                       "adjacency matrices" = data[[3]],
                       "parameter grid" = est[[2]])

        class(output) <- "simBgms"

      } else{

        summary <- summarize(est = est[[1]],
                             level = "Discrete",
                             repetitions  = repetitions,
                             no_observations = no_observations,
                             no_variables = no_variables,
                             no_categories = no_categories,
                             interaction_scale = interaction_scale,
                             density = density,
                             adj_mat = data[[3]],
                             average = average)

        output <- list("summarized results" = summary,
                       "estimated models" = est[[1]],
                       "simulated data" = data[[1]],
                       "mple parameters" = data[[2]],
                       "adjacency matrices" = data[[3]],
                       "parameter grid" = est[[2]])

        class(output) <- "simBgms"


      }

      return(output)

    } else {

      data <- generate_data(level = "Gaussian",
                            repetitions  = repetitions ,
                            no_observations = no_observations,
                            no_variables = no_variables,
                            density = density)

      est <- estimate_models(level = "Gaussian",
                             repetitions  = repetitions ,
                             data = data[[1]],
                             no_observations = no_observations,
                             no_variables = no_variables,
                             iter = 1e4,
                             df.prior = df.prior,
                             algorithm = algorithm,
                             density = density,
                             no_cores = no_cores,
                             inclusion_probability = inclusion_probability,
                             burnin = burnin,
                             save = save)


      if(save == TRUE){
        warning("The samples from the MCMC algorithm are saved from all the iterations,
                 therefore, the results will not be summarized. If you wish obtain
                 summarized results, set save = FALSE.")

        output <- list("estimated models" = est[[1]],
                       "simulated data" = data[[1]],
                       "precision matrices" = data[[2]],
                       "adjacency matrices" = data[[3]],
                       "parameter grid" = est[[2]])
      } else{

        summary <- summarize(est = est[[1]],
                             level = "Gaussian",
                             repetitions  = repetitions,
                             no_observations = no_observations,
                             no_variables = no_variables,
                             no_categories = no_categories,
                             interaction_scale = interaction_scale,
                             density = density,
                             adj_mat = data[[3]],
                             average = average)

        output <- list("summarized results" = summary,
                       "estimated models" = est[[1]],
                       "simulated data" = data[[1]],
                       "precision matrices" = data[[2]],
                       "adjacency matrices" = data[[3]],
                       "parameter grid" = est[[2]])

        class(output) <- "simBgms"

      }
      return(output)
    }
  }
}
