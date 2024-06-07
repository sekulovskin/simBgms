#' Generate Data from a Markov Random Field Model
#'
#' This function generates data from a Markov Random Field (MRF) model, which can either be Gaussian or Discrete. For the Discrete case, the function can generate data from an Item Response Theory (IRT) model if specified. The generated data can be used for various simulation studies.
#'
#' @param level The type of data to be simulated. It can be either "Gaussian" or "Discrete".
#' @param variable_type The type of ordinal variables to be simulated if \code{level = "Discrete"}? This can be a single character string specifying the variable type of all p variables at once. Currently, the \code{R} package \code{bgms} supports \code{“ordinal”} and \code{“blume-capel”} variable types. Binary variables are automatically treated as \code{“ordinal”}. Defaults to \code{variable_type = "ordinal"}. For more details see the documentation of the \code{R} package \code{bgms}.
#' @param reference_category An integer vector of length \code{no_variables} specifying the \code{reference_category} that is used in case \code{variable_type = 'blume-capel'}. Can be any integer value between 0 and no_categories (or \code{no_categories[i]}).
#' @param repetitions   The number of repetitions in the simulation.
#' @param no_observations The sample size of the simulated datasets.
#' @param no_variables The number of variables in the simulated data.
#' @param no_categories The number of ordinal categories (only relevant if level = "Discrete").
#' @param induce_sparsity If TRUE then it induces sparsity in the network structure based on the values provided in the \code{density} argument.
#' @param density The density/sparsity of the adjacency matrix used to generate the data sets (only relevant if \code{induce_sparsity = TRUE}).
#' @param irt If TRUE, an IRT model is used to generate the maximum pseudolikelihood parameters used to generate the simulated data sets (defaults to TRUE).
#' @param dataset  In case level = "Discrete" and irt = FALSE, the user can provide a dataset which will be used to generate the maximum pseudolikelihood parameters. Note that the number of variables in the provided data set needs to be equal to or larger than the highest number provided in the argument \code{no_variables}. Also, the number of ordinal categories in the provided data set must be larger than or equal to the number of ordinal categories provided in \code{no_categories}.
#' @param bc_alpha,bc_beta  bc_alpha the linear contribution of the Blume-Capel model andbc_beta the quadratic contribution. Defaults to 0.5 and 0.5, respectively.
#'
#' @return A list containing the generated data, maximum pseudolikelihood estimates (MPLEs) that were used to generate the data sets, and the graph adjacency matrices based on the provided network density through the \code{density} argument.
#'
#' @details The parameters no_observations, no_variables, and no_categories can be specified as vectors to generate data for multiple scenarios in one function call.
#'
#' @export
#'
#' @examples
#' # Generate Ordinal data sets using parameters from an IRT model
#' generate_data(level = "Discrete",
#'               variable_type = "ordinal",
#'               repetitions  = 3,
#'               no_observations = 200,
#'               no_variables = 15,
#'               no_categories = 4,
#'               induce_sparsity = TRUE,
#'               density = 0.1,
#'               irt = TRUE)


generate_data <- function(level = c("Gaussian", "Discrete"),
                          variable_type = "ordinal",
                          reference_category,
                          repetitions,
                          no_observations,
                          no_variables,
                          no_categories,
                          induce_sparsity = TRUE,
                          density = 0.5,
                          irt = TRUE,
                          dataset,
                          bc_alpha,
                          bc_beta){


  if(induce_sparsity == TRUE){

    sparsity_matrix <- sparsity(density, no_variables)
    if(level == "Discrete") {
      if(irt == TRUE) {
        gen_data <- list()
        for(i in 1:length(no_categories)) {
          gen_data[[i]] <- list()
            gen_data[[i]] <- generate_IRT_data(max(no_variables), 1e3,
                                                    no_categories[i])
        }
      } else {
        if (ncol(dataset) < max(no_variables)) {
          stop("The number of variables provided in the data set used to
               generate MPLE parameters must be larger than or equal to the
               largest number of variables (no_variables) in your simulated data.")
        } else {
          gen_data <- list()
          for(i in 1:length(no_categories)) {
            gen_data[[i]] <- list()
              gen_data[[i]] <- dataset

          }
        }
      }
      MPLEs <- list()
      for(h in 1:length(no_categories)) {
        MPLEs[[h]] <- list()
        for (i in 1:length(no_variables)) {
          MPLEs[[h]][[i]] <- list()
          for (j in 1:length(density)) {
            tries <- 0
            while (tries < 10) {
              tries <- tries + 1
              tryCatch({
                MPLEs[[h]][[i]][[j]] <- mple_G(x = gen_data[[h]][, c(1:no_variables[i])],
                                               G = sparsity_matrix[[i]][[j]])
              }, error = function(e) {
                if (tries == 10) {
                  stop("Failed to compute MPLEs after 10 attempts.
                       Please supply different value(s) for the density of the
                       network or try re-running the example with a different seed.")
                } else {
                  sparsity_matrix <- sparsity(density, no_variables)
                }
              })
            }
          }
        }
      }

      data <- generate_omrf_data(repetitions  = repetitions,
                                 variable_type = variable_type,
                                 reference_category = reference_category,
                                 no_observations = no_observations,
                                 no_variables = no_variables,
                                 no_categories = no_categories, iter = 1e3,
                                 density = density,
                                 MPLEs = MPLEs,
                                 bc_alpha = bc_alpha,
                                 bc_beta = bc_beta)

      return(list("simulated data" = data[[1]],
                  "mple parameters" = data[[2]],
                  "adjacency matrices" = sparsity_matrix))
    } else {
      data <- generate_ggm_data(repetitions  = repetitions ,
                                no_observations = no_observations,
                                no_variables = no_variables,
                                density = density,
                                adj_matrix = sparsity_matrix)

      return(list("simulated data" = data[[1]],
                  "precision matrices" = data[[2]],
                  "adjacency matrices" = sparsity_matrix))

    }
  } else {
    density <- 1
    sparsity_matrix <- NULL
    if(level == "Discrete") {
      if(irt == TRUE) {
        gen_data <- list()
        for(i in 1:length(no_categories)) {
          gen_data[[i]] <- list()
          gen_data[[i]] <- generate_IRT_data(max(no_variables), 1e3,
                                             no_categories[i])
        }
      } else {
        if (ncol(dataset) < max(no_variables)) {
          stop("The number of variables provided in the data set used to
               generate MPLE parameters must be larger than or equal to the largest
               number of variables (no_variables) in your simulated data.")
        } else {
          gen_data <- list()
          for(i in 1:length(no_categories)) {
            gen_data[[i]] <- list()
            gen_data[[i]] <- dataset

          }
        }
      }
      MPLEs <- list()
      for(h in 1:length(no_categories)) {
        MPLEs[[h]] <- list()
        for (i in 1:length(no_variables)) {
          MPLEs[[h]][[i]] <- list()
          for (j in 1:length(density)) {
            tries <- 0
            while (tries < 10) {
              tries <- tries + 1
              tryCatch({
                MPLEs[[h]][[i]][[j]] <- mple(x = gen_data[[h]][, c(1:no_variables[i])])
              }, error = function(e) {
                if (tries == 10) {
                  stop("Failed to compute MPLEs after 10 attempts.
                       Please supply different value(s) for the density of the
                       network or try re-running the example with a different seed.")
                } else {
                  sparsity_matrix <- sparsity(density, no_variables)
                }
              })
            }
          }
        }
      }

      data <- generate_omrf_data(repetitions  = repetitions,
                                 variable_type = variable_type,
                                 reference_category = reference_category,
                                 no_observations = no_observations,
                                 no_variables = no_variables,
                                 no_categories = no_categories, iter = 1e3,
                                 density = density,
                                 MPLEs = MPLEs,
                                 bc_alpha = bc_alpha,
                                 bc_beta = bc_beta)

      return(list("simulated data" = data[[1]],
                  "mple parameters" = data[[2]],
                  "adjacency matrices" = sparsity_matrix))
    } else {
      data <- generate_ggm_data(repetitions  = repetitions ,
                                no_observations = no_observations,
                                no_variables = no_variables,
                                density = density,
                                adj_matrix = sparsity_matrix)

      return(list("simulated data" = data[[1]],
                  "precision matrices" = data[[2]],
                  "adjacency matrices" = sparsity_matrix))

    }

  }
}
