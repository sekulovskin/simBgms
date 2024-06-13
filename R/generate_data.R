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
#' @param ordinal_sim_type The mechanism for obtaining the parameters for simulating ordinal data. It can be \code{irt}, so that maximum pseudolikelihood estimates are obtained from an IRT data set (this is the recommended option if you would like to control for different sparsity levels of the graph), \code{dataset}, so that the user provides a data set to obtain the maximum pseudolikelihood estimates, or \code{parameters}, so that the user provides the thresholds and interactions to simulate the data (specified through the \code{thrseholds} and \code{interactions} arguments). Defaults to \code{irt}.
#' @param induce_sparsity If TRUE then it induces sparsity in the network structure based on the values provided in the \code{density} argument. Only relevant if \code{ordinal_sim_type = 'irt'}. Defaults to TRUE.
#' @param density The density/sparsity of the adjacency matrix used to generate the data sets (only relevant if \code{induce_sparsity = TRUE}).
#' @param dataset  In case \code{level = "Discrete"} and \code{ordinal_sim_type = 'dataset'}, the user can provide a dataset which will be used to generate the maximum pseudolikelihood parameters. Note that the number of variables in the provided data set needs to be equal to or larger than the highest number provided in the argument \code{no_variables}. Also, the number of ordinal categories in the provided data set must be larger than or equal to the number of ordinal categories provided in \code{no_categories}.
#' @param thresholds In case \code{level = "Discrete"} and \code{ordinal_sim_type = 'parameters'}, the user can provide a matrix containing the thresholds for the ordinal variables. The matrix should have the same number of rows as the maximum number of variables and the same number of columns as the maximum number of categories minus one.
#' @param interactions In case \code{level = "Discrete"} and \code{ordinal_sim_type = 'parameters'}, the user can provide a matrix containing the interaction parameters for the ordinal variables. The matrix should be symmetric with zeros on the diagonal.
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
#'               ordinal_sim_type = "irt",
#'               density = 0.1)
#


generate_data <- function(level = c("Gaussian", "Discrete"),
                          variable_type = "ordinal",
                          reference_category,
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
                          bc_alpha = 0.5,
                          bc_beta = 0.5){
  if(level == "Discrete"){ # for discrete variables

    if(ordinal_sim_type == "irt"){
      if(induce_sparsity == TRUE){
        sparsity_matrix <- sparsity(density, no_variables)
        gen_data <- list()
        for(i in 1:length(no_categories)) {
          gen_data[[i]] <- list()
          gen_data[[i]] <- generate_IRT_data(max(no_variables), 1e3,
                                             no_categories[i])
        }
      } else {
        density <- 1
        sparsity_matrix <- NULL
        gen_data <- list()
        for(i in 1:length(no_categories)) {
          gen_data[[i]] <- list()
          gen_data[[i]] <- generate_IRT_data(max(no_variables), 1e3,
                                             no_categories[i])
        }
      }

      # estimate the MLE parameters


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
    } else if (ordinal_sim_type == "dataset"){# end of irt
      induce_sparsity <- FALSE
      density <- 1
      sparsity_matrix <- NULL
      warning("Since you will be providing a dataset or your own, the
            induce_sparsity is automatically set to FALSE.")

        if (ncol(dataset) < max(no_variables)) {
          stop("The number of variables in the data set used to generate the
                MPLE parameter must be greater than or equal to the largest
               number in your no_variables vector.")
        } else {
          gen_data <- list()
          for(i in 1:length(no_categories)) {
            gen_data[[i]] <- list()
            gen_data[[i]] <- dataset
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
                       Please provide another data set (make sure
                       all variables in the supplied data set are
                       are either inary or ordinal.")
                  }
                })
              }
            }
          }
        }
       } else {
         induce_sparsity <- FALSE
         density <- 1
         sparsity_matrix <- NULL
        warning("Because you have specified values for the thresholds and
                interactions, sparsity is not applied to the network structure.
                Note that the number of variables in the interaction matrix
                (i.e., the number of rows and columns in the matrix) should
                be less than or equal to the largest number in the no_variables
                vector. Then, the passed parameters will be equal to the number
                of variables specified in the no_variables vector. If you've also
                specified variable_type = 'blume-capel', then the thresholds will
                be constructed based on what you've passed to the 'bc_alpha' and
                'bc_beta' arguments and the values in the argument  'thesholds' will
                be ignored.")


        if(!isSymmetric(interactions, mode = "l")) {
          stop("The matrix containing the provided interaction parameters
           should be symmetric with zeros on the diagonal.")
        }


        if(ncol(thresholds) != (max(no_categories) - 1) & nrow(thresholds) != max(no_variables)) {
          stop("The matrix containing the thresholds should have the same number of rows as the maximum number of variables
           and the same number of columns as the maximum number of categories minus one.")
        }

        MPLEs <- list()
        for(h in 1:length(no_categories)) {
          MPLEs[[h]] <- list()
          for (i in 1:length(no_variables)) {
            MPLEs[[h]][[i]] <- list()
            for (j in 1:length(density)) {
              MPLEs[[h]][[i]][[j]] <- list()
              MPLEs[[h]][[i]][[j]][["interactions"]] <- interactions[1:no_variables[i], 1:no_variables[i]]
              MPLEs[[h]][[i]][[j]][["thresholds"]] <- thresholds[1:no_variables[i], 1:(no_categories[h] - 1)]
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



  } else { #end of discrete variables

    if(induce_sparsity == TRUE){
      sparsity_matrix <- sparsity(density, no_variables)

      data <- generate_ggm_data(repetitions  = repetitions ,
                                no_observations = no_observations,
                                no_variables = no_variables,
                                density = density,
                                adj_matrix = sparsity_matrix)

      return(list("simulated data" = data[[1]],
                  "precision matrices" = data[[2]],
                  "adjacency matrices" = sparsity_matrix))


    } else{
      density <- 1
      sparsity_matrix <- NULL

      data <- generate_ggm_data(repetitions  = repetitions ,
                                no_observations = no_observations,
                                no_variables = no_variables,
                                density = density,
                                adj_matrix = sparsity_matrix)

      return(list("simulated data" = data[[1]],
                  "precision matrices" = data[[2]],
                  "adjacency matrices" = sparsity_matrix))
    }



  } # end of continuous variables


} # end of function
