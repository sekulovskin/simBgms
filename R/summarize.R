# Function that summarizes the output from estimate_models functiom

summarize <- function(est,
                     # average = TRUE
                      level = c("Gaussian", "Discrete"),
                      repetitions,
                      interaction_scale,
                      repetitionsetitions,
                      no_observations,
                      no_variables,
                      no_categories,
                      density = density){


  # Count the instances with error messages
  error_count <- sum(sapply(est, function(x) is.character(x) && !is.list(x)))
  if(error_count/length(est) != 0){
    warning("It seems that the models for some data sets have failed
            during the estimation process. In the post-processing, these
            non-converged models are replaced with models from the same
            simulation settings that did converge. To identify the failed models,
            inspect the list of estimated models in the output and compare
            their indices with those in the parameter grid.
            If a significant number of models have failed to converge,
            this might cause errors in the post-processing.
            If you encounter such issues, consider adjusting the
            data generating process. If the problem
            persists, please contact the package maintainer.")


    # Find indices of failed models (character elements) and succeeded models (list elements)
    indices_failed <- which(sapply(est, is.character))
    indices_succeeded <- which(sapply(est, is.list))

    # Iterate over the failed models
    for (failed_index in indices_failed) {
      # Calculate the range of indices for successful models based on the failed index
      range_start <- ((failed_index - 1) %/% repetitions) * repetitions + 1
      range_end <- range_start + 99

      # Find the appropriate successful index within the range
      successful_indices <- which(indices_succeeded >= range_start & indices_succeeded <= range_end)

      # If there are no successful indices within the range, skip this failed index
      if (length(successful_indices) == 0) next

      # Choose a random successful index within the range
      selected_successful_index <- sample(successful_indices, 1)

      # repetitionslace the failed model with the successful one
      est[[failed_index]] <- est[[indices_succeeded[selected_successful_index]]]
    }

  }


  # take out the inclusion indicators and parameters  ------------------------

  gammas <- list()
  thetas <- list()
  for(i in 1:length(est)){
    gammas[[i]] <- est[[i]][1]
    thetas[[i]] <- est[[i]][2]
  }

  # take the upper triangular part of the matrix ----------------------------

  gammas_1 <- list()
  thetas_1 <- list()
  for(i in 1:length(est)){
    gammas_1[[i]] <- gammas[[i]][[1]][upper.tri(gammas[[i]][[1]])]
    thetas_1[[i]] <- thetas[[i]][[1]][upper.tri(thetas[[i]][[1]])]
  }


  if(level == "Discrete"){ # For discrete networks


    param_grid <- expand.grid(
      s = interaction_scale,
      k = no_categories,
      l = density,
      i = no_observations,
      j = no_variables
    )

   # Separate the matrices w.r.t. to repetitions and store them in a list

   step <- length(no_observations)*length(no_variables)*length(no_categories)*length(density)*length(interaction_scale)
   index <- rep(1:step, each = repetitions)

   gammas_mat <- t(sapply(gammas_1, '[', 1:max(sapply(gammas_1, length))))
   thetas_mat <- t(sapply(thetas_1, '[', 1:max(sapply(thetas_1, length))))

   # bind with the index value
   gammas_mat <- cbind(gammas_mat, index)
   thetas_mat <- cbind(thetas_mat, index)

   # Separate the matrices w.r.t. to repetitions and store them in a list
   gammas_comb <- lapply(by(gammas_mat,gammas_mat[,max(no_variables)*(max(no_variables)-1)/2+1],identity),as.matrix)
   thetas_comb <- lapply(by(thetas_mat,thetas_mat[,max(no_variables)*(max(no_variables)-1)/2+1],identity),as.matrix)

  # The if statement for the averageing should be somewhere here ---------------

   # Take the average
   mean_gammas <- lapply(gammas_comb, colMeans)
   mean_thetas <- lapply(thetas_comb, colMeans)

   # Put the elements in a matrix
   gammas <- matrix(0, nrow = length(mean_gammas), ncol = ncol(gammas_mat))
   thetas <- matrix(0, nrow = length(mean_thetas), ncol = ncol(thetas_mat))

   for(i in 1:length(mean_gammas)){
     gammas[i, ] <- mean_gammas[[i]]
     thetas[i, ] <- mean_thetas[[i]]
   }

   # take out the index variable
   gammas <- gammas[, -ncol(gammas)]
   thetas <- thetas[, -ncol(thetas)]


   gammas <- cbind(param_grid, gammas)
   thetas <- cbind(param_grid, thetas)

    # Sort by setting ----------------------------------------------------------------

    # Split w.r.t. no_observations
    gammas_N <- split(gammas, gammas$i)
    thetas_N <- split(thetas, thetas$i)

    # Split w.r.t. no_variables
    gammas_p <- list()
    thetas_p <- list()
    for(i in 1:length(no_observations)){
      gammas_p[[i]] <- split(gammas_N[[i]],gammas_N[[i]]$j)
      thetas_p[[i]] <- split(thetas_N[[i]],thetas_N[[i]]$j)
    }


    # Split w.r.t. number of cateogires
    gammas_spar <- vector(mode = "list", length = length(no_observations))
    thetas_spar <- vector(mode = "list", length = length(no_observations))
    for (i in seq_along(gammas_spar)) {
      gammas_spar[[i]] <- vector(mode = "list", length = length(no_variables))
      thetas_spar[[i]] <- vector(mode = "list", length = length(no_variables))
      for (j in seq_along(gammas_spar[[i]])) {
        gammas_spar[[i]][[j]] <- vector(mode = "list", length = length(no_categories))
        thetas_spar[[i]][[j]] <- vector(mode = "list", length = length(no_categories))
        for (k in seq_along(gammas_spar[[i]][[j]])) {
          gammas_spar[[i]][[j]][[k]] <- matrix(0)
          thetas_spar[[i]][[j]][[k]] <- matrix(0)
        }
      }
    }


    for(i in 1:length(thetas_N)){
      for(j in 1:length(thetas_p[[i]])){
        gammas_spar[[i]][[j]] <- split(gammas_p[[i]][[j]], gammas_p[[i]][[j]]$k)
      }
    }

    for(i in 1:length(thetas_N)){
      for(j in 1:length(thetas_p[[i]])){
        thetas_spar[[i]][[j]] <- split(thetas_p[[i]][[j]], thetas_p[[i]][[j]]$k)
      }
    }


    # Split w.r.t. density
    # storage space
    gammas_cat <- vector(mode = "list", length = length(no_observations))
    thetas_cat <- vector(mode = "list", length = length(no_observations))
    for (i in seq_along(  gammas_cat)) {
      gammas_cat[[i]] <- vector(mode = "list", length = length(no_variables))
      thetas_cat[[i]] <- vector(mode = "list", length = length(no_variables))
      for (j in seq_along(  gammas_cat[[i]])) {
        gammas_cat[[i]][[j]] <- vector(mode = "list", length = length(no_categories))
        thetas_cat[[i]][[j]] <- vector(mode = "list", length = length(no_categories))
        for (k in seq_along(gammas_cat[[i]][[j]])) {
          gammas_cat [[i]][[j]][[k]] <- vector(mode = "list", length = length(density))
          thetas_cat [[i]][[j]][[k]] <- vector(mode = "list", length = length(density))
          for (l in seq_along(gammas_cat[[i]][[j]][[k]])) {
            gammas_cat[[i]][[j]][[k]][[l]] <- matrix(0)
            thetas_cat[[i]][[j]][[k]][[l]] <- matrix(0)

          }
        }
      }
    }



    for(i in 1:length(gammas_N)){
      for(j in 1:length(gammas_p[[i]])){
        for (k in seq_along(gammas_spar[[i]][[j]])) {
          gammas_cat[[i]][[j]][[k]] <- split(gammas_spar[[i]][[j]][[k]], gammas_spar[[i]][[j]][[k]]$l)
          thetas_cat[[i]][[j]][[k]] <- split(thetas_spar[[i]][[j]][[k]], thetas_spar[[i]][[j]][[k]]$l)
        }
      }
    }

    # rename
    gammas <- gammas_cat
    thetas <- thetas_cat

    # take out the parameter grids and the NAs
    for(i in 1:length(no_observations)){
      for(j in 1:length(no_variables)){
        for(k in 1:length(no_categories)){
          for(l in 1:length(density)){
            gammas[[i]][[j]][[k]][[l]] <- gammas[[i]][[j]][[k]][[l]][, c(1:(((no_variables[j]*(no_variables[j]-1))/2) + 5))]
            thetas[[i]][[j]][[k]][[l]] <- thetas[[i]][[j]][[k]][[l]][, c(1:(((no_variables[j]*(no_variables[j]-1))/2) + 5))]
          }
        }
      }
    }




    # Reshape
  suppressWarnings({
    for(i in 1:length(no_observations)){
      for(j in 1:length(no_variables)){
        for(k in 1:length(no_categories)){
          for(l in 1:length(density)){
            gammas[[i]][[j]][[k]][[l]] <- cbind(gammas[[i]][[j]][[k]][[l]][, 1:5], stack(as.data.frame(gammas[[i]][[j]][[k]][[l]][, 6:ncol(gammas[[i]][[j]][[k]][[l]])])))
            thetas[[i]][[j]][[k]][[l]] <- cbind(thetas[[i]][[j]][[k]][[l]][, 1:5], stack(as.data.frame(thetas[[i]][[j]][[k]][[l]][, 6:ncol(thetas[[i]][[j]][[k]][[l]])])))
          }
        }
      }
    }
  })

    # Save everything into one data frame

  combined_df_list_gammas <- list()
  combined_df_list_thetas <- list()

  # Iterate over the nested list
  for(i in 1:length(no_observations)){
    for(j in 1:length(no_variables)){
      for(k in 1:length(no_categories)){
        for(l in 1:length(density)){
          # Access each data frame within the nested list
          df_gammas <- gammas[[i]][[j]][[k]][[l]]
          df_thetas <- thetas[[i]][[j]][[k]][[l]]
          # Append the data frames to the respective lists
          combined_df_list_gammas <- append(combined_df_list_gammas, list(df_gammas))
          combined_df_list_thetas <- append(combined_df_list_thetas, list(df_thetas))
        }
      }
    }
  }

    # Combine all data frames into a single data frame
    full <- do.call(rbind, combined_df_list_gammas)
    full_thetas <- do.call(rbind, combined_df_list_thetas)
    # take out the index value
    full <- full[, -ncol(full)]
    # take out only the thetas
    full <- cbind(full, full_thetas$values)
    # add names
    names(full) <- c("interaction_scale", "no_categories", "density",
                     "no_observations", "no_variables",
                     "posterior_inclusion_probability",
                     "mean_posterior_edge_weight")

    } else{ # For Gaussian networks

      param_grid <- expand.grid(
        k = density,
        i = no_observations,
        j = no_variables
      )


      # Separate the matrices w.r.t. to repetitions and store them in a list

      step <- length(no_observations)*length(no_variables)*length(density)
      index <- rep(1:step, each = repetitions)

      gammas_mat <- t(sapply(gammas_1, '[', 1:max(sapply(gammas_1, length))))
      thetas_mat <- t(sapply(thetas_1, '[', 1:max(sapply(thetas_1, length))))

      # bind with the index value
      gammas_mat <- cbind(gammas_mat, index)
      thetas_mat <- cbind(thetas_mat, index)

      # Separate the matrices w.r.t. to repetitions and store them in a list
      gammas_comb <- lapply(by(gammas_mat,gammas_mat[,max(no_variables)*(max(no_variables)-1)/2+1],identity),as.matrix)
      thetas_comb <- lapply(by(thetas_mat,thetas_mat[,max(no_variables)*(max(no_variables)-1)/2+1],identity),as.matrix)

      # The if statement for the averageing should be somewhere here ---------------

      # Take the average
      mean_gammas <- lapply(gammas_comb, colMeans)
      mean_thetas <- lapply(thetas_comb, colMeans)

      # Put the elements in a matrix
      gammas <- matrix(0, nrow = length(mean_gammas), ncol = ncol(gammas_mat))
      thetas <- matrix(0, nrow = length(mean_thetas), ncol = ncol(thetas_mat))

      for(i in 1:length(mean_gammas)){
        gammas[i, ] <- mean_gammas[[i]]
        thetas[i, ] <- mean_thetas[[i]]
      }

      # take out the index variable
      gammas <- gammas[, -ncol(gammas)]
      thetas <- thetas[, -ncol(thetas)]


      gammas <- cbind(param_grid, gammas)
      thetas <- cbind(param_grid, thetas)

      # Sort by setting ----------------------------------------------------------------

      # Split w.r.t. no_observations
      gammas_N <- split(gammas, gammas$i)
      thetas_N <- split(thetas, thetas$i)

      # Split w.r.t. no_variables
      gammas_p <- list()
      thetas_p <- list()
      for(i in 1:length(no_observations)){
        gammas_p[[i]] <- split(gammas_N[[i]],gammas_N[[i]]$j)
        thetas_p[[i]] <- split(thetas_N[[i]],thetas_N[[i]]$j)
      }


      # Split w.r.t. number the density
      gammas_spar <- vector(mode = "list", length = length(no_observations))
      thetas_spar <- vector(mode = "list", length = length(no_observations))
      for (i in seq_along(gammas_spar)) {
        gammas_spar[[i]] <- vector(mode = "list", length = length(no_variables))
        thetas_spar[[i]] <- vector(mode = "list", length = length(no_variables))
        for (j in seq_along(gammas_spar[[i]])) {
          gammas_spar[[i]][[j]] <- vector(mode = "list", length = length(density))
          thetas_spar[[i]][[j]] <- vector(mode = "list", length = length(density))
          for (k in seq_along(gammas_spar[[i]][[j]])) {
            gammas_spar[[i]][[j]][[k]] <- matrix(0)
            thetas_spar[[i]][[j]][[k]] <- matrix(0)
          }
        }
      }


      for(i in 1:length(thetas_N)){
        for(j in 1:length(thetas_p[[i]])){
          gammas_spar[[i]][[j]] <- split(gammas_p[[i]][[j]], gammas_p[[i]][[j]]$k)
        }
      }

      for(i in 1:length(thetas_N)){
        for(j in 1:length(thetas_p[[i]])){
          thetas_spar[[i]][[j]] <- split(thetas_p[[i]][[j]], thetas_p[[i]][[j]]$k)
        }
      }


      # rename
      gammas <- gammas_spar
      thetas <- thetas_spar

      # take out the parameter grids and the NAs
      for(i in 1:length(no_observations)){
        for(j in 1:length(no_variables)){
            for(k in 1:length(density)){
              gammas[[i]][[j]][[k]] <- gammas[[i]][[j]][[k]][, c(1:(((no_variables[j]*(no_variables[j]-1))/2) + 3))]
              thetas[[i]][[j]][[k]] <- thetas[[i]][[j]][[k]][, c(1:(((no_variables[j]*(no_variables[j]-1))/2) + 3))]
            }
          }
        }




      # Reshape
      suppressWarnings({
        for(i in 1:length(no_observations)){
          for(j in 1:length(no_variables)){
            for(k in 1:length(density)){
                gammas[[i]][[j]][[k]] <- cbind(gammas[[i]][[j]][[k]][, 1:3], stack(as.data.frame(gammas[[i]][[j]][[k]][, 4:ncol(gammas[[i]][[j]][[k]])])))
                thetas[[i]][[j]][[k]] <- cbind(thetas[[i]][[j]][[k]][, 1:3], stack(as.data.frame(thetas[[i]][[j]][[k]][, 4:ncol(thetas[[i]][[j]][[k]])])))
              }
            }
          }
      })

      # Save everything into one data frame

      combined_df_list_gammas <- list()
      combined_df_list_thetas <- list()

      # Iterate over the nested list
      for(i in 1:length(no_observations)){
        for(j in 1:length(no_variables)){
          for(k in 1:length(density)){
              # Access each data frame within the nested list
              df_gammas <- gammas[[i]][[j]][[k]]
              df_thetas <- thetas[[i]][[j]][[k]]
              # Append the data frames to the respective lists
              combined_df_list_gammas <- append(combined_df_list_gammas, list(df_gammas))
              combined_df_list_thetas <- append(combined_df_list_thetas, list(df_thetas))
          }
        }
      }

      # Combine all data frames into a single data frame
      full <- do.call(rbind, combined_df_list_gammas)
      full_thetas <- do.call(rbind, combined_df_list_thetas)
      # take out the index value
      full <- full[, -ncol(full)]
      # take out only the thetas
      full <- cbind(full, full_thetas$values)
      # add names
      names(full) <- c("density",
                       "no_observations", "no_variables",
                       "posterior_inclusion_probability",
                       "mean_posterior_partial_correlation")


      }# END OF else STATEMENT

  return(full)

} # END OF FUNCTION



