# runction to compute the partial correlation matrix manually
precision_to_partial_correlation <- function(precision_matrix) {
  # Ensure the input is a matrix
  if (!is.matrix(precision_matrix)) {
    stop("Input must be a matrix.")
  }

  # Get the dimensions of the precision matrix
  p <- nrow(precision_matrix)

  # Initialize the partial correlation matrix
  partial_corr_matrix <- matrix(0, nrow = p, ncol = p)

  # Compute the partial correlations
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) {
        # Diagonal elements are 1
        partial_corr_matrix[i, j] <- 1
      } else {
        # Off-diagonal elements use the formula
        partial_corr_matrix[i, j] <- -precision_matrix[i, j] / sqrt(precision_matrix[i, i] * precision_matrix[j, j])
      }
    }
  }

  return(partial_corr_matrix)
}

# function to combine the parameter grid with the estimates
combine_with_param_grid <- function(matrix_list, param_grid) {
  result_list <- list()

  for (i in 1:length(matrix_list)) {
    # Get the current matrix
    current_matrix <- matrix_list[[i]]

    # Get the number of rows in the current matrix
    n_i <- nrow(current_matrix)

    # Duplicate the corresponding row in param_grid n_i times
    duplicated_params <- param_grid[rep(i, n_i), ]

    # Combine the duplicated params with the matrix
    combined_matrix <- cbind(duplicated_params, current_matrix)

    # Add the combined matrix to the result list
    result_list[[i]] <- combined_matrix
  }

  return(result_list)
}


# Function to add lower triangular elements as a new column
add_lower_tri_elements <- function(combined_list, adj_mat) {
  for (i in 1:length(combined_list)) {
    # Extract no_variables and density from the current combined_list element
    no_vars <- combined_list[[i]]$no_variables[1]
    dens <- combined_list[[i]]$density[1]

    # Construct the key for the nested list
    no_vars_key <- paste0("no_variables", no_vars)
    dens_key <- paste0("density", dens)

    # Check if the adjacency matrix exists
    if (!is.null(adj_mat[[no_vars_key]][[dens_key]])) {
      # Extract the corresponding adjacency matrix
      adj_matrix <- adj_mat[[no_vars_key]][[dens_key]]

      # Extract the lower triangular elements of the adjacency matrix
      lower_tri_elements <- adj_matrix[lower.tri(adj_matrix)]

      # Get the number of rows in the current data frame
      n_rows <- nrow(combined_list[[i]])

      # Check if the number of lower triangular elements matches the number of rows
      if (length(lower_tri_elements) == n_rows) {
        # Add the lower triangular elements as a new column
        combined_list[[i]]$lower_tri <- lower_tri_elements
      } else if (length(lower_tri_elements) > n_rows) {
        # Truncate the lower triangular elements to match the number of rows
        combined_list[[i]]$lower_tri <- lower_tri_elements[1:n_rows]
      } else {
        # Pad the lower triangular elements with NA to match the number of rows
        combined_list[[i]]$lower_tri <- c(lower_tri_elements, rep(NA, n_rows - length(lower_tri_elements)))
      }
    } else {
      # If the adjacency matrix does not exist, add a column of NA values
      combined_list[[i]]$lower_tri <- rep(NA, nrow(combined_list[[i]]))
    }
  }

  return(combined_list)
}

# the main function

summarize <- function(estimates,
                      level = c("Gaussian", "Discrete"),
                      repetitions,
                      interaction_scale,
                      no_observations,
                      no_variables,
                      no_categories,
                      density,
                      adj_mat,
                      average = TRUE) {  # New argument to control averaging

  # Convert estimates to a list if not already
  est <- as.list(estimates)

  # Count the instances with error messages
  error_count <- sum(sapply(est, function(x) is.character(x) && !is.list(x)))
  if (error_count / length(est) != 0) {
    warning("There are models that did not converge,
            so in post-processing they are replaced
            with models from the same simulation settings
            that did converge. To check which models have failed,
            check the list of estimated models
            in the output and compare the index with
            the one in the parameter grid.")

    # Find indices of failed models (character elements) and succeeded models (list elements)
    indices_failed <- which(sapply(est, is.character))
    indices_succeeded <- which(sapply(est, is.list))

    # Iterate over the failed models
    for (failed_index in indices_failed) {
      # Calculate the range of indices for successful models based on the failed index
      range_start <- ((failed_index - 1) %/% repetitions) * repetitions + 1
      range_end <- range_start + repetitions - 1

      # Find the appropriate successful index within the range
      successful_indices <- which(indices_succeeded >= range_start & indices_succeeded <= range_end)

      # If there are no successful indices within the range, skip this failed index
      if (length(successful_indices) == 0) next

      # Choose a random successful index within the range
      selected_successful_index <- sample(successful_indices, 1)

      # Replace the failed model with the successful one
      est[[failed_index]] <- est[[indices_succeeded[selected_successful_index]]]
    }
  }

  # Create the parameter grid and extract elements

  if (level == "Gaussian") {
    param_grid <- expand.grid(
      k = density,
      i = no_observations,
      j = no_variables
    )
      gammas <- lapply(est, `[[`, "p_links")
      thetas <- lapply(est, `[[`, "K_hat")
      thetas <- lapply(thetas, precision_to_partial_correlation)

  } else {
    param_grid <- expand.grid(
      rep = 1:repetitions,
      s = interaction_scale,
      k = no_categories,
      l = density,
      i = no_observations,
      j = no_variables
    )
    # Extract inclusion indicators and parameters
    gammas <- lapply(est, `[[`, "indicator")
    thetas <- lapply(est, `[[`, "interactions")
  }

  # Extract the lower triangular elements
  gammas_1 <- lapply(gammas, function(x) x[lower.tri(x)])
  thetas_1 <- lapply(thetas, function(x) x[lower.tri(x)])

  # Combine into a list of data frames
  combined <- lapply(1:length(est), function(i) {
    data.frame(
      estimate = thetas_1[[i]],
      inclusion = gammas_1[[i]]
    )
  })

  # Combine with parameter grid
  combined_list <- combine_with_param_grid(combined, param_grid)

  # Add column names
  if (level == "Gaussian") {
    for (i in 1:length(combined_list)) {
      colnames(combined_list[[i]]) <- c("density", "no_observations", "no_variables",
                                        "estimate", "inclusion")
    }
  } else {
    for (i in 1:length(combined_list)) {
      colnames(combined_list[[i]]) <- c(
        "repetition", "interaction_scale", "no_categories", "density",
        "no_observations", "no_variables", "estimate", "inclusion"
      )
    }
  }
    # Add the true inclusions (lower triangular elements of adjacency matrices)
    combined_list <- add_lower_tri_elements(combined_list, adj_mat)

    # Unlist into a single data frame
    combined_df <- dplyr::bind_rows(combined_list)

    # Average the results if requested
    if (average) {
      # Group by the specified columns and calculate the mean for each group
      combined_df <- aggregate(
        cbind(estimate, inclusion, lower_tri) ~ interaction_scale + no_categories + density + no_observations + no_variables,
        data = combined_df,
        FUN = function(x) mean(x, na.rm = TRUE)
      )
    }
    return(combined_df)
  }


