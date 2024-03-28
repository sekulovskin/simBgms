#' @importFrom stats rbinom rnorm runif
#'
# Generate data from an IRT model ==============================================

generate_IRT_data <- function(no_variables,
                              no_observations,
                              no_categories) {
  result_list <- list()

  alpha <- t(apply(matrix(runif(no_variables * no_categories, .3, 1), no_variables),
                   1, cumsum))
  alpha <- -(alpha - mean(alpha))
  alpha <- alpha + rnorm(no_variables)
  alpha <- alpha * -1

  beta <- mean(alpha)

  steps <- cbind(beta, t(apply(alpha, 1, function(row) row - beta)))

  theta <- rnorm(no_observations)

  y <- matrix(NA, nrow = length(theta), ncol = no_variables)

  for (k in 1:length(theta)) {
    for (i in 1:no_variables) {
      measure <- 0
      q <- numeric(no_categories)
      q[1] <- 1

      for (j in 2:no_categories) {
        measure <- measure + theta[k] - steps[i, 1] - steps[i, j]
        q[j] <- q[(j-1)] + exp(measure)
      }

      U <- runif(1, 0, 1)
      U <- U * q[no_categories]

      for (j in 1:no_categories) {
        if (U <= q[j]) {
          y[k, i] <- (j-1)
          break
        }
      }
    }
  }

  return(y)
}


# sparsity matrix ==============================================================

sparsity <- function(density,
                     no_variables){
  sparsity <- list()

  for (p_val in no_variables) {
    p_sublist <- list()

    for (spar_val in density) {
      sparsity_matrix <- matrix(0, nrow = p_val, ncol = p_val)

      for (i in 1:(p_val - 1)) {
        for (j in (i + 1):p_val) {
          sparsity_matrix[i, j] <- rbinom(1, 1, spar_val)
          sparsity_matrix[j, i] <- sparsity_matrix[i, j]
        }
      }

      p_sublist[[paste0("density", spar_val)]] <- sparsity_matrix
    }

    sparsity[[paste0("no_variables", p_val)]] <- p_sublist
  }

  sparsity
}


# Generate GGM data ============================================================

generate_ggm_data <- function(repetitions,
                              no_observations,
                              no_variables,
                              density,
                              adj_matrix) {

  # Storage

  data <- vector(mode = "list", length = repetitions )
  for(h in 1:repetitions ) {
    data[[h]] <- vector(mode = "list", length = length(no_observations))
    for(i in 1:length(no_observations)) {
      data[[h]][[i]] <- vector(mode = "list", length = length(no_variables))
      for(j in 1:length(no_variables)) {
        data[[h]][[i]][[j]] <- vector(mode = "list", length = length(density))
        for(k in 1:length(density)) {
          data[[h]][[i]][[j]][[k]] <- list()
        }
      }
    }
  }

  K <- vector(mode = "list", length = length(no_variables))
  for(j in 1:length(no_variables)) {
    K[[j]] <- vector(mode = "list", length = length(density))
    for(k in 1:length(density)) {
      K[[j]][[k]] <- list()
    }
  }



  # generate precision matrices according to the required nework density

  for(j in 1:length(no_variables)){
    for(k in 1:length(density)){
      K[[j]][[k]] <- BDgraph::rgwish(n = 1, adj = adj_matrix[[j]][[k]])
    }
  }

  # simulate data
  for(h in 1:repetitions ){
    for(i in 1:length(no_observations)){
      for(j in 1:length(no_variables)){
        for(k in 1:length(density)){
          data[[h]][[i]][[j]][[k]] <- BDgraph::bdgraph.sim(p = no_variables[j],
                                                           graph = "fixed",
                                                           n = no_observations[i],
                                                           K = K[[j]][[k]])
        }
      }
    }
  }


  return(list(data, K))
}



#Generate ordinal data =========================================================

generate_omrf_data <- function(repetitions,
                               variable_type,
                               reference_category,
                               no_observations,
                               no_variables,
                               no_categories,
                               iter = 1e3,
                               density,
                               MPLEs,
                               bc_alpha,
                               bc_beta){

  # Storage

  data <- vector(mode = "list", length = repetitions )
  for(h in 1:repetitions ) {
    data[[h]] <- vector(mode = "list", length = length(no_observations))
    for(i in 1:length(no_observations)) {
      data[[h]][[i]] <- vector(mode = "list", length = length(no_variables))
      for(j in 1:length(no_variables)) {
        data[[h]][[i]][[j]] <- vector(mode = "list", length = length(density))
        for(k in 1:length(density)) {
          data[[h]][[i]][[j]][[k]] <- vector(mode = "list", length = length(no_categories))
          for(l in 1:length(no_categories)) {
            data[[h]][[i]][[j]][[k]][[l]] <- list()
          }
        }
      }
    }
  }


  if(variable_type == "ordinal"){
    # simulate data
    for(h in 1:repetitions ){
      for(i in 1:length(no_observations)){
        for(j in 1:length(no_variables)){
          for(k in 1:length(density)){
            for(l in 1:length(no_categories)){
              data[[h]][[i]][[j]][[k]][[l]] <-  bgms::mrfSampler(no_states = no_observations[i],
                                                                 no_variables = no_variables[j],
                                                                 no_categories = rep(no_categories[l] - 1, no_variables[j]),
                                                                 interactions = MPLEs[[l]][[j]][[k]][["interactions"]],
                                                                 thresholds = as.matrix(MPLEs[[l]][[j]][[k]][["thresholds"]]),
                                                                 variable_type = variable_type,
                                                                 iter = iter)
            }

          }

        }

      }

    }
  } else{

    # simulate data
    for(h in 1:repetitions ){
      for(i in 1:length(no_observations)){
        for(j in 1:length(no_variables)){
          for(k in 1:length(density)){
            for(l in 1:length(no_categories)){
              MPLEs[[l]][[j]][[k]][["thresholds"]] <- cbind(rep(bc_alpha, no_variables[j]), rep(bc_beta, no_variables[j]))
              data[[h]][[i]][[j]][[k]][[l]] <-  bgms::mrfSampler(no_states = no_observations[i],
                                                                 no_variables = no_variables[j],
                                                                 no_categories = rep(no_categories[l] - 1, no_variables[j]),
                                                                 interactions = MPLEs[[l]][[j]][[k]][["interactions"]],
                                                                 thresholds = as.matrix(MPLEs[[l]][[j]][[k]][["thresholds"]]),
                                                                 variable_type = variable_type,
                                                                 reference_category = reference_category,
                                                                 iter = iter)
            }

          }

        }

      }

    }

  }


  return(list(data, MPLEs))
}
