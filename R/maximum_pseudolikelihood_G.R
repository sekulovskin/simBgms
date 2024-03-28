# Maximum Pseudolikelihood estimation for a Markov Random Field model for
# ordinal variables with the possibility to induce sparsity by providing
# an adjacency matrix

mple_G = function(x,
                  convergence_criterion = sqrt(.Machine$double.eps),
                  maximum_iterations = 1e3,
                  thresholds,
                  interactions,
                  G) {
  #Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix"))
    stop("The input x is supposed to be a matrix.")

  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  #Format the data input -------------------------------------------------------
  data = reformat_data(x = x)
  x = data$x
  no_categories = data$no_categories
  no_nodes = ncol(x)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_thresholds = sum(no_categories)
  no_parameters = no_thresholds + no_interactions

  #Check NR input --------------------------------------------------------------
  if(convergence_criterion <= 0)
    stop("Parameter ``convergence_criterion'' needs to be positive.")
  if(maximum_iterations <= 0 ||
     abs(maximum_iterations - round(maximum_iterations)) >
     sqrt(.Machine$double.eps))
    stop("Parameter ``maximum_iterations'' needs to be a positive integer.")

  # Starting values -----------------------------------------------------------
  if(!hasArg("thresholds")) {
    thresholds = matrix(0,
                        nrow = no_nodes,
                        ncol = max(no_categories))
  }
  if(!hasArg("interactions")) {
    interactions = matrix(0,
                          nrow = no_nodes,
                          ncol = no_nodes)
  }

  # Newton-Raphson ------------------------------------------------------------
  log_pl = log_pseudolikelihood(interactions * G,
                                thresholds,
                                observations = x,
                                no_categories)

  row = rep(1, no_nodes * (no_nodes - 1) / 2)
  cntr = 0
  for(s in 1:(no_nodes-1)) {
    for(t in (s+1):no_nodes) {
      cntr = cntr + 1
      if(G[s, t] == 0) {
        row[cntr] = 0
      }
    }
  }

  hessian = matrix(data = NA,
                   nrow = no_parameters - length(which(row == 0)),
                   ncol = no_parameters - length(which(row == 0)))

  gradient = matrix(data = NA,
                    nrow = 1,
                    ncol = no_parameters - length(which(row == 0)))

  for(iteration in 1:maximum_iterations) {
    old_log_pl = log_pl

    #Compute gradient vector (first order derivatives) ------------------------
    gradient[1:no_thresholds] =
      gradient_thresholds_pseudolikelihood(interactions = G * interactions,
                                           thresholds = thresholds,
                                           observations = x,
                                           no_categories)
    tmp =
      gradient_interactions_pseudolikelihood(interactions = G * interactions,
                                             thresholds = thresholds,
                                             observations = x,
                                             no_categories)
    tmp = tmp[-which(row == 0)]

    gradient[-c(1:no_thresholds)] = tmp
    # Compute Hessian matrix (second order partial derivatives) ---------------
    hessian[1:no_thresholds, 1:no_thresholds] =
      hessian_thresholds_pseudolikelihood(interactions = G * interactions,
                                          thresholds = thresholds,
                                          observations = x,
                                          no_categories)

    tmp =  hessian_interactions_pseudolikelihood(interactions = G * interactions,
                                                 thresholds = thresholds,
                                                 observations = x,
                                                 no_categories)
    tmp = tmp[-which(row == 0),]
    tmp = tmp[, -which(row == 0)]




    hessian[-(1:no_thresholds), -(1:no_thresholds)] = tmp


    tmp = hessian_crossparameters(interactions = G * interactions,
                                  thresholds = thresholds,
                                  observations = x,
                                  no_categories)
    tmp = tmp[-which(row == 0),]

    hessian[-(1:no_thresholds), 1:no_thresholds] = tmp

    hessian[1:no_thresholds, -(1:no_thresholds)] =
      t(hessian[-(1:no_thresholds), 1:no_thresholds])

    # Update parameter values (Newton-Raphson step) ---------------------------
    Delta = gradient %*% solve(hessian)
    if(any(is.nan(Delta)) || any(is.infinite(Delta)))
      stop("log_pseudolikelihood optimization failed. Please check the data. If the data checks out, please try different starting values.")

    cntr = 0
    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        cntr = cntr + 1
        thresholds[node, category] = thresholds[node, category] - Delta[cntr]
      }
    }
    for(node in 1:(no_nodes - 1)) {
      for(node_2 in (node + 1):no_nodes) {
        if(G[node, node_2] != 0) {
          cntr = cntr + 1
          interactions[node, node_2] = interactions[node, node_2] - Delta[cntr]
          interactions[node_2, node] = interactions[node, node_2]
        }
      }
    }

    # Check for convergence ---------------------------------------------------
    log_pl = log_pseudolikelihood(interactions * G,
                                  thresholds,
                                  observations = x,
                                  no_categories)

    if(abs(log_pl - old_log_pl) < convergence_criterion)
      break
  }

  if(abs(log_pl - old_log_pl) >= convergence_criterion &&
     iteration == maximum_iterations)
    warning(paste("The optimization procedure did not convergence in",
                  maximum_iterations, "iterations.",
                  sep = " "),
            call. = FALSE)

  colnames(interactions) = paste0("node ", 1:no_nodes)
  rownames(interactions) = paste0("node ", 1:no_nodes)
  colnames(thresholds) = paste0("category ", 1:max(no_categories))
  rownames(thresholds) = paste0("node ", 1:no_nodes)

  return(list(interactions = interactions, thresholds = thresholds))
}
