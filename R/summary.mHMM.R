#' @export
#'
summary.mHMM <- function(object, ...){
  input   <- object$input
  dep_labels <- input$dep_labels
  n_subj  <- input$n_subj
  burn_in <- input$burn_in
  J       <- input$J
  m       <- input$m
  q_emiss <- input$q_emiss
  n_dep   <- input$n_dep
  xx <- input$xx
  gamma_pop <- matrix(round(apply(object$gamma_prob_bar[((burn_in + 1): J),], 2, median),3), byrow = TRUE, ncol = m, nrow = m)
  colnames(gamma_pop) <- paste("To state", 1:m)
  rownames(gamma_pop) <- paste("From state", 1:m)
  cat("State transition probability matrix","\n",  "(at the group level):", "\n", "\n")
  print(gamma_pop)
  cat("\n")

  if(!is.null(dim(object$gamma_cov_bar))){
    #mean_beta <- colMeans(object$gamma_cov_bar[-(1:burn_in), ])
    median_beta <- matrix(round(apply(object$gamma_cov_bar[-(1:burn_in), ], 2, median), 3), byrow = TRUE, ncol = m, nrow = 1)
    gamma_cov <- t(rbind(-median_beta, median_beta))
    colnames(gamma_cov) <- paste("To state", 1:m)
    rownames(gamma_cov) <- paste("From state", 1:m)
    cat("Covariate contribution to transition probabilities", "\n")
    print(gamma_cov)
  }
  cat("\n", "\n")

  cat("Emission distribution for each of the dependent variables","\n",  "(at the group level):", "\n", "\n")
  EM_pop <- vector("list", n_dep)
  names(EM_pop) <- dep_labels
  for(i in 1:n_dep){
    EM_pop[[i]] <- matrix(round(apply(object$emiss_prob_bar[[i]][((burn_in + 1): J),], 2, median),3), byrow = TRUE, ncol = q_emiss[i], nrow = m)
    colnames(EM_pop[[i]]) <- paste("Category", 1:q_emiss[i])
    rownames(EM_pop[[i]]) <- paste("State", 1:m)
  }
  print(EM_pop)
  cat("\n")

mean_theta <- vector("list", n_dep)
names(mean_theta) <- dep_labels
for(i in 1:n_dep){
  mean_theta[[i]] <- matrix(round(apply(object$emiss_cov_bar[[i]][((burn_in + 1): J),], 2, median), 3), byrow = TRUE, ncol = q_emiss[i]-1, nrow = m)
  colnames(mean_theta[[i]]) <- paste("Category", 2:q_emiss[i])
  rownames(mean_theta[[i]]) <- paste("State", 1:m)
}
if(!is.null(dim(object$emiss_cov_bar))){
  cat("Covariate contribution to emission distributions", "\n")
  print(mean_theta)
  cat("\n")
}
}
