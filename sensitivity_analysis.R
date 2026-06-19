sensitivity_analysis_ <- function(n, d, U, dataset, betas){
  SA_effects <- vector("numeric", length=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  main_effect_r <- U[rowSums(U > 0) == 1,]
  # SA
  main_effect_betas <- betas[rowSums(U > 0) == 1]
  two_factor_effect_r <- U[rowSums(U > 0) == 2,]
  
  largest_two_factor = rep(0,d)
  largest_two_factor_beta <- -1
  two_factor_effect_betas <- betas[rowSums(U > 0) == 2]
  total_var <- sum(betas[-1]^2)
  
  for (j in 1:d){
    SA_effects[j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
  }
  
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      # TODO FIX
      curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                          two_factor_effect_r[,k] > 0))]^2) / total_var
      if (curr > largest_two_factor_beta){
        largest_two_factor_beta = curr
        largest_two_factor = rep(0,d)
        largest_two_factor[j] = 1
        largest_two_factor[k] = 1
      }
    }
  }
  
  # Total percent variation explained by the main effects.
  SA_effects[d+1] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 1]]^2) / total_var
  # Total percent variation explained by the two-factor effects.
  SA_effects[d+2] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 2]]^2) / total_var
  # Total percent variation explained by the three-factor effects.
  SA_effects[d+3] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 3]]^2) / total_var
  # variation of the main effects by variable.
  for (i in 1:d){
    print(paste(predictor_names[i], format(mean(SA_effects[i]), digits=3, scientific=TRUE)))
  }
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(largest_two_factor_beta, digits=3, scientific=TRUE)))
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[d+1]), digits=3, scientific=TRUE)))
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[d+2]), digits=3, scientific=TRUE)))
  print(paste("Total percent variation explained by the three-factor effects:",
              format(mean(SA_effects[d+3]), digits=3, scientific=TRUE)))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[d+1]) +  mean(SA_effects[d+2]) > 0.99)
    eff_dim = 2
  else
    eff_dim = 3
  
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
}

sensitivity_analysis_P <- function(n, d, U, dataset, betas, ISE_esti, runtime){
  SA_effects <- vector("numeric", length=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  main_effect_r <- U[rowSums(U > 0) == 1,]
  # SA
  main_effect_betas <- betas[rowSums(U > 0) == 1]
  two_factor_effect_r <- U[rowSums(U > 0) == 2,]
  
  largest_two_factor = rep(0,d)
  largest_two_factor_beta <- -1
  two_factor_effect_betas <- betas[rowSums(U > 0) == 2]
  total_var <- sum(betas[-1]^2)
  
  for (j in 1:d){
    SA_effects[j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
  }
  
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      # TODO FIX
      curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                          two_factor_effect_r[,k] > 0))]^2) / total_var
      if (curr > largest_two_factor_beta){
        largest_two_factor_beta = curr
        largest_two_factor = rep(0,d)
        largest_two_factor[j] = 1
        largest_two_factor[k] = 1
      }
    }
  }
  
  # Total percent variation explained by the main effects.
  SA_effects[d+1] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 1]]^2) / total_var
  # Total percent variation explained by the two-factor effects.
  SA_effects[d+2] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 2]]^2) / total_var
  
  # sorted main effects.
  sort_ind =  sort(SA_effects[1:d],decreasing=T, index.return=T)$ix[1:13]
  
  # variation of the main effects by variable.
  for (i in sort_ind){
    print(paste(predictor_names[i], format(SA_effects[i], digits=3, scientific=TRUE)))
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(largest_two_factor_beta, digits=3, scientific=TRUE)))
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[d+1]), digits=3, scientific=TRUE)))
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[d+2]), digits=3, scientific=TRUE)))
  
  print(paste("Interpretable Set Size:" , nrow(U)))
  
  # Effective Dimension
  #eff_dim <- 0
  
  #if (mean(SA_effects[d+1]) > 0.99)
  #  eff_dim = 1
  #else if (mean(SA_effects[d+1]) +  mean(SA_effects[d+2]) > 0.99)
  #  eff_dim = 2
  #else
  #  eff_dim = 3
  
  #print(paste("The Effective Dimension in the superposition sense of 99%:",
  #            eff_dim))
  
  print(paste("ISE estimate:", ISE_esti$ISE_mean))
  print(paste("ISE sd:", ISE_esti$ISE_sd))
  print(paste("Runtime:", runtime))
  print(paste("Efficiency ISE:", 1 / (ISE_esti$ISE_mean * runtime)))
}

sensitivity_analysis_sobol <- function(n, d, U, dataset, betas){
  SA_effects <- vector("numeric", length=d+4)
  predictor_names = colnames(dataset)[-(d+1)]
  main_effect_r <- U[rowSums(U > 0) == 1,]
  # SA
  main_effect_betas <- betas[rowSums(U > 0) == 1]
  two_factor_effect_r <- U[rowSums(U > 0) == 2,]
  
  largest_two_factor = rep(0,d)
  largest_two_factor_beta <- -1
  two_factor_effect_betas <- betas[rowSums(U > 0) == 2]
  total_var <- sum(betas[-1]^2)
  
  for (j in 1:d){
    SA_effects[j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
  }
  
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      # TODO FIX
      curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                          two_factor_effect_r[,k] > 0))]^2) / total_var
      if (curr > largest_two_factor_beta){
        largest_two_factor_beta = curr
        largest_two_factor = rep(0,d)
        largest_two_factor[j] = 1
        largest_two_factor[k] = 1
      }
    }
  }
  
  # Total percent variation explained by the main effects.
  SA_effects[d+1] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 1]]^2) / total_var
  # Total percent variation explained by the two-factor effects.
  SA_effects[d+2] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 2]]^2) / total_var
  # Total percent variation explained by the three-factor effects.
  SA_effects[d+3] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 3]]^2) / total_var
  # Total percent variation explained by the four-factor effects.
  SA_effects[d+4] = sum(betas[(1:nrow(U))[rowSums(U > 0) == 4]]^2) / total_var
  # variation of the main effects by variable.
  for (i in 1:d){
    print(paste(predictor_names[i], mean(SA_effects[i])))
  }
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, largest_two_factor_beta))
  print(paste("Total percent variation explained by the main effects:",
              mean(SA_effects[d+1])))
  print(paste("Total percent variation explained by the two-factor effects:",
              mean(SA_effects[d+2])))
  print(paste("Total percent variation explained by the three-factor effects:",
              mean(SA_effects[d+3])))
  print(paste("Total percent variation explained by the four-factor effects:",
              mean(SA_effects[d+4])))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[d+1]) +  mean(SA_effects[d+2]) > 0.99)
    eff_dim = 2
  else if (mean(SA_effects[d+1]) +  mean(SA_effects[d+2]) + mean(SA_effects[d+3]) > 0.99)
    eff_dim = 3
  else
    eff_dim = 4
  
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
}

glob_sen_ind_lower <- function(n, d, f)
{
  X <- matrix(nrow=n, ncol=d, runif(n*d))
  Y <- sobol(2^17,d,randomize = "digital.shift")
  X_prime <- matrix(nrow=n, ncol=d, runif(n*d))
  
  tau_lower <- vector("numeric", length=d)
  
  mu_sqr <- mean(f(Y)^2) - var(f(Y))
  
  for (i in 1:d){
    X_copy = X_prime
    X_copy[,i] = X[,i]
    
    tau_lower[i] = mean(f(X)*f(X_copy)) - mu_sqr
  }
  tau_lower
}

sensitivity_analysis <- function(n, d, U, Q, dataset, betas, LOF_fit){
  print(format(mean(LOF_fit[n,]), digits=3, scientific=TRUE))
  print(format(sd(LOF_fit[n,]) / sqrt(Q), digits=3, scientific=TRUE))
  SA_effects <- matrix(0, nrow=Q, ncol=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  main_effect_r <- U[rowSums(U > 0) == 1,]
  two_factor_effect_r <- U[rowSums(U > 0) == 2,]
  
  largest_two_factor = rep(0,d)
  largest_two_factor_beta <- rep(-1,Q)
  
  # SA over each Q
  for (i in 1:Q){
    main_effect_betas <- betas[,i][rowSums(U > 0) == 1]
    two_factor_effect_betas <- betas[,i][rowSums(U > 0) == 2]
    total_var <- sum(betas[,i][-1]^2)
    
    for (j in 1:d){
      SA_effects[i,j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
    }
    
    for (j in 1:(d-1)){
      for (k in (j+1):d){
        # TODO FIX
        curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                             two_factor_effect_r[,k] > 0))]^2) / total_var
        if (curr > largest_two_factor_beta[i]){
          largest_two_factor_beta[i] = curr
          largest_two_factor = rep(0,d)
          largest_two_factor[j] = 1
          largest_two_factor[k] = 1
        }
      }
    }
    
    # Total percent variation explained by the main effects.
    SA_effects[i,d+1] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 1]]^2) / total_var
    # Total percent variation explained by the two-factor effects.
    SA_effects[i,d+2] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 2]]^2) / total_var
    # Total percent variation explained by the three-factor effects.
    SA_effects[i,d+3] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 3]]^2) / total_var
  }
  # variation of the main effects by variable.
  for (i in 1:d){
    print(paste(predictor_names[i], format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "(",
                format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(mean(largest_two_factor_beta), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  print(paste("Total percent variation explained by the three-factor effects:",
              format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[,d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) > 0.99)
    eff_dim = 2
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) + mean(SA_effects[,d+3]) > 0.99)
    eff_dim = 3
  else
    eff_dim = 4
  
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
  
  print(paste("The approximate convergence rate:", lm(log(apply(LOF_fit, MARGIN=1, mean)) ~ log(1:n))$coefficients[[2]]))
}

sensitivity_analysis_out <- function(n, d, U, Q, dataset, betas, LOF_fit, LOF_mean, LOF_sd, var_rate, runtime){
  output = list()
  output[[1]] = paste0(format(mean(LOF_mean), digits=3, scientific=TRUE), "}_{(",format(LOF_sd, digits=3, scientific=TRUE), "})")
  output[[2]] = paste0(format(mean(LOF_fit[n,]), digits=3, scientific=TRUE), "}_{(",format(sd(LOF_fit[n,]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(format(mean(LOF_mean), digits=3, scientific=TRUE))
  print(format(LOF_sd, digits=3, scientific=TRUE))
  print(format(mean(LOF_fit[n,]), digits=3, scientific=TRUE))
  print(format(sd(LOF_fit[n,]) / sqrt(Q), digits=3, scientific=TRUE))
  SA_effects <- matrix(0, nrow=Q, ncol=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  main_effect_r <- U[rowSums(U > 0) == 1,]
  two_factor_effect_r <- U[rowSums(U > 0) == 2,]
  
  largest_two_factor = rep(0,d)
  largest_two_factor_beta <- rep(-1,Q)
  
  # SA over each Q
  for (i in 1:Q){
    main_effect_betas <- betas[,i][rowSums(U > 0) == 1]
    two_factor_effect_betas <- betas[,i][rowSums(U > 0) == 2]
    total_var <- sum(betas[,i][-1]^2)
    
    for (j in 1:d){
      SA_effects[i,j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
    }
    
    for (j in 1:(d-1)){
      for (k in (j+1):d){
        # TODO FIX
        curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                            two_factor_effect_r[,k] > 0))]^2) / total_var
        if (curr > largest_two_factor_beta[i]){
          largest_two_factor_beta[i] = curr
          largest_two_factor = rep(0,d)
          largest_two_factor[j] = 1
          largest_two_factor[k] = 1
        }
      }
    }
    
    # Total percent variation explained by the main effects.
    SA_effects[i,d+1] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 1]]^2) / total_var
    # Total percent variation explained by the two-factor effects.
    SA_effects[i,d+2] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 2]]^2) / total_var
    # Total percent variation explained by the three-factor effects.
    SA_effects[i,d+3] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 3]]^2) / total_var
  }
  # variation of the main effects by variable.
  for (i in 1:d){
    output[[i+2]] = paste0(format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), "})")
    print(paste(predictor_names[i], format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "(",
                format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  output[[d + 3]] = two_factor_name
  output[[d + 4]] = paste0(format(mean(largest_two_factor_beta), digits=3, scientific=TRUE), "}_{(",format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(mean(largest_two_factor_beta), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  output[[d + 5]] = paste0(format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 6]] = paste0(format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 7]] =  paste0(format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the three-factor effects:",
              format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[,d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) > 0.99)
    eff_dim = 2
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) + mean(SA_effects[,d+3]) > 0.99)
    eff_dim = 3
  else
    eff_dim = 4
  
  output[[d + 8]] = eff_dim
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
  
  output[[d + 9]] = format(lm(log(apply(LOF_fit, MARGIN=1, mean)) ~ log(1:n))$coefficients[[2]], digits=3, scientific=FALSE)
  output[[d + 10]] = format(var_rate, digits=3, scientific=FALSE)
  print(paste("The approximate convergence rate:", lm(log(apply(LOF_fit, MARGIN=1, mean)) ~ log(1:n))$coefficients[[2]]))
  print(paste("The approximate var convergence rate:", var_rate))
  
  output[[d + 11]] = format(runtime, digits=5, scientific=FALSE)
  output[[d + 12]] = format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)
  output[[d + 13]] = format(1 / (mean(LOF_fit[n,]) * runtime), digits=3, scientific=TRUE)
  print(paste("The runtime in seconds:", runtime))
  print(paste("The Effciency of LOF:",format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)))
  print(paste("The Effciency of Beta Error:",format(1 / (mean(LOF_fit[n,]) * runtime), digits=3, scientific=TRUE)))
  output_vec = unlist(output)
  edit_out = paste0("$", gsub("e", "\\\\times 10^{",gsub("-0", "-", output_vec)), "}$")
  edit_out[d + 3] = output_vec[d + 3]
  edit_out[c(d + 8,d + 9,  d + 10, d + 11)] = gsub("}", "", edit_out[c(d + 8,d + 9,  d + 10, d + 11)])
  edit_out
}

sensitivity_analysis_QR_ABS <- function(n, d, U, Q, dataset, betas, LOF_mean, LOF_sd, runtime){
  output = list()
  output[[1]] = paste0(format(mean(LOF_mean), digits=3, scientific=TRUE), "}_{(",format(LOF_sd, digits=3, scientific=TRUE), "})")
  print(format(mean(LOF_mean), digits=3, scientific=TRUE))
  print(format(LOF_sd, digits=3, scientific=TRUE))
  
  SA_effects <- matrix(0, nrow=Q, ncol=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  main_effect_r <- U[rowSums(U > 0) == 1,]
  two_factor_effect_r <- U[rowSums(U > 0) == 2,]
  
  largest_two_factor = rep(0,d)
  largest_two_factor_beta <- rep(-1,Q)
  
  # SA over each Q
  for (i in 1:Q){
    main_effect_betas <- betas[,i][rowSums(U > 0) == 1]
    two_factor_effect_betas <- betas[,i][rowSums(U > 0) == 2]
    total_var <- sum(betas[,i][-1]^2)
    
    for (j in 1:d){
      SA_effects[i,j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
    }
    
    for (j in 1:(d-1)){
      for (k in (j+1):d){
        # TODO FIX
        curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                            two_factor_effect_r[,k] > 0))]^2) / total_var
        if (curr > largest_two_factor_beta[i]){
          largest_two_factor_beta[i] = curr
          largest_two_factor = rep(0,d)
          largest_two_factor[j] = 1
          largest_two_factor[k] = 1
        }
      }
    }
    
    # Total percent variation explained by the main effects.
    SA_effects[i,d+1] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 1]]^2) / total_var
    # Total percent variation explained by the two-factor effects.
    SA_effects[i,d+2] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 2]]^2) / total_var
    # Total percent variation explained by the three-factor effects.
    SA_effects[i,d+3] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 3]]^2) / total_var
  }
  # variation of the main effects by variable.
  for (i in 1:d){
    output[[i+1]] = paste0(format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), "})")
    print(paste(predictor_names[i], format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "(",
                format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  output[[d + 2]] = two_factor_name
  output[[d + 3]] = paste0(format(mean(largest_two_factor_beta), digits=3, scientific=TRUE), "}_{(",format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(mean(largest_two_factor_beta), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  output[[d + 4]] = paste0(format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 5]] = paste0(format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 6]] =  paste0(format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the three-factor effects:",
              format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[,d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) > 0.99)
    eff_dim = 2
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) + mean(SA_effects[,d+3]) > 0.99)
    eff_dim = 3
  else
    eff_dim = 4
  
  output[[d + 7]] = eff_dim
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
  
  output[[d + 8]] = format(runtime, digits=5, scientific=FALSE)
  output[[d + 9]] = format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)
  print(paste("The runtime in seconds:", runtime))
  print(paste("The Effciency of LOF:",format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)))
  output_vec = unlist(output)
  edit_out = paste0("$", gsub("e", "\\\\times 10^{",gsub("-0", "-", output_vec)), "}$")
  edit_out[d + 2] = output_vec[d + 2]
  edit_out[c(d + 7,d + 8,  d + 9)] = gsub("}", "", edit_out[c(d + 7,d + 8,  d + 9)])
  edit_out
}

sensitivity_analysis_ABS_out <- function(n, d, U, Q, dataset, betas, LOF_fit, LOF_mean, LOF_sd, runtime){
  output = list()
  output[[1]] = paste0(format(mean(LOF_mean), digits=3, scientific=TRUE), "}_{(",format(LOF_sd, digits=3, scientific=TRUE), "})")
  output[[2]] = paste0(format(mean(LOF_fit[n,]), digits=3, scientific=TRUE), "}_{(",format(sd(LOF_fit[n,]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(format(mean(LOF_mean), digits=3, scientific=TRUE))
  print(format(LOF_sd, digits=3, scientific=TRUE))
  print(format(mean(LOF_fit[n,]), digits=3, scientific=TRUE))
  print(format(sd(LOF_fit[n,]) / sqrt(Q), digits=3, scientific=TRUE))
  SA_effects <- matrix(0, nrow=Q, ncol=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  
  inter_avg_size = 0
  
  # SA over each Q
  
  largest_two_factor_beta <- rep(0,Q)
  
  for (i in 1:Q){
    inter_avg_size = inter_avg_size + nrow(U[[i]])
    main_effect_r <- U[[i]][rowSums(U[[i]] > 0) == 1,]
    two_factor_effect_r <- U[[i]][rowSums(U[[i]] > 0) == 2,]
    
    largest_two_factor = rep(0,d)
    
    main_effect_betas <- betas[[i]][rowSums(U[[i]] > 0) == 1]
    two_factor_effect_betas <- betas[[i]][rowSums(U[[i]] > 0) == 2]
    total_var <- sum(betas[[i]][-1]^2)
    
    for (j in 1:d){
      SA_effects[i,j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
    }
    
    for (j in 1:(d-1)){
      for (k in (j+1):d){
        # TODO FIX
        curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                            two_factor_effect_r[,k] > 0))]^2) / total_var
        if (curr >= largest_two_factor_beta[i]){
          largest_two_factor_beta[i] = curr
          largest_two_factor = rep(0,d)
          largest_two_factor[j] = 1
          largest_two_factor[k] = 1
        }
      }
    }
    
    # Total percent variation explained by the main effects.
    SA_effects[i,d+1] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 1]]^2) / total_var
    # Total percent variation explained by the two-factor effects.
    SA_effects[i,d+2] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 2]]^2) / total_var
    # Total percent variation explained by the three-factor effects.
    SA_effects[i,d+3] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 3]]^2) / total_var
  }
  # variation of the main effects by variable.
  for (i in 1:d){
    output[[i+2]] = paste0(format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), "})")
    print(paste(predictor_names[i], format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "(",
                format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  output[[d + 3]] = two_factor_name
  output[[d + 4]] = paste0(format(mean(largest_two_factor_beta), digits=3, scientific=TRUE), "}_{(",format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(mean(largest_two_factor_beta), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  output[[d + 5]] = paste0(format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 6]] = paste0(format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 7]] =  paste0(format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the three-factor effects:",
              format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[,d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) > 0.99)
    eff_dim = 2
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) + mean(SA_effects[,d+3]) > 0.99)
    eff_dim = 3
  else
    eff_dim = 4
  
  output[[d + 8]] = eff_dim
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
  
  output[[d + 9]] = format(lm(log(apply(LOF_fit, MARGIN=1, mean)) ~ log(1:n))$coefficients[[2]], digits=3, scientific=FALSE)
  output[[d + 10]] = round(inter_avg_size / Q)
  print(paste("The approximate convergence rate:", lm(log(apply(LOF_fit, MARGIN=1, mean)) ~ log(1:n))$coefficients[[2]]))
  print(paste("The average interpretable set size:", output[[d + 10]]))
  
  output[[d + 11]] = format(runtime, digits=5, scientific=FALSE)
  output[[d + 12]] = format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)
  output[[d + 13]] = format(1 / (mean(LOF_fit[n,]) * runtime), digits=3, scientific=TRUE)
  print(paste("The runtime in seconds:", runtime))
  print(paste("The Effciency of LOF:",format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)))
  print(paste("The Effciency of Beta Error:",format(1 / (mean(LOF_fit[n,]) * runtime), digits=3, scientific=TRUE)))
  output_vec = unlist(output)
  edit_out = paste0("$", gsub("e", "\\\\times 10^{",gsub("-0", "-", output_vec)), "}$")
  edit_out[d + 3] = output_vec[d + 3]
  edit_out[d + 10] = output_vec[d + 10]
  edit_out[c(d + 8,d + 9,  d + 11)] = gsub("}", "", edit_out[c(d + 8,d + 9,  d + 11)])
  edit_out
}

sensitivity_analysis_QR_ABS <- function(n, d, U, Q, dataset, betas, LOF_mean, LOF_sd, runtime){
  output = list()
  output[[1]] = paste0(format(mean(LOF_mean), digits=3, scientific=TRUE), "}_{(",format(LOF_sd, digits=3, scientific=TRUE), "})")
  print(format(mean(LOF_mean), digits=3, scientific=TRUE))
  print(format(LOF_sd, digits=3, scientific=TRUE))
  
  SA_effects <- matrix(0, nrow=Q, ncol=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  main_effect_r <- U[rowSums(U > 0) == 1,]
  two_factor_effect_r <- U[rowSums(U > 0) == 2,]
  
  largest_two_factor = rep(0,d)
  largest_two_factor_beta <- rep(-1,Q)
  
  # SA over each Q
  for (i in 1:Q){
    main_effect_betas <- betas[,i][rowSums(U > 0) == 1]
    two_factor_effect_betas <- betas[,i][rowSums(U > 0) == 2]
    total_var <- sum(betas[,i][-1]^2)
    
    for (j in 1:d){
      SA_effects[i,j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
    }
    
    for (j in 1:(d-1)){
      for (k in (j+1):d){
        # TODO FIX
        curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                            two_factor_effect_r[,k] > 0))]^2) / total_var
        if (curr > largest_two_factor_beta[i]){
          largest_two_factor_beta[i] = curr
          largest_two_factor = rep(0,d)
          largest_two_factor[j] = 1
          largest_two_factor[k] = 1
        }
      }
    }
    
    # Total percent variation explained by the main effects.
    SA_effects[i,d+1] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 1]]^2) / total_var
    # Total percent variation explained by the two-factor effects.
    SA_effects[i,d+2] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 2]]^2) / total_var
    # Total percent variation explained by the three-factor effects.
    SA_effects[i,d+3] = sum(betas[,i][(1:nrow(U))[rowSums(U > 0) == 3]]^2) / total_var
  }
  # variation of the main effects by variable.
  for (i in 1:d){
    output[[i+1]] = paste0(format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), "})")
    print(paste(predictor_names[i], format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "(",
                format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  output[[d + 2]] = two_factor_name
  output[[d + 3]] = paste0(format(mean(largest_two_factor_beta), digits=3, scientific=TRUE), "}_{(",format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(mean(largest_two_factor_beta), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  output[[d + 4]] = paste0(format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 5]] = paste0(format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 6]] =  paste0(format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the three-factor effects:",
              format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[,d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) > 0.99)
    eff_dim = 2
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) + mean(SA_effects[,d+3]) > 0.99)
    eff_dim = 3
  else
    eff_dim = 4
  
  output[[d + 7]] = eff_dim
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
  output[[d+8]] = nrow(U)
  
  output[[d + 9]] = format(runtime, digits=5, scientific=FALSE)
  output[[d + 10]] = format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)
  print(paste("The runtime in seconds:", runtime))
  print(paste("The Effciency of LOF:",format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)))
  output_vec = unlist(output)
  edit_out = paste0("$", gsub("e", "\\\\times 10^{",gsub("-0", "-", output_vec)), "}$")
  edit_out[d + 2] = output_vec[d + 2]
  edit_out[c(d + 7,d + 8,  d + 9)] = gsub("}", "", edit_out[c(d + 7,d + 8,  d + 9)])
  edit_out
}

sensitivity_analysis_ABS_out <- function(n, d, U, Q, dataset, betas, LOF_fit, LOF_mean, LOF_sd, runtime){
  output = list()
  output[[1]] = paste0(format(mean(LOF_mean), digits=3, scientific=TRUE), "}_{(",format(LOF_sd, digits=3, scientific=TRUE), "})")
  output[[2]] = paste0(format(mean(LOF_fit[n,]), digits=3, scientific=TRUE), "}_{(",format(sd(LOF_fit[n,]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(format(mean(LOF_mean), digits=3, scientific=TRUE))
  print(format(LOF_sd, digits=3, scientific=TRUE))
  print(format(mean(LOF_fit[n,]), digits=3, scientific=TRUE))
  print(format(sd(LOF_fit[n,]) / sqrt(Q), digits=3, scientific=TRUE))
  SA_effects <- matrix(0, nrow=Q, ncol=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  
  inter_avg_size = 0
  
  # SA over each Q
  
  largest_two_factor_beta <- rep(0,Q)
  
  for (i in 1:Q){
    inter_avg_size = inter_avg_size + nrow(U[[i]])
    main_effect_r <- U[[i]][rowSums(U[[i]] > 0) == 1,]
    two_factor_effect_r <- U[[i]][rowSums(U[[i]] > 0) == 2,]
    
    largest_two_factor = rep(0,d)
    
    main_effect_betas <- betas[[i]][rowSums(U[[i]] > 0) == 1]
    two_factor_effect_betas <- betas[[i]][rowSums(U[[i]] > 0) == 2]
    total_var <- sum(betas[[i]][-1]^2)
    
    for (j in 1:d){
      SA_effects[i,j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
    }
    
    for (j in 1:(d-1)){
      for (k in (j+1):d){
        # TODO FIX
        curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                            two_factor_effect_r[,k] > 0))]^2) / total_var
        if (curr >= largest_two_factor_beta[i]){
          largest_two_factor_beta[i] = curr
          largest_two_factor = rep(0,d)
          largest_two_factor[j] = 1
          largest_two_factor[k] = 1
        }
      }
    }
    
    # Total percent variation explained by the main effects.
    SA_effects[i,d+1] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 1]]^2) / total_var
    # Total percent variation explained by the two-factor effects.
    SA_effects[i,d+2] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 2]]^2) / total_var
    # Total percent variation explained by the three-factor effects.
    SA_effects[i,d+3] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 3]]^2) / total_var
  }
  # variation of the main effects by variable.
  for (i in 1:d){
    output[[i+2]] = paste0(format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), "})")
    print(paste(predictor_names[i], format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "(",
                format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  output[[d + 3]] = two_factor_name
  output[[d + 4]] = paste0(format(mean(largest_two_factor_beta), digits=3, scientific=TRUE), "}_{(",format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(mean(largest_two_factor_beta), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  output[[d + 5]] = paste0(format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 6]] = paste0(format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 7]] =  paste0(format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the three-factor effects:",
              format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[,d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) > 0.99)
    eff_dim = 2
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) + mean(SA_effects[,d+3]) > 0.99)
    eff_dim = 3
  else
    eff_dim = 4
  
  output[[d + 8]] = eff_dim
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
  
  output[[d + 9]] = format(lm(log(apply(LOF_fit, MARGIN=1, mean)) ~ log(1:n))$coefficients[[2]], digits=3, scientific=FALSE)
  output[[d + 10]] = round(inter_avg_size / Q)
  print(paste("The approximate convergence rate:", lm(log(apply(LOF_fit, MARGIN=1, mean)) ~ log(1:n))$coefficients[[2]]))
  print(paste("The average interpretable set size:", output[[d + 10]]))
  
  output[[d + 11]] = format(runtime, digits=5, scientific=FALSE)
  output[[d + 12]] = format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)
  output[[d + 13]] = format(1 / (mean(LOF_fit[n,]) * runtime), digits=3, scientific=TRUE)
  print(paste("The runtime in seconds:", runtime))
  print(paste("The Effciency of LOF:",format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)))
  print(paste("The Effciency of Beta Error:",format(1 / (mean(LOF_fit[n,]) * runtime), digits=3, scientific=TRUE)))
  output_vec = unlist(output)
  edit_out = paste0("$", gsub("e", "\\\\times 10^{",gsub("-0", "-", output_vec)), "}$")
  edit_out[d + 3] = output_vec[d + 3]
  edit_out[d + 10] = output_vec[d + 10]
  edit_out[c(d + 8,d + 9,  d + 11)] = gsub("}", "", edit_out[c(d + 8,d + 9,  d + 11)])
  edit_out
}

sensitivity_analysis_MARQR_ABS <- function(n, d, U, Q, dataset, betas, LOF_mean, LOF_sd, runtime){
  output = list()
  output[[1]] = paste0(format(mean(LOF_mean), digits=3, scientific=TRUE), "}_{(",format(LOF_sd, digits=3, scientific=TRUE), "})")
  print(format(mean(LOF_mean), digits=3, scientific=TRUE))
  print(format(LOF_sd, digits=3, scientific=TRUE))
  SA_effects <- matrix(0, nrow=Q, ncol=d+3)
  predictor_names = colnames(dataset)[-(d+1)]
  
  inter_avg_size = 0
  
  # SA over each Q
  
  largest_two_factor_beta <- rep(0,Q)
  
  for (i in 1:Q){
    inter_avg_size = inter_avg_size + nrow(U[[i]])
    main_effect_r <- U[[i]][rowSums(U[[i]] > 0) == 1,]
    two_factor_effect_r <- U[[i]][rowSums(U[[i]] > 0) == 2,]
    
    largest_two_factor = rep(0,d)
    
    main_effect_betas <- betas[[i]][rowSums(U[[i]] > 0) == 1]
    two_factor_effect_betas <- betas[[i]][rowSums(U[[i]] > 0) == 2]
    total_var <- sum(betas[[i]][-1]^2)
    
    for (j in 1:d){
      SA_effects[i,j] = sum(main_effect_betas[main_effect_r[,j] > 0]^2) / total_var
    }
    
    for (j in 1:(d-1)){
      for (k in (j+1):d){
        # TODO FIX
        curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                            two_factor_effect_r[,k] > 0))]^2) / total_var
        if (curr >= largest_two_factor_beta[i]){
          largest_two_factor_beta[i] = curr
          largest_two_factor = rep(0,d)
          largest_two_factor[j] = 1
          largest_two_factor[k] = 1
        }
      }
    }
    
    # Total percent variation explained by the main effects.
    SA_effects[i,d+1] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 1]]^2) / total_var
    # Total percent variation explained by the two-factor effects.
    SA_effects[i,d+2] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 2]]^2) / total_var
    # Total percent variation explained by the three-factor effects.
    SA_effects[i,d+3] = sum(betas[[i]][(1:nrow(U[[i]]))[rowSums(U[[i]] > 0) == 3]]^2) / total_var
  }
  # variation of the main effects by variable.
  for (i in 1:d){
    output[[i+1]] = paste0(format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), "})")
    print(paste(predictor_names[i], format(mean(SA_effects[,i]), digits=3, scientific=TRUE), "(",
                format(sd(SA_effects[,i]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  output[[d + 2]] = two_factor_name
  output[[d + 3]] = paste0(format(mean(largest_two_factor_beta), digits=3, scientific=TRUE), "}_{(",format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(mean(largest_two_factor_beta), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  output[[d + 4]] = paste0(format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the main effects:",
              format(mean(SA_effects[,d+1]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+1]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 5]] = paste0(format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the two-factor effects:",
              format(mean(SA_effects[,d+2]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+2]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  output[[d + 6]] =  paste0(format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE), "}_{(",format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), "})")
  print(paste("Total percent variation explained by the three-factor effects:",
              format(mean(SA_effects[,d+3]), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(SA_effects[,d+3]) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  # Effective Dimension
  eff_dim <- 0
  
  if (mean(SA_effects[,d+1]) > 0.99)
    eff_dim = 1
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) > 0.99)
    eff_dim = 2
  else if (mean(SA_effects[,d+1]) +  mean(SA_effects[,d+2]) + mean(SA_effects[,d+3]) > 0.99)
    eff_dim = 3
  else
    eff_dim = 4
  
  output[[d + 7]] = eff_dim
  print(paste("The Effective Dimension in the superposition sense of 99%:",
              eff_dim))
  
  output[[d + 8]] = round(inter_avg_size / Q)
  print(paste("The average interpretable set size:", output[[d + 8]]))
  
  output[[d + 9]] = format(runtime, digits=5, scientific=FALSE)
  output[[d + 10]] = format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)
  print(paste("The runtime in seconds:", runtime))
  print(paste("The Effciency of LOF:",format(1 / (mean(LOF_mean) * runtime), digits=3, scientific=TRUE)))
  output_vec = unlist(output)
  edit_out = paste0("$", gsub("e", "\\\\times 10^{",gsub("-0", "-", output_vec)), "}$")
  edit_out[d + 2] = output_vec[d + 2]
  edit_out[d + 8] = output_vec[d + 8]
  edit_out[c(d + 7,  d + 9)] = gsub("}", "", edit_out[c(d + 7,d + 9)])
  edit_out
}

lte_ABS_out <- function(d, U, Q, dataset, betas){
  predictor_names = colnames(dataset)[-(d+1)]
  
  largest_two_factor_beta <- rep(0,Q)
  largest_two_factor = rep(0,d)
  
  # SA over each Q
  for (i in 1:Q){
    two_factor_effect_r <- U[[i]][rowSums(U[[i]] > 0) == 2,]
    
    two_factor_effect_betas <- betas[[i]][rowSums(U[[i]] > 0) == 2]
    total_var <- sum(betas[[i]][-1]^2)
    
    for (j in 1:(d-1)){
      for (k in (j+1):d){
        # TODO FIX
        curr = sum(two_factor_effect_betas[Reduce("&", list(two_factor_effect_r[,j] > 0,  
                                                            two_factor_effect_r[,k] > 0))]^2) / total_var
        if (curr >= largest_two_factor_beta[i]){
          largest_two_factor_beta[i] = curr
          largest_two_factor = rep(0,d)
          largest_two_factor[j] = 1
          largest_two_factor[k] = 1
        }
      }
    }
  }
  
  two_factor_name = paste(predictor_names[largest_two_factor==1], collapse = ' ')
  
  print(paste("The percent variation explained by the largest 2 factor effect is:",
              two_factor_name, format(mean(largest_two_factor_beta), digits=3, scientific=TRUE)))
  print(paste("(", format(sd(largest_two_factor_beta) / sqrt(Q), digits=3, scientific=TRUE), ")"))
  
  print(largest_two_factor_beta)
}