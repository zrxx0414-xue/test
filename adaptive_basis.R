# Teacher source: main.R; exact origins and interface changes are documented in README.md.
# fit_adaptive_basis uses quantile_est through eigenfun_discrete; no quantiles are estimated.

fit_adaptive_basis <- function(data, kernel_bw = "adaptive", dim_basis = "adaptive"){
  m <- length(data) - 1
  sizes <- numeric(m+1)
  combined_data <- NULL
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)

  if (m < 1) stop("Error: number of samples should be at least 2!")
  if ((dim_basis != "adaptive") & (dim_basis > m)) stop("Error: dimension of the basis function should be smaller than the number of samples!")

  unique_data <- list()
  for (i in 1:(m+1)){
    unique_pt <- unique(data[[i]])
    num_unique <- length(unique_pt)

    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2)
    for (j in 1:num_unique){
      unique_data[[i]][j,1] <- unique_pt[j]
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j])
    }
  }
  if (identical(kernel_bw, "silverman")) {
    kernel_bw <- vapply(data, function(x) {
      0.9 * min(sd(x), IQR(x) / 1.34) * length(x)^(-1/5)
    }, numeric(1))
  }
  if (identical(kernel_bw, "adaptive")) kernel_bw <- adaptive_bw(data = data, ref_model = "normal")
  if (length(kernel_bw) == 1L) kernel_bw <- rep(kernel_bw, m+1)
  if (!is.numeric(kernel_bw) || length(kernel_bw) != m+1 || any(!is.finite(kernel_bw) | kernel_bw <= 0))
    stop("kernel_bw must be adaptive, a positive number, or one positive number per group.")
  if (dim_basis == "adaptive") dim_basis <- adaptive_dim(data = data, kernel_bw = kernel_bw)

  M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1))

  ff_tmp_seq <- list()
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))
  for (i in 1:length(combined_data)){
    x_M <- combined_data[i]
    log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2]))
    for (j in 1:(m+1)){
      ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M
    }
  }

  f_sum_M <- numeric(length(combined_data))
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j]
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]
  }
  f_avg_M <- f_sum_M/(m+1)
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M

  for (i_1 in 1:(m+1)){
    for (i_2 in i_1:(m+1)){
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]
      M_mat[i_1, i_2] <- mean(ff_tmp_mat)
      M_mat[i_2, i_1] <- M_mat[i_1, i_2]
    }
  }
  eigen_sys <- eigen(M_mat)

  eigenfun_unique_data <- matrix(NA, nrow = length(unique(combined_data)), ncol = dim_basis)
  for (ii in 1:length(unique(combined_data))){
    x <- unique(combined_data)[ii]
    x_loc <- which(combined_data == x)[1]
    log_DR <- numeric(m+1)
    for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc]

    fpca_basis_fn <- numeric(dim_basis)
    for (jj in 1:dim_basis) fpca_basis_fn[jj] <- (eigen_sys$values[jj]^(-1))*sum((eigen_sys$vectors[(1:(m+1)),jj])*log_DR)
    eigenfun_unique_data[ii,] <- fpca_basis_fn
  }

  eigenfun_discrete <- function(x){
    loc_unique_data <- which(x == unique(combined_data))
    return(eigenfun_unique_data[loc_unique_data,])
  }

  return(list(d = dim_basis, fun = eigenfun_discrete))
}

adaptive_bw <- function(data, ref_model){
  m <- length(data) - 1
  sizes <- numeric(m+1)
  combined_data <- NULL
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)

  bw_RoT <- numeric(m+1)
  for (i in 1:(m+1)) bw_RoT[i] <- 0.9*(sizes[i]^(-1/5))*min(sd(data[[i]]), IQR(data[[i]])/1.34)

  unique_data <- list()
  for (i in 1:(m+1)){
    unique_pt <- unique(data[[i]])
    num_unique <- length(unique_pt)

    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2)
    for (j in 1:num_unique){
      unique_data[[i]][j,1] <- unique_pt[j]
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j])
    }
  }

  M_mat_ref <- matrix(NA, nrow = (m+1), ncol = (m+1))

  if (ref_model == "normal"){
    ff_tmp_seq <- list()
    for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))
    for (i in 1:length(combined_data)){
      x_M <- combined_data[i]
      log_base_M <- log(dnorm(x = x_M, mean = median(data[[1]]), sd = min(sd(data[[1]]), IQR(data[[1]])/1.34)))
      for (j in 1:(m+1)){
        ff_tmp_seq[[j]][i] <- log(dnorm(x = x_M, mean = median(data[[j]]), sd = min(sd(data[[j]]), IQR(data[[j]])/1.34))) - log_base_M
      }
    }
  }

  if (ref_model == "gamma"){
    ff_tmp_seq <- list()
    for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))
    for (i in 1:length(combined_data)){
      x_M <- combined_data[i]
      log_base_M <- log(dgamma(x = x_M, shape = mean(data[[1]])^2/var(data[[1]]), scale = var(data[[1]])/mean(data[[1]])))
      for (j in 1:(m+1)){
        ff_tmp_seq[[j]][i] <- log(dgamma(x = x_M, shape = mean(data[[j]])^2/var(data[[j]]), scale = var(data[[j]])/mean(data[[j]]))) - log_base_M
      }
    }
  }

  f_sum_M <- numeric(length(combined_data))
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j]
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]
  }
  f_avg_M <- f_sum_M/(m+1)
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M

  for (i_1 in 1:(m+1)){
    for (i_2 in i_1:(m+1)){
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]
      M_mat_ref[i_1, i_2] <- mean(ff_tmp_mat)
      M_mat_ref[i_2, i_1] <- M_mat_ref[i_1, i_2]
    }
  }
  eigen_sys_ref <- eigen(M_mat_ref)

  eigenfun_unique_data_ref <- matrix(NA, nrow = length(unique(combined_data)), ncol = 2)
  for (ii in 1:length(unique(combined_data))){
    x <- unique(combined_data)[ii]
    x_loc <- which(combined_data == x)[1]
    log_DR <- numeric(m+1)
    for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc]

    fpca_basis_fn <- numeric(2)
    for (jj in 1:2) fpca_basis_fn[jj] <- (eigen_sys_ref$values[jj]^(-1))*sum((eigen_sys_ref$vectors[(1:(m+1)),jj])*log_DR)
    eigenfun_unique_data_ref[ii,] <- fpca_basis_fn
  }

  obj_fun <- function(bw_factor){
    if ((min(bw_factor) < 1) | (max(bw_factor) > 5)){
      return(Inf)
    } else{

      kernel_bw <- bw_factor*bw_RoT
      M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1))

      ff_tmp_seq <- list()
      for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))
      for (i in 1:length(combined_data)){
        x_M <- combined_data[i]
        log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2]))
        for (j in 1:(m+1)){
          ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M
        }
      }

      f_sum_M <- numeric(length(combined_data))
      mean_log_DR <- numeric(m+1)
      for (j in 1:(m+1)){
        mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
        ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j]
        f_sum_M <- f_sum_M + ff_tmp_seq[[j]]
      }
      f_avg_M <- f_sum_M/(m+1)
      for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M

      for (i_1 in 1:(m+1)){
        for (i_2 in i_1:(m+1)){
          ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]
          M_mat[i_1, i_2] <- mean(ff_tmp_mat)
          M_mat[i_2, i_1] <- M_mat[i_1, i_2]
        }
      }
      eigen_sys <- eigen(M_mat)

      eigenfun_unique_data <- matrix(NA, nrow = length(unique(combined_data)), ncol = 2)
      for (ii in 1:length(unique(combined_data))){
        x <- unique(combined_data)[ii]
        x_loc <- which(combined_data == x)[1]
        log_DR <- numeric(m+1)
        for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc]

        fpca_basis_fn <- numeric(2)
        for (jj in 1:2) fpca_basis_fn[jj] <- (eigen_sys$values[jj]^(-1))*sum((eigen_sys$vectors[(1:(m+1)),jj])*log_DR)
        eigenfun_unique_data[ii,] <- fpca_basis_fn
      }

      return(sum(colMeans((eigenfun_unique_data - eigenfun_unique_data_ref)^2)))
    }
  }

  opt <- optimize(f = obj_fun, interval = c(1, 5))
  adaptive_bw_factor <- opt$minimum

  return(adaptive_bw_factor*bw_RoT)
}

library(drmdel)

adaptive_dim <- function(data, kernel_bw){
  m <- length(data) - 1
  sizes <- numeric(m+1)
  combined_data <- NULL
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)

  unique_data <- list()
  for (i in 1:(m+1)){
    unique_pt <- unique(data[[i]])
    num_unique <- length(unique_pt)

    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2)
    for (j in 1:num_unique){
      unique_data[[i]][j,1] <- unique_pt[j]
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j])
    }
  }

  M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1))

  ff_tmp_seq <- list()
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))
  for (i in 1:length(combined_data)){
    x_M <- combined_data[i]
    log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2]))
    for (j in 1:(m+1)){
      ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M
    }
  }

  f_sum_M <- numeric(length(combined_data))
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j]
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]
  }
  f_avg_M <- f_sum_M/(m+1)
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M

  for (i_1 in 1:(m+1)){
    for (i_2 in i_1:(m+1)){
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]
      M_mat[i_1, i_2] <- mean(ff_tmp_mat)
      M_mat[i_2, i_1] <- M_mat[i_1, i_2]
    }
  }
  eigen_sys <- eigen(M_mat)

  eigenvalues <- eigen_sys$values
  eigenvalues_ratio <- numeric(m)
  for (j in 1:m) eigenvalues_ratio[j] <- sum(eigenvalues[1:j])/sum(eigenvalues)

  eigenfun_unique_data <- matrix(NA, nrow = length(unique(combined_data)), ncol = m)
  for (ii in 1:length(unique(combined_data))){
    x <- unique(combined_data)[ii]
    x_loc <- which(combined_data == x)[1]
    log_DR <- numeric(m+1)
    for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc]

    fpca_basis_fn <- numeric(m)
    for (jj in 1:m) fpca_basis_fn[jj] <- (eigen_sys$values[jj]^(-1))*sum((eigen_sys$vectors[(1:(m+1)),jj])*log_DR)
    eigenfun_unique_data[ii,] <- fpca_basis_fn
  }

  eigenfun_discrete <- function(x, nFPC){
    loc_unique_data <- which(x == unique(combined_data))
    return(eigenfun_unique_data[loc_unique_data,(1:nFPC)])
  }

  BIC_candidates <- numeric(m)
  for (j in 1:m){
    drmfit_FPCA <- drmdel(x=combined_data, n_samples=sizes, basis_func=function(x) return(eigenfun_discrete(x, nFPC=j)))
    likelihood_FPCA <- -drmfit_FPCA$negldl - n*log(n)
    BIC_candidates[j] <- -2*likelihood_FPCA + j*m*log(n)
  }

  nFPC_95 <- min(m, which(eigenvalues_ratio >= 0.95)[1])

  nFPC_BIC <- which.min(BIC_candidates)

  nFPC_adaptive <- max(nFPC_BIC, nFPC_95)

  return(nFPC_adaptive)
}
