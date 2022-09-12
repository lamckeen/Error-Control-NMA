# All functions for NMA analysis/simulation


# New Functions -----------------------------------------------------------

## Function to get design/error matrices and implements frequentist NMA for arm level data
## Returns matrices, estimates, and parameters for the contrast of interest in this simulation

nma_arm_analysis <- function(MTC_data){
  
  main_trts= setdiff(unique(c(MTC_data$Arm.1,MTC_data$Arm.2,
                              MTC_data$Arm.3)),"")
  baseline = names(which(table(MTC_data$Arm.2) == max(table(MTC_data$Arm.2)))) #NAC
  non_baseline_trts = setdiff(main_trts, baseline)
  
  data <- MTC_data
  Y <- list()
  S <- list()
  X <- list()
  
  for(i in 1:nrow(data)){
    
    n_arms <- data[i,"Number.of.arms"]
    
    if(n_arms == 2){
      
      p1 <- data[i,"Number.of.Event.in.arm.1"]/data[i,"Total.number.in.arm.1"]
      p2 <- data[i,"Number.of.Event.in.arm.2"]/data[i,"Total.number.in.arm.2"]
      n1 <- data[i,"Total.number.in.arm.1"]
      n2 <- data[i,"Total.number.in.arm.2"]
      
      lor = log((p1/(1-p1))/(p2/(1-p2)))
      var =  (1/(n1*p1*(1-p1)))+(1/(n2*p2*(1-p2)))
        
      y <- lor
      s <- var
      
      x <- rep(0, length(non_baseline_trts))
      baseline_arm <- data[i, "Arm.2"]
      if(baseline_arm == baseline){
        x[which(non_baseline_trts == data[i, "Arm.1"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.2"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.1"])
        x[trt1] <- -1
        x[trt2] <- 1
      }
      
      Y[[i]] <- y
      S[[i]] <- s
      X[[i]] <- x
    }
    
    if(n_arms == 3){
      
      p1 <- data[i,"Number.of.Event.in.arm.1"]/data[i,"Total.number.in.arm.1"]
      p2 <- data[i,"Number.of.Event.in.arm.2"]/data[i,"Total.number.in.arm.2"]
      p3 <- data[i,"Number.of.Event.in.arm.3"]/data[i,"Total.number.in.arm.3"]
      n1 <- data[i,"Total.number.in.arm.1"]
      n2 <- data[i,"Total.number.in.arm.2"]
      n3 <- data[i,"Total.number.in.arm.3"]
      
      lor1 = log((p1/(1-p1))/(p2/(1-p2)))
      var1 =  (1/(n1*p1*(1-p1)))+(1/(n2*p2*(1-p2)))
      
      lor3 = log((p3/(1-p3))/(p2/(1-p2)))
      var3 =  (1/(n3*p3*(1-p3)))+(1/(n2*p2*(1-p2)))
      
      y1 <- lor1
      s1 <- var1
      
      y2 <- lor3
      s2 <- var3
      
      v <- (1/(n2*p2*(1-p2)))
      
      x1 <- rep(0, length(non_baseline_trts))
      x2 <- rep(0, length(non_baseline_trts))
      
      baseline_arm <- data[i, "Arm.2"]
      if(baseline_arm == baseline){
        x1[which(non_baseline_trts == data[i, "Arm.1"])] <- 1
        x2[which(non_baseline_trts == data[i, "Arm.3"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.2"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.1"])
        trt3 <- which(non_baseline_trts == data[i, "Arm.3"])
        x1[trt1] <- -1
        x1[trt2] <- 1
        x2[trt1] <- -1
        x2[trt3] <- 1
      }
      
      Y[[i]] <- c(y1, y2)
      S[[i]] <- matrix(data = c(s1, v, v, s2), nrow = 2, byrow = TRUE)
      X[[i]] <- rbind(x1, x2)
    }
    
  }  
  
  X_final <- do.call(rbind, X)
  
  Y_final <- do.call(c, Y)
  
  S_final <- do.call(adiag, S)
  
  ## Analysis
  
  mu_hat <- solve(t(X_final)%*%solve(S_final)%*%X_final) %*% t(X_final)%*%solve(S_final)%*% Y_final
  # names(mu_hat) <- non_baseline_trts
  v_hat <- solve(t(X_final)%*%solve(S_final)%*%X_final)
  
  contrast_vec = c(rep(0,7),1,0,0,0,-1) #Current indirect of interest
  contrast_mu = contrast_vec %*% mu_hat
  contrast_var = t(contrast_vec) %*% v_hat %*% contrast_vec
  contrast_z = contrast_mu/sqrt(contrast_var)
  contrast_pval = 2*pnorm(-abs(contrast_z))
  
  output <- list(Y = Y_final, S = S_final, X = X_final,
                 main_trts = main_trts, baseline = baseline,
                 non_baseline_trts = non_baseline_trts,
                 mu_hat = mu_hat,
                 v_hat = v_hat,
                 c_mu = contrast_mu,
                 c_var = contrast_var,
                 c_z = contrast_z,
                 c_pval = contrast_pval,
                 df_nw = as.data.frame(data) #also need to return original data to be used later
                 )
  return(output)
}

nma_new_trial <- function(iter,sample_size,c_1,c_2){
  
  x= iter
  n = sample_size/2
  new_event = rbinom(2,n, risk_AB)
  new_event[new_event == n] = new_event-.5
  new_event[new_event == 0] = .5
  
  # New data frame
  df_new = x$df_nw %>% add_row(Study.number = 99,
                                Number.of.Event.in.arm.1 = new_event[1],
                                Number.of.Event.in.arm.2 = new_event[2],
                                Total.number.in.arm.1 = n,
                                Total.number.in.arm.2 = n,
                                Number.of.arms = 2,
                                Total = sample_size,
                                Arm.1 = "Ceftiofur hydrochloride",
                                Arm.2 = "Trimethoprim",
                                Arm.3 = "")
  
  p_new_1 = new_event[1]/n
  p_new_2 = new_event[2]/n
  lor_new = log((p_new_1/(1-p_new_1))/(p_new_2/(1-p_new_2)))
  s_new = (1/(n*p_new_1*(1-p_new_1)))+(1/(n*p_new_2*(1-p_new_2)))
  
  # Without network analysis (glm)
  z_without = lor_new/sqrt(s_new)
  p_val_without = 2*pnorm(-abs(z_without))
  p_val_without_ind = ifelse(p_val_without < .05,1,0)
  
  # New X
  x_new = c(rep(0,7),1,0,0,0,-1)
  x_all = rbind(x$X,x_new)
  
  # New Y
  y_all = c(x$Y,lor_new)
  
  # New S
  s_all = rbind(x$S,c(rep(0,nrow(x$S)))) %>%
    cbind(c(rep(0,nrow(.))))
  
  s_all[nrow(s_all),nrow(s_all)] <- s_new
  
  # New NMA estimation
  mu_hat <- solve(t(x_all)%*%solve(s_all)%*%x_all) %*% t(x_all)%*%solve(s_all)%*% y_all
  v_hat <- solve(t(x_all)%*%solve(s_all)%*%x_all)
  
  contrast_vec = c(rep(0,7),1,0,0,0,-1) #Current indirect of interest
  contrast_mu = contrast_vec %*% mu_hat
  contrast_var = t(contrast_vec) %*% v_hat %*% contrast_vec
  contrast_z = contrast_mu/sqrt(contrast_var)
  contrast_pval = 2*pnorm(-abs(contrast_z))
  p_val_ind = ifelse(contrast_pval < .05,1,0)
  
  # Sequential NMA
  p_val_seq_ind = ifelse(abs(x$c_z) > c_1 | abs(contrast_z) > c_2, 1,0)
  
  return(list(c_mu = contrast_mu,
           c_var = contrast_var,
           c_pval = contrast_pval,
           power = p_val_ind,
           power_glm = p_val_without_ind,
           power_seq = p_val_seq_ind,
           df_new = df_new,
           decision_p = x$c_pval))
}

# Fangshu's Functions -----------------------------------------------------

# it seems like this didn't use at this stage
# log(k1/k2)=lor
lor2prob <- function(p1, lor){
  p2 <- p1/(p1 + exp(lor)*(1-p1))
  return(p2)
}

# rearrange the data to the long type: one arm one study per row
wide2long <- function(MTCdata){
  N <- max(MTCdata$Number.of.arms)
  study <- rep(MTCdata$Study.number,N)
  t <- NULL
  n <- NULL
  r <- NULL
  for(i in c(2,1,3)){
    r <- c(r, eval(parse(text = paste0("MTCdata$Number.of.Event.in.arm.",i, sep = ""))))
    n <- c(n, eval(parse(text = paste0("MTCdata$Total.number.in.arm.",i, sep = ""))))
    t <- c(t, eval(parse(text = paste0("MTCdata$Arm.",i, sep = ""))))
  }
  res <- data.frame(id = study, t = t, r = r, n = n)
  res <- res %>% dplyr::filter(!is.na(n)) %>% arrange(id)
  res
}



# the type of the new_s is the same as what wide2long() generate above.

new_study <- function(p_baseline, p_trt2, sigma, pig_alloc = c(50,50)){
  r_vec <- rbinom(2, size = pig_alloc, prob = c(p_baseline, p_trt2))
  new_s <- data.frame(id = rep(100,2), t = c("Trimethoprim", "Ceftiofur hydrochloride"),
                      r = r_vec, n = pig_alloc)
  return(new_s)
}

# Simulation

### with previous network
bio_equal <- function(p_baseline, p_trt2, sigma, re_sigma = 0, pig_alloc, data_prev){
  #res_TE <- numeric(1)
  #res_seTE <- numeric(1)
  #power <- numeric(length(p_trt2_vector))
  new_s <- new_study(p_baseline, p_trt2, sigma = sigma, pig_alloc = pig_alloc)
  data_final <- rbind(data_prev, new_s)
  BRD_new <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = data_final, allstudies = T, sm = "OR")
  nma_res <- netmeta::netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new,sm="OR",
                              comb.fixed =F,comb.random = T, tol.multiarm.se = 0.01)
  #z <- true_lor/se
  #power <- pnorm(z-qnorm(0.975))+pnorm(-z-qnorm(0.975))
  lor_hat <- nma_res$TE.fixed[12,1]
  se <- nma_res$seTE.fixed[12,1]
  power <- ifelse(abs(lor_hat)/se > qnorm(0.975),1,0)
  
  return(matrix(c(power,lor_hat,se),ncol = 3))
}

### non-inferiority and superiority without previous network
bioeq_single_s <- function(p_baseline, p_trt2, pig_alloc){
  #power <- numeric(1)
  
  # random part in the simulation, rbinom
  BRD_s <- new_study(p_baseline, p_trt2, sigma = 0, pig_alloc = pig_alloc)
  BRD_new_s <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_s, allstudies = T, sm = "OR")
  # TE: Estimate of treatment effect (log odds ratio, mean difference)
  # seTE: S.E. of TE
  # sm: summary measure,
  # comb.fixed: whether a fixed effects (common effects) network meta-analysis should be conducted.
  # comb.random: whether a random effects network meta-analysis should be conducted.
  nma_s <- netmeta::netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new_s,sm="OR",
                            comb.fixed =F,comb.random = T, tol.multiarm.se = 0.01)
  #z <- true_lor/se
  #power <- pnorm(z-qnorm(0.975))+pnorm(-z-qnorm(0.975))
  # lor_TILD_2_CEFTP
  lor_hat <- nma_s$TE.fixed[2,1]
  se <- nma_s$seTE.fixed[2,1]
  power <- ifelse(abs(lor_hat)/se > qnorm(0.975),1,0)
  return(matrix(c(power,lor_hat,se),ncol = 3))
}

con_fn <- function(n, beta1, beta2,s){
  n1 <- n[1]
  n2 <- n[2]
  n1 + n2 - s
}

var_fn <- function(n, beta1, beta2, s){
  mu_1 <- exp(beta1)/(1+exp(beta1))^2
  mu_2 <- exp(beta1 + beta2)/(1+exp(beta1 + beta2))^2
  1/(mu_1 * n[1]) + 1/(mu_2 * n[2])
}


Withprev_formula <- function(lor_old_hat, lor_new_hat, sigma_old_hat, sigma_new_hat){
  X <- matrix(c(1,1),nrow = 2)
  sigma <- diag(c(sigma_old_hat^2, sigma_new_hat^2))
  y <- matrix(c(lor_old_hat, lor_new_hat),nrow = 2)
  lor_hat <- solve(t(X) %*% solve(sigma) %*% X) %*% t(X) %*% solve(sigma) %*% y
  se <- sqrt(solve(t(X) %*% solve(sigma) %*% X))
  power <- ifelse(abs(lor_hat)/se > qnorm(0.975),1,0)
  return(matrix(c(power,lor_hat,se),ncol = 3))
}

Single_study <- function(p_baseline, p_trt2, pig_alloc){
  #power <- numeric(1)
  
  # random part in the simulation, rbinom
  BRD_s <- new_study(p_baseline, p_trt2, sigma = 0, pig_alloc = pig_alloc)
  BRD_new_s <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_s, allstudies = T, sm = "OR")
  # TE: Estimate of treatment effect (log odds ratio, mean difference)
  # seTE: S.E. of TE
  # sm: summary measure,
  # comb.fixed: whether a fixed effects (common effects) network meta-analysis should be conducted.
  # comb.random: whether a random effects network meta-analysis should be conducted.
  nma_s <- netmeta::netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new_s,sm="OR",
                            comb.fixed =F,comb.random = T, tol.multiarm.se = 0.01)
  
  lor_hat <- nma_s$TE.fixed[2,1]
  se <- nma_s$seTE.fixed[2,1]
  
  return(c(lor_hat,se))
}


