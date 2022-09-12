# Resampling under the null 
# All sample sizes (using 100 test data sets)

library(tidyverse)
library(readxl)
library(gemtc)
library(igraph)
library(magic)
library(xtable)
library(msm)
library(netmeta)
library(tictoc)
source("Code/Functions_Project_2.R")


# Simulation (Null) -------------------------------------------------------

library(doParallel)
myCluster <- makeCluster(10)
registerDoParallel(myCluster)
res_null = foreach(i=c(50,100,150,200), .combine='list', .multicombine=TRUE) %dopar% {
  
  library(tidyverse)
  library(gemtc)
  library(igraph)
  library(magic)
  library(xtable)
  library(msm)
  library(netmeta)
  
  sample_size=i
  load(paste0("Data/New_Trial_Data_n",sample_size,".RData"))
  
  trial_1_cond = Filter(function(x){
    p = as.numeric(x[["decision_p"]])
    p < .1 & p > .05},
    trial_all_sim)
  
  trial_1_cond_100 = trial_1_cond[c(1:100)]
  
  resample_sim = lapply(trial_1_cond_100,function(x){
    
    data_new = x$df_new[,-c(24:26)]
    
    # Existing network info
    obs_pval_exg = as.numeric(x$decision_p)
      
    # NMA on new
    nma_new = nma_arm_analysis(data_new)
    
    # Observed test statistic
    obs_test = as.numeric(nma_new$c_z)
    obs_pval = as.numeric(nma_new$c_pval)
    
    # Get risks
    
    BRD_long = wide2long(data_new)
    risk_NAC <- BRD_long %>% 
      filter(t == "No active control") %>% 
      summarise(p = sum(r)/sum(n)) %>% 
      pull
    
    risk_all = c(risk_NAC*exp(nma_new$mu_hat)/(risk_NAC*exp(nma_new$mu_hat) + 1 - risk_NAC))
    names(risk_all) = nma_new$non_baseline_trts
    
    risk_table = data.frame(p = risk_all,
                            trt = names(risk_all)) %>%
      rbind(c(risk_NAC,"No active control")) %>%
      mutate(p = as.numeric(p))
    
    rownames(risk_table) = seq(1,13)
    
    # Common risk
    n_A = sample_size/2
    n_B = sample_size/2
    p_A_all = risk_table[risk_table$trt=="Ceftiofur hydrochloride","p"]
    p_B_all = risk_table[risk_table$trt=="Trimethoprim","p"]
    p_AB_pooled = p_A_all*(180+n_A)/(180+n_A+233+n_B)+p_B_all*(233+n_B)/(180+n_A+233+n_B)
    risk_table[risk_table$trt=="Trimethoprim","p"] <- p_AB_pooled
    risk_table[risk_table$trt=="Ceftiofur hydrochloride","p"] <- p_AB_pooled
    
    # Add risk info
    
    BRD_extend <- data_new
    risk_1 <- risk_table
    colnames(risk_1) <- c("p1","Arm.1")
    risk_2 <- risk_table
    colnames(risk_2) <- c("p2","Arm.2")
    risk_3<- risk_table
    colnames(risk_3) <- c("p3","Arm.3")
    
    BRD_extend <- BRD_extend %>%
      merge(risk_1, by="Arm.1",all.x=T) %>%
      merge(risk_2, by="Arm.2",all.x=T) %>%
      merge(risk_3, by="Arm.3",all.x=T) %>%
      arrange(Study.number)
    
    # Re-sample data
    
    n_iter <- 10000 #10000
    
    sim_test = sapply(1:n_iter, function(x){
      
      # Replacing number of events in a trial from simulated binomial
      BRD_extend$Number.of.Event.in.arm.1=rbinom(nrow(BRD_extend),
                                                 size =BRD_extend$Total.number.in.arm.1,
                                                 prob = BRD_extend$p1)
      BRD_extend$Number.of.Event.in.arm.2=rbinom(nrow(BRD_extend),
                                                 size =BRD_extend$Total.number.in.arm.2,
                                                 prob = BRD_extend$p2)
      BRD_extend$Number.of.Event.in.arm.3=rbinom(nrow(BRD_extend),
                                                 size =BRD_extend$Total.number.in.arm.3,
                                                 prob = BRD_extend$p3)
      
      # If number of events is 0 or equal to n, we need to correct
      
      BRD_extend[BRD_extend$Number.of.Event.in.arm.1 == BRD_extend$Total.number.in.arm.1,"Number.of.Event.in.arm.1"] <-
        BRD_extend[BRD_extend$Number.of.Event.in.arm.1 == BRD_extend$Total.number.in.arm.1,"Number.of.Event.in.arm.1"] - .5
      BRD_extend[BRD_extend$Number.of.Event.in.arm.2 == BRD_extend$Total.number.in.arm.2,"Number.of.Event.in.arm.2"] <-
        BRD_extend[BRD_extend$Number.of.Event.in.arm.2 == BRD_extend$Total.number.in.arm.2,"Number.of.Event.in.arm.2"] - .5
      # BRD_extend[BRD_extend$Number.of.Event.in.arm.3 == BRD_extend$Total.number.in.arm.3,"Number.of.Event.in.arm.3"] <-
      #   BRD_extend[BRD_extend$Number.of.Event.in.arm.3 == BRD_extend$Total.number.in.arm.3,"Number.of.Event.in.arm.3"] - .5
      #
      
      BRD_extend[BRD_extend$Number.of.Event.in.arm.1 == 0,"Number.of.Event.in.arm.1"] <- .5
      BRD_extend[BRD_extend$Number.of.Event.in.arm.2 == 0,"Number.of.Event.in.arm.2"] <- .5
      # BRD_extend[BRD_extend$Number.of.Event.in.arm.3 == 0,"Number.of.Event.in.arm.3"] <- .5
      
      # New NMA analysis
      nma_brd_new = nma_arm_analysis(BRD_extend)
      return(c(nma_brd_new$c_z,nma_brd_new$c_pval))
    },simplify=F)
    
    # Filter the data to meet the criteria
    df = do.call(rbind,sim_test) %>% as.data.frame() %>%
      filter(V2 < .1 & V2 > .05)
    
    # Return re-sampled p-value
    
    p_1 = sum(abs(df$V1) >= abs(obs_test)) / nrow(df)
    
    # Other criteria of filtering the data and corresponding p-value
    
    df_2 = do.call(rbind,sim_test) %>% as.data.frame() %>%
      filter(V2 < obs_pval_exg)
    p_2 = sum(abs(df_2$V1) >= abs(obs_test)) / nrow(df_2)
    
    return(c(p_1 = p_1,
             p_2 = p_2))
  })
  
  df_final = do.call(rbind,resample_sim) %>%
    as.data.frame()
  
  return(c(i,
         sum(df_final$p_1 < .05)/nrow(df_final),
         sum(df_final$p_2 < .05)/nrow(df_final)))
}
stopCluster(myCluster)


# Saving Results ----------------------------------------------------------

df_final_all = do.call(rbind, res_null)
colnames(df_final_all) = c("n", 
                       "resample_error_1",
                       "resample_error_2")
xtable(df_final_all %>% 
         as.data.frame(),
       digits = c(0,0,2,2))

write.csv(df_final_all,file="Data/sim_null_05_1_100.csv")