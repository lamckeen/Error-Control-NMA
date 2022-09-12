# Resampling procedure: Example
# For filtering criteria under null of p<.1

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


# Sample Data -------------------------------------------------------------

## Original Nework data (from Fangshu)
nw_dat = read.csv("Data/BRD_Arm_Network.csv", stringsAsFactors = F)
nw_dat$Number.of.Event.in.arm.2[nw_dat$Study.number == 2] <- 
  nw_dat$Number.of.Event.in.arm.2[nw_dat$Study.number == 2] + .5
nw_dat$Total.number.in.arm.2[nw_dat$Study.number == 2] <- 
  nw_dat$Total.number.in.arm.2[nw_dat$Study.number == 2] + 1

## Example data set with added new study between treatments A and B
sample_size = 50 #25 in arm 1, 25 in arm 2
test_dat = nw_dat %>% add_row(Study.number = 99,
                              Number.of.Event.in.arm.1 = 10,
                              Number.of.Event.in.arm.2 = 18,
                              Total.number.in.arm.1 = 25,
                              Total.number.in.arm.2 = 25,
                              Number.of.arms = 2,
                              Total = 50,
                              Arm.1 = "Ceftiofur hydrochloride",
                              Arm.2 = "Trimethoprim",
                              Arm.3 = "")

## NMA analysis on new data set (using my own function - 100x quicker than netmeta)
nma_brd_matrix_test = nma_arm_analysis(test_dat)
nma_brd_matrix_test$c_pval
nma_brd_matrix_test$c_mu

## Using package netmeta (verifying my function)
BRD_long <- wide2long(test_dat)
BRD_pair <- netmeta::pairwise(treat = t, 
                              event = r, 
                              n = n, 
                              studlab = id, 
                              data = BRD_long, 
                              allstudies = T, 
                              sm = "OR",
                              reference.group = "No active control")
nma_old <- netmeta(TE,
                   seTE,
                   treat1,
                   treat2,
                   studlab,
                   data=BRD_pair,
                   sm="OR",
                   fixed = T,
                   random = F,
                   reference.group = "No active control" )
nma_old
nma_old$pval.fixed[12,1]

# Resampling Procedure ----------------------------------------------------

resample_v3 = function(nw_dat_new,n_iter){
  
  data_new = nw_dat_new

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
  
  # Estimate common risk under null
  n_A = sample_size/2
  n_B = sample_size/2
  p_A_all = risk_table[risk_table$trt=="Ceftiofur hydrochloride","p"]
  p_B_all = risk_table[risk_table$trt=="Trimethoprim","p"]
  p_AB_pooled = p_A_all*(180+n_A)/(180+n_A+233+n_B)+p_B_all*(233+n_B)/(180+n_A+233+n_B) #180,233 from original network
  risk_table[risk_table$trt=="Trimethoprim","p"] <- p_AB_pooled
  risk_table[risk_table$trt=="Ceftiofur hydrochloride","p"] <- p_AB_pooled
  
  # Add risk info to df
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
  
  # Re-sample data under the null
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
    
    # If number of events is 0 or equal to n, we need to correct (add .5/1)
    BRD_extend = BRD_extend %>%
      mutate(Number.of.Event.in.arm.1 = ifelse(Number.of.Event.in.arm.1 == 0 | Number.of.Event.in.arm.1 == Total.number.in.arm.1,
                                               Number.of.Event.in.arm.1+.5, Number.of.Event.in.arm.1),
             Total.number.in.arm.1 = ifelse(Number.of.Event.in.arm.1 == .5 | Number.of.Event.in.arm.1 == Total.number.in.arm.1+.5,
                                            Total.number.in.arm.1+1, Total.number.in.arm.1)) %>%
      mutate(Number.of.Event.in.arm.2 = ifelse(Number.of.Event.in.arm.2 == 0 | Number.of.Event.in.arm.2 == Total.number.in.arm.2,
                                               Number.of.Event.in.arm.2+.5, Number.of.Event.in.arm.2),
             Total.number.in.arm.2 = ifelse(Number.of.Event.in.arm.2 == 0 | Number.of.Event.in.arm.2 == Total.number.in.arm.2+.5,
                                            Total.number.in.arm.2+1,Total.number.in.arm.2)) %>%
      mutate(Number.of.Event.in.arm.3 = ifelse(Number.of.Event.in.arm.3 == 0+.5 | Number.of.Event.in.arm.3 == Total.number.in.arm.3,
                                               Number.of.Event.in.arm.3+.5,Number.of.Event.in.arm.3),
             Total.number.in.arm.3 = ifelse(Number.of.Event.in.arm.3 == 0+.5 | Number.of.Event.in.arm.3 == Total.number.in.arm.3+.5,
                                            Total.number.in.arm.3+1,Total.number.in.arm.3))
    
    # New NMA analysis
    nma_brd_new = nma_arm_analysis(BRD_extend)
    
    # Return test statistic and p-value
    return(c(nma_brd_new$c_z,nma_brd_new$c_pval))
    
  },simplify=F)
  
  # Filter the data to meet the criteria
  df = do.call(rbind,sim_test) %>% as.data.frame() %>%
    filter(V2 < .1)
    # filter(V2 < .1 & V2 > .05)
  
  # Return re-sampled p-value
  p_1 = sum(abs(df$V1) >= abs(obs_test)) / nrow(df)
  return(p_1)
}

# Testing Resampling on Sample Data ---------------------------------------

resample_v3(test_dat, n_iter=10000)
