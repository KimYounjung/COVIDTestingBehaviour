# log-likelihood function
log_post = function(par, pop, fitted, alpha.dom, delta.dom, omicron.dom, covid, vacc.dom, term, event, lockdown, capacity, lfd, pos.age.region, neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3) {
  
  # true_pos: true positives
  # false_pos: false positives
  true_pos  = pop*fitted*inv.logit((I_age %*% par[c(1:4)]) + (I_region_3 %*% par[c(9:16)]) + par[17]*alpha.dom + par[19]*delta.dom + par[21]*omicron.dom + (I_age_3 %*% par[c(23:25)])*vacc.dom + (I_age_2 %*% par[c(29:32)])*covid/1e2 + par[33]*capacity/1e5 + par[34]*lockdown + par[35]*I_age_1*term + par[36]*event + par[38]*lfd)*TP 
  false_pos = pop*(1-fitted)*inv.logit((I_age %*% par[c(5:8)]) + (I_region_3 %*% par[c(9:16)]) + par[18]*alpha.dom + par[20]*delta.dom + par[22]*omicron.dom + (I_age_3 %*% par[c(26:28)])*vacc.dom + (I_age_2 %*% par[c(29:32)])*covid/1e2 + par[33]*capacity/1e5 + par[34]*lockdown + par[35]*I_age_1*term + par[37]*event + par[39]*lfd)*FP 
  
  # true_neg: true negatives
  # false_neg: false negatives
  true_neg  = pop*(1-fitted)*inv.logit((I_age %*% par[c(5:8)]) + (I_region_3 %*% par[c(9:16)]) + par[18]*alpha.dom + par[20]*delta.dom + par[22]*omicron.dom + (I_age_3 %*% par[c(26:28)])*vacc.dom + (I_age_2 %*% par[c(29:32)])*covid/1e2 + par[33]*capacity/1e5 + par[34]*lockdown + par[35]*I_age_1*term + par[37]*event + par[39]*lfd)*TN
  false_neg = pop*fitted*inv.logit((I_age %*% par[c(1:4)]) + (I_region_3 %*% par[c(9:16)]) + par[17]*alpha.dom + par[19]*delta.dom + par[21]*omicron.dom + (I_age_3 %*% par[c(23:25)])*vacc.dom + (I_age_2 %*% par[c(29:32)])*covid/1e2 + par[33]*capacity/1e5 + par[34]*lockdown + par[35]*I_age_1*term + par[36]*event + par[38]*lfd)*FN 
  
  pos_mean = true_pos + false_pos # no. of people determined as infected = true_pos + false_pos 
  neg_mean = true_neg + false_neg # no. of people determined as non-infected = true_neg + false_neg 

  # no. of people determined as infected or non-infected follows 
  # the negative binomial distribution with mean of pos_mean (NHS testing data) 
  # and dispersion parameters set separately for the two groups. 
  LLik_pos = sum(dnbinom(x = pos.age.region, size = par[40], mu = pos_mean, log = T))
  LLik_neg = sum(dnbinom(x = neg.age.region, size = par[41], mu = neg_mean, log = T))
  
  return(sum(LLik_pos, LLik_neg))
  
}

no.par.1 = 41         # total no. of parameters 
no.par.2 = no.par.1-2
no.par.3 = no.par.1+1
no.dis.1 = no.par.1-1
no.dis.2 = no.par.1

library(truncnorm)
library(boot)

for(c in 1:4) {
  
  # No. of iterations for checking the acceptance rate     
  check = 200  
  
  # Set random initial parameter and step size values 
  par = c(runif(no.par.2, min = -5, max = 5), runif(2, min = 1e-10, max = 5))
  par_sd = c(rep(1, no.par.1))
  
  iter = 1e5          # the number of iterations  
  TP = 0.950          # test sensitivity 
  TN = 0.999          # test specificity
  FN = 0.050          # 1 - test sensitivity
  FP = 0.001          # 1 - test specificity
  
  # A data frame to store accepted parameter values 
  MC_P = as.data.frame(matrix(nrow = iter, ncol = no.par.3)) 
  
  # Other parameters for the MH mcmc algorithm
  att = rep(0,no.par.1); acc = rep(0,no.par.1); Prob = rep(0,no.par.1)
  
  # MH mcmc process
  for(i in 1:iter) {
    
    # For parameters excluding dispersion parameters 
    for(j in 1:no.par.2) {
      
      att[j] = att[j] + 1  
      
      can = par 
      can[j] = rnorm(1, mean = can[j], sd = par_sd[j]) # Normal distributions are used for proposing parameter values
      
      # Post_current_par: posterior log-likelihood based on parameter values in the preceding iteration 
      # Post_updated_par: posterior log-likelihood based on newly proposed parameter values
      Post_current_par = log_post(par, data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
      Post_updated_par = log_post(can, data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
      
      # Obtain the ratio of posterior probabilities 
      Prob[j] = exp(Post_updated_par - Post_current_par)
      
      # Determine whether or not to select newly proposedd parameter values 
      if(runif(1) <= min(1,Prob[j])) {
        
        par = can
        acc[j] = acc[j] + 1
        
      }}
    
    # For dispersion parameters 
    for(j in no.dis.1:no.dis.2) {
      
      att[j] = att[j] + 1  
      
      can = par 
      can[j] = truncnorm::rtruncnorm(1, mean = can[j], sd = par_sd[j], a = 1e-10, b = Inf) # Truncated normal distributions are used for proposing parameter values
      
      Post_current_par = log_post(par, data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
      Post_updated_par = log_post(can, data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
      
      # Obtain the ratio of posterior probabilities, accounting for the asymmetrical property of truncated normal distributions 
      Prob[j] = exp(Post_updated_par - Post_current_par)*
        truncnorm::dtruncnorm(par[j], mean = can[j], sd = par_sd[j], a = 1e-10, b = Inf)/
        truncnorm::dtruncnorm(can[j], mean = par[j], sd = par_sd[j], a = 1e-10, b = Inf)
      
      if(runif(1) <= min(1,Prob[j])) {
        
        par = can
        acc[j] = acc[j] + 1
        
      }}
    
    # Store updated parameter values and the log-likelihood based on those values 
    MC_P[i,c(1:no.par.1)] = par 
    MC_P[i,no.par.3]  = log_post(par, data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
    
    # Tune the acceptance rate (20-30% limits)
    for(j in 1:no.par.1){
      
      if(att[j] == check){
        
        print(paste0("Can sd of ", round(par_sd[j], 3),
                     " for par[",j,"] gave acc rate ",acc[j]/att[j], " at ",i," iterations")) 
        
        if(acc[j]/att[j] < 0.2){ par_sd[j] = par_sd[j]*0.5 }
        if(acc[j]/att[j] > 0.3){ par_sd[j] = par_sd[j]*1.5 }
        
        acc[j] = att[j] = 0 
        
      }}
  } 
  write.csv(MC_P, file = paste(c, "Model_370.csv", sep =""),  row.names = FALSE)
}


