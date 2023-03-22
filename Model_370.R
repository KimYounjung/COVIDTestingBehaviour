setwd("C:/Users/yk329/OneDrive - University of Sussex/Documents/GitHub/COVIDSurveillance")


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

for(c in 1:1) {
  
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
    
    # Real-time trace-plots 
    # t = i %% 1001 
    # if(t == 0) {
    #   graphics.off()
    #   par(mfrow=c(3,3))
    #   for(k in 1:no.par.3) {
    #     plot(MC_P[,k][c(1:which(is.na(MC_P[,k]))[1]-1)], type = "l", ylab = "", xlab = MC_P[,k][which(is.na(MC_P[,k]))[1]-1])
    #     title(main= paste("par", k, sep = ""))
    #   }}
  } 
  write.csv(MC_P, file = paste(c, "Model_370.csv", sep =""),  row.names = FALSE)
}


################################################################################

setwd("C:/Users/yk329/OneDrive - University of Sussex/Documents/GitHub/COVIDSurveillance")

library(scales)
library(boot)
library(coda)
source("DBDA2E-utilities.R")

MC = vector(mode = "list", length = 2)

MC[[1]] = read.csv(file = "1Model_370.csv", header = TRUE)
MC[[2]] = read.csv(file = "2Model_370.csv", header = TRUE)

tail(MC[[1]])

for(i in 1:2) {
  MC[[i]][,33] = MC[[i]][,33]*(max(data$capacity)-min(data$capacity))/1e5
  MC[[i]][,c(29:32)] = MC[[i]][,c(29:32)]*(max(data$covid)-min(data$covid))/1e2
}

burn = 1e4+1
iter = 1e5
post_mcmc = vector(mode = "list", length = 2)
for(i in 1:2) {
  post_mcmc[[i]] = coda::mcmc(MC[[i]][burn:iter,])}
post_mcmc_list = coda::mcmc.list(post_mcmc)

coda::effectiveSize(post_mcmc_list)
coda::autocorr.plot(post_mcmc_list, lag.max = 150)
gelman = coda::gelman.diag(post_mcmc_list)
gelman$psrf
summary(post_mcmc_list)

library(hBayesDM)
MC = rbind(MC[[1]][burn:iter,], MC[[2]][burn:iter,])

median = vector(length = ncol(MC))
for(i in 1:ncol(MC)) {
  median[i] = exp(median(MC[,i]))}

pos.1 = data.frame(var = colnames(I_age),
                   median = median[c(1:4)],
                   t(exp(HDInterval::hdi(post_mcmc_list[,c(1:4)], credMass = 0.95))))
pos.2 = data.frame(var = c("alpha",
                           "delta",
                           "omicron",
                           "age group 2", "age group 3", "age group 4", 
                           "event", "lfd"),
                   median = median[c(17,19,21,23:25,36,38)],
                   t(exp(HDInterval::hdi(post_mcmc_list[,c(17,19,21,23:25,36,38)], credMass = 0.95))))

neg.1 = data.frame(var = colnames(I_age),
                   median = median[c(5:8)],
                   t(exp(HDInterval::hdi(post_mcmc_list[,c(5:8)], credMass = 0.95))))
neg.2 = data.frame(var = c("alpha",
                           "delta",
                           "omicron",
                           "age group 2", "age group 3", "age group 4", 
                           "event","lfd"),
                   median = median[c(18,20,22,26:28,37,39)],
                   t(exp(HDInterval::hdi(post_mcmc_list[,c(18,20,22,26:28,37,39)], credMass = 0.95))))

uni.1 = data.frame(var = c("age group 1", "age group 2", "age group 3", "age group 4",
                           "capacity",
                           "lockdown",
                           "term"),
                   median = median[c(29:32,33,34,35)],
                   t(exp(HDInterval::hdi(post_mcmc_list[,c(29:32,33,34,35)], credMass = 0.95))))
uni.2 = data.frame(var = colnames(I_region_3),
                   median = median[c(9:16)],
                   t(exp(HDInterval::hdi(post_mcmc_list[,c(9:16)], credMass = 0.95))))

pos.odds = rbind(pos.1, pos.2)
neg.odds = rbind(neg.1, neg.2)
uni.odds = rbind(uni.1, uni.2)

pos.odds$type = "pos"
neg.odds$type = "neg"
uni.odds$type = "uni"

pos.odds = pos.odds[-1,]
neg.odds = neg.odds[-1,]

odds = rbind(pos.odds, neg.odds, uni.odds)

uni.odds$y = rev(seq(10, 300, 20))
pos.odds$y = rev(seq(310, 520, 20))
neg.odds$y = rev(seq(310, 520, 20))

pos.odds$y = pos.odds$y + 2
neg.odds$y = neg.odds$y - 2



library(scales)
par(mar = c(2, 15, 2, 2)) 
plot(neg.odds$median, neg.odds$y, pch = 15, cex = 1.7, 
     col = alpha("darkgreen", 1), 
     xlim = c(min(pos.odds$lower, neg.odds$lower, uni.odds$lower), 
              max(pos.odds$upper, neg.odds$upper, uni.odds$upper)),
     ylim = c(20, 505),
     yaxt="n", xaxt="n", ylab = "", log = "x")
rect(xleft = 1e-10, ybottom = 440, xright = 1000, ytop = 460,
     col = alpha("gray90", 0.6), border = NA)
rect(xleft = 1e-10, ybottom = 400, xright = 1000, ytop = 420,
     col = alpha("gray90", 0.6), border = NA)
rect(xleft = 1e-10, ybottom = 320, xright = 1000, ytop = 340,
     col = alpha("gray90", 0.6), border = NA)
rect(xleft = 1e-10, ybottom = 220, xright = 1000, ytop = 300,
     col = alpha("gray90", 0.6), border = NA)
rect(xleft = 1e-10, ybottom = 180, xright = 1000, ytop = 200,
     col = alpha("gray90", 0.6), border = NA)
rect(xleft = 1e-10, ybottom = 0, xright = 1000, ytop = 160,
     col = alpha("gray90", 0.6), border = NA)
abline(h = seq(0,500,20), col = "gray50", lty = 1, lwd = 0.1)
abline(v = 1, col = "black", lty = 2, lwd = 1.5)
points(neg.odds$median, neg.odds$y, pch = 15, cex = 1.7, col = alpha("darkgreen", 1))
arrows(neg.odds$lower, neg.odds$y,
       neg.odds$upper, neg.odds$y,
       length = 0, code = 3, col = alpha("darkgreen", 1), lwd = 1.5) 
points(pos.odds$median, pos.odds$y, pch = 17, cex = 1.7, col = alpha("red", 1))
arrows(pos.odds$lower, pos.odds$y,
       pos.odds$upper, pos.odds$y,
       length = 0, code = 3, col = alpha("red", 1), lwd = 1.5) 
points(uni.odds$median, uni.odds$y, pch = 16, cex = 1.7, col = alpha("black", 1))
arrows(uni.odds$lower, uni.odds$y,
       uni.odds$upper, uni.odds$y,
       length = 0, code = 3, col = alpha("black", 1), lwd = 1.5) 
axis(side = 2, at=rev(seq(10, 520, 20)), 
     labels=c(neg.odds$var,uni.odds$var),
     las = 2, cex.axis = 1.3)
axis(side = 1, at=c(0.25, 0.5, 1e0, 2, 4, 8, 16, 32), 
     labels=c("0.25", "0.5", "1", "2", "4", "8", "16", "32"), cex.axis = 1.3)
title(main= "Odds ratio for taking a PCR test") 


pos.odds$CI = paste(round(pos.odds$lower,2), "-", round(pos.odds$upper,2), sep="")
neg.odds$CI = paste(round(neg.odds$lower,2), "-", round(neg.odds$upper,2), sep="")
uni.odds$CI = paste(round(uni.odds$lower,2), "-", round(uni.odds$upper,2), sep="")

head(pos.odds)
write.csv(pos.odds, file = "pos.OR from M370.csv", row.names = FALSE)
write.csv(neg.odds, file = "neg.OR from M370.csv", row.names = FALSE)
write.csv(uni.odds, file = "uni.OR from M370.csv", row.names = FALSE)

################################################################################

MC = vector(mode = "list", length = 3)

MC[[1]] = read.csv(file = "1Model_200_NB2.csv", header = TRUE)
MC[[2]] = read.csv(file = "3Model_200_NB2.csv", header = TRUE)
MC[[3]] = read.csv(file = "4Model_200_NB2.csv", header = TRUE)

for(i in 1:3) {
  MC[[i]][,c(27:34)] = MC[[i]][,c(27:34)]*(max(data$covid)-min(data$covid))/1e2
  MC[[i]][,c(35,36)] = MC[[i]][,c(35,36)]*0.7250503 # max vaccinated proportion across regions and age groups 
  MC[[i]][,c(37:44)] = MC[[i]][,c(37:44)]/max(data$alpha)
  MC[[i]][,c(45:52)] = MC[[i]][,c(45:52)]/max(data$delta)
  MC[[i]][,c(53:60)] = MC[[i]][,c(53:60)]/max(data$omicron)
}
burn = 1e4+1
iter = 1e5
MC = rbind(MC[[1]][burn:iter,], MC[[2]][burn:iter,], MC[[3]][burn:iter,])

pos.age.group = vector(mode = "list", length = 4)
neg.age.group = vector(mode = "list", length = 4)

for(i in 1:4) {
  pos.age.group[[i]] = vector(mode = "list", length = 4)
  neg.age.group[[i]] = vector(mode = "list", length = 4)
}


### wild ### 
pos.age.group[[1]][[1]] = inv.logit(MC[,1]) 
pos.age.group[[1]][[2]] = inv.logit(MC[,1] + MC[,10]) 
pos.age.group[[1]][[3]] = inv.logit(MC[,1] + MC[,11]) 
pos.age.group[[1]][[4]] = inv.logit(MC[,1] + MC[,12]) 
neg.age.group[[1]][[1]] = inv.logit(MC[,13]) 
neg.age.group[[1]][[2]] = inv.logit(MC[,13] + MC[,22]) 
neg.age.group[[1]][[3]] = inv.logit(MC[,13] + MC[,23]) 
neg.age.group[[1]][[4]] = inv.logit(MC[,13] + MC[,24]) 

### alpha ### 
pos.age.group[[2]][[1]] = inv.logit(MC[,1] + MC[,37]) 
pos.age.group[[2]][[2]] = inv.logit(MC[,1] + MC[,10] + MC[,38]) 
pos.age.group[[2]][[3]] = inv.logit(MC[,1] + MC[,11] + MC[,39]) 
pos.age.group[[2]][[4]] = inv.logit(MC[,1] + MC[,12] + MC[,40]) 
neg.age.group[[2]][[1]] = inv.logit(MC[,13] + MC[,41]) 
neg.age.group[[2]][[2]] = inv.logit(MC[,13] + MC[,22] + MC[,42]) 
neg.age.group[[2]][[3]] = inv.logit(MC[,13] + MC[,23] + MC[,43]) 
neg.age.group[[2]][[4]] = inv.logit(MC[,13] + MC[,24] + MC[,44]) 

### delta ### 
pos.age.group[[3]][[1]] = inv.logit(MC[,1] + MC[,45]) 
pos.age.group[[3]][[2]] = inv.logit(MC[,1] + MC[,10] + MC[,46]) 
pos.age.group[[3]][[3]] = inv.logit(MC[,1] + MC[,11] + MC[,47]) 
pos.age.group[[3]][[4]] = inv.logit(MC[,1] + MC[,12] + MC[,48]) 
neg.age.group[[3]][[1]] = inv.logit(MC[,13] + MC[,49]) 
neg.age.group[[3]][[2]] = inv.logit(MC[,13] + MC[,22] + MC[,50]) 
neg.age.group[[3]][[3]] = inv.logit(MC[,13] + MC[,23] + MC[,51]) 
neg.age.group[[3]][[4]] = inv.logit(MC[,13] + MC[,24] + MC[,52]) 

### omicron ### 
pos.age.group[[4]][[1]] = inv.logit(MC[,1] + MC[,53]) 
pos.age.group[[4]][[2]] = inv.logit(MC[,1] + MC[,10] + MC[,54]) 
pos.age.group[[4]][[3]] = inv.logit(MC[,1] + MC[,11] + MC[,55]) 
pos.age.group[[4]][[4]] = inv.logit(MC[,1] + MC[,12] + MC[,56]) 
neg.age.group[[4]][[1]] = inv.logit(MC[,13] + MC[,57]) 
neg.age.group[[4]][[2]] = inv.logit(MC[,13] + MC[,22] + MC[,58]) 
neg.age.group[[4]][[3]] = inv.logit(MC[,13] + MC[,23] + MC[,59]) 
neg.age.group[[4]][[4]] = inv.logit(MC[,13] + MC[,24] + MC[,60]) 

pos.ratio.1 = vector(mode = "list", length = 4)
pos.ratio.2 = vector(mode = "list", length = 4)
pos.ratio.3 = vector(mode = "list", length = 4)
pos.ratio.4 = vector(mode = "list", length = 4)
neg.ratio.1 = vector(mode = "list", length = 4)
neg.ratio.2 = vector(mode = "list", length = 4)
neg.ratio.3 = vector(mode = "list", length = 4)
neg.ratio.4 = vector(mode = "list", length = 4)

for(i in 1:4) {
  pos.ratio.1[[i]] = vector(mode = "list", length = 4)
  pos.ratio.2[[i]] = vector(mode = "list", length = 4)
  pos.ratio.3[[i]] = vector(mode = "list", length = 4)
  pos.ratio.4[[i]] = vector(mode = "list", length = 4)
  neg.ratio.1[[i]] = vector(mode = "list", length = 4)
  neg.ratio.2[[i]] = vector(mode = "list", length = 4)
  neg.ratio.3[[i]] = vector(mode = "list", length = 4)
  neg.ratio.4[[i]] = vector(mode = "list", length = 4)
}

for(i in 1:4) {
  for(j in 1:4) {
    pos.ratio.1[[i]][[j]] = pos.age.group[[i]][[1]]/pos.age.group[[j]][[1]]
    pos.ratio.2[[i]][[j]] = pos.age.group[[i]][[2]]/pos.age.group[[j]][[2]]
    pos.ratio.3[[i]][[j]] = pos.age.group[[i]][[3]]/pos.age.group[[j]][[3]]
    pos.ratio.4[[i]][[j]] = pos.age.group[[i]][[4]]/pos.age.group[[j]][[4]]
    neg.ratio.1[[i]][[j]] = neg.age.group[[i]][[1]]/neg.age.group[[j]][[1]]
    neg.ratio.2[[i]][[j]] = neg.age.group[[i]][[2]]/neg.age.group[[j]][[2]]
    neg.ratio.3[[i]][[j]] = neg.age.group[[i]][[3]]/neg.age.group[[j]][[3]]
    neg.ratio.4[[i]][[j]] = neg.age.group[[i]][[4]]/neg.age.group[[j]][[4]]
    }}


variant = c("pre-alpha", "alpha", "delta", "omicron")

library(scales)
par(mfrow=c(3,1))
pos.delta.alpha = data.frame(x=c(1,4,7,10),
                             y=c(median(pos.ratio.1[[3]][[2]]),median(pos.ratio.2[[3]][[2]]),median(pos.ratio.3[[3]][[2]]),median(pos.ratio.4[[3]][[2]])),
                             y.l = c(HDInterval::hdi(pos.ratio.1[[3]][[2]],credMass = 0.95)[1],
                                     HDInterval::hdi(pos.ratio.2[[3]][[2]],credMass = 0.95)[1],
                                     HDInterval::hdi(pos.ratio.3[[3]][[2]],credMass = 0.95)[1],
                                     HDInterval::hdi(pos.ratio.4[[3]][[2]],credMass = 0.95)[1]),
                             y.h = c(HDInterval::hdi(pos.ratio.1[[3]][[2]],credMass = 0.95)[2],
                                     HDInterval::hdi(pos.ratio.2[[3]][[2]],credMass = 0.95)[2],
                                     HDInterval::hdi(pos.ratio.3[[3]][[2]],credMass = 0.95)[2],
                                     HDInterval::hdi(pos.ratio.4[[3]][[2]],credMass = 0.95)[2]))

neg.delta.alpha = data.frame(x=c(2,5,8,11),
                             y=c(median(neg.ratio.1[[3]][[2]]),median(neg.ratio.2[[3]][[2]]),median(neg.ratio.3[[3]][[2]]),median(neg.ratio.4[[3]][[2]])),
                             y.l = c(HDInterval::hdi(neg.ratio.1[[3]][[2]],credMass = 0.95)[1],
                                     HDInterval::hdi(neg.ratio.2[[3]][[2]],credMass = 0.95)[1],
                                     HDInterval::hdi(neg.ratio.3[[3]][[2]],credMass = 0.95)[1],
                                     HDInterval::hdi(neg.ratio.4[[3]][[2]],credMass = 0.95)[1]),
                             y.h = c(HDInterval::hdi(neg.ratio.1[[3]][[2]],credMass = 0.95)[2],
                                     HDInterval::hdi(neg.ratio.2[[3]][[2]],credMass = 0.95)[2],
                                     HDInterval::hdi(neg.ratio.3[[3]][[2]],credMass = 0.95)[2],
                                     HDInterval::hdi(neg.ratio.4[[3]][[2]],credMass = 0.95)[2]))

plot(pos.delta.alpha$x, pos.delta.alpha$y, 
     pch = 16, col = "red",
     xaxt="n", yaxt="n", ylab = "", xlab = "",
     xlim = c(0,12), ylim = c(0,max(pos.delta.alpha$y.h)), cex = 2)
abline(h = 1, lty = 5, lwd = 1.5)
points(neg.delta.alpha$x, neg.delta.alpha$y, pch = 15, col = "darkgreen", cex = 2)
points(pos.delta.alpha$x, pos.delta.alpha$y, pch = 16, col = "red", cex = 2)
arrows(pos.delta.alpha$x, pos.delta.alpha$y.l,
       pos.delta.alpha$x, pos.delta.alpha$y.h,
       length = 0, code = 3, col = alpha("red", 1), lwd = 2) 
arrows(neg.delta.alpha$x, neg.delta.alpha$y.l,
       neg.delta.alpha$x, neg.delta.alpha$y.h,
       length = 0, code = 3, col = alpha("darkgreen", 1), lwd = 2) 
axis(side = 1, at=c(1.5, 4.5, 7.5, 10.5), 
     labels=c("age group 1", "age group 2", "age group 3", "age group 4"), cex.axis = 2)
axis(side = 2, at=c(0,1,2,3,4,5,6,7,8,9), 
     labels=c(0,1,2,3,4,5,6,7,8,9), cex.axis = 2)
title(main= paste("A. ", variant[3], " over ", variant[2], sep = ""), cex = 2) 

pos.delta.omicron = data.frame(x=c(1,4,7,10),
                               y=c(median(pos.ratio.1[[3]][[4]]),median(pos.ratio.2[[3]][[4]]),median(pos.ratio.3[[3]][[4]]),median(pos.ratio.4[[3]][[4]])),
                               y.l = c(HDInterval::hdi(pos.ratio.1[[3]][[4]],credMass = 0.95)[1],
                                       HDInterval::hdi(pos.ratio.2[[3]][[4]],credMass = 0.95)[1],
                                       HDInterval::hdi(pos.ratio.3[[3]][[4]],credMass = 0.95)[1],
                                       HDInterval::hdi(pos.ratio.4[[3]][[4]],credMass = 0.95)[1]),
                               y.h = c(HDInterval::hdi(pos.ratio.1[[3]][[4]],credMass = 0.95)[2],
                                       HDInterval::hdi(pos.ratio.2[[3]][[4]],credMass = 0.95)[2],
                                       HDInterval::hdi(pos.ratio.3[[3]][[4]],credMass = 0.95)[2],
                                       HDInterval::hdi(pos.ratio.4[[3]][[4]],credMass = 0.95)[2]))

neg.delta.omicron = data.frame(x=c(2,5,8,11),
                               y=c(median(neg.ratio.1[[3]][[4]]),median(neg.ratio.2[[3]][[4]]),median(neg.ratio.3[[3]][[4]]),median(neg.ratio.4[[3]][[4]])),
                               y.l = c(HDInterval::hdi(neg.ratio.1[[3]][[4]],credMass = 0.95)[1],
                                       HDInterval::hdi(neg.ratio.2[[3]][[4]],credMass = 0.95)[1],
                                       HDInterval::hdi(neg.ratio.3[[3]][[4]],credMass = 0.95)[1],
                                       HDInterval::hdi(neg.ratio.4[[3]][[4]],credMass = 0.95)[1]),
                               y.h = c(HDInterval::hdi(neg.ratio.1[[3]][[4]],credMass = 0.95)[2],
                                       HDInterval::hdi(neg.ratio.2[[3]][[4]],credMass = 0.95)[2],
                                       HDInterval::hdi(neg.ratio.3[[3]][[4]],credMass = 0.95)[2],
                                       HDInterval::hdi(neg.ratio.4[[3]][[4]],credMass = 0.95)[2]))
plot(pos.delta.omicron$x, pos.delta.omicron$y, 
     pch = 16, col = "red",
     xaxt="n", yaxt="n", ylab = "", xlab = "",
     xlim = c(0,12), ylim = c(0,max(pos.delta.omicron$y.h)), cex = 2)
abline(h = 1, lty = 5, lwd = 1.5)
points(neg.delta.omicron$x, neg.delta.omicron$y, pch = 15, col = "darkgreen", cex = 2)
points(pos.delta.omicron$x, pos.delta.omicron$y, pch = 16, col = "red", cex = 2)
arrows(pos.delta.omicron$x, pos.delta.omicron$y.l,
       pos.delta.omicron$x, pos.delta.omicron$y.h,
       length = 0, code = 3, col = alpha("red", 1), lwd = 2) 
arrows(neg.delta.omicron$x, neg.delta.omicron$y.l,
       neg.delta.omicron$x, neg.delta.omicron$y.h,
       length = 0, code = 3, col = alpha("darkgreen", 1), lwd = 2) 
axis(side = 1, at=c(1.5, 4.5, 7.5, 10.5), 
     labels=c("age group 1", "age group 2", "age group 3", "age group 4"), cex.axis = 2)
axis(side = 2, at=c(0,1,2,3,4,5,6,7), 
     labels=c(0,1,2,3,4,5,6,7), cex.axis = 2)
title(main= paste("B. ", variant[3], " over ", variant[4], sep = ""), cex = 2) 

pos.omicron.alpha = data.frame(x=c(1,4,7,10),
                               y=c(median(pos.ratio.1[[4]][[2]]),median(pos.ratio.2[[4]][[2]]),median(pos.ratio.3[[4]][[2]]),median(pos.ratio.4[[4]][[2]])),
                               y.l = c(HDInterval::hdi(pos.ratio.1[[4]][[2]],credMass = 0.95)[1],
                                       HDInterval::hdi(pos.ratio.2[[4]][[2]],credMass = 0.95)[1],
                                       HDInterval::hdi(pos.ratio.3[[4]][[2]],credMass = 0.95)[1],
                                       HDInterval::hdi(pos.ratio.4[[4]][[2]],credMass = 0.95)[1]),
                               y.h = c(HDInterval::hdi(pos.ratio.1[[4]][[2]],credMass = 0.95)[2],
                                       HDInterval::hdi(pos.ratio.2[[4]][[2]],credMass = 0.95)[2],
                                       HDInterval::hdi(pos.ratio.3[[4]][[2]],credMass = 0.95)[2],
                                       HDInterval::hdi(pos.ratio.4[[4]][[2]],credMass = 0.95)[2]))

neg.omicron.alpha = data.frame(x=c(2,5,8,11),
                               y=c(median(neg.ratio.1[[4]][[2]]),median(neg.ratio.2[[4]][[2]]),median(neg.ratio.3[[4]][[2]]),median(neg.ratio.4[[4]][[2]])),
                               y.l = c(HDInterval::hdi(neg.ratio.1[[4]][[2]],credMass = 0.95)[1],
                                       HDInterval::hdi(neg.ratio.2[[4]][[2]],credMass = 0.95)[1],
                                       HDInterval::hdi(neg.ratio.3[[4]][[2]],credMass = 0.95)[1],
                                       HDInterval::hdi(neg.ratio.4[[4]][[2]],credMass = 0.95)[1]),
                               y.h = c(HDInterval::hdi(neg.ratio.1[[4]][[2]],credMass = 0.95)[2],
                                       HDInterval::hdi(neg.ratio.2[[4]][[2]],credMass = 0.95)[2],
                                       HDInterval::hdi(neg.ratio.3[[4]][[2]],credMass = 0.95)[2],
                                       HDInterval::hdi(neg.ratio.4[[4]][[2]],credMass = 0.95)[2]))
plot(pos.omicron.alpha$x, pos.omicron.alpha$y, 
     pch = 16, col = "red",
     xaxt="n", yaxt="n", ylab = "", xlab = "",
     xlim = c(0,12), ylim = c(0,max(pos.omicron.alpha$y.h)), cex = 2)
abline(h = 1, lty = 5, lwd = 1.5)
points(neg.omicron.alpha$x, neg.omicron.alpha$y, pch = 15, col = "darkgreen", cex = 2)
points(pos.omicron.alpha$x, pos.omicron.alpha$y, pch = 16, col = "red", cex = 2)
arrows(pos.omicron.alpha$x, pos.omicron.alpha$y.l,
       pos.omicron.alpha$x, pos.omicron.alpha$y.h,
       length = 0, code = 3, col = alpha("red", 1), lwd = 1.5) 
arrows(neg.omicron.alpha$x, neg.omicron.alpha$y.l,
       neg.omicron.alpha$x, neg.omicron.alpha$y.h,
       length = 0, code = 3, col = alpha("darkgreen", 1), lwd = 2) 
axis(side = 1, at=c(1.5, 4.5, 7.5, 10.5), 
     labels=c("age group 1", "age group 2", "age group 3", "age group 4"), cex.axis = 1.5)
axis(side = 2, at=c(0,1,2,3,4,5,6,7), 
     labels=c(0,1,2,3,4,5,6,7), cex.axis = 2)
title(main= paste("C. ", variant[4], " over ", variant[2], sep = ""), cex = 2) 

################################################################################

MC = vector(mode = "list", length = 2)

MC[[1]] = read.csv(file = "1Model_370.csv", header = TRUE)
MC[[2]] = read.csv(file = "2Model_370.csv", header = TRUE)

MC = rbind(MC[[1]][burn:iter,], MC[[2]][burn:iter,])

log_post_pos = function(par, pop, fitted, alpha.dom, delta.dom, omicron.dom, covid, vacc.dom, term, event, lockdown, capacity, lfd, pos.age.region, neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3) {
  
  true_pos  = pop*fitted*inv.logit((I_age %*% par[c(1:4)]) + (I_region_3 %*% par[c(9:16)]) + par[17]*alpha.dom + par[19]*delta.dom + par[21]*omicron.dom + (I_age_3 %*% par[c(23:25)])*vacc.dom + (I_age_2 %*% par[c(29:32)])*covid/1e2 + par[33]*capacity/1e5 + par[34]*lockdown + par[35]*I_age_1*term + par[36]*event + par[38]*lfd)*TP 
  false_pos = pop*(1-fitted)*inv.logit((I_age %*% par[c(5:8)]) + (I_region_3 %*% par[c(9:16)]) + par[18]*alpha.dom + par[20]*delta.dom + par[22]*omicron.dom + (I_age_3 %*% par[c(26:28)])*vacc.dom + (I_age_2 %*% par[c(29:32)])*covid/1e2 + par[33]*capacity/1e5 + par[34]*lockdown + par[35]*I_age_1*term + par[37]*event + par[39]*lfd)*FP 
  
  LLik_pos = sum(dnbinom(x = pos.age.region, size = par[no.dis.1], mu = true_pos + false_pos, log = T))
  return(sum(LLik_pos))}

log_post_neg = function(par, pop, fitted, alpha.dom, delta.dom, omicron.dom, covid, vacc.dom, term, event, lockdown, capacity, lfd, pos.age.region, neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3) {
  
  true_neg  = pop*(1-fitted)*inv.logit((I_age %*% par[c(5:8)]) + (I_region_3 %*% par[c(9:16)]) + par[18]*alpha.dom + par[20]*delta.dom + par[22]*omicron.dom + (I_age_3 %*% par[c(26:28)])*vacc.dom + (I_age_2 %*% par[c(29:32)])*covid/1e2 + par[33]*capacity/1e5 + par[34]*lockdown + par[35]*I_age_1*term + par[37]*event + par[39]*lfd)*TN
  false_neg = pop*fitted*inv.logit((I_age %*% par[c(1:4)]) + (I_region_3 %*% par[c(9:16)]) + par[17]*alpha.dom + par[19]*delta.dom + par[21]*omicron.dom + (I_age_3 %*% par[c(23:25)])*vacc.dom + (I_age_2 %*% par[c(29:32)])*covid/1e2 + par[33]*capacity/1e5 + par[34]*lockdown + par[35]*I_age_1*term + par[36]*event + par[38]*lfd)*FN 
  
  LLik_neg = sum(dnbinom(x = neg.age.region, size = par[no.dis.1], mu = true_neg + false_neg, log = T))
  return(sum(LLik_neg))}

par_m = apply(MC, 2, median)
logL_median_theta = log_post(par_m[-length(par_m)], data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
logL_median_theta_pos = log_post_pos(par_m[-length(par_m)], data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
logL_median_theta_neg = log_post_neg(par_m[-length(par_m)], data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)

nSim = 3000
logL = vector(length = nSim)
logL_pos = vector(length = nSim)
logL_neg = vector(length = nSim)

for(i in 1:nSim) {
  
  k = sample(1:nrow(MC), 1, replace = FALSE)
  par_r = as.numeric(MC[k,])
  logL[i] = log_post(par_r[-length(par_r)], data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
  logL_pos[i] = log_post_pos(par_r[-length(par_r)], data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
  logL_neg[i] = log_post_neg(par_r[-length(par_r)], data$pop, data$fitted, data$alpha.dom, data$delta.dom, data$omicron.dom, data$covid, data$vacc.dom, data$term, data$event, data$lockdown, data$capacity, data$lfd, data$pos.age.region, data$neg.age.region, TP, TN, FP, FN, I_age, I_age_1, I_age_2, I_age_3, I_region_3)
}

dev1 = -2*median(logL)
dev1_pos = -2*median(logL_pos)
dev1_neg = -2*median(logL_neg)

dev2 = -2*logL_median_theta
dev2_pos = -2*logL_median_theta_pos
dev2_neg = -2*logL_median_theta_neg

Pd = dev1-dev2
Pd_pos = dev1_pos-dev2_pos
Pd_neg = dev1_neg-dev2_neg

DIC = Pd+dev1
DIC_pos = Pd_pos+dev1_pos
DIC_neg = Pd_neg+dev1_neg

DIC
DIC_pos
DIC_neg

################################################################################

MC = vector(mode = "list", length = 2)

burn = 1e4+1
iter = 1e5
MC[[1]] = read.csv(file = "1Model_300.csv", header = TRUE)
MC[[2]] = read.csv(file = "2Model_300.csv", header = TRUE)

MC = rbind(MC[[1]][burn:iter,], MC[[2]][burn:iter,])

for(i in 1:nrow(data)) {
  
  true_pos  = data$pop[i]*data$fitted[i]*inv.logit(rowSums(mapply("%*%", I_age[i,], MC[, c(1:4)])) + rowSums(mapply("%*%", I_region_3[i,], MC[, c(9:16)])) + as.numeric(MC[, 17])*data$alpha.dom[i] + as.numeric(MC[, 19])*data$delta.dom[i] + as.numeric(MC[, 21])*data$omicron.dom[i] + rowSums(mapply("%*%", I_age_3[i,], MC[, c(23:25)]))*(data$vacc.dom[i]) + rowSums(mapply("%*%", I_age_2[i,], MC[, c(29:32)]))*(data$covid[i]/1e2) + as.numeric(MC[, 33])*(data$capacity[i]/1e5) + as.numeric(MC[, 34])*(data$lockdown[i]) + as.numeric(MC[, 35])*(I_age_1[i])*(data$term[i]) + as.numeric(MC[, 36])*(data$event[i]))*TP
  false_neg = data$pop[i]*data$fitted[i]*inv.logit(rowSums(mapply("%*%", I_age[i,], MC[, c(1:4)])) + rowSums(mapply("%*%", I_region_3[i,], MC[, c(9:16)])) + as.numeric(MC[, 17])*data$alpha.dom[i] + as.numeric(MC[, 19])*data$delta.dom[i] + as.numeric(MC[, 21])*data$omicron.dom[i] + rowSums(mapply("%*%", I_age_3[i,], MC[, c(23:25)]))*(data$vacc.dom[i]) + rowSums(mapply("%*%", I_age_2[i,], MC[, c(29:32)]))*(data$covid[i]/1e2) + as.numeric(MC[, 33])*(data$capacity[i]/1e5) + as.numeric(MC[, 34])*(data$lockdown[i]) + as.numeric(MC[, 35])*(I_age_1[i])*(data$term[i]) + as.numeric(MC[, 36])*(data$event[i]))*FN 
  
  false_pos = data$pop[i]*(1-data$fitted[i])*inv.logit(rowSums(mapply("%*%", I_age[i,], MC[, c(5:8)])) + rowSums(mapply("%*%", I_region_3[i,], MC[, c(9:16)])) + as.numeric(MC[, 18])*data$alpha.dom[i] + as.numeric(MC[, 20])*data$delta.dom[i] + as.numeric(MC[, 22])*data$omicron.dom[i] + rowSums(mapply("%*%", I_age_3[i,], MC[, c(26:28)]))*(data$vacc.dom[i]) + rowSums(mapply("%*%", I_age_2[i,], MC[, c(29:32)]))*(data$covid[i]/1e2) + as.numeric(MC[, 33])*(data$capacity[i]/1e5) + as.numeric(MC[, 34])*(data$lockdown[i]) + as.numeric(MC[, 35])*(I_age_1[i])*(data$term[i]) + as.numeric(MC[, 37])*(data$event[i]))*FP
  true_neg  = data$pop[i]*(1-data$fitted[i])*inv.logit(rowSums(mapply("%*%", I_age[i,], MC[, c(5:8)])) + rowSums(mapply("%*%", I_region_3[i,], MC[, c(9:16)])) + as.numeric(MC[, 18])*data$alpha.dom[i] + as.numeric(MC[, 20])*data$delta.dom[i] + as.numeric(MC[, 22])*data$omicron.dom[i] + rowSums(mapply("%*%", I_age_3[i,], MC[, c(26:28)]))*(data$vacc.dom[i]) + rowSums(mapply("%*%", I_age_2[i,], MC[, c(29:32)]))*(data$covid[i]/1e2) + as.numeric(MC[, 33])*(data$capacity[i]/1e5) + as.numeric(MC[, 34])*(data$lockdown[i]) + as.numeric(MC[, 35])*(I_age_1[i])*(data$term[i]) + as.numeric(MC[, 37])*(data$event[i]))*TN
  
  data$post_pos_mean[i] = mean(rnbinom(n=nrow(MC), size=as.numeric(MC[, 38]), mu = true_pos + false_pos))
  data$post_neg_mean[i] = mean(rnbinom(n=nrow(MC), size=as.numeric(MC[, 39]), mu = true_neg + false_neg))
  
}

write.csv(data, file = "Model_300_PPC.csv", row.names = FALSE)

################################################################################
setwd("C:/Users/yk329/OneDrive - University of Sussex/Documents/GitHub/COVIDSurveillance")
data = read.csv(file = "Model_300_PPC.csv", header = TRUE)

head(data)

data$date = as.Date(data$date)
colour = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3")

missing = data.frame(missing.1 = as.Date(c("2021-07-19")),
                     missing.2 = as.Date(c("2021-09-05")))
lockdown1st = data.frame(lockdown1st.1 = as.Date(c("2020-03-23")),
                         lockdown1st.2 = as.Date(c("2020-06-23")))
lockdown2nd = data.frame(lockdown2nd.1 = as.Date(c("2020-11-05")),
                         lockdown2nd.2 = as.Date(c("2020-12-02")))
lockdown3rd = data.frame(lockdown3rd.1 = as.Date(c("2021-01-06")),
                         lockdown3rd.2 = as.Date(c("2021-06-14")))
term.1 = as.Date(c("2020-09-03","2020-10-24")) # First autumn term 
term.2 = as.Date(c("2020-10-29","2020-12-17")) # Second autumn term
term.3 = as.Date(c("2021-03-08","2021-04-08")) # Pupils began returning to schools from 8 March 2021
term.4 = as.Date(c("2021-04-25","2021-05-29")) # First summer term   
term.5 = as.Date(c("2021-06-02","2021-07-21")) # Second summer term 
term.6 = as.Date(c("2021-09-03","2021-10-24")) # First autumn term 
term.7 = as.Date(c("2021-10-29","2021-12-17")) # Second autumn term

#omicron arrival: https://www.gov.uk/government/news/first-uk-cases-of-omicron-variant-identified

library(scales)
data_new = vector(mode = "list", length = length(region.list))
for(i in 1:length(region.list)) {
  data_new[[i]] = vector(mode = "list", length = length(age.class.list))}

for(i in 1:length(region.list)) {
  par(mfrow=c(3,4), oma=c(3,3,3,3), mar=c(2,2,4,2))
  data_region = data[data$region == region.list[i],]
  for(j in 1:length(age.class.list)) {
    data_new[[i]][[j]] = data[data$region == region.list[i] & data$age.class == age.class.list[j],]
    plot(x=data_new[[i]][[j]]$date, y=data_new[[i]][[j]]$pos.age.region, 
         type = "l", col = alpha("black", 0.7), ylim = c(0, max(data_region$post_pos_mean, data_region$pos.age.region)),
         xlab="", ylab="No. of PCR positives", xaxt='n', lwd = 1)
    rect(lockdown1st$lockdown1st.1, -1e10, lockdown1st$lockdown1st.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(lockdown2nd$lockdown2nd.1, -1e10, lockdown2nd$lockdown2nd.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(lockdown3rd$lockdown3rd.1, -1e10, lockdown3rd$lockdown3rd.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(term.1[1], -1e10, term.1[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.2[1], -1e10, term.2[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.3[1], -1e10, term.3[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.4[1], -1e10, term.4[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.5[1], -1e10, term.5[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.6[1], -1e10, term.6[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.7[1], -1e10, term.7[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    lines(x=data_new[[i]][[j]]$date, y=data_new[[i]][[j]]$pos.age.region, 
         col = alpha("black", 0.7), ylim = c(0, max(data_region$post_pos_mean, data_region$pos.age.region)),
         xlab="", ylab="", xaxt='n', lwd = 1)
    title(main=age.class.list[j]) 
    axis.Date(1, at=pretty(range(data_new[[i]][[j]]$date), 10), format='%b %Y', las=2, 
              cex.axis=0.8)
    polygon(c(data_new[[i]][[j]]$date, rev(data_new[[i]][[j]]$date)), 
            c(data_new[[i]][[j]]$pos.age.region, rev(data_new[[i]][[j]]$post_pos_mean)), 
            col = alpha("blue", 0.1), border=NA)
    lines(data_new[[i]][[j]]$date, data_new[[i]][[j]]$post_pos_mean, lwd = 2, col = "blue", lty = 1)}

  for(j in 1:length(age.class.list)) {
    
    plot(x=data_new[[i]][[j]]$date, y=data_new[[i]][[j]]$neg.age.region, 
         type = "l", col = "black", ylim = c(0, max(data_region$post_neg_mean, data_region$neg.age.region)),
         xlab="", ylab="No. of PCR negatives", xaxt='n', lwd = 1)
    #rect(missing$missing.1, -1e10, missing$missing.2, 1e10,
    #     col = alpha("gray90", 0.6), border = NA)
    rect(lockdown1st$lockdown1st.1, -1e10, lockdown1st$lockdown1st.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(lockdown2nd$lockdown2nd.1, -1e10, lockdown2nd$lockdown2nd.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(lockdown3rd$lockdown3rd.1, -1e10, lockdown3rd$lockdown3rd.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(term.1[1], -1e10, term.1[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.2[1], -1e10, term.2[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.3[1], -1e10, term.3[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.4[1], -1e10, term.4[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.5[1], -1e10, term.5[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.6[1], -1e10, term.6[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.7[1], -1e10, term.7[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    lines(x=data_new[[i]][[j]]$date, y=data_new[[i]][[j]]$neg.age.region, 
         col = "black", ylim = c(0, max(data_region$post_neg_mean, data_region$neg.age.region)),
         xlab="", ylab="No. of PCR negatives", xaxt='n', lwd = 1)
    axis.Date(1, at=pretty(range(data_new[[i]][[j]]$date), 10), format='%b %Y', las=2, 
              cex.axis=0.8)
    polygon(c(data_new[[i]][[j]]$date, rev(data_new[[i]][[j]]$date)), 
            c(data_new[[i]][[j]]$neg.age.region, rev(data_new[[i]][[j]]$post_neg_mean)), 
            col = alpha("red", 0.1), border=NA)
    lines(data_new[[i]][[j]]$date, data_new[[i]][[j]]$post_neg_mean, lwd = 2, col = "red", lty = 1)}
  
  
  for(j in 1:length(age.class.list)) {
    
    plot(x=data_new[[i]][[j]]$date, y=data_new[[i]][[j]]$fitted, 
         type = "l", col = alpha("black", 0.7), ylim = c(0, max(data_region$upper, data_region$pre_fitted_median)),
         xlab="", ylab="REACT prevalence", xaxt='n', lwd = 1)
    rect(lockdown1st$lockdown1st.1, -1e10, lockdown1st$lockdown1st.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(lockdown2nd$lockdown2nd.1, -1e10, lockdown2nd$lockdown2nd.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(lockdown3rd$lockdown3rd.1, -1e10, lockdown3rd$lockdown3rd.2, 1e10,
         col = alpha("gray", 0.3), border = NA)
    rect(term.1[1], -1e10, term.1[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.2[1], -1e10, term.2[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.3[1], -1e10, term.3[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.4[1], -1e10, term.4[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.5[1], -1e10, term.5[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.6[1], -1e10, term.6[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    rect(term.7[1], -1e10, term.7[2], 1e10, col = alpha("yellow", 0.5), border = NA)
    polygon(c(data_new[[i]][[j]]$date, rev(data_new[[i]][[j]]$date)), 
            c(data_new[[i]][[j]]$lower, rev(data_new[[i]][[j]]$upper)), col = "#D7E5F0", border = NA)
    lines(data_new[[i]][[j]]$date, data_new[[i]][[j]]$upper,  type = "l", col = "#D7E5F0")
    lines(data_new[[i]][[j]]$date, data_new[[i]][[j]]$lower,  type = "l", col = "#D7E5F0")
    lines(data_new[[i]][[j]]$date, data_new[[i]][[j]]$fitted,  type = "l", col = "black", lwd = 1)
    axis.Date(1, at=pretty(range(data_new[[i]][[j]]$date), 10), format='%b %Y', las=2, 
              cex.axis=0.8)
  }
  mtext(region.list[i], side = 3, line = 0, outer = TRUE)
}

