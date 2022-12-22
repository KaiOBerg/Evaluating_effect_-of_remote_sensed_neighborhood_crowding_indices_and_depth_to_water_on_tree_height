###Note: Bayesian models need a long time for calculation especially when using a lot of sampling iterations and multiple monte-carlo chains.
###Split data set into species if you have multiple computers and run them simultaneously, or be sure to include saving commands on a regular
###basis throughout the script.
###Time for calculation needed approximately 12 - 16 h depending on engine.



#load and attach specific add on packages
library(tidyverse)
library(brms)
library(rstanarm)
library(readr)
library(dplyr)
library(loo)
library(clipr)

#delete working environment
rm(list = ls())


#####1. Get data set ready####

#import used data set for all species
df_model_man <- read_csv("~/Studium_Tue/Bachelorarbeit Kanada/Arbeit/Daten Arbeitscomputer/thesis/files for submission/data set/comp_ind_dtw_mean_seg_1806.csv")
df_model_man <- subset(df_model_man, select = -c(...1))  #delete unnecessary columns
df_model_man$id <- as.factor(df_model_man$id)            #change id to factor
df_model_man$Species <- as.factor(df_model_man$Species)  #change species to factor
df_model_man$Site <- as.factor(df_model_man$Site)        #change site to factor
df_model_man <- df_model_man %>%                         #rename column names
  rename(
    hg_19 = crw__19,
    hg_21 = crw__21
  )

#create data set for each species
WS_man <- df_model_man[df_model_man$Species=="WS", ]   #white spruce
LP_man <- df_model_man[df_model_man$Species=="LP", ]   #lodgepole pine
TA_man <- df_model_man[df_model_man$Species=="TA", ]   #trembling aspen



#####2. Create models####


#####2.1 White spruce####


#model for white spruce using initial height as only predictor
WS_height <- brm(hg_gr ~ log(hg_19) + (log(hg_19) | Site),  #model formula
                 data   = WS_man,     #white spruce data set
                 warmup = 2000,       #discarded warmup iterations
                 iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                 chains = 4,          #define number of monte carlo chains
                 prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                 init  = "random",    #set init to random
                 save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                 control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                 cores  = 4)          #define number of computer cores to calculate model

#model for white spruce only exculuding neighborhood crowding
WS_height_dtw <- brm(hg_gr ~ log(hg_19) * dtw_values + (log(hg_19) + dtw_values | Site),  #model formula
                     data   = WS_man,     #white spruce data set
                     warmup = 2000,       #discarded warmup iterations
                     iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                     chains = 4,          #define number of monte carlo chains
                     prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                     init  = "random",    #set init to random
                     save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                     control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                     cores  = 4)          #define number of computer cores to calculate model

#model for white spruce using mean canopy height 
WS_mch_model_man <- brm(hg_gr ~ log(hg_19) + mch + dtw_values + dtw_values:mch + log(hg_19):dtw_values + log(hg_19):mch + (log(hg_19) + dtw_values + mch | Site),  #model formula
                        data   = WS_man,     #white spruce data set
                        warmup = 2000,       #discarded warmup iterations
                        iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                        chains = 4,          #define number of monte carlo chains
                        prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                        init  = "random",    #set init to random
                        save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                        control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                        cores  = 4)          #define number of computer cores to calculate model

#model for white spruce using canopy cover 
WS_cc_model_man <- brm(hg_gr ~ log(hg_19) + cc + dtw_values + dtw_values:cc + log(hg_19):dtw_values + log(hg_19):cc + (log(hg_19) + dtw_values + cc | Site),  #model formula 
                       data   = WS_man,     #white spruce data set
                       warmup = 2000,       #discarded warmup iterations
                       iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                       chains = 4,          #define number of monte carlo chains
                       prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                       init  = "random",    #set init to random
                       save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                       control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                       cores  = 4)          #define number of computer cores to calculate model

#model for white spruce using canopy above
WS_ca_model_man <- brm(hg_gr ~ log(hg_19) + ca + dtw_values + dtw_values:ca + log(hg_19):dtw_values + log(hg_19):ca + (log(hg_19) + dtw_values + ca | Site),  #model formula
                       data   = WS_man,     #white spruce data set
                       warmup = 2000,       #discarded warmup iterations
                       iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                       chains = 4,          #define number of monte carlo chains
                       prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                       init  = "random",    #set init to random
                       save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                       control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                       cores  = 4)          #define number of computer cores to calculate model

#model for white spruce using canopy angle
WS_can_ang_model_man <- brm(hg_gr ~ log(hg_19) + canopy_angle + dtw_values + dtw_values:canopy_angle + log(hg_19):dtw_values + log(hg_19):canopy_angle + (log(hg_19) + dtw_values + canopy_angle | Site),  #model formula  
                            data   = WS_man,     #white spruce data set
                            warmup = 2000,       #discarded warmup iterations
                            iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                            chains = 4,          #define number of monte carlo chains
                            prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                            init  = "random",    #set init to random
                            save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                            control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                            cores  = 4)          #define number of computer cores to calculate model




#####2.2 Lodgepole pine####

#model for lodgepole pine using initial height as only predictor
LP_height <- brm(hg_gr ~ log(hg_19) + (log(hg_19) | Site),  #model formula
                 data   = LP_man,     #lodgepole pine data set
                 warmup = 2000,       #discarded warmup iterations
                 iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                 chains = 4,          #define number of monte carlo chains
                 prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                 init  = "random",    #set init to random
                 save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                 control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                 cores  = 4)          #define number of computer cores to calculate model

#model for lodgepole pine only exculuding neighborhood crowding
LP_height_dtw <- brm(hg_gr ~ log(hg_19) * dtw_values + (log(hg_19) + dtw_values | Site),  #model formula
                     data   = LP_man,     #lodgepole pine data set
                     warmup = 2000,       #discarded warmup iterations
                     iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                     chains = 4,          #define number of monte carlo chains
                     prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                     init  = "random",    #set init to random
                     save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                     control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                     cores  = 4)          #define number of computer cores to calculate model

#model for lodgepole pine using mean canopy height 
LP_mch_model_man <- brm(hg_gr ~ log(hg_19) + mch + dtw_values + dtw_values:mch + log(hg_19):dtw_values + log(hg_19):mch + (log(hg_19) + dtw_values + mch | Site),  #model formula
                        data   = LP_man,     #lodgepole pine data set
                        warmup = 2000,       #discarded warmup iterations
                        iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                        chains = 4,          #define number of monte carlo chains
                        prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                        init  = "random",    #set init to random
                        save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                        control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                        cores  = 4)          #define number of computer cores to calculate model

#model for lodgepole pine using canopy cover 
LP_cc_model_man <- brm(hg_gr ~ log(hg_19) + cc + dtw_values + dtw_values:cc + log(hg_19):dtw_values + log(hg_19):cc + (log(hg_19) + dtw_values + cc | Site),  #model formula 
                       data   = LP_man,     #lodgepole pine data set
                       warmup = 2000,       #discarded warmup iterations
                       iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                       chains = 4,          #define number of monte carlo chains
                       prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                       init  = "random",    #set init to random
                       save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                       control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                       cores  = 4)          #define number of computer cores to calculate model

#model for lodgepole pine using canopy above
LP_ca_model_man <- brm(hg_gr ~ log(hg_19) + ca + dtw_values + dtw_values:ca + log(hg_19):dtw_values + log(hg_19):ca + (log(hg_19) + dtw_values + ca | Site),  #model formula
                       data   = LP_man,     #lodgepole pine data set
                       warmup = 2000,       #discarded warmup iterations
                       iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                       chains = 4,          #define number of monte carlo chains
                       prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                       init  = "random",    #set init to random
                       save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                       control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                       cores  = 4)          #define number of computer cores to calculate model

#model for lodgepole pine using canopy angle
LP_can_ang_model_man <- brm(hg_gr ~ log(hg_19) + canopy_angle + dtw_values + dtw_values:canopy_angle + log(hg_19):dtw_values + log(hg_19):canopy_angle + (log(hg_19) + dtw_values + canopy_angle | Site),  #model formula  
                            data   = LP_man,     #lodgepole pine data set
                            warmup = 2000,       #discarded warmup iterations
                            iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                            chains = 4,          #define number of monte carlo chains
                            prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                            init  = "random",    #set init to random
                            save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                            control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                            cores  = 4)          #define number of computer cores to calculate model



#####2.3 Trembling aspen####


#model for trembling aspen using initial height as only predictor
TA_height <- brm(hg_gr ~ log(hg_19) + (log(hg_19) | Site),  #model formula
                 data   = TA_man,     #trembling aspen data set
                 warmup = 2000,       #discarded warmup iterations
                 iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                 chains = 4,          #define number of monte carlo chains
                 prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                 init  = "random",    #set init to random
                 save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                 control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                 cores  = 4)          #define number of computer cores to calculate model

#model for trembling aspen only exculuding neighborhood crowding
TA_height_dtw <- brm(hg_gr ~ log(hg_19) * dtw_values + (log(hg_19) + dtw_values | Site),  #model formula
                     data   = TA_man,     #trembling aspen data set
                     warmup = 2000,       #discarded warmup iterations
                     iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                     chains = 4,          #define number of monte carlo chains
                     prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                     init  = "random",    #set init to random
                     save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                     control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                     cores  = 4)          #define number of computer cores to calculate model

#model for trembling aspen using mean canopy height 
TA_mch_model_man <- brm(hg_gr ~ log(hg_19) + mch + dtw_values + dtw_values:mch + log(hg_19):dtw_values + log(hg_19):mch + (log(hg_19) + dtw_values + mch | Site),  #model formula
                        data   = TA_man,     #trembling aspen data set
                        warmup = 2000,       #discarded warmup iterations
                        iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                        chains = 4,          #define number of monte carlo chains
                        prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                        init  = "random",    #set init to random
                        save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calculating r2/loo
                        control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                        cores  = 4)          #define number of computer cores to calculate model

#model for trembling aspen using canopy cover 
TA_cc_model_man <- brm(hg_gr ~ log(hg_19) + cc + dtw_values + dtw_values:cc + log(hg_19):dtw_values + log(hg_19):cc + (log(hg_19) + dtw_values + cc | Site),  #model formula 
                       data   = TA_man,     #trembling aspen data set
                       warmup = 2000,       #discarded warmup iterations
                       iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                       chains = 4,          #define number of monte carlo chains
                       prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                       init  = "random",    #set init to random
                       save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                       control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                       cores  = 4)          #define number of computer cores to calculate model

#model for trembling aspen using canopy above
TA_ca_model_man <- brm(hg_gr ~ log(hg_19) + ca + dtw_values + dtw_values:ca + log(hg_19):dtw_values + log(hg_19):ca + (log(hg_19) + dtw_values + ca | Site),  #model formula
                       data   = TA_man,     #trembling aspen data set
                       warmup = 2000,       #discarded warmup iterations
                       iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                       chains = 4,          #define number of monte carlo chains
                       prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                       init  = "random",    #set init to random
                       save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                       control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                       cores  = 4)          #define number of computer cores to calculate model

#model for trembling aspen using canopy angle
TA_can_ang_model_man <- brm(hg_gr ~ log(hg_19) + canopy_angle + dtw_values + dtw_values:canopy_angle + log(hg_19):dtw_values + log(hg_19):canopy_angle + (log(hg_19) + dtw_values + canopy_angle | Site),  #model formula  
                            data   = TA_man,     #trembling aspen data set
                            warmup = 2000,       #discarded warmup iterations
                            iter   = 6000,       #total sampling iterations (increase when warning about Bulk_ess, etc.)
                            chains = 4,          #define number of monte carlo chains
                            prior = prior(normal(0, 1), class = "b" ),   #add uninformative prior
                            init  = "random",    #set init to random
                            save_pars = save_pars(all = TRUE),           #save pars to avoid warnings when calcualting r2/loo
                            control = list(max_treedepth = 15, adapt_delta = 0.99),   #increase maximum tree depth or adapt_delta when model not converges
                            cores  = 4)          #define number of computer cores to calculate model



#####3. Evaluate model performance####


#####3.1 Create model summaries####

#####3.1.1 White spruce####
WS_height_sum <- summary(WS_height)
WS_height_dtw_sum <- summary(WS_height_dtw)
WS_mch_man_sum <- summary(WS_mch_model_man)
WS_cc_man_sum <- summary(WS_cc_model_man)
WS_ca_man_sum <- summary(WS_ca_model_man)
WS_can_ang_man_sum <- summary(WS_can_ang_model_man)

#####3.1.2 Lodgepole pine####
LP_height_sum <- summary(LP_height)
LP_height_dtw_sum <- summary(LP_height_dtw)
LP_mch_man_sum <- summary(LP_mch_model_man)
LP_cc_man_sum <- summary(LP_cc_model_man)
LP_ca_man_sum <- summary(LP_ca_model_man)
LP_can_ang_man_sum <- summary(LP_can_ang_model_man)

#####3.1.3 Trembling aspen####
TA_height_sum <- summary(TA_height)
TA_height_dtw_sum <- summary(TA_height_dtw)
TA_mch_man_sum <- summary(TA_mch_model_man)
TA_cc_man_sum <- summary(TA_cc_model_man)
TA_ca_man_sum <- summary(TA_ca_model_man)
TA_can_ang_man_sum <- summary(TA_can_ang_model_man)



#####3.2 Calculate Bayes R2####

#define function to calculate bayes r2 
bayes_R2_res <- function(fit) {
  y <- rstanarm::get_y(fit)
  ypred <- rstanarm::posterior_epred(fit)
  if (family(fit)$family == "binomial" && NCOL(y) == 2) {
    trials <- rowSums(y)
    y <- y[, 1]
    ypred <- ypred %*% diag(trials)
  }
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- apply(ypred, 1, var)
  var_e <- apply(e, 1, var)
  var_ypred / (var_ypred + var_e)
}



#####3.2.1 White spruce####
#calculate r2 for all models
WS_height_man <- (round(median(bayes_R2_res(WS_height)), 4))
WS_height_dtw_man <- (round(median(bayes_R2_res(WS_height_dtw)), 4))
WS_mch_man <- (round(median(bayes_R2_res(WS_mch_model_man)), 4))
WS_cc_man <- (round(median(bayes_R2_res(WS_cc_model_man)), 4))
WS_ca_man <- (round(median(bayes_R2_res(WS_ca_model_man)), 4))
WS_can_ang_man <- (round(median(bayes_R2_res(WS_can_ang_model_man)), 4))
#show bayes r2 values for white spruce in one data frame
models_r2_WS <- data.frame(WS_height_man,WS_height_dtw_man,WS_mch_man,WS_cc_man,WS_ca_man,
                           WS_can_ang_man)

#####3.2.2 Lodgepole pine####
#calculate r2 for all models
LP_height_man <- (round(median(bayes_R2_res(LP_height)), 4))
LP_height_dtw_man <- (round(median(bayes_R2_res(LP_height_dtw)), 4))
LP_mch_man <- (round(median(bayes_R2_res(LP_mch_model_man)), 4))
LP_cc_man <- (round(median(bayes_R2_res(LP_cc_model_man)), 4))
LP_ca_man <- (round(median(bayes_R2_res(LP_ca_model_man)), 4))
LP_can_ang_man <- (round(median(bayes_R2_res(LP_can_ang_model_man)), 4))
#show bayes r2 values for lodgepole pine in one data frame
models_r2 <- data.frame(LP_height_man,LP_height_dtw_man,LP_mch_man,LP_cc_man,LP_ca_man,
                        LP_can_ang_man)

#####3.2.3 Trembling aspen####
#calculate r2 for all models
TA_height_man <- (round(median(bayes_R2_res(TA_height)), 4))
TA_height_dtw_man <- (round(median(bayes_R2_res(TA_height_dtw)), 4))
TA_mch_man <- (round(median(bayes_R2_res(TA_mch_model_man)), 4))
TA_cc_man <- (round(median(bayes_R2_res(TA_cc_model_man)), 4))
TA_ca_man <- (round(median(bayes_R2_res(TA_ca_model_man)), 4))
TA_can_ang_man <- (round(median(bayes_R2_res(TA_can_ang_model_man)), 4))
#show bayes r2 values for trembling aspen in one data frame
models_r2 <- data.frame(TA_height_man,TA_height_dtw_man,TA_mch_man,TA_cc_man,TA_ca_man,
                        TA_can_ang_man)



#####3.2 Calculate leave-one-out cross-validation for model selection####


#####3.3.1 White spruce####
#calculate leave one out cross validation 
WS_height_loo <- loo::loo(WS_height, moment_match = TRUE)
WS_height_dtw_loo <- loo::loo(WS_height_dtw, moment_match = TRUE)
WS_mch_man_loo <- loo::loo(WS_mch_model_man, moment_match = TRUE)
WS_cc_man_loo <- loo::loo(WS_cc_model_man, moment_match = TRUE)
WS_ca_man_loo <- loo::loo(WS_ca_model_man, moment_match = TRUE)
WS_can_ang_man_loo <- loo::loo(WS_can_ang_model_man, moment_match = TRUE)
#compare results of leave one out cross validation for all white spruce models
WS_loo_comp_man <- as.data.frame(loo::loo_compare(WS_height_loo,WS_height_dtw_loo,WS_mch_man_loo,WS_cc_man_loo,WS_ca_man_loo,
                                                  WS_can_ang_man_loo))


#####3.3.2 Lodgepole pine####
#calculate leave one out cross validation 
LP_height_loo <- loo::loo(LP_height, moment_match = TRUE)
LP_height_dtw_loo <- loo::loo(LP_height_dtw, moment_match = TRUE)
LP_mch_man_loo <- loo::loo(LP_mch_model_man, moment_match = TRUE)
LP_cc_man_loo <- loo::loo(LP_cc_model_man, moment_match = TRUE)
LP_ca_man_loo <- loo::loo(LP_ca_model_man, moment_match = TRUE)
LP_can_ang_man_loo <- loo::loo(LP_can_ang_model_man, moment_match = TRUE)
#compare results of leave one out cross validation for all lodgepole pine models
LP_loo_comp_man <- as.data.frame(loo::loo_compare(LP_height_loo,LP_height_dtw_loo,LP_mch_man_loo,LP_cc_man_loo,LP_ca_man_loo,
                                                  LP_can_ang_man_loo))


#####3.3.3 Trembling aspen####
#calculate leave one out cross validation 
TA_height_loo <- loo::loo(TA_height, moment_match = TRUE)
TA_height_dtw_loo <- loo::loo(TA_height_dtw, moment_match = TRUE)
TA_mch_man_loo <- loo::loo(TA_mch_model_man, moment_match = TRUE)
TA_cc_man_loo <- loo::loo(TA_cc_model_man, moment_match = TRUE)
TA_ca_man_loo <- loo::loo(TA_ca_model_man, moment_match = TRUE)
TA_can_ang_man_loo <- loo::loo(TA_can_ang_model_man, moment_match = TRUE)
#compare results of leave one out cross validation for all trembling aspen models
TA_loo_comp_man <- as.data.frame(loo::loo_compare(TA_height_loo,TA_height_dtw_loo,TA_mch_man_loo,TA_cc_man_loo,TA_ca_man_loo,
                                                  TA_can_ang_man_loo))

#set working directory
setwd("H:/r_workspace_model_comparison/")
#save r workspace (Important!)
save.image(file = "model_comaprison_final.RData")


####ceck model assumptions####

WS_fit_mch <- as.data.frame(fitted(WS_mch_model_man,WS_man))
plot(WS_fit_mch$Estimate,WS_man$hg_gr)
abline(0,1)

LP_fit_cc <- as.data.frame(fitted(LP_cc_model_man,LP_man))
plot(LP_fit_cc$Estimate,LP_man$hg_gr)
abline(0,1)

TA_fit_mch <- as.data.frame(fitted(TA_mch_model_man,TA_man))
plot(TA_fit_mch$Estimate,TA_man$hg_gr)
abline(0,1)



