library(ggplot2)
library(GridLMM)
library(foreach)
library(doParallel)
library(Matrix)
library(data.table)
library(dplyr)
library(sommer)


load('sample_size_simulations_setup.workspace')
get_chol_V = function(h2s,Ks) {
  n = nrow(Ks[[1]])
  V = (1-sum(h2s))*diag(1,n)
  for(i in 1:length(h2s)) V = V + h2s[i]*Ks[[i]]
  chol_V = chol_c(V)
  return(chol_V)
}
get_p = function(y,X_cov,X,chol_V) {
  X_list = list(as.matrix(X))
  X_indices = 1:ncol(X_list[[1]])
  p = 2
  r1 = GridLMM_SS_matrix(matrix(y,nc=1),chol_V,X_cov,X_list,X_indices,rep(0,3))
  r1 = get_LL(r1,X_cov,X_list,X_indices,n,1,T,T,F)
  r1$nl10p = -pf(r1$F_hats[3,],1,n - p - 1,lower.tail=F,log.p = TRUE)/log(10)
  return(r1)
}

sample = as.numeric(commandArgs(t=T)[1])
run_ID = as.numeric(commandArgs(t=T)[2])
if(is.na(sample)) sample = 1
if(is.na(run_ID)) run_ID = 1
var_bs = c(0,10^seq(-3,log10(.15),length=15))
nSNPs = 10

try(dir.create(sprintf('Run_%d_%d',sample,run_ID)))
setwd(sprintf('Run_%d_%d',sample,run_ID))
LDAK_path = '../misc_software/ldak5' # set path to LDAK program

# registerDoParallel(my_detectCores())
results = c()
var_b = max(var_bs)
for(sample in 1:3) {
  print(paste('Sample',sample))
for(i in 1:nSNPs) {
  print(paste('SNP',i))
    n = ns[sample]
    K_list = K_lists[[sample]]
    e = t(chol_Rs[[sample]]) %*% rnorm(n)
    e = e/sd(e)
    x = Xs[[sample]][,nSNPs*(run_ID-1)+i]
    xb = x
    xb = xb/sd(xb) * sqrt(var_b/var(e))[1]
    
    e = e * sqrt(1-var_b)
    data = data.frame(y = xb + e,x=x,ID=1:n)
    
    X_cov = X_covs[[sample]]
    
    null_time0 = system.time(optim_null <- get_h2_LDAK(data$y,X_cov,K_list,LDAK_path,maxiter = 0)) # no iterations
    null_time = system.time(optim_null <- get_h2_LDAK(data$y,X_cov,K_list,LDAK_path))
    
    chol_time = system.time(chol_V <- get_chol_V(optim_null,Ks[[sample]]))
    ps = 10^(0:5)
    p_times = sapply(ps,function(p) {
      X = matrix(x,nr=n,nc=p)
      system.time(get_p(data$y,X_cov,X,chol_V))
    })
    
    
    grid_size = nrow(get_h2s_ball(optim_null,0.001))
    
    results = rbind(results,data.frame(n=n,
                                       p = ps,
                                       grid_size = grid_size,
                                       null_time = null_time[3],
                                       null_time0 = null_time0[3],
                                       chol_time = chol_time[3],
                                       p_times[3,])
    )
  }
}
  
  setwd('..')
  saveRDS(results,file = 'Results/Sample_size_v3/timings.rds')

  
sim_results = readRDS('Results/Sample_size_v3/compiled_results.rds')  
registerDoParallel(my_detectCores())
sim_results$n_steps_0.01 = foreach(i = seq_len(nrow(sim_results)),.combine = 'c') %dopar% {
  if(i %% 100 == 0) print(i)
  null_h2s = sim_results[i,grep('null_h2.',colnames(sim_results),fixed=T)]
  names(null_h2s) = sub('null_h2.','',names(null_h2s))
  full_h2s = sim_results[i,grep('full_h2.',colnames(sim_results),fixed=T)]
  names(full_h2s) = sub('full_h2.','',names(full_h2s))
  
  h2s_ball = get_h2s_ball(null_h2s,0.01) # get ball around null_h2s
  tested_h2s = h2s_ball
  step = 1
  while(TRUE) {
    step = step + 1
    current_h2s = tested_h2s[which.min(rowSums(abs(sweep(tested_h2s,2,unlist(full_h2s),'-')))),,drop=FALSE]  # find closest h2 to full_h2s
    new_h2s = get_h2s_to_test(current_h2s,tested_h2s,h2_step = 0.01,F,F)
    if(nrow(new_h2s) == 0) break
    tested_h2s = rbind(tested_h2s,new_h2s)
  }
  nrow(tested_h2s)
}

mean_n_steps = aggregate(n_steps_0.01~n+var_b,sim_results,FUN=mean)
mean_n_steps_optim = subset(mean_n_steps,var_b == 0.15)[,3]
names(mean_n_steps_optim) = c(362,1814,9070)

mean_n_steps_optim = c(362 =  mean_n_steps[43,3], # var_b == 0.1048
                       1814 = mean_n_steps[35,3], # var_b == 0.035838589
                       9070 = mean_n_steps[21,3]  # var_b == 0.005986534
                    )
  
results$LDAK_time2 = pmax(0,results$LDAK_time)
results$Exact_LMM = results$p * results$LDAK_time2
results$null_LMM = results$LDAK_time2 + results$p_times.3... + results$chol_time
results$`Grid_LMM(0.1)` = results$LDAK_time2 + ncol(setup_Grid(1:3,.1))*(results$p_times.3... + results$chol_time)
results$`Grid_LMM_fast(0.01)` = results$LDAK_time2 + results$grid_size*(results$p_times.3... + results$chol_time) + 0*results$chol_time * mean_n_steps_optim[as.character(results$n)]
results_tall = melt(results,id.vars = c('n','p'),measure.vars = c('Exact_LMM','null_LMM','Grid_LMM(0.1)','Grid_LMM_fast(0.01)'))
results_tall_mean = aggregate(value~n+p+variable,results_tall,FUN=mean)
results_tall_mean$minutes = results_tall_mean$value/60

ggplot(results_tall_mean,aes(x=p,y=minutes)) + geom_line(aes(group = variable,color=variable)) + 
  facet_wrap(~n,scales = 'free') + scale_x_log10() + scale_y_log10() + coord_cartesian(ylim = range(results_tall_mean$minutes)) + 
  xlim('number of markers')

