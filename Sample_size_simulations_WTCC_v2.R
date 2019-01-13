library(ggplot2)
library(GridLMM)
library(foreach)
library(doParallel)
library(Matrix)
library(data.table)
library(dplyr)
library(sommer)


# # prep data
# 
# LDAK_path = 'misc_software/ldak5' # set path to LDAK program
# 
# data(mice,package='BGLR')
# 
# trait_data = droplevels(data.frame(ID = mice.pheno$SUBJECT.NAME,GENDER = mice.pheno$GENDER, cage = mice.pheno$cage, y = mice.pheno$EndNormalBW))
# trait_data$ID2 = trait_data$ID3 = trait_data$ID
# 
# trait_X = mice.X
# 
# # download map (Table S5, convert to csv): https://journals.plos.org/plosbiology/article/file?type=supplementary&id=info:doi/10.1371/journal.pbio.0040395.st005
# map = fread('Data/WTCC/Mice_map_Plos_T5.csv',data.table = F)
# colnames(trait_X) = substr(colnames(trait_X),1,nchar(colnames(trait_X))-2)
# map = subset(map,chromosome != 'X' & marker %in% colnames(trait_X))
# trait_X = trait_X[,match(subset(map,chromosome != 'X')$marker,colnames(trait_X))]
# 
# trait_X = trait_X[,apply(trait_X,2,var)>0] # remove non-variable markers
# trait_A = A.mat(trait_X-1)
# 
# trait_E = E.mat(trait_X-1)
# cage_K = tcrossprod(model.matrix(~0+factor(cage),trait_data))
# rownames(cage_K) = colnames(cage_K) = rownames(trait_X)
# 
# K_list = c(
#   A = prep_LDAK_Kinship(trait_A,'K_A',LDAK_path),
#   E = prep_LDAK_Kinship(trait_E,'K_E',LDAK_path)
#   # cage = prep_LDAK_Kinship(cage_K,'K_cage',LDAK_path)
# )
# 
# X_cov = model.matrix(~GENDER,trait_data)
# base_h2s = get_h2_LDAK(trait_data$y,X_cov,K_list[1:2],LDAK_path)
# 
# 
# n = nrow(trait_A)
# 
# ns = c(362,n,9070)
# 
# Ks = list(
#   small = list(
#     A = trait_A[1:ns[1],1:ns[1]],
#     E = trait_E[1:ns[1],1:ns[1]]
#     # ,cage = cage_K[1:ns[1],1:ns[1]]
#   ),
#   medium = list(
#     A = trait_A,
#     E = trait_E
#     # ,cage = cage_K
#   ),
#   large = list(
#     A = as.matrix(do.call(bdiag,lapply(1:5,function(x) trait_A))),
#     E = as.matrix(do.call(bdiag,lapply(1:5,function(x) trait_E)))
#     # ,cage = as.matrix(do.call(bdiag,lapply(1:5,function(x) cage_K)))
#   )
# )
# 
# for(i in seq_along(Ks)) {
#   Ks[[i]] = lapply(Ks[[i]],function(x) x/mean(diag(x)))
# }
# 
# for(i in seq_along(Ks[[3]])) rownames(Ks[[3]][[i]]) = 1:nrow(Ks[[3]][[i]])
# 
# X_covs = list(
#   X_cov[1:ns[1],,drop=FALSE],
#   X_cov,
#   as.matrix(do.call(rbind,lapply(1:5,function(x) X_cov)))
# )
# 
# chol_Rs = list(
#   chol(base_h2s[1]*Ks[[1]]$A + base_h2s[2]*Ks[[1]]$E  + (1-sum(base_h2s))*diag(1,ns[1])), #+ base_h2s[3]*Ks[[1]]$cage
#   chol(base_h2s[1]*Ks[[2]]$A + base_h2s[2]*Ks[[2]]$E  + (1-sum(base_h2s))*diag(1,ns[2])) #+ base_h2s[3]*Ks[[1]]$cage
# )
# chol_Rs[[3]] = as.matrix(do.call(bdiag,lapply(1:5,function(x) chol_Rs[[2]])))
# 
# Xs = list(
#   trait_X[1:ns[1],],
#   trait_X,
#   do.call(rbind,lapply(1:5,function(x) trait_X))
# )
# 
# try(dir.create('Run_1_1'))
# setwd('Run_1_1')
# LDAK_path = '../misc_software/ldak5' # set path to LDAK program
# K_lists = lapply(seq_along(ns),function(i) {
#   K_list = c(
#     A = prep_LDAK_Kinship(Ks[[i]]$A,sprintf('../K_A_%d',i),LDAK_path),
#     E = prep_LDAK_Kinship(Ks[[i]]$E,sprintf('../K_E_%d',i),LDAK_path)
#     # ,cage = prep_LDAK_Kinship(Ks[[i]]$cage,sprintf('../K_cage_%d',i),LDAK_path)
#   )
# })
# setwd('..')
# 
# 
# save.image(file = 'sample_size_simulations_setup.workspace')

print(rnorm(3))

load('sample_size_simulations_setup.workspace')
get_p = function(y,X_cov,X,h2s,Ks) {
  n = length(y)
  V = (1-sum(h2s))*diag(1,n)
  for(i in 1:length(h2s)) V = V + h2s[i]*Ks[[i]]
  chol_V = chol_c(V)
  X_list = list(matrix(X,nc=1))
  X_indices = 1
  p = 2
  r1 = GridLMM_SS_matrix(matrix(y,nc=1),chol_V,X_cov,X_list,X_indices,rep(0,3))
  r1 = get_LL(r1,X_cov,X_list,X_indices,n,1,T,T,F)
  r1$nl10p = -pf(r1$F_hats[3],1,n - p - 1,lower.tail=F,log.p = TRUE)/log(10)
  return(r1)
}

sample = as.numeric(commandArgs(t=T)[1])
run_ID = as.numeric(commandArgs(t=T)[2])
if(is.na(sample)) sample = 1
if(is.na(run_ID)) run_ID = 1
var_bs = c(0,10^seq(-3,log10(.15),length=15))
nSNPs = 1

try(dir.create(sprintf('Run_%d_%d',sample,run_ID)))
setwd(sprintf('Run_%d_%d',sample,run_ID))
LDAK_path = '../misc_software/ldak5' # set path to LDAK program

# run simulations 

# registerDoParallel(my_detectCores())
results = foreach(i = 1:nSNPs,.combine='rbind') %do% {
      foreach(var_b = var_bs,.combine = 'rbind') %do% {
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
        
        optim_null = get_h2_LDAK(data$y,X_cov,K_list,LDAK_path)
        optim_full = get_h2_LDAK(data$y,cbind(X_cov,data$x),K_list,LDAK_path)
        
        res_null = get_p(data$y,X_cov,x,optim_null,Ks[[sample]])
        res_full = get_p(data$y,X_cov,x,optim_full,Ks[[sample]])
        
        grid = do.call(expand.grid,lapply(seq_along(optim_full),function(x) c(-1,1))) * 0.005
        h2s_ball_0.01 = sweep(grid,2,optim_full,'+')
        res_ball_0.01 = apply(h2s_ball_0.01,1,function(h2s) get_p(data$y,X_cov,x,h2s,Ks[[sample]]))
        res_Grid0.01 = res_ball_0.01[[order(sapply(res_ball_0.01,function(x) x$REML),decreasing = T)[1]]]
        
        grid = do.call(expand.grid,lapply(seq_along(optim_full),function(x) c(-1,1))) * 0.05
        h2s_ball_0.1 = sweep(grid,2,optim_full,'+')
        res_ball_0.1 = apply(h2s_ball_0.1,1,function(h2s) get_p(data$y,X_cov,x,h2s,Ks[[sample]]))
        res_Grid0.1 = res_ball_0.1[[order(sapply(res_ball_0.1,function(x) x$REML),decreasing = T)[1]]]
        
        data.frame(sample_id=sample,n=n,var_b=var_b,i=i,
                   var_b_act = var(xb)/var(data$y),#conc_x = conc_x,
                   null_h2 = optim_null,
                   full_h2 = optim_full,
                   null_REML = res_null$REML,
                   full_REML = res_full$REML,
                   grid0.01_REML = res_Grid0.01$REML,
                   grid0.1_REML = res_Grid0.1$REML,
                   null_F = res_null$F_hats[3],
                   full_F = res_full$F_hats[3],
                   null_beta = res_null$beta_hats[3],
                   full_beta = res_full$beta_hats[3],
                   # null_beta2 = var(qtx*res_null$beta_hats[2]),
                   # full_beta2 = var(qtx*res_full$beta_hats[2]),
                   null_l10p = res_null$nl10p,
                   full_l10p = res_full$nl10p,
                   grid0.01_l10p = res_Grid0.01$nl10p,
                   grid0.1_l10p = res_Grid0.1$nl10p
        )
      }
}
setwd('..')

write.csv(results,file = sprintf('Results/Sample_size_v2/sample_size_sims_WTCC_%d_%d.csv',sample,run_ID),row.names=F)

# library(foreach)
# results = foreach(file = list.files(path='Results/Sample_size_v2',pattern = 'csv',full.names = T),.combine='rbind') %do% {
#   read.csv(file)
# }
# saveRDS(results,file = 'Results/Sample_size_v2/compiled_results.rds')
# 
# 
# results = readRDS('compiled_results.rds')
# results$sample_id = factor(results$sample_id,levels = c(1,2,3),labels = c('s1','m1','l1'))
# results = filter(results,sample_id != 'xl1')
# 
# library(reshape2)
# library(cowplot)
# library(dplyr)
# results_tall = melt(results[,c('n','var_b','null_l10p','full_l10p','grid0.1_l10p','grid0.01_l10p')],
#                     id.vars = c('n','var_b'))
# results_tall$model = factor(results_tall$variable,
#                             levels = c('full_l10p','null_l10p','grid0.1_l10p','grid0.01_l10p'),
#                             labels = c('Exact-LMM','null-LMM','Grid-LMM(0.1)','Grid-LMM(0.01)'))
# 
# p1 = ggplot(results_tall,aes(x=var_b)) +
#   geom_hline(yintercept = c(0,8)) +
#   geom_smooth(aes(y=value,group=model,color=model),se=F) +
#   facet_wrap(~n,scales='free') +  expand_limits(y=0) +
#   xlab('Percent of variance from X') + ylab(expression(-log[10](p)));p1
# legend <- get_legend(p1)
# p1 = p1 + theme(legend.position = 'none')
# 
# sub_results_tall = filter(results_tall,!(n == 1814 & var_b > 0.05) & !(n == 9070 & var_b > 0.01) & !(n == 90700 & var_b > 0.001))
# # sub_results_tall_mean = aggregate(value ~ n+var_b+variable+model,sub_results_tall,FUN=mean)
# # p2 = ggplot(sub_results_tall_mean,aes(x=var_b)) +
# #   geom_hline(yintercept = c(0,8)) +
# #   geom_point(aes(y=value,color = model)) + 
# #   # stat_smooth(aes(y=value,group=model,color=model),se=F,method='loess',geom='point') +
# #   facet_wrap(~n,scales='free') + theme(legend.position = 'none') + coord_cartesian(ylim=c(0, 10)) +
# #   xlab('Percent of variance from X') + ylab(expression(-log[10](p)));p2
# p2 = ggplot(sub_results_tall,aes(x=var_b)) +
#   geom_hline(yintercept = c(0,8)) +
#   geom_smooth(aes(y=value,group=model,color=model),se=F,method='loess') +
#   facet_wrap(~n,scales='free') + theme(legend.position = 'none') + coord_cartesian(ylim=c(0, 10)) +
#   xlab('Percent of variance from X') + ylab(expression(-log[10](p)));p2
# ss_p_values = plot_grid(p1,legend,p2,NULL,nrow=2,rel_widths = c(6,1),labels = c('a','','b',''))
# save_plot('Figures/Sample_size_p_values.pdf',ss_p_values,base_aspect_ratio = 1.2,base_height = 10)
# 
# comp_results_tall = melt(results[,c('n','var_b','null_l10p','full_l10p','grid0.1_l10p','grid0.01_l10p')],
#                          id.vars = c('n','var_b','full_l10p'))
# comp_results_tall$model = factor(comp_results_tall$variable,
#                                  levels = c('full_l10p','null_l10p','grid0.1_l10p','grid0.01_l10p'),
#                                  labels = c('Exact-LMM','null-LMM','Grid-LMM(0.1)','Grid-LMM(0.01)'))
# p3 = ggplot(comp_results_tall) +
#   geom_point(aes(x=jitter(var_b,factor = 50),y=abs(full_l10p-value),color = model),size=.5) +
#   geom_smooth(aes(x=var_b,y=abs(full_l10p-value),color = model),se=F,size=1) +
#   geom_hline(yintercept = 0) +
#   scale_color_discrete(drop=FALSE) +
#   theme(legend.position = 'none') +
#   xlab('Percent of variance from X') +
#   ylim(c(0,max(abs(comp_results_tall$full_l10p - comp_results_tall$value)))) +
#   ylab(expression(atop(paste('abs(',Delta,-log[10],'p)'),'relative to Exact-LMM'))) +
#   facet_wrap(~n,scales = 'free');p3
# 
# ss_plots = plot_grid(p1,NULL,p2,legend,p3,NULL,nrow=3,rel_widths = c(6,1),labels = c('a','','b','','c',''))
# save_plot('Figures/Sample_size_plots.pdf',ss_plots,base_aspect_ratio = 1.6,base_height = 8)
# 
