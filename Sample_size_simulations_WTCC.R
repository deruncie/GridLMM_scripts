library(ggplot2)
library(GridLMM)
library(foreach)
library(doParallel)
library(Matrix)
library(data.table)
library(sommer)

# load the WTCC genotypes and parse for autosomal markers
data(mice,package='BGLR')
# mouse_data=readRDS('Data/WTCC/Mouse_data.rds')
X = mice.X

colnames(X) = sapply(colnames(X),function(x) strsplit(x,'_')[[1]][1])
map = fread('Data/WTCC/Mice_map_Plos_T5.csv')
map = subset(map,marker %in% colnames(X) & chromosome != 'X')
X = X[,map$marker]
# X = sweep(X,2,colMeans(X),'-')
n = nrow(X)

# make Kinship matrices
A = A.mat(X-1)
E = E.mat(X-1)
Cage = tcrossprod(model.matrix(~0+factor(mice.pheno$cage)))
rownames(A) = rownames(E) = rownames(Cage) = rownames(X)

LDAK_path = "misc_software/ldak5"
K_list = c(
  A = prep_LDAK_Kinship(A,'K_A',LDAK_path),
  E = prep_LDAK_Kinship(E,'K_E',LDAK_path),
  cage = prep_LDAK_Kinship(Cage,'K_cage',LDAK_path)
)

X_cov = model.matrix(~GENDER,mice.pheno)
h2s_null = get_h2_LDAK(mice.pheno$EndNormalBW,X_cov,K_list,LDAK_path)

error = h2s_null[1,'E'] * E + h2s_null[1,'cage'] * Cage + (1-sum(h2s_null)) * diag(1,n)
chol_error = chol(error)

chol_inv_error = solve(chol_error)

K = t(chol_inv_error) %*% A %*% chol_inv_error
sK = svd(K)
rot_matrix = t(sK$u) %*% t(chol_inv_error)
qtX =  rot_matrix %*% (X-1)
qt1_n = rot_matrix %*% matrix(1,n,1)


# simulation parameters

h2 = h2s_null[1,'A']
var_bs = c(0,10^seq(-3,-1,length=12))


# function for rapidly finding REML solution for mixed model with two diagonal random effects
fit_REML = function(h2,data,chol_I=NULL,X = T,return_result = F) {
  n = nrow(data)
  d = data$d
  qty = data$qty
  qt1 = data$qt1
  qtx = data$qtx
  h2s = c(h2,1-h2)
  diag_R = h2s[1]*d + h2s[2]
  if(is.null(chol_I)) {
    chol_Vi_R = Diagonal(n,sqrt(diag_R))
    chol_Vi_R = as(as(chol_Vi_R,'CsparseMatrix'),'dgCMatrix')
  } else{
    chol_Vi_R = chol_I
    diag(chol_Vi_R) = sqrt(diag_R)
  }
  X_indices = c(1)
  X_cov = matrix(qt1,n,1)
  if(X) {
    p = 1
    X_list = list(matrix(qtx,nc=1))
  } else{
    p = 0
    X_list = list()
  }
  r1 = GridLMM_SS_matrix(matrix(qty,nc=1),chol_Vi_R,X_cov,X_list,X_indices,rep(0,2))
  r1 = get_LL(r1,X_cov,X_list,p,n,X_indices,T,T,F)
  if(return_result){
    r1$nl10p = -pf(r1$F_hats[2],1,n - p - 1,lower.tail=F,log.p = TRUE)/log(10)
    return(r1)
  }
  return(-r1$REML[1])
}

samples = list(
  s1 = 1:(n/5),  # first ~n/4 eigenvalues
  #s2 = floor(seq(1,n,length = 500)), # 500 evenly sampled eigenvalues
  m1 = 1:n,   # original data
  l1 = rep(1:n,each = 5), # 5n samples with original eigenvalues each 5x
  xl1 = rep(1:n,each = 50) # 50n samples with original eigenvalues each 50x
  # l2 = sort(c(1:n,rep(floor(n/2),4*n))), # 5n samples with a median eigenvalue repeated
  # xl1 = rep(1:n,each = 50), # 50n samples with original eigenvalues each 5x
  # xl2 = sort(c(1:n,rep(floor(n/2),49*n))) # 50n samples with a median eigenvalue repeated
)



# run simulations

registerDoParallel(my_detectCores())
results = foreach(i = sample(1:ncol(qtX),1000),.combine='rbind') %dopar% {
  foreach(sample_id = names(samples),.combine='rbind') %do% {
    sample = samples[[sample_id]]
    n = length(sample)
    d = sK$d[sample]
    d = d / sum(d) * n
    qtx = qtX[sample,i]
    qt1 = qt1_n[sample]
    chol_I = Diagonal(n,1)
    chol_I = as(as(chol_I,'CsparseMatrix'),'dgCMatrix')
    # foreach(h2 = h2s,.combine = 'rbind') %do% {
      foreach(var_b = var_bs,.combine = 'rbind') %do% {
        e = sqrt(h2) * sqrt(d) * rnorm(n) + sqrt(1-h2) * rnorm(n)
        e = e / sd(e)
        # qtx = qtx* sample(c(-1,1),n,replace=T)
        qtxb = qtx 
        qtxb = qtxb/sd(qtxb) * sqrt(var_b/var(e)) 
        # e = e * sqrt(1-var_b)
        
        cov_be = cov(qtxb,e)
        coef = -sqrt(1-var_b + cov_be^2)+cov_be
        coef = 1
        data = data.frame(qty = qtxb + e*coef, qt1 = qt1,qtx = qtx,d=d,ID = 1:n)
        
        optim_null = optimize(fit_REML,interval = c(0,.9999),data=data,chol_I = chol_I,X=FALSE,return_result=F)
        optim_full = optimize(fit_REML,interval = c(0,.9999),data=data,chol_I = chol_I,X=TRUE,return_result=F)
        
        res_null = fit_REML(optim_null$minimum,data=data,chol_I = chol_I,X=TRUE,return_result = T)
        res_full = fit_REML(optim_full$minimum,data=data,chol_I = chol_I,X=TRUE,return_result = T)
        
        res_np0.005 = fit_REML(optim_full$minimum+0.005,data=data,chol_I = chol_I,X=TRUE,return_result = T)
        res_nm0.005 = fit_REML(optim_full$minimum-0.005,data=data,chol_I = chol_I,X=TRUE,return_result = T)
        res_Grid0.01 = res_np0.005
        if(!is.na(res_nm0.005$REML + res_np0.005$REML) && res_nm0.005$REML > res_np0.005$REML) res_Grid0.01 = res_nm0.005

        res_np0.05 = fit_REML(optim_full$minimum+0.05,data=data,chol_I = chol_I,X=TRUE,return_result = T)
        res_nm0.05 = fit_REML(optim_full$minimum-0.05,data=data,chol_I = chol_I,X=TRUE,return_result = T)
        res_Grid0.1 = res_np0.05
        if(!is.na(res_nm0.05$REML + res_np0.05$REML) && res_nm0.05$REML > res_np0.05$REML) res_Grid0.1 = res_nm0.05
        
                
        data.frame(sample_id=sample_id,n=n,h2=h2,var_b=var_b,i=i,
                   var_b_act = var(qtxb)/var(data$qty),#conc_x = conc_x,
                   null_h2 = optim_null$minimum,
                   full_h2 = optim_full$minimum,
                   null_REML = res_null$REML,
                   full_REML = res_full$REML,
                   grid0.01_REML = res_Grid0.01$REML,
                   grid0.1_REML = res_Grid0.1$REML,
                   null_F = res_null$F_hats[2],
                   full_F = res_full$F_hats[2],
                   null_beta = res_null$beta_hats[2],
                   full_beta = res_full$beta_hats[2],
                   null_beta2 = var(qtx*res_null$beta_hats[2]),
                   full_beta2 = var(qtx*res_full$beta_hats[2]),
                   null_l10p = res_null$nl10p,
                   full_l10p = res_full$nl10p,
                   grid0.01_l10p = res_Grid0.01$nl10p,
                   grid0.1_l10p = res_Grid0.1$nl10p
        )
      }
    # }
  }
}

write.csv(results,file = 'Results/Sample_size/sample_size_sims_WTCC.csv',row.names=F)

# library(data.table)
# results = fread('Results/Sample_size/sample_size_sims_WTCC.csv',data.table=F)
# # results$sample_id = factor(results$sample_id,levels = c('s1','s2','m1','l1','l2','xl1','xl2'))
# results$sample_id = factor(results$sample_id,levels = c('s1','m1','l1'))
# 
# ggplot(results,aes(x=factor(sample_id),y=null_h2 - h2)) + 
#   geom_boxplot(aes(group = interaction(sample_id,var_b),color=var_b)) +
#   facet_grid(h2~.)
# ggplot(results,aes(x=factor(n),y=full_beta-null_beta)) + geom_boxplot(aes(group = interaction(n,var_b),color=var_b)) +
#   facet_grid(h2~.)
# ggplot(results,aes(x=factor(n),y=null_beta2 - var_b)) + geom_boxplot(aes(group = interaction(n,var_b),color=var_b)) +
#   facet_grid(h2~.)
# 
# # results$Grid_l10p = results$np0.01_l10p
# # results$Grid_l10p[results$np0.01_l10p < results$nm0.01_l10p] = results$nm0.01_l10p
# 
# results_tall = melt(results[,c('n','h2','var_b','null_l10p','full_l10p','grid_l10p')],
#                     id.vars = c('n','h2','var_b'))
# results_tall$model = factor(results_tall$variable,
#                                levels = c('full_l10p','null_l10p','grid_l10p'),
#                                labels = c('Exact-LMM','Null-LMM','Grid-LMM(0.01)'))
# 
# p1 = ggplot(subset(results_tall,h2 == 0.1),aes(x=var_b)) + 
#   geom_hline(yintercept = c(0,8)) +
#   geom_smooth(aes(y=value,group=model,color=model),se=F) + 
#   facet_wrap(~n,scales='free') +  expand_limits(y=0) + 
#   xlab('Percent of variance from X') + ylab(expression(-log[10](p)));p1
# legend <- get_legend(p1)
# p1 = p1 + theme(legend.position = 'none') 
# 
# sub_results_tall = filter(results_tall,!(n == 1814 & var_b > 0.05) & !(n == 9070 & var_b > 0.01))
# p2 = ggplot(subset(sub_results_tall,h2 == 0.1),aes(x=var_b)) + 
#   geom_hline(yintercept = c(0,8)) +
#   geom_smooth(aes(y=value,group=model,color=model),se=F,method='loess') + 
#   facet_wrap(~n,scales='free') + theme(legend.position = 'none') + #ylim(c(0,11)) + 
#   xlab('Percent of variance from X') + ylab(expression(-log[10](p)));p2
# ss_p_values = plot_grid(p1,legend,p2,NULL,nrow=2,rel_widths = c(6,1),labels = c('a','','b',''))
# save_plot('Figures/Sample_size_p_values.pdf',ss_p_values,base_aspect_ratio = 1.2,base_height = 10)
# 
# 
# 
# cols = c('Exact-LMM' = 'Black','Null-LMM' = 'Blue','Grid-LMM(0.01)' = 'Orange')
# p1 = ggplot(subset(results,h2 == 0.1),aes(x=(var_b))) + 
#   # geom_point(aes(color = var_b)) +
#   # geom_point(aes(y=full_l10p),color='Black') +
#   # geom_point(aes(x=var_b/.9,y=null_l10p),color='Blue') +
#   # geom_point(aes(x=var_b/1.1,y=grid_l10p),color='Orange') +
#   geom_smooth(aes(y=full_l10p),se=F,size=2,color='Exact-LMM') + 
#   geom_smooth(aes(y=null_l10p),se=F,size=2,color='Null-LMM') + 
#   geom_smooth(aes(y=grid_l10p),se=F,size=2,color='Grid-LMM(0.01)') + 
#   facet_wrap(~n,scales = 'free') + 
#   geom_hline(yintercept = 8) +
#   xlab('Percent of variance from X') + ylab(expression(-log[10](p))) +
#   scale_color_manual(name='Model',values = cols);p1
# 
# sub_results = filter(results,!(n == 1814 & var_b > 0.05) & !(n == 9070 & var_b > 0.01))
# p2= ggplot(subset(sub_results,h2 == 0.1),aes(x=(var_b))) + 
#   # geom_point(aes(color = var_b)) +
#   # geom_point(aes(y=full_l10p),color='Black') +
#   # geom_point(aes(x=var_b/.9,y=null_l10p),color='Blue') +
#   # geom_point(aes(x=var_b/1.1,y=grid_l10p),color='Orange') +
#   geom_smooth(aes(y=full_l10p),se=F,size=1,color='Black',method='loess') + 
#   geom_smooth(aes(y=null_l10p),se=F,size=1,color='Blue',method='loess') +
#   geom_smooth(aes(y=grid_l10p),se=F,size=1,color='Orange',method='loess') +
#   facet_wrap(~n,scales = 'free_x') + 
#   geom_hline(yintercept = 8) +
#   xlab('Percent of variance from X') + ylab(expression(-log[10](p)))
# 
# p2
# ss_p_values = plot_grid(p1,p2,nrow = 2,labels = c('a','b'),rel_widths = c(1,1));ss_p_values
# save_plot('Figures/Sample_size_p_values.pdf',ss_p_values,base_aspect_ratio = 1.2,base_height = 10)
# 
# 
# ggplot(subset(results,h2 == 0.1 & var_b>-10),aes(x=(var_b))) + 
#   geom_point(aes(x=var_b*0.97,y=full_l10p-null_l10p),color='Red') + 
#   geom_point(aes(x=var_b*1.03,y=full_l10p-grid_l10p),color='Blue') + 
#   # geom_smooth(aes(y=full_l10p-null_l10p),color='Red') + 
#   # geom_smooth(aes(y=full_l10p-grid_l10p),color='Blue') + 
#   facet_wrap(~n,scales = 'free_x')# + ylim(c(0,20)) #+ scale_x_log10() #+ scale_y_log10()
# 
# ggplot(subset(results,h2 == 0.1),aes(x=(var_b),y=full_h2-null_h2)) + 
#   geom_smooth() + 
#   facet_wrap(~n,scales = 'free') 
# 
# 
# ggplot(subset(results,h2 == 0.1),aes(x=null_l10p,y=full_l10p)) + 
#   geom_point(aes(color=var_b)) +
#   facet_grid(n~.) + geom_abline(slope=1,intercept=0)
# 
# ggplot(subset(results,h2 == 0.1),aes(x=null_l10p,y=full_l10p-null_l10p)) + 
#   geom_point(aes(color=var_b)) +
#   facet_grid(n~.) + geom_abline(slope=0,intercept=0)
# 
# ggplot(subset(results,h2 == 0.1),aes(x=null_h2,y=full_l10p-null_l10p)) + 
#   geom_point(aes(color=var_b)) +
#   facet_grid(n~.) + geom_abline(slope=0,intercept=0)
# 
# 
