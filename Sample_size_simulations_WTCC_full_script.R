
library(GridLMM)
library(foreach)
library(doParallel)
library(dplyr)
library(reshape2)
library(cowplot)
library(Matrix)
library(data.table)
library(sommer)


#-----------------------------------------------------------------#
# prep data
#-----------------------------------------------------------------#

LDAK_path = 'misc_software/ldak5' # set path to LDAK program

data(mice,package='BGLR')

trait_data = droplevels(data.frame(ID = mice.pheno$SUBJECT.NAME,GENDER = mice.pheno$GENDER, cage = mice.pheno$cage, y = mice.pheno$EndNormalBW))
trait_data$ID2 = trait_data$ID3 = trait_data$ID

trait_X = mice.X

# download map (Table S5, convert to csv): https://journals.plos.org/plosbiology/article/file?type=supplementary&id=info:doi/10.1371/journal.pbio.0040395.st005
map = fread('Data/WTCC/Mice_map_Plos_T5.csv',data.table = F)
colnames(trait_X) = substr(colnames(trait_X),1,nchar(colnames(trait_X))-2)
map = subset(map,chromosome != 'X' & marker %in% colnames(trait_X))
trait_X = trait_X[,match(subset(map,chromosome != 'X')$marker,colnames(trait_X))]

trait_X = trait_X[,apply(trait_X,2,var)>0] # remove non-variable markers

trait_A = A.mat(trait_X-1)

trait_E = E.mat(trait_X-1)
cage_K = tcrossprod(model.matrix(~0+factor(cage),trait_data))
rownames(cage_K) = colnames(cage_K) = rownames(trait_X)

K_list = c(
  A = prep_LDAK_Kinship(trait_A,'K_A',LDAK_path),
  E = prep_LDAK_Kinship(trait_E,'K_E',LDAK_path)
  ,cage = prep_LDAK_Kinship(cage_K,'K_cage',LDAK_path)
)

X_cov = model.matrix(~GENDER,trait_data)
base_h2s = get_h2_LDAK(trait_data$y,X_cov,K_list[1:3],LDAK_path)



#-----------------------------------------------------------------#
# Set up simulated populations
#-----------------------------------------------------------------#

trait_X = trait_X[,sample(1:ncol(trait_X),300)] # select 300 random markers

n = nrow(trait_A)

ns = c(362,n,9070)

Ks = list(
  small = list(
    A = trait_A[1:ns[1],1:ns[1]],
    E = trait_E[1:ns[1],1:ns[1]]
    ,cage = cage_K[1:ns[1],1:ns[1]]
  ),
  medium = list(
    A = trait_A,
    E = trait_E
    ,cage = cage_K
  ),
  large = list(
    A = as.matrix(do.call(bdiag,lapply(1:5,function(x) trait_A))),
    E = as.matrix(do.call(bdiag,lapply(1:5,function(x) trait_E)))
    ,cage = as.matrix(do.call(bdiag,lapply(1:5,function(x) cage_K)))
  )
)

for(i in seq_along(Ks)) {
  Ks[[i]] = lapply(Ks[[i]],function(x) x/mean(diag(x)))
}

for(i in seq_along(Ks[[3]])) rownames(Ks[[3]][[i]]) = 1:nrow(Ks[[3]][[i]])

X_covs = list(
  X_cov[1:ns[1],,drop=FALSE],
  X_cov,
  as.matrix(do.call(rbind,lapply(1:5,function(x) X_cov)))
)

chol_Rs = list(
  chol(base_h2s[1]*Ks[[1]]$A + base_h2s[2]*Ks[[1]]$E + base_h2s[3]*Ks[[1]]$cage + (1-sum(base_h2s))*diag(1,ns[1])),
  chol(base_h2s[1]*Ks[[2]]$A + base_h2s[2]*Ks[[2]]$E + base_h2s[3]*Ks[[2]]$cage + (1-sum(base_h2s))*diag(1,ns[2]))
)
chol_Rs[[3]] = as.matrix(do.call(bdiag,lapply(1:5,function(x) chol_Rs[[2]])))

Xs = list(
  trait_X[1:ns[1],],
  trait_X,
  do.call(rbind,lapply(1:5,function(x) trait_X))
)





#-----------------------------------------------------------------#
# Run simulations
#-----------------------------------------------------------------#
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


# setup up LDAK
K_lists = lapply(seq_along(ns),function(i) {
  K_list = c(
    A = prep_LDAK_Kinship(Ks[[i]]$A,sprintf('K_A_%d',i),LDAK_path),
    E = prep_LDAK_Kinship(Ks[[i]]$E,sprintf('K_E_%d',i),LDAK_path)
    ,cage = prep_LDAK_Kinship(Ks[[i]]$cage,sprintf('K_cage_%d',i),LDAK_path)
  )
})

nSNPs = ncol(Xs[[1]])

# registerDoParallel(my_detectCores())
results = foreach(sample = 1:3,.combine = 'rbind') %do% 
      foreach(i = 1:nSNPs,.combine='rbind') %do% {
      foreach(var_b = var_bs,.combine = 'rbind') %do% {
        n = ns[sample]
        K_list = K_lists[[sample]]
        e = t(chol_Rs[[sample]]) %*% rnorm(n)
        e = e/sd(e)
        x = Xs[[sample]][,i]
        xb = x
        xb = xb/sd(xb) * sqrt(var_b/var(e))[1]
        
        e = e * sqrt(1-var_b)
        data = data.frame(y = xb + e,x=x,ID=1:n)
        
        X_cov = X_covs[[sample]]
        
        null_time = system.time(optim_null <- get_h2_LDAK(data$y,X_cov,K_list,LDAK_path))  # also record the time for this calculation
        optim_full = get_h2_LDAK(data$y,cbind(X_cov,data$x),K_list,LDAK_path)
        
        # measure times for GridLMM / null-LMM parts
        chol_time = system.time(chol_V <- get_chol_V(optim_null,Ks[[sample]])) 
        ps = 10^5
        p_times = sapply(ps,function(p) {
          X = matrix(x,nr=n,nc=p)
          system.time(get_p(data$y,X_cov,X,chol_V))
        })
        
        grid_size = nrow(get_h2s_ball(optim_null,0.001))
        
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
                   var_b_act = var(xb)/var(data$y),
                   null_time = null_time,
                   chol_time = chol_time,
                   p_time = p_times[1],
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

saveRDS(results,file = 'Results/compiled_results.rds')



