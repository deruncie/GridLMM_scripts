# This set of power simulations will be different from Zhou and Stephens 2012
# Instead of adding to the existing phenotype, I will generate de-novo phenotypes with different relative weightings of the two variance components
# In each simulation, I will then add a random SNP effect of known proprotion of variance, and then evaluate the power with the following methods
# 1) LMM - 1 RE
# 2) LMM - 2 RE
# 5) EMMAX: 1 RE
# 6) EMMAX: 2 RE
# 7) Grid-LMM: step = 0.1
# 8) Grid-LMM: step = 0.01

run = commandArgs(t=T)[1]
if(is.na(run)) run = 1

# use the sommer package for calculations (only REML)
library(lme4qtl)
library(lme4)
library(foreach)
library(doParallel)
library(GxEMMA)
library(reshape2)
library(sommer)

nearest_h2 = function(h2,all_h2s){
  dist = rowSums(sweep(all_h2s,2,h2,'-')^2)
  all_h2s[order(dist)[1],]
}

h2_divisions = 10
all_h2s = as.matrix(expand.grid(seq(0,1,length=h2_divisions+1),seq(0,1,length=h2_divisions+1)))#+1/(2*h2_divisions)
all_h2s = all_h2s[rowSums(all_h2s)<1,]
all_h2s10 = all_h2s

h2_divisions = 20
all_h2s = as.matrix(expand.grid(seq(0,1,length=h2_divisions+1),seq(0,1,length=h2_divisions+1)))#+1/(2*h2_divisions)
all_h2s = all_h2s[rowSums(all_h2s)<1,]
all_h2s20 = all_h2s

h2_divisions = 100
all_h2s = as.matrix(expand.grid(seq(0,1,length=h2_divisions+1),seq(0,1,length=h2_divisions+1)))#+1/(2*h2_divisions)
all_h2s = all_h2s[rowSums(all_h2s)<1,]
all_h2s100 = all_h2s


# Load the dataset

dataset = readRDS('Data/dataset.rds')

phen = dataset$phenotypes
map = dataset$map
X = dataset$X
rm(dataset)
gc()

trait_1 = '5_FT10'
trait_2 = '45_8W GH FT'
# trait_2 = '57_FT Field'

phen = na.omit(phen[,c('ecotype_id',trait_1,trait_2)])

# standardize traits
phen[,-1] = apply(phen[,-1],2,function(x) (x-mean(x))/sd(x))

phen_tall = melt(phen,id.vars = 'ecotype_id')
phen_tall$ecotype_id2 = phen_tall$ecotype_id
phen_tall$ID = as.character(1:nrow(phen_tall))
phen$plasticity = phen[,3]-phen[,2]
X = X[phen$ecotype_id,]
X = X[,apply(X,2,var)>0]

n_full = nrow(phen_tall)
n_strain = nrow(phen)

K = A.mat(X-1)


options(contrasts = c('contr.sum','contr.poly'))
Z_eco = model.matrix(~0+ecotype_id,phen_tall);colnames(Z_eco) = sub('ecotype_id','',colnames(Z_eco))
Z_trt = model.matrix(~1+variable,phen_tall)[,2]
Z_plasticity = model.matrix(~0+ecotype_id,phen)

K_g = Z_eco %*% K[colnames(Z_eco),colnames(Z_eco)] %*% t(Z_eco)
K_gxe = (Z_trt*Z_eco) %*% K[colnames(Z_eco),colnames(Z_eco)] %*% t(Z_trt*Z_eco)

chol_K = chol(K[colnames(Z_eco),colnames(Z_eco)])

ZK_g = Z_eco %*% t(chol_K)
ZK_gxe = (Z_trt*ZK_g)
ZK_gxe_plasticity = Z_plasticity %*% t(chol_K)


sK = svd(K)
U = sK$u
Ut = t(U)
d = diag(sK$d)
rownames(d) = colnames(d) = rownames(X)

sK_g = svd(K_g)
d_g = diag(sK_g$d)
sK_gxe = svd(K_gxe)
d_gxe = diag(sK_gxe$d)
rownames(d_g) = colnames(d_g) = rownames(d_gxe) = colnames(d_gxe) = phen_tall$ID

h2s = as.matrix(expand.grid(G = c(0,0.4,.8), GxE = c(0,0.4,.8)))
h2s = h2s[rowSums(h2s) < 1,]

prop_vars = c(0,0.01,0.025,0.05,0.1,0.15,0.2)
nReps = 4
colnames(X) = 1:ncol(X)
X = X[,sample(1:ncol(X),nReps)]


get_h2_lmer = function(model){
  vars = as.data.frame(VarCorr(model))$vcov
  vars[-length(vars)]/sum(vars)
}

get_h2_sommer = function(model){
  vars = unlist(model$var.comp)
  vars[-length(vars)]/sum(vars)
}
get_F = function(h2s,Vs,y,x,X_cov){
  n = nrow(Vs[[1]])
  V = diag(1-sum(h2s),n)
  for(i in 1:length(h2s)) V = V + h2s[i]*Vs[[i]]
  chol_V = chol(V)
  V_log_det = 2*sum(log(diag(chol_V)))
  calc_LL(matrix(y,dimnames = list(c(),1)),X_cov = matrix(X_cov,n),X_list = lapply(1:ncol(x),function(i) x[,i,drop=FALSE]),h2s = 1,chol_Vi = chol_V,V_log_det = V_log_det,
          inv_prior_X = rep(0,4),downdate_Xs = NULL,n_SNPs_downdated_RRM = NULL,REML = T,BF = F)
}

registerDoParallel(my_detectCores())
results = foreach(h2 = iter(h2s,by='row'),.combine = 'rbind') %do% {
  foreach(prop_var = prop_vars,.combine = 'rbind') %do% {
    print(c(h2,prop_var=prop_var))
    foreach(rep = 1:nReps,.combine = 'rbind') %dopar% {
      # snp = sample(1:ncol(X),1)
      snp = colnames(X)[rep]
      x = X[,snp]
      e1 = rnorm(n_strain);e1=e1/sd(e1)
      e2 = rnorm(n_strain);e2=e2/sd(e2)
      e3 = rnorm(n_full);e3=e3/sd(e3)
      e = sqrt(h2[[1]])*ZK_g %*% e1 + sqrt(h2[[2]]) * ZK_gxe %*% e2 + sqrt(1-sum(h2)) * e3
      e = e/sd(e)
      g = Z_trt * x
      g = g/sd(g)
      phen_tall$x = rep(x,2)
      phen_tall$y = sqrt(prop_var)*g + sqrt(1-prop_var)*e
      phen$plasticity = tapply(phen_tall$y,phen_tall$ecotype_id,function(x) x[2]-x[1])[phen$ecotype_id]
      phen$x = x
      
      
      X_tall = model.matrix(~x*variable,phen_tall)
      X_plasticity = model.matrix(~x,phen)
      
      # 1RE model - ecotype_id
      lm1 = mmer2(y~variable+x+variable:x,random = ~g(ecotype_id),data = phen_tall,G = list(ecotype_id = K),iters = 100,silent = T)
      lm1_F = get_F(get_h2_sommer(lm1),list(K_g),phen_tall$y,X_tall[,c(2,4)],X_tall[,c(1,3)])$F.2
      
      # 1RE model - ecotype_id x env
      lm1b = mmer2(y~variable+x+variable:x,random = ~g(ID),data = phen_tall,G = list(ID = K_gxe),iters = 100,silent = T)
      lm1b_F = get_F(get_h2_sommer(lm1b),list(K_gxe),phen_tall$y,X_tall[,c(2,4)],X_tall[,c(1,3)])$F.2

      # # 2RE model
      lm2 = mmer2(y~variable+x+variable:x,random = ~g(ecotype_id) + g(ID),data = phen_tall,G = list(ecotype_id = K, ID = K_gxe),iters = 100,silent = T)
      lm2_0 = mmer2(y~variable,random = ~g(ecotype_id) + g(ID),data = phen_tall,G = list(ecotype_id = K, ID = K_gxe),iters = 100,silent = T)
      lm2_F = get_F(get_h2_sommer(lm2),list(K_g,K_gxe),phen_tall$y,X_tall[,c(2,4)],X_tall[,c(1,3)])$F.2

      # 1RE plasticity
      lm_plasticity = mmer2(plasticity~x,random=~g(ecotype_id),data=phen,G = list(ecotype_id = K[phen$ecotype_id,phen$ecotype_id]),iters = 100,silent = T)
      lm_plasticity_F = get_F(get_h2_sommer(lm_plasticity),list(K[phen$ecotype_id,phen$ecotype_id]),phen$plasticity,X_plasticity[,2,drop=FALSE],X_plasticity[,1,drop=FALSE])$F

      # EMMAX-2
      EMMAX2_F = get_F(get_h2_sommer(lm2_0),list(K_g,K_gxe),phen_tall$y,X_tall[,c(2,4)],X_tall[,c(1,3)])$F.2

      # Grid-LMM - 10 steps
      Grid_10_F = get_F(nearest_h2(get_h2_sommer(lm2),all_h2s10),list(K_g,K_gxe),phen_tall$y,X_tall[,c(2,4)],X_tall[,c(1,3)])$F.2

      # Grid-LMM - 20 steps
      Grid_20_F = get_F(nearest_h2(get_h2_sommer(lm2),all_h2s20),list(K_g,K_gxe),phen_tall$y,X_tall[,c(2,4)],X_tall[,c(1,3)])$F.2

      # Grid-LMM - 100 steps
      Grid_100_F = get_F(nearest_h2(get_h2_sommer(lm2),all_h2s100),list(K_g,K_gxe),phen_tall$y,X_tall[,c(2,4)],X_tall[,c(1,3)])$F.2
            
      data.frame(h2,prop_var,rep,snp
                 ,get_h2_sommer(lm1),get_h2_sommer(lm1b)
                 ,p_1 = -log10(pf(lm1_F,1,n_full - 4,lower.tail=F))
                 ,p_1b = -log10(pf(lm1b_F,1,n_full - 4,lower.tail=F))
                 ,p_2 = -log10(pf(lm2_F,1,n_full - 4,lower.tail=F))
                 ,p_3 = -log10(pf(lm_plasticity_F,1,n_strain - 2,lower.tail=F))
                 ,p_6 = -log10(pf(EMMAX2_F,1,n_full - 4,lower.tail=F))
                 ,p_7 = -log10(pf(Grid_10_F,1,n_full - 4,lower.tail=F))
                 ,p_8 = -log10(pf(Grid_20_F,1,n_full - 4,lower.tail=F))
                 ,p_9 = -log10(pf(Grid_100_F,1,n_full - 4,lower.tail=F))
                 )
      
    }
  }
}

saveRDS(results,file = sprintf('Results/Atwell_power_Wald/set1_%s.rds',run))

# collect all results
results = foreach(file = list.files(path='Results/Atwell_power_Wald',pattern = 'set1',full.names = T),.combine = 'rbind') %do% {
  readRDS(file)
}
saveRDS(results,file = 'Results/Atwell_power_Wald/compiled_results.rds')
