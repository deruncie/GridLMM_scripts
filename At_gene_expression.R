library(GridLMM)
library(data.table)
library(foreach)
library(doParallel)

library(cowplot)
library(rstan)
# download file from: http://signal.salk.edu/1001.php
exp = fread('Data/At_1001/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv',data.table=F,check.names = F)

rownames(exp) = exp[,1]
exp = exp[,-1]
ACC_ID = colnames(exp)
ACC_ID = sub('X','',ACC_ID)
colnames(exp) = ACC_ID
exp = as.matrix(exp)
exp = exp[rowMeans(exp)>=10,]

# write.table(cbind(ACC_ID,ACC_ID),file = 'ACC_ID.txt',row.names=F,quote=F,col.names = F)

# load kinship matrix from 1001 genomes
# download file from: http://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5
K_1001_accessions = fread('Data/At_1001/1001_accessions.csv')
K_1001 = fread('Data/At_1001/1001_kinship.csv')
K_1001 = as.matrix(K_1001[,-1])
rownames(K_1001) = colnames(K_1001) = K_1001_accessions$x

exp = exp[,colnames(exp) %in% K_1001_accessions$x]

K = K_1001[colnames(exp),colnames(exp)]

normalize_K = function(K){
  # centers K and then normalizes to have mean diagonal = 1
  n = nrow(K)
  J = matrix(1,n,1)
  M = diag(1,n) - J %*% solve(t(J) %*% J) %*% t(J)
  centered_K = M %*% K %*% t(M)
  centered_K = centered_K/mean(diag(centered_K))
  rownames(centered_K) = rownames(K)
  colnames(centered_K) = colnames(K)
  centered_K
}
K = normalize_K(K)
write.csv(K,file = 'Data/At_1001/AT_GE_K.csv',row.names=T,col.names = T)

# normalize data

colData = data.frame(ACC_ID = colnames(exp))
design = model.matrix(~1,colData)

dds = DESeqDataSetFromMatrix(countData = exp,colData = colData,design=~1)
dds_vst = varianceStabilizingTransformation(dds,blind=T)
vst_matrix = assay(dds_vst)
fwrite(as.data.frame(vst_matrix),file = 'Data/At_1001/At_GE_VST.csv',row.names = T)

K = read.csv('Data/At_1001/AT_GE_K.csv',check.names = F,row.names=1)
K = as.matrix(K)

vst_matrix = fread('Data/At_1001/At_GE_VST.csv',data.table=F,header = T)
rownames(vst_matrix) = vst_matrix[,1]
vst_matrix = vst_matrix[,-1]
vst_matrix = as.matrix(vst_matrix)

data = data.frame(ACC_ID = rownames(K))

K_E = K*K
K = K/mean(diag(K))
K_E = K_E/mean(diag(K_E))
data$ID2 = data$ACC_ID

# prep GridLMM
data$y = vst_matrix[1,]
start = Sys.time()
m_null = GridLMM_ML(y~(1|ACC_ID)+(1|ID2),data,relmat = list(ACC_ID = K,ID2 = K_E),tolerance = 0.01,diagonalize = F)

h2_grid = t(setup_Grid(names(m_null$setup$RE_setup),h2_step = .05))

X_cov = model.matrix(~1,data)
b = ncol(X_cov)+1
n = nrow(data)
X_indices_h2 = cbind(1:nrow(vst_matrix),0,0)
registerDoParallel(10)
res_h2 = foreach(h2 = iter(h2_grid,by='row')) %dopar% {
  print(h2)
  colnames(h2) = NULL
  chol_Vi_R = make_chol_V_setup(m_null$setup,h2)$chol_V
  SSs = GridLMM_SS_matrix(t(vst_matrix),chol_Vi_R,X_cov,list(),numeric(),rep(0,2))
  SSs = get_LL(SSs,X_cov,list(),numeric(),n,1,T,T,T)
  res_i = data.frame(Trait = 1,X_ID = 1:nrow(SSs$beta_hats), REML_logLik = SSs$REML[1,], ID.REML = h2, 
                     LP = SSs$log_posterior_factor[,1],stringsAsFactors = F)
  return(res_i)
}
saveRDS(res_h2,file = 'At_1001_eqtl_GridLMM.rds')

# results for grid when h2[2] == 0

res_LP = res_h2

res_h2 = compile_results(res_h2)
Sys.time() - start

h2_grid1 = h2_grid[h2_grid[,2] == 0,]
res_LP_1 = res_LP[h2_grid[,2] == 0]
res_h2_1 = compile_results(res_LP_1)

# approximate student-t priors for variance components
# translate these to discrete-grid prior on h2s
prior_draws = matrix(rt(3*1e5,3)*10,ncol=3)^2
prior_h2s = prior_draws/rowSums(prior_draws)
prior_h2s = data.frame(ACC_ID = prior_h2s[,1],ID2 = prior_h2s[,2])

# make 2d histogram
# grid boundaries are 0.025, 0.075,... up to 0.925, and then 1
breaks = c(0,seq(0.05/2,1-0.05/2,by=0.05),1)[-21]
index.x = cut(prior_h2s$ACC_ID,breaks,include.lowest = T)
index.y = cut(prior_h2s$ID2,breaks,include.lowest = T)
counts = tapply(1:nrow(prior_h2s),list(index.x,index.y),length)
counts[is.na(counts)] = 0
prior_hist = list(x=seq(0,1,by=.05),y=seq(0,1,by=.05),counts = counts)
# prior_density = data.frame(expand.grid(ACC_ID = seq(0,1,by=.05),ID2=seq(0,1,by=.05)))
# prior_density$density = c(m)
prior_density = data.frame(h2_grid)
for(i in 1:nrow(prior_density)){
  ix = which.min(abs(prior_density$ACC_ID[i] - prior_hist$x))
  iy = which.min(abs(prior_density$ID2[i] - prior_hist$y))
  prior_density$density[i] = prior_hist$counts[20*(iy-1)+ix]
}
prior_density = subset(prior_density,rowSums(prior_density[,1:2])<=1)
prior_density$density = prior_density$density/sum(prior_density$density,na.rm=T)
ggplot(prior_density,aes(x=ACC_ID,y=ID2)) + geom_point(aes(size = density)) +
  stat_density_2d(data=prior_h2s)
hist(prior_h2s[,1],breaks=seq(0,1,length=21),prob=T)
mp = marginal_posterior(prior_density[,2],t(log(prior_density[,3,drop=FALSE])));plot(as.numeric(colnames(mp)),mp,ylim = c(0,max(mp,na.rm=T)))
# ggplot(prior_h2s,aes(x=ACC_ID,y=ID2)) + stat_density_2d()

# prior_density = subset(prior_density,rowSums(prior_density[,1:2])< 1)



# calculate posteriors
log_sum_log = function(x) {
  mx = max(x)
  mx + log(sum(exp(x-mx)))
}
normalize_posterior = function(lp_mat) {
  lp_mat = t(apply(lp_mat,1,function(x) {
    x - log_sum_log(x)
  }))
}
marginal_posterior = function(x_value,lp_mat){
  values = sort(unique(x_value))
  res = sapply(values,function(x) {
    rowSums(exp(lp_mat[,x_value == x,drop=FALSE]))
  })
  if(is.null(dim(res))) {
    res = matrix(res,nr=1)
  }
  colnames(res) = values
  res
}

lp_mat_uniform = sapply(res_LP,function(i) i$LP)
lp_mat_t = sweep(lp_mat_uniform,2,log(prior_density$density),'+')
lp_mat_uniform = normalize_posterior(lp_mat_uniform)
lp_mat_t = normalize_posterior(lp_mat_t)

h2_post = as.data.frame(h2_grid)

lp_add = marginal_posterior(h2_grid[,1],lp_mat_uniform)
pm_add = rowSums(sweep(lp_add,2,as.numeric(colnames(lp_add)),'*'))
lp_epi = marginal_posterior(h2_grid[,2],lp_mat_uniform)
pm_epi = rowSums(sweep(lp_epi,2,as.numeric(colnames(lp_epi)),'*'))


lp_add_t = marginal_posterior(h2_grid[,1],lp_mat_t)
pm_add_t = rowSums(sweep(lp_add_t,2,as.numeric(colnames(lp_add_t)),'*'))
lp_epi_t = marginal_posterior(h2_grid[,2],lp_mat_t)
pm_epi_t = rowSums(sweep(lp_epi_t,2,as.numeric(colnames(lp_epi_t)),'*'))


# example genes:
i1 = which(abs(res_h2$ID.REML.1 - 0.05)<0.1 & abs(res_h2$ID.REML.2 - 0.75)<0.05)[1];
i2 = which(abs(res_h2$ID.REML.1 - 0.55)<0.1 & abs(res_h2$ID.REML.2 - 0.05)<0.05)[1];
i1 = 78
i2 = 32


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)#parallel::detectCores())

const = 1e-1
K1 = K*(1-const) + diag(const,nrow(K));K1[lower.tri(K1)] = t(K1)[lower.tri(K1)]
K_E1 = K_E*(1-const) + diag(const,nrow(K));K_E[lower.tri(K_E)] = t(K_E)[lower.tri(K_E)]


sK1 = svd(K1)
Ut = t(sK1$u)
d_sqrt = sqrt(sK1$d)
D = diag(d_sqrt^2);rownames(D) = colnames(D) = rownames(K1)
UtK2U = Ut %*% K_E %*% t(Ut)
rownames(UtK2U) = colnames(UtK2U) = rownames(K_E)
cUtK2U = t(chol(UtK2U))

d = data.frame(y = vst_matrix[i1,])
d$Uty = Ut %*% d$y
d$Ut1 = Ut %*% matrix(1,n)

stan_fit1 = stan('MM_twoRE_diagonal1.stan', 
                 data =list(
                   N = n,
                   Y = c(d$Uty),
                   b=1,
                   X = matrix(d$Ut1,n),
                   D_sqrt_cov_1 = d_sqrt,
                   N2 = n,
                   Lcov_2 = cUtK2U
                 ), 
                 control = list(adapt_delta = 0.99,max_treedepth = 15),
                 iter=10000,
                 warmup = 5000,
                 thin = 20,
                 chains = 4,
                 cores = 4,
                 pars=c('beta','sigma','sds')
)
saveRDS(stan_fit1,file = 'i1_v3.rds')
stan_fit1 = readRDS('i1_v3.rds')

d = data.frame(y = vst_matrix[i2,])
d$Uty = Ut %*% d$y
d$Ut1 = Ut %*% matrix(1,n)

stan_fit2 = stan('MM_twoRE_diagonal1.stan', 
                 data =list(
                   N = n,
                   Y = c(d$Uty),
                   b=1,
                   X = matrix(d$Ut1,n),
                   D_sqrt_cov_1 = d_sqrt,
                   N2 = n,
                   Lcov_2 = cUtK2U
                 ), 
                 control = list(adapt_delta = 0.999,max_treedepth = 15),
                 iter=10000,
                 warmup = 5000,
                 thin = 20,
                 chains = 4,
                 cores = 4,
                 pars=c('beta','sigma','sds')
)
saveRDS(stan_fit2,file = 'i2_v3.rds')
stan_fit2 = readRDS('i2_v3.rds')

stan_fit1
stan_fit2
sum(get_elapsed_time(stan_fit1))/60/60 # hours
sum(get_elapsed_time(stan_fit2))/60/60


samples = do.call(cbind,extract(stan_fit1,c('sigma','sds')))
samples = samples^2;samples = (samples/rowSums(samples))[,2:3]
h2s_bm1 = data.frame(ACC_ID = samples[,1],ID2 = samples[,2]) * (1-const)
samples = do.call(cbind,extract(stan_fit2,c('sigma','sds')))
samples = samples^2;samples = (samples/rowSums(samples))[,2:3]
h2s_bm2 = data.frame(ACC_ID = samples[,1],ID2 = samples[,2]) * (1-const)

# plots
p1 = ggplot(res_h2,aes(x=ID.REML.1,y=ID.REML.2)) + geom_jitter(size = .1,alpha = 0.1) + 
  geom_point(data = res_h2[c(i1,i2),,drop=FALSE],aes(x=ID.REML.1,y=ID.REML.2),col='red',size=3) + 
  xlab('K_A') + ylab('K_E');p1
h2_post$posterior = exp(lp_mat_uniform[i1,])
p2 = ggplot(h2_post,aes(x=ACC_ID,y=ID2)) + xlim(c(0,1)) + ylim(c(0,1)) + 
  stat_density_2d(data = h2s_bm1,contour = T,size=.5,color='grey30',bins = 20) + 
  geom_point(aes(size = posterior)) + scale_size_continuous(range = c(0.0001,4)) + 
  geom_point(data = res_h2[i1,,drop=FALSE],aes(x=ID.REML.1,y=ID.REML.2),col='red',size=1) + 
  geom_point(data = data.frame(ACC_ID = pm_add[i1],ID2 = pm_epi[i1]),col='red',size=3,shape=3,stroke=2) +  
  geom_point(data = data.frame(ACC_ID = mean(h2s_bm1[,1]),ID2 = mean(h2s_bm1[,2])),col='blue',size=3,shape=3,stroke=2) + 
  xlab('K_A') + ylab('K_E') + ggtitle(rownames(vst_matrix)[i1]) + theme(legend.position = 'none');p2
h2_post$posterior = exp(lp_mat_uniform[i2,])
p3 = ggplot(h2_post,aes(x=ACC_ID,y=ID2)) + xlim(c(0,1)) + ylim(c(0,1)) + 
  stat_density_2d(data = h2s_bm2,contour = T,size=.5,color='grey30',bins = 20) + 
  geom_point(aes(size = posterior)) + scale_size_continuous(range = c(0.0001,3)) + 
  geom_point(data = res_h2[i2,,drop=FALSE],aes(x=ID.REML.1,y=ID.REML.2),col='red',size=1) + 
  geom_point(data = data.frame(ACC_ID = pm_add[i2],ID2 = pm_epi[i2]),col='red',size=3,shape=3,stroke=2) +  
  geom_point(data = data.frame(ACC_ID = mean(h2s_bm2[,1]),ID2 = mean(h2s_bm2[,2])),col='blue',size=3,shape=3,stroke=2) + 
  xlab('K_A') + ylab('K_E') + ggtitle(rownames(vst_matrix)[i2]) + theme(legend.position = 'none');p3
p = plot_grid(p1,p2,p3,nrow = 1,labels = 'auto');p
save_plot(plot=p,base_aspect_ratio = 3,filename = 'AT_GE_uniform_prior.pdf',base_height = 4)

h2_A_comparison = data.frame(M1 = res_h2_1$ID.REML.1,M2 = res_h2$ID.REML.1)
p4 = ggplot(h2_A_comparison,aes(x=M2,y=M1)) + geom_jitter(size = 0.1,alpha = 0.1) + xlab('K_A with K_E') + ylab('K_A without K_E') +
  geom_abline(slope=1,intercept=0);p4
save_plot(plot=p4,filename= 'AT_GE_model_comparison_KA.pdf',base_height=6)


h2_post$posterior = exp(lp_mat_t[i1,])
p2b = ggplot(h2_post,aes(x=ACC_ID,y=ID2)) + xlim(c(0,1)) + ylim(c(0,1)) + 
  stat_density_2d(data = h2s_bm1,contour = T,size=.5,color='grey30',bins = 20) + 
  geom_point(aes(size = posterior)) + scale_size_continuous(range = c(0.0001,4)) + 
  geom_point(data = res_h2[i1,,drop=FALSE],aes(x=ID.REML.1,y=ID.REML.2),col='red',size=1) + 
  geom_point(data = data.frame(ACC_ID = pm_add_t[i1],ID2 = pm_epi_t[i1]),col='red',size=3,shape=3,stroke=2) +  
  geom_point(data = data.frame(ACC_ID = mean(h2s_bm1[,1]),ID2 = mean(h2s_bm1[,2])),col='blue',size=3,shape=3,stroke=2) + 
  xlab('K_A') + ylab('K_E') + ggtitle(rownames(vst_matrix)[i1]) + theme(legend.position = 'none');p2b
h2_post$posterior = exp(lp_mat_t[i2,])
p3b = ggplot(h2_post,aes(x=ACC_ID,y=ID2)) + xlim(c(0,1)) + ylim(c(0,1)) + 
  stat_density_2d(data = h2s_bm2,contour = T,size=.5,color='grey30',bins = 20) + 
  geom_point(aes(size = posterior)) + scale_size_continuous(range = c(0.0001,3)) + 
  geom_point(data = res_h2[i2,,drop=FALSE],aes(x=ID.REML.1,y=ID.REML.2),col='red',size=1) + 
  geom_point(data = data.frame(ACC_ID = pm_add_t[i2],ID2 = pm_epi_t[i2]),col='red',size=3,shape=3,stroke=2) +  
  geom_point(data = data.frame(ACC_ID = mean(h2s_bm2[,1]),ID2 = mean(h2s_bm2[,2])),col='blue',size=3,shape=3,stroke=2) + 
  xlab('K_A') + ylab('K_E') + ggtitle(rownames(vst_matrix)[i2]) + theme(legend.position = 'none');p3b
p = plot_grid(p2b,p3b,nrow = 1,labels = 'auto');p
save_plot(plot=p,base_aspect_ratio = 2,filename = 'AT_GE_t_prior.pdf',base_height = 4)


# different priors
# student-t
prior_draws = matrix(rt(3*1e5,3)*10,ncol=3)^2
prior_h2s = prior_draws/rowSums(prior_draws)
prior_h2s = data.frame(K_A = prior_h2s[,1],K_E = prior_h2s[,2])
p1=ggplot(prior_h2s,aes(x=K_A,y=K_E)) + ggtitle(expression(paste(sigma,'~half-t(3,0,10)'))) + 
  stat_density_2d(aes(fill = stat(density)),geom='raster',contour=F) + 
  theme(legend.position = 'none') # + stat_density_2d(bins=30)
p2=ggplot(prior_h2s,aes(x=K_A)) + stat_density(geom='line') + xlim(c(0,1))

prior_draws = matrix(1/rgamma(3*1e5,shape=2,rate=1),ncol=3)
prior_h2s = prior_draws/rowSums(prior_draws)
prior_h2s = data.frame(K_A = prior_h2s[,1],K_E = prior_h2s[,2])
p3=ggplot(prior_h2s,aes(x=K_A,y=K_E)) + ggtitle(expression(paste(sigma^2,'~invGamma(2,1)'))) + 
  stat_density_2d(aes(fill = stat(density)),geom='raster',contour=F) + 
  theme(legend.position = 'none') # + stat_density_2d(bins=30)
p4=ggplot(prior_h2s,aes(x=K_A)) + stat_density(geom='line') + xlim(c(0,1))

library(gtools)
prior_h2s = data.frame(rdirichlet(1e5,rep(1,3)))
colnames(prior_h2s) = c('K_A','K_E','e')
p5=ggplot(prior_h2s,aes(x=K_A,y=K_E)) + ggtitle(expression(paste(h^2,'~uniform(0,1)'))) + 
  stat_density_2d(aes(fill = stat(density)),geom='raster',contour=F) + 
  theme(legend.position = 'none') # + stat_density_2d(bins=30)
p6=ggplot(prior_h2s,aes(x=K_A)) + stat_density(geom='line') + xlim(c(0,1))

p = plot_grid(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,labels = 'auto');p
save_plot(plot=p,filename='AT_priors_comparison.pdf',base_aspect_ratio = 2/3,base_height=8)


# LDAK


LDAK_path = file.path('./../misc_software/ldak5') # change path to your LDAK program
folder = paste0('G',0)
try(dir.create(folder))
setwd(folder)
K_list = c(K_A = prep_LDAK_Kinship(K,'../K_A',LDAK_path),K_E = prep_LDAK_Kinship(K_E,'../K_E',LDAK_path))
setwd('..')
system(sprintf('rm -rf %s',folder))

r = foreach(i = 1:nrow(vst_matrix),.combine='rbind') %dopar% {
  folder = paste0('G',i)
  try(dir.create(folder))
  setwd(folder)
  data$y = vst_matrix[i,]
  LDAK_call = prep_h2_LDAK(data$y,matrix(1,nrow(data)),K_list,LDAK_path,1000)
  system(LDAK_call)
  time_i = (Sys.time() - start)
  total_time = total_time + time_i
  res = fread('h2_LDAK.progress')
  nIter = nrow(res)
  times = rbind(times,data.frame(Time = time_i,nIter = nIter))
  h2_start = get_LDAK_result(K_list)
  setwd('..')
  system(sprintf('rm -rf %s',folder))
  h2_start
}
saveRDS(r,file = 'At_1001_eQTL_LDAK_h2s.rds')


## estimate I/O


r = matrix(0,nrow(vst_matrix),2)
total_time = 0
times = c()
for(i in 1:nrow(vst_matrix)){
  
  start = Sys.time()
  data$y = vst_matrix[i,]
  LDAK_call = prep_h2_LDAK(data$y,matrix(1,nrow(data)),K_list,LDAK_path,1000)
  system(LDAK_call)
  time_i = (Sys.time() - start)
  total_time = total_time + time_i
  res = fread('h2_LDAK.progress')
  nIter = nrow(res)
  times = rbind(times,data.frame(Time = time_i,nIter = nIter))
  h2_start = get_LDAK_result(K_list)
  r[i,] = h2_start
}
saveRDS(times,file = 'At_1001_eQTL_LDAK_times.rds')
saveRDS(r,file = 'At_1001_eQTL_LDAK_h2s.rds')
const = 0 # what is this?

ldak_io_test = foreach(i = 1:3) %do% {
  times = time_LDAK(vst_matrix[i,],matrix(1,nrow(data)),K_list,LDAK_path)
}
ldak_io_test_results = sapply(ldak_io_test,function(x) coef(lm(time~n_iter,x)))
ldak_io_test_results[1,]/ldak_io_test_results[2,]
ldak_io_test_results[1,]/(ldak_io_test_results[1,]+8*ldak_io_test_results[2,])  # estimate 6

h2_A = data.frame(GridLMM=res_h2$ID.REML.1, LDAK = r[,1]*(1-const))
h2_E = data.frame(GridLMM=res_h2$ID.REML.2, LDAK = r[,2]*(1-const))

p5 = ggplot(h2_A,aes(x=LDAK,y=GridLMM)) + geom_point(alpha = 0.1) + geom_abline(slope=1,intercept=0) + xlim(c(0,1)) + ylim(c(0,1)) + ggtitle('K_A');p5
p6 = ggplot(h2_E,aes(x=LDAK,y=GridLMM)) + geom_point(alpha = 0.1) + geom_abline(slope=1,intercept=0) + xlim(c(0,1)) + ylim(c(0,1)) + ggtitle('K_E');p6
p = plot_grid(p5,p6,labels = 'auto')
save_plot(p,base_aspect_ratio = 2,filename = 'AT_GE_LDAK_comparison.pdf')

