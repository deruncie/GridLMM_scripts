# script to run GridLMM on a single trait from Atwell data

library(GridLMM)
library(data.table)
library(sommer)
library(snpStats)
library(foreach)
library(doParallel)
registerDoParallel(my_detectCores())

dataset = readRDS('Data/dataset.rds') # create with prep_Arabidopsis_data.R
log_traits = fread('Datapheno_transforms.csv',data.table = F)

phen = dataset$phenotypes
map = dataset$map
X = dataset$X

trait = as.numeric(commandArgs(t=T)[1])
if(is.na(trait)) trait = 7
trait_name = colnames(phen)[trait]
print(trait)
print(trait_name)
phen$y = phen[[trait_name]]
if(log_traits$Log[log_traits$Trait == trait_name] == 1) phen$y = log(phen$y)

phen = subset(phen,!is.na(y))
X = X[phen$ecotype_id,]
X = X[,apply(X,2,var)>0]

# use sommer to make K
K = A.mat(X - 1)

# run GEMMA
gemma_dir = sprintf('GEMMA_%03d',trait)
dir.create(gemma_dir)
write.plink(file.base = sprintf('%s/At',gemma_dir),snps = as(X,'SnpMatrix'),id = phen$ecotype_id,phenotype = phen$y)
write.table(K,file = sprintf('%s/At_K.txt',gemma_dir),row.names=F,col.names=F)
# system(sprintf('./gemma.macosx -bfile %s/At -k %s/At_K.txt -lmm 1 -maf 0 -r2 1 -outdir ./%s -o At_test',gemma_dir,gemma_dir,gemma_dir)) # -region 100 -lmax 1e4
start = Sys.time()
system(sprintf('~/software/gemma/bin/gemma -bfile %s/At -k %s/At_K.txt -lmm 1 -maf 0 -r2 1 -outdir ./%s -o At_test',gemma_dir,gemma_dir,gemma_dir))
gemma_time = Sys.time() - start

results_GEMMA = fread(sprintf('%s/At_test.assoc.txt',gemma_dir),data.table = F)
results_GEMMA$p_value_REML = results_GEMMA$p_wald
system(sprintf('rm -rf %s',gemma_dir))

# Run Grid-LMM - 0.1
start = Sys.time()
results_0.1 = GridLMM_GWAS(y~(1|ecotype_id),~1,~0,data = phen,X = X,X_ID = 'ecotype_id',relmat = list(ecotype_id = K),X_map = map,h2_divisions = 10,REML = T,ML=F,mc.cores=1)
Grid_0.1_time = Sys.time() - start


# Run Grid-LMM - 0.01
start = Sys.time()
results_0.01_full = GridLMM_GWAS(y~(1|ecotype_id),~1,~0,data = phen,X = X,X_ID = 'ecotype_id',relmat = list(ecotype_id = K),X_map = map,h2_divisions = 100,REML = T,ML=F,mc.cores=1)
Grid_0.01_full_time = Sys.time() - start

# find ML
start = Sys.time()
null_model = GridLMM_ML(y~(1|ecotype_id),data = phen,relmat = list(ecotype_id = K),mc.cores=1)
(h2_start = get_current_h2s(null_model$results,'ecotype_id',F,T))
null_time = Sys.time() - start

# Grid-LMM - 0.01 adaptive
start = Sys.time()
results_0.01 = GridLMM_GWAS_fast(y~(1|ecotype_id),~1,~0,data = phen,X = X,X_ID = 'ecotype_id',relmat = list(ecotype_id = K),X_map = map,
                                V_setup = null_model$setup,h2_start = h2_start,h2_step = 0.01,method = 'REML',mc.cores = 1)
Grid_0.01_time = Sys.time() - start

# Grid-LMM - EMMAX
start = Sys.time()
results_EMMAX = GridLMM_GWAS_fast(y~(1|ecotype_id),~1,~0,data = phen,X = X,X_ID = 'ecotype_id',relmat = list(ecotype_id = K),X_map = map,
                                 V_setup = null_model$setup,h2_start = h2_start,h2_step = 1,method = 'REML',mc.cores = 1)
EMMAX_time = Sys.time() - start

results = list(
  results_GEMMA = results_GEMMA,
  results_0.01 = results_0.01$results,
  results_0.01_full = results_0.01_full$results,
  results_0.1 = results_0.1$results,
  results_EMMAX = results_EMMAX$results,
  times = c(gemma_time,null_time,Grid_0.1_time,Grid_0.01_time,EMMAX_time)
)
saveRDS(results,file = sprintf('Results/Atwell_oneTrait/results_%03db.rds',trait))

# script to collect all traits into a summary file
# errors = foreach(trait = 3:109,.combine = 'rbind') %dopar% {
#   if(sprintf('results_%03d.rds',trait) %in% list.files(path='Results/Atwell_oneTrait')) {
#     results = readRDS(sprintf('Results/Atwell_oneTrait/results_%03d.rds',trait))
#     results_wide = do.call(cbind,lapply(results,function(x) x$p_value_REML))
#     # i = match(results$results_GEMMA$rs,results$results_0.01$X_ID)
#     # results_wide = cbind(results_GEMMA=results$results_GEMMA$p_value_REML,do.call(cbind,lapply(results[-1],function(x) x$p_value_REML[i])))
#     # results_wide = data.frame(-log10(results_wide))
#     ylim = c(0,max(results_wide,na.rm=T))
#     png(sprintf('Results/Atwell_oneTrait/plot_results_%03d.png',trait,colnames(phen)[trait]))
#     plot(results_wide$results_GEMMA,results_wide$results_EMMAX,ylim=ylim,xlim=ylim,cex=0.5,pch=21,bg=1,col=NULL,main = colnames(phen)[trait]);abline(0,1)
#     points(results_wide$results_GEMMA,results_wide$results_0.1,cex=0.5,pch=21,bg=2,col=NULL)
#     points(results_wide$results_GEMMA,results_wide$results_0.01,cex=0.5,pch=21,bg=3,col=NULL)
#     dev.off()
#     c(cor(results_wide,use='p')[1,-1],sqrt(colMeans((results_wide[,-1]-results_wide[,1])^2,na.rm=T)))
#   }
# }
# write.csv(errors,file = 'Results/Atwell_oneTrait/errors.csv',row.names=F)
# 
