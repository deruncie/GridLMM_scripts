
library(GridLMM)
library(data.table)
library(sommer)
library(snpStats)
library(reshape2)
library(foreach)

LDAK_path = 'misc_software/ldak5' # set path to LDAK program

dataset = readRDS('Data/Atwell/dataset.rds') # create with prep_Arabidopsis_data.R
log_traits = fread('Data/Atwell/pheno_transforms.csv',data.table = F)

phen = dataset$phenotypes
map = dataset$map
X = dataset$X

trait_1 = '5_FT10'
trait_2 = '57_FT Field'

phen = na.omit(phen[,c('ecotype_id',trait_1,trait_2)])

# standardize traits
phen[,-1] = apply(phen[,-1],2,function(x) (x-mean(x))/sd(x))

phen_tall = melt(phen,id.vars = 'ecotype_id')
phen_tall$ecotype_id2 = phen_tall$ecotype_id
phen$plasticity = phen[,3]-phen[,2]
X = X[phen$ecotype_id,]
X = X[,apply(X,2,var)>0]

K = A.mat(X-1)

options(contrasts = c('contr.sum','contr.poly'))
null_model = GridLMM_ML(value ~ variable + (1+variable|ecotype_id) + (1|ecotype_id2),phen_tall,relmat = list(ecotype_id = K),tolerance = 10)
RE_setup = null_model$setup$RE_setup
full_RE_names = names(RE_setup)
names(RE_setup) = c('G','GxE','accession')
Ks = lapply(names(RE_setup),function(x) {
  K = as.matrix(with(RE_setup[[x]], Z %*% K %*% t(Z)))
})
names(Ks) = names(RE_setup)
try(dir.create('Atwell_working'))
K_list = sapply(names(RE_setup),function(x) {
  prep_LDAK_Kinship(Ks[[x]],paste0('Atwell_working/',x),LDAK_path)
  paste0('Atwell_working/',x)
})


X_cov = model.matrix(~variable,phen_tall)
get_h2_LDAK(phen_tall$value,X_cov,K_list,LDAK_path,weights = sapply(Ks,function(x) mean(diag(x))))


method_times = c()
all_results = data.frame(X_ID = colnames(X))

# ---------- LDAK ----------- #
# force the number of iterations of the optimization algorithm to change to estimate I/O limits to time
ldak_io_test = foreach(i = 1:3) %do% {
  times = time_LDAK(phen_tall$value,X_cov,K_list,LDAK_path)
}
ldak_io_test_results = sapply(ldak_io_test,function(x) coef(lm(time~n_iter,x)))
ldak_io_test_results[1,]/ldak_io_test_results[2,]
ldak_io_test_results[1,]/(ldak_io_test_results[1,]+11*ldak_io_test_results[2,])

# 
ldak_times = foreach(i = 1:10,.combine = 'c') %do% {
  start = Sys.time()
  h2_start = get_h2_LDAK(phen_tall$value,X_cov,K_list,LDAK_path,weights = sapply(Ks,function(x) mean(diag(x))))
  Sys.time() - start  # 3 minutes on my laptop
}
method_times$LDAK_AEcage = ldak_times



## GEMMA-plasticity
gemma_dir = sprintf('GEMMA_plasticity_%s',make.names(trait_2))
dir.create(gemma_dir)
write.plink(file.base = sprintf('%s/At',gemma_dir),snps = as(X,'SnpMatrix'),id = phen$ecotype_id,phenotype = phen$plasticity)
write.table(K,file = sprintf('%s/At_K.txt',gemma_dir),row.names=F,col.names=F)
start = Sys.time()
system(sprintf('./gemma.macosx -bfile %s/At -k %s/At_K.txt -lmm 1 -maf 0 -r2 1 -outdir ./%s -o At_test',gemma_dir,gemma_dir,gemma_dir)) # -region 100 -lmax 1e4
# system(sprintf('~/software/gemma/bin/gemma -bfile %s/At -k %s/At_K.txt -lmm 1 -maf 0 -r2 1 -outdir ./%s -o At_test',gemma_dir,gemma_dir,gemma_dir))
(method_times$GEMMA = Sys.time() - start)

results_GEMMA = fread(sprintf('%s/At_test.assoc.txt',gemma_dir),data.table = F)
results_GEMMA$p_value_REML = results_GEMMA$p_wald
all_results$GEMMA = results_GEMMA$p_value_REML
system(sprintf('rm -rf %s',gemma_dir))


## Grid-LMM no K_GxE
start = Sys.time()
# h2_start = get_h2_LDAK(phen_tall$value,X_cov,K_list[c(1,3)],LDAK_path,weights = sapply(Ks[c(1,3)],function(x) mean(diag(x))))
# colnames(h2_start) = c('ecotype_id','ecotype_id2')
Grid_GxE_noK_GxE = GridLMM_GWAS_fast(value ~ variable + (1|ecotype_id) + (1|ecotype_id2),~variable,~1,data=phen_tall,X = X,X_ID = 'ecotype_id',
                                     # h2_start = h2_start,
                                    relmat = list(ecotype_id = K),mc.cores = 1,max_steps = 100,h2_step = 0.01,method = 'REML')
(method_times$Grid_GxE_noK_GxE_time = Sys.time() - start)
all_results$Grid_GxE_noK_GxE = Grid_GxE_noK_GxE$results$p_value_REML.2


## Grid-LMM with K_GxE
start = Sys.time()
h2_start = get_h2_LDAK(phen_tall$value,X_cov,K_list,LDAK_path,weights = sapply(Ks,function(x) mean(diag(x))))
colnames(h2_start) = full_RE_names
Grid_GxE_GxE = GridLMM_GWAS_fast(value ~ variable + (1+variable|ecotype_id) + (1|ecotype_id2),~variable,~1,data=phen_tall,X = X,X_ID = 'ecotype_id',
                                 h2_start = h2_start,
                                 relmat = list(ecotype_id = K),mc.cores = 1,max_steps = 100,h2_step = 0.01,method = 'REML')
(method_times$Grid_GxE_GxE_time = Sys.time() - start)
all_results$Grid_GxE_GxE = Grid_GxE_GxE$results$p_value_REML.2


## Sul with K_GxE Sul
start = Sys.time()
h2_start = get_h2_LDAK(phen_tall$value,X_cov,K_list,LDAK_path,weights = sapply(Ks,function(x) mean(diag(x))))
colnames(h2_start) = full_RE_names
Sul_GxE_GxE = GridLMM_GWAS_fast(value ~ variable + (1+variable|ecotype_id) + (1|ecotype_id2),~variable,~1,data=phen_tall,X = X,X_ID = 'ecotype_id',
                                h2_start = h2_start,
                                relmat = list(ecotype_id = K),mc.cores = 1,max_steps = 100,h2_step = 10.01,method = 'REML')
(method_times$Sul_GxE_GxE_time = Sys.time() - start)
all_results$Sul_GxE_GxE = Sul_GxE_GxE$results$p_value_REML.2


saveRDS(method_times,file = 'Atwell_working/timings.rds')
saveRDS(all_results,file = 'Atwell_working/results.rds')
results_list = list(
  results_GEMMA = results_GEMMA,
  Grid_GxE_noK_GxE = Grid_GxE_noK_GxE,
  Grid_GxE_GxE = Grid_GxE_GxE,
  Sul_GxE_GxE = Sul_GxE_GxE
)
saveRDS(results_list,file = 'Atwell_working/all_results.rds')

