# using pheno data from BGLR
library(GridLMM)
library(data.table)
library(snpStats)
library(foreach)
library(doParallel)
library(sommer)
library(lme4)

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
  E = prep_LDAK_Kinship(trait_E,'K_E',LDAK_path),
  cage = prep_LDAK_Kinship(cage_K,'K_cage',LDAK_path)
)

X_cov = model.matrix(~GENDER,trait_data)
get_h2_LDAK(trait_data$y,X_cov,K_list[1:2],LDAK_path)


method_times = c()
all_results = data.frame(X_ID = colnames(trait_X))

# ---------- LDAK ----------- #
# force the number of iterations of the optimization algorithm to change to estimate I/O limits to time
ldak_io_test = foreach(i = 1:3) %do% {
  times = time_LDAK(trait_data$y,cbind(X_cov,trait_X[,i]),K_list[c(1,2,3)],LDAK_path)
}
ldak_io_test_results = sapply(ldak_io_test,function(x) coef(lm(time~n_iter,x)))
# 

ldak_times = foreach(i = 1:10,.combine = 'c') %do% {
  start = Sys.time()
  h2_start = get_h2_LDAK(trait_data$y,cbind(X_cov,trait_X[,i]),K_list[c(1,2,3)],LDAK_path)
  Sys.time() - start  # 3 minutes on my laptop
}
method_times$LDAK_AEcage = ldak_times

# ---------- GEMMA ----------- #
gemma_dir = 'GEMMA'
dir.create(gemma_dir)
write.plink(file.base = sprintf('%s/WTCC',gemma_dir),snps = as(trait_X,'SnpMatrix'),id = trait_data$ID,phenotype = trait_data$y)
write.table(trait_A,file = sprintf('%s/WTCC_A.txt',gemma_dir),row.names=F,col.names=F)
cov = model.matrix(~1+GENDER,trait_data)
write.table(cov,file = sprintf('%s/cov.txt',gemma_dir),row.names=F,col.names=F)
start = Sys.time()
system(sprintf('../gemma.macosx -bfile %s/WTCC -k %s/WTCC_A.txt -c %s/cov.txt -lmm 1 -maf 0 -r2 1 -outdir ./%s -o WTCC_test',gemma_dir,gemma_dir,gemma_dir,gemma_dir)) # -region 100 -lmax 1e4
# system(sprintf('~/software/gemma/bin/gemma -bfile %s/WTCC -k %s/WTCC_A.txt -lmm 1 -maf 0 -r2 1 -outdir ./%s -o WTCC_test',gemma_dir,gemma_dir,gemma_dir))
gemma_time = Sys.time() - start
method_times$GEMMA = gemma_time

results_GEMMA = fread(sprintf('%s/WTCC_test.assoc.txt',gemma_dir),data.table = F)
results_GEMMA$p_value_REML = results_GEMMA$p_wald
all_results$GEMMA = results_GEMMA$p_value_REML



# ---------- Grid-A ----------- #
start = Sys.time()
h2_start = get_h2_LDAK(trait_data$y,X_cov,K_list[1],LDAK_path)
colnames(h2_start)[1] = 'ID'
grid_A = GridLMM_GWAS_fast(y~GENDER + (1|ID),~1,~0,data=trait_data,X = trait_X,X_ID = 'ID',relmat = list(ID = trait_A), h2_start = h2_start,#diagonalize = F,
                                clusterType = 'mclapply',mc.cores = 1,max_steps = 100,h2_step = 0.01,method = 'REML')
(grid_A_time = Sys.time() - start)
method_times$grid_A = grid_A_time
all_results$grid_A = grid_A$results$p_value_REML

# ---------- Grid-AEcage ----------- #
start = Sys.time()
h2_start = get_h2_LDAK(trait_data$y,X_cov,K_list[c(1,2,3)],LDAK_path)
colnames(h2_start) = c('ID','ID2','ID3')
grid_AEcage = GridLMM_GWAS_fast(y~GENDER + (1|ID) + (1|ID2) + (1|ID3),~1,~0,data=trait_data,X = trait_X,X_ID = 'ID',relmat = list(ID = trait_A,ID2 = trait_E,ID3 = cage_K), 
                               h2_start = h2_start,
                               clusterType = 'mclapply',mc.cores = 1,max_steps = 100,h2_step = 0.01,method = 'REML')
(grid_AEcage_time = Sys.time() - start)
method_times$grid_AEcage = grid_AEcage_time
all_results$grid_AEcage = grid_AEcage$results$p_value_REML


# ---------- Grid-A-residuals ----------- #
trait_data$y_resid = resid(lmer(y~GENDER + (1|cage),trait_data))
start = Sys.time()
h2_start = get_h2_LDAK(trait_data$y_resid,X_cov,K_list[1],LDAK_path)
colnames(h2_start)[1] = 'ID'
grid_A_resid = GridLMM_GWAS_fast(y_resid~GENDER + (1|ID),~1,~0,data=trait_data,X = trait_X,X_ID = 'ID',relmat = list(ID = trait_A), h2_start = h2_start,#diagonalize = F,
                           clusterType = 'mclapply',mc.cores = 1,max_steps = 100,h2_step = 0.01,method = 'REML')
(grid_A_resid_time = Sys.time() - start)
method_times$grid_A_resid = grid_A_resid_time
all_results$grid_A_resid = grid_A_resid$results$p_value_REML


# ---------- Grid-AEcage-BF ----------- #
start = Sys.time()
null_model = GridLMM_posterior(y~GENDER + (1|ID)+(1|ID2) + (1|ID3),data = trait_data,relmat = list(ID = trait_A,ID2 = trait_E,ID3 = cage_K),save_V_folder = 'V_folder',h2_divisions = 3)
V_setup = null_model$V_setup
colSums(null_model$h2s_results*null_model$h2s_results$posterior)
(h2_step = min(diff(unique(null_model$h2s_results[,1]))))
h2_step
grid_AEcage_BF = GridLMM_GWAS_fast(y~GENDER + (1|ID)+(1|ID2) + (1|ID3),~1,~0,data=trait_data,X = trait_X,X_ID = 'ID',relmat = list(ID = trait_A,ID2 = trait_E,ID3 = cage_K),
                                     h2_start = null_model$h2s_results[order(null_model$h2s_results$posterior,decreasing = T)[1:12],1:3],V_setup = V_setup,
                                     inv_prior_X = 1,
                                     clusterType = 'mclapply',mc.cores = 1,max_steps = 100,h2_step = h2_step,method = 'BF')
(grid_AEcage_BF_time = Sys.time() - start)
method_times$grid_AEcage_BF = grid_AEcage_BF_time
all_results$grid_AEcage_BF = grid_AEcage_BF$results$BF*log10(exp(1))

saveRDS(method_times,file = 'WTCC_timings.rds')
saveRDS(all_results,file = 'WTCC_results.rds')
results_list = list(
  results_GEMMA = results_GEMMA,
  grid_A = grid_A,
  grid_AE = grid_AE,
  grid_Acage = grid_Acage,
  grid_AEcage = grid_AEcage,
  grid_A_resid = grid_A_resid,
  grid_AEcage_pylmm = grid_AEcage_pylmm,
  grid_AEcage_BF = grid_AEcage_BF
)
saveRDS(results_list,file = 'WTCC_all_results.rds')

# calculate genomic control values:
sapply(results_list,function(x) median(x$results$F.1)/qf(0.5,1,nrow(trait_data)-2))

