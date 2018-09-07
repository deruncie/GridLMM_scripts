library(data.table)
# genotype data: https://github.com/Gregor-Mendel-Institute/atpolydb/blob/master/250k_snp_data/call_method_32.tar.gz
# phenotype data: https://github.com/Gregor-Mendel-Institute/atpolydb/blob/master/miscellaneous_data/phenotype_published_raw.tsv

# extract genotype data into a matrix
X = as.data.frame(fread('call_method_32/Network/Data/250k/db/dataset/call_method_32.b',skip=1,header = T))
map = X[,1:2]
X = X[,-c(1:2)]
rownames(X) = paste(map$Chromosome,map$Positions,sep='_')
X = apply(as.matrix(X),1,function(x) as.numeric(as.factor(x))) - 1
X = apply(X,2,function(x) {if(mean(x)<0.5) {return(x)};return(1-x)})

phen = fread('phenotype_published_raw.tsv')

dataset = list(
  X = X,
  map = map,
  phenotypes = phen
)
saveRDS(dataset,file = 'Data/dataset.rds')