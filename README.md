
# Scripts to replicate analyses in GridLMM paper.

### Necessary datasets:

## Figure 1: Arabidopsis phenotypes

[genotype data](https://github.com/Gregor-Mendel-Institute/atpolydb/blob/master/250k_snp_data/call_method_32.tar.gz)

[phenotype data](https://github.com/Gregor-Mendel-Institute/atpolydb/blob/master/miscellaneous_data/phenotype_published_raw.tsv)

[Atwell_one_trait.R](Atwell_one_trait.R)

[Atwell_GxE.R](Atwell_GxE.R)

[Power_simulations_Atwell_Wald.R](Power_simulations_Atwell_Wald.R)

## Figure 2: WTCC body weight

Data is provided in the **BGLR** package

[Genetic map](https://journals.plos.org/plosbiology/article/file?type=supplementary&id=info:doi/10.1371/journal.pbio.0040395.st005)

[WTCC_BW.R](WTCC_BW.R)

## Figure 3: Arabidopsis gene expression

[Expression data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80744)

[Kinship matrix](http://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5)

[At_gene_expression.R](At_gene_expression.R)

## External software used

[GEMMA](https://github.com/genetics-statistics/GEMMA)

[LDAK](http://dougspeed.com/ldak/)