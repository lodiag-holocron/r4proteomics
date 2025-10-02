# Functional Analysis, Longitudinal & Public Data {#day4}

## Objectives

-   Interpret DE results biologically\
-   Explore longitudinal trajectories\
-   Work with public datasets

## Module 1: Functional Enrichment & GSEA

Introduce enrichment (GO, KEGG, Reactome).\
GSEA principles and visualization.


``` r
# example with clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
# Suppose `de_genes` is a vector of gene IDs with statistics
ego <- enrichGO(de_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")
dotplot(ego)
```


``` r
library(lme4)
#> Loading required package: Matrix
df_long <- data.frame(
  sample = rep(1:10, each = 3),
  time = rep(c(0,1,2), times = 10),
  intensity = rnorm(30)
)
m <- lmer(intensity ~ time + (1 + time | sample), data = df_long)
#> boundary (singular) fit: see help('isSingular')
summary(m)
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: intensity ~ time + (1 + time | sample)
#>    Data: df_long
#> 
#> REML criterion at convergence: 84.8
#> 
#> Scaled residuals: 
#>      Min       1Q   Median       3Q      Max 
#> -2.54442 -0.47295  0.00704  0.75598  1.33428 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev. Corr 
#>  sample   (Intercept) 0.23219  0.4819        
#>           time        0.06152  0.2480   -1.00
#>  Residual             0.88093  0.9386        
#> Number of obs: 30, groups:  sample, 10
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  -0.3296     0.3109  -1.060
#> time          0.2320     0.2241   1.036
#> 
#> Correlation of Fixed Effects:
#>      (Intr)
#> time -0.804
#> optimizer (nloptwrap) convergence code: 0 (OK)
#> boundary (singular) fit: see help('isSingular')
```

### Exercise

On a longitudinal proteomics dataset, plot trajectories per protein or per cluster, fit simple linear/mixed model.

## Module 2: Public Data Integration

Demonstrate downloading, cleaning, integrating a public proteomics or expression dataset (e.g. from PRIDE, GEO). Run pipeline: QC → normalization → DE → enrichment.

### Exercise

Pick a public dataset, import into R, preprocess, analyze, interpret results.
