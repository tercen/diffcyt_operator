# diffcyt

##### Description

Statistical methods for differential discovery analyses in high-dimensional cytometry data.

##### Usage

Input projection|.
---|---
`y-axis`        | numeric, cell count
`row`           | factor, cluster IDs 
`column`        | factor, sample IDs 
`colors`        | factor, group IDs (fixed effect)
`labels`        | factor, optional, patient / batch IDs (random effect)

Input parameters|.
---|---
`method`          | statistical method to be used (any of `DA_edgeR` or `DA_GLMM` for Differential Abundance, or `DS_limma` or `DS_LMM` for Differential State)
`reference.index` | Index of the reference category to be used. Default is 1, meaning that the first color specified in the crosstab will be used as a reference. In case more than 2 colors are present in the data, each of them will be compared to the reference group.

Output relations|.
---|---
`logFC`     | log fold change
`logCPM`    | log of counts per million
`LR`        | likelihood ratio
`p_val`     | p-value
`p_adj`     | adjusted p-value
`group_1`     | First group in the comparison
`group_2`     | Second group in the comparison

##### References

[diffcyt Bioconductor package](https://www.bioconductor.org/packages/release/bioc/html/diffcyt.html)

> Weber, L.M., Nowicka, M., Soneson, C. et al. diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. Commun Biol 2, 183 (2019). 

[Publication](https://doi.org/10.1038/s42003-019-0415-5)

