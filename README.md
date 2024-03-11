# NanoCMSer

An R package designed for subtyping FF / FFPE patient samples profiled by the NanoString platform. It's versatile and can also handle RNAseq and microarray datasets for CMS classification.

An example data set is provided to explain how to work with `NanoCMSer` R package and `NanoCMSer` function.

To install the `NanoCMSer` package use the following code:

```         
devtools::install_github("atorang/NanoCMSer")
```

To run the demo below, you need to first load the package and the dataset:

```         
library("NanoCMSer")
load("exprs.test")
```

Now you can utilize an elastic-net model for data classification. The function takes five arguments:

1.  `data`: A numeric matrix or data frame containing expression levels. It's recommended to use raw count data. Rows represent genes, and columns represent samples.

2.  `sample_type`: A character string specifying the type of samples. Accepted values are `tumorFF` for fresh frozen patient samples, `tumorFFPE` for formalin-fixed paraffin-embedded patient samples, and `models` for human in vitro models including cell lines, primary cultures, and organoids.

3.  `gene_names`: A character string specifying the gene annotation used in the data. Accepted values are `ensembl`, `symbol`, and `entrez`. The default value is `ensembl`, which is recommended to mitigate the risk of missing values due to gene symbol updates. The current version encompasses all previous versions of gene symbols up to 2024.

4.  `perform_log2`: A logical value determining whether data needs log2-transformation (`TRUE`) or if the data is already log2-transformed (`FALSE`). The default is `perform_log2 = TRUE`.

5.  `impute`: A logical value. If `impute = TRUE`, missing genes (up to 10% of the genes utilized in the classifiers) will be imputed using a trained linear regression model. If `impute = FALSE`, the classifier will generate an error message in the event of missing genes. The default is `impute = TRUE`.

```         
res <- NanoCMSer(data=exprs.test, 
                 sample_type="tumorFFPE", 
                 perform_log2=TRUE, 
                 gene_names="ensembl",
                 impute=TRUE)
```
