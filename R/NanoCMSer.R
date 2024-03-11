#' NanoCMSer
#'
#' \code{NanoCMSer} is a classifier developed for the NanoString platform to predict
#'  Consensus Molecular Subtypes (CMS) in colorectal cancer.
#'
#' @param data a numeric matrix or data.frame with expression levels. Rows are genes and columns are samples.
#' We recommend using raw count data.
#' @param perform_log2 if \code{perform_log2 = TRUE}, the data will be log2-transformed.
#' It is a logical determining if data needs log2-transformation or if the data is already log2-transformed.
#' Default is \code{perform_log2 = TRUE}.
#' @param gene_names a character string specifying the gene annotation used in \code{data}. Accepted values are
#' \code{"ensembl"}, \code{"symbol"} and \code{"entrez"}. Default value is \code{"ensembl"}.
#' We advise using Ensembl IDs to mitigate the risk of missing values caused by gene symbol updates.
#' The current version encompasses all previous versions of gene symbols up to 2024.
#' @param sample_type a character string. Accepted values are: \code{"tumorFF"}; fresh frozen patient samples,
#' \code{"tumorFFPE"}; formalin-fixed paraffin-embedded patient samples,
#' \code{"models"}; human in vitro models including cell lines, primary cultures and organoids.
#' @param impute if \code{impute = TRUE}, for missing genes, up to 10% of the genes utilized in the classifiers will be imputed
#' utilizing a trained linear regression model. If \code{impute = FALSE}, the classifier will generate an error message in the event of missing genes.
#' Default is \code{impute = TRUE}.
#' @return An object of class "\code{data.frame}" containing the classification results.
#' Columns CMS1-4 show the chance of being each subtype, based on an elastic-net model.
#' nearestCMS indicates subtype with highest chance, regardless of the cutoff.
#' predictedCMS, shows subtype of samples that confidently predicted (cuttoff = 0.6).
#' @export

#' @examples
#' #classify a predefined set
#' data("exprs.test")
#' NanoCMSer(data=exprs.test, sample_type="tumorFFPE", perform_log2=TRUE,
#' gene_names="ensembl",impute=TRUE)


NanoCMSer <- function(data,sample_type=NA, perform_log2=T,gene_names="ensembl",impute=T){
  requireNamespace("glmnet", quietly = TRUE)

  #validate inputs
  validate_sample_type(sample_type)
  validate_gene_names(gene_names)
  data=validate_data_na(data)
  data=perform_log2_transformation(data, perform_log2)

  #select models
  if(sample_type=="tumorFF"){
    exprs_model=exprs.FF
    genes_model=genes.FF
    cms_model=cms.FF
    Ens_model=Ens.FF
    model=model.FF
    imputer=imputer.FF
  }else if(sample_type=="tumorFFPE"){
    exprs_model=exprs.FFPE
    genes_model=genes.FFPE
    cms_model=cms.FFPE
    Ens_model=Ens.FFPE
    model=model.FFPE
    imputer=imputer.FFPE
  }else if(sample_type=="models"){
    exprs_model=exprs.CL
    genes_model=genes.CL
    cms_model=cms.CL
    Ens_model=Ens.CL
    model=model.CL
    imputer=imputer.CL
  }
  colnames(Ens_model)=c("symbol","ensembl","entrez","other")

  #update possible old gene symbol
  data=handle_old_gene_symbol(data,Ens_model,gene_names,exprs_model)
  data_genes=intersect(Ens_model[,gene_names],rownames(data))

  #check for percentage of missing genes and imputing option
  check_impute(data_genes, Ens_model, gene_names, exprs_model, impute)

  #order genes and change the row names
  data=data[Ens_model[,gene_names],]
  rownames(data)=rownames(exprs_model)

  #impute missing values
  if(sum(is.na(data))>0){data=impute_missing_values(data, imputer)}

  for (i in 1:ncol(data)) {
    comb=preprocessCore::normalize.quantiles(as.matrix(data.frame(exprs_model,data[,i])))
    data[,i]=comb[,ncol(comb)]
  }
  prediction=stats::predict(model, newx=t(data), interval ="prediction",type="response")
  cls=stats::predict(model, newx=t(data), interval ="prediction",type="class")
  Max=apply(prediction, 1, max)
  Predicted=cls
  Predicted[Max<0.6]<-NA
  results=data.frame(prediction,NearestCMS=cls,CMS=Predicted)
  colnames(results)=c("CMS1","CMS2","CMS3","CMS4","nearestCMS","predictedCMS")
  return(results)
}


# functions-------------------------
# Function to validate sample type argument
validate_sample_type <- function(sample_type) {
  if (!sample_type %in% c("tumorFF", "tumorFFPE", "models")) {
    stop("Set 'sample_type' argument from one of these values: 'tumorFF', 'tumorFFPE', 'models'.\nFor more information run '?NanoCMSer'")
  }
}

# Function to validate data
validate_data_na <- function(data) {
  data <- data.frame(data)
  if (sum(is.na(data)) > 0) {
    stop("There are missing values (NA) in the data")
  }
  for(i in 1:ncol(data)) {data[,i]=as.numeric(data[,i])}
  return(data)
}

# Function to validate gene names
validate_gene_names <- function(gene_names) {
  if(!gene_names%in%c("ensembl", "symbol", "entrez")){
    stop("'gene_names' can accept one of the following values:\n
         'ensembl', 'entrez', 'symbol'")
  }
}

# Function to perform log2 transformation
perform_log2_transformation <- function(data, perform_log2) {
  if (perform_log2) {
    warning("Log2-transformation is applied to the data. If the input data is already log2-transformed, set 'perform_log2 = FALSE'")
    data <- log2(data + 1)
  } else if (max(data) > 30) {
    data <- log2(data + 1)
    warning("Log2-transformation is applied to the data because values > 30 are observed")
  }
  return(data)
}

# Function to handle old gene ids
handle_old_gene_symbol<-function(data,Ens_model,gene_names,exprs_model){
  data_genes=intersect(Ens_model[,gene_names],rownames(data))
  if((gene_names=="symbol") & ((length(data_genes)<nrow(exprs_model)))){
    symb=unlist(apply(Ens_model[!is.na(Ens_model$other),], 1,
                      function(x){strsplit(x[4],', ')}))
    x=intersect(rownames(data),symb)
    if(length(x)>0){
      for(i in x){
        rownames(data)[which(rownames(data)==i)]=
          gsub('\\..*','',names(symb)[symb==i])
      }
    }
    data_genes=intersect(Ens_model[,gene_names],rownames(data))
  }
  return(data)
}

# Function to check missing genes
check_impute <- function(data_genes, Ens_model, gene_names, exprs_model, impute) {
  if(length(data_genes)<0.9*nrow(exprs_model)){
    stop(paste0("There is not enough overlapping genes (Maximum imputing power is 10% of genes),
         or rownames of data are not set to the right value\nSetting 'gene_names' argument to 'ensembl', 'entrez' or 'symbol' may solve the problem!"
                ,"Missing genes:\n",paste0(setdiff(Ens_model[,gene_names],data_genes), collapse = " ")))
  }
  if((length(data_genes)<nrow(exprs_model)) & (impute)){
    warning(paste0("There are some missing genes,
            the classification is less reliable with missing genes:","\n",
                   paste0( setdiff(Ens_model[,gene_names],data_genes),collapse = " ")))
  }else if((length(data_genes)<nrow(exprs_model)) & ((!impute))){
    stop(paste0("There are missing genes in the data,
                set 'impute=TRUE' to impute missing genes!"
                ,"Missing genes:\n",paste0(setdiff(Ens_model[,gene_names],data_genes), collapse = " ")))
  }
}

# Function to impute missing values
impute_missing_values <- function(data, imputer) {
  data[is.na(data)]=-1
  imputResult=stats::predict(imputer, newx = t(data), s = "lambda.min")
  imputResult=t(imputResult[, , 1])
  id=which(data[,1]==-1)
  data[id,]=imputResult[id,]
  return(data)
}
