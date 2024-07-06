library(MatrixEQTL)

#analysis_mode='ewas'
#analysis_mode='kp'
analysis_mode='cis_trans_meth'
#analysis_mode='cis_trans_metab'

# from https://slowkow.com/notes/ggplot2-qqplot/
#
#' Create a quantile-quantile plot with ggplot2.
#'
#' Assumptions:
#'   - Expected P values are uniformly distributed.
#'   - Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @param ps Vector of p-values.
#' @param ci Size of the confidence interval, 95% by default.
#' @return A ggplot2 plot.
#' @examples
#' library(ggplot2)
#' gg_qqplot(runif(1e2)) + theme_grey(base_size = 24)
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}


# begin analysis here


## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste("genotypes_t.tsv", sep="");
SNP_file_name


# Gene expression (methylation) file name
if(analysis_mode=='ewas' || analysis_mode=='cis_trans_metab'){
  if(analysis_mode=='ewas'){
    expression_file_name = paste("methylation_t.tsv", sep="");
  }else if(analysis_mode=='cis_trans_metab'){
    expression_file_name = paste("metabolites_t.tsv", sep="");
  }
}else if(analysis_mode=='kp'){
  expression_file_name = paste("methylation_nosnp_probes_t.tsv", sep="");
}else if(analysis_mode=='cis_trans_meth'){
  snps_location_file_name = paste("../snp_map.txt", sep="");
  gene_location_file_name = paste("../gene_map.txt", sep="")
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  expression_file_name = paste("methylation_nosnp_probes_t.tsv", sep="");
}


# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste("covariates_t.tsv", sep="");
covariates_file_name

# Output file name
if(analysis_mode=='cis_trans_meth'){
  output_file_name_cis = tempfile();
  output_file_name_tra = tempfile();
}else{
  output_file_name = tempfile();
}

# Only associations significant at this level will be saved
if(analysis_mode=='ewas'|| analysis_mode=='cis_trans_metab'){
  pvOutputThreshold = 1e-8;
}else if(analysis_mode=='kp'){
  pvOutputThreshold = 1e-6;
}else if(analysis_mode=='cis_trans_meth' ){
  pvOutputThreshold_cis = 2e-6;
  #pvOutputThreshold_cis = 1;
  pvOutputThreshold_tra = 1e-6;
  cisDist = 1e6
}

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}
cvrt

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);
snps

## Load gene methylation data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);
gene


## Run the analysis

if(analysis_mode=='cis_trans_meth'){
  me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
}else{
  me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
}


if(analysis_mode=='cis_trans_meth'){
  unlink(output_file_name_cis);
  unlink(output_file_name_tra);
}else{
  unlink(output_file_name);
}

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Number of Tests:', me$all$ntests, '\n');
cat('Detected meQTLs:', me$all$neqtls, '\n');
cat('Detected meQTLs (FDR<0.05):', sum(me$all$eqtls$FDR<0.05), '\n');
cat('Detected meQTLs (FDR<0.001):', sum(me$all$eqtls$FDR<0.001), '\n');

show(me$all$eqtls[1:25,-3])

names(me)

## Plot the histogram of all p-values
plot(me)

save(me, file=paste("me",analysis_mode,"RData",sep="."))
if(analysis_mode=='cis_trans_meth' ){
  write.table(me$cis$eqtls,file='cis_eqtls.raw',quote=F,row.names=F,sep='\t')
  write.table(me$trans$eqtls,file='trans_eqtls.raw',quote=F,row.names=F,sep='\t')
}else if(analysis_mode=='kp'|| analysis_mode=='cis_trans_metab'){
  write.table(me$all$eqtls,file='all_eqtls.raw',quote=F,row.names=F,sep='\t')
}
#ps<-me$all$eqtls$pvalue

#library(ggplot2)

#gg_qqplot(ps) + theme_bw(base_size = 8) + theme( axis.ticks = element_line(size = 0.5), panel.grid = element_blank() 
# panel.grid = element_line(size = 0.5, color = "grey80")
#) 
#ggsave(filename='qqplot_all_pvalues.png',plot=last_plot(),device='png',path = NULL, scale = 1, width = 1000, height = 1000, units = "px", dpi = 300, limitsize = TRUE, bg = NULL)
