#!/usr/bin/env Rscript
#This script was provided by Daniel Straub (QBIC)

# benchmark_taxonomy.R
#
# The script expects the following arguments: expected_taxa_file.tsv, calculated_taxa_file.tsv

# WARNING: taxa have to be unique at each level! I.e. Streptococcus;longus & Herbespirillum;longus will be counted as one species (longus), optimally: Streptococcus;Streptococcus_longus & Herbespirillum;Herbespirillum_longus !

# Get params and files from the command line
args            <- commandArgs(trailingOnly=TRUE)
expFILE         <- args[1] # expected taxa*
resFILE         <- args[2] # calculated taxa*
taxaseparator   <- args[3] # string that separates taxonomic level, default "; "
mergedtaxasep   <- args[4] # string that separates taxonomic level, default "\\|" (i.e. "|")
fbeta           <- args[5] # Fbeta weight, default 2; 1=precision and recall equally weighted (=F1 score), 2=weighs recall higher than precision, 0.5=weighs recall lower than precision.
prefix          <- args[6] # prefix string for output files
#*: first column with taxonomic strings, following none or many columns (=samples) with numeric values (=abundance), only presence/absence are used here

# Read params, use defaults if not provided
if ( is.na(taxaseparator) ) { taxaseparator <- "; " }
if ( is.na(mergedtaxasep) ) { mergedtaxasep <- "\\|" }
if ( is.na(fbeta) )         { fbeta <- 2 }
if ( is.na(prefix) )        { prefix <- "" }
fbeta <- as.numeric(fbeta)

# Input

# Read stuff
exp = read.table( expFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
res = read.table( resFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)

colnames(exp)[1] <- "tax"
colnames(res)[1] <- "tax"

# function to produce statistics per taxonomic level
get_stats <- function(s_exp,s_res,df,taxaseparator,sample) {
	# ignore here how many taxa levels, just do 10 max for now
	for (i in 1:10) {
		i_exp <- unique( sapply(strsplit(s_exp, taxaseparator, fixed = TRUE), `[`, i) )
		i_res <- sapply(strsplit(s_res, taxaseparator, fixed = TRUE), `[`, i)
		i_res <- unique( unlist( strsplit(i_res, mergedtaxasep, fixed = FALSE) ) )
		
		# drop empty entries (NA and "")
		i_exp <- unlist( lapply(i_exp, function(z){ z[!is.na(z) & z != ""]}) )
		i_res <- unlist( lapply(i_res, function(z){ z[!is.na(z) & z != ""]}) )
	
		# if none are expected, skip!
		if( length(i_exp)==0 ) {
			print( paste("No expected taxa in sample",sample,"level",i, "- skipping") ); next
		} else {
			print( paste(length(i_exp),"expected taxa in sample",sample,"level",i) )
		}
		
		# stats
		TP <- intersect(i_exp, i_res)
		FN <- setdiff(i_exp, i_res)
		FP <- setdiff(i_res, i_exp)
		
		# save
		types <- c(
			"detected",
			"expected",
			"TP",
			"FN",
			"FP",
			"recall",
			"precision",
			"F1",
			"Fbeta",
			"TP_taxa",
			"FN_taxa",
			"FP_taxa")
		values <- c(
			length(i_res),
			length(i_exp),
			length(TP),
			length(FN),
			length(FP),
			length(TP)/length(i_exp),
			length(TP)/length(i_res),
			( 2 * length(TP)/length(i_res) * length(TP)/length(i_exp) ) / ( length(TP)/length(i_res) + length(TP)/length(i_exp) ),
			( (1+fbeta^2) * length(TP)/length(i_res) * length(TP)/length(i_exp) ) / ( (fbeta^2) * length(TP)/length(i_res) + length(TP)/length(i_exp) ),
			paste(TP, collapse=','),
			paste(FN, collapse=','),
			paste(FP, collapse=',') )
		ids <- rep(sample, length(types))
		levels <- rep(i, length(types))
		df_append <- data.frame(
			sample=ids,
			level=levels,
			type= types,
			value= values,  
			stringsAsFactors=FALSE)
		df <- rbind( df, df_append )
	}
	return(df)
}

# Output dataframe
df <- data.frame(
	sample=character(),
	level=character(),
	type=character(),
	value=character(),  
	stringsAsFactors=FALSE) 

if( ncol(exp)>=2 && ncol(res)>=2 ) {
	# Find common sample IDs
	exp_samples <- colnames(exp)[2:ncol(exp)]
	res_samples <- colnames(res)[2:ncol(res)]
	samples <- intersect(exp_samples, res_samples)

	if( length(samples)>0 ) { 
		print(paste("Found",length(samples),"shared expected and measured samples:"))
		print(paste(samples,collapse=","))
		# for each shared sample
		for (sample in samples) {

			# keep only columns with taxa and sample abundances
			keep_cols <- c("tax",sample)
			s_exp <- subset(exp, select = keep_cols)
			s_res <- subset(res, select = keep_cols)

			# keep only detected taxa (abundance > 0)
			s_exp = s_exp[s_exp[,2] > 0,]
			s_res = s_res[s_res[,2] > 0,]

			# for each taxonomic level
			df <- get_stats(s_exp$tax, s_res$tax, df, taxaseparator, sample)
		}
	} else {
		# if there isnt any shared sample
		print("Found no shared expected and measured samples, compared all taxonomic strings.")
		print(paste( "Expected samples:", paste(exp_samples,collapse=",")))
		print(paste( "Measured samples:", paste(res_samples,collapse=",")))
		df <- get_stats(exp$tax, res$tax, df, taxaseparator, "unknown")
	}
} else {
	print("Found no abundance information in at least one file, compared all taxonomic strings.")
	df <- get_stats(exp$tax, res$tax, df, taxaseparator, "unknown")
}

# Write detailed output
outfile <- paste0(prefix,"benchmarking_stats_long.tsv")
print(paste("write",outfile))
write.table(df, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Make & write summary output
df_sum <- data.frame(level=character(), type=character(), mean=character(), stringsAsFactors=FALSE) 
outfile <- paste0(prefix,"benchmarking_stats_mean.tsv")
print("Warnings are expected at that point:")
for (type in unique(df$type)) {
	for (level in unique(df$level)) {
		df_subset <- df[df$type == type & df$level == level, ]
		MEAN <- mean( as.numeric(df_subset$value) )
		df_sum[nrow(df_sum) + 1,] <- list(level, type, MEAN)
	}
}
print(paste("write",outfile))
write.table(df_sum, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Plot Values
df_subset <- subset(df, type == "recall")
df_subset$value <- as.numeric(df_subset$value)
svg(paste0(prefix,"benchmark_recall.svg"), height = 4, width = 4)
boxplot(value~level, data=df_subset, xlab="Taxonomic level", ylab="Recall", ylim = c(0, 1))
invisible(dev.off())

df_subset <- subset(df, type == "precision")
df_subset$value <- as.numeric(df_subset$value)
svg(paste0(prefix,"benchmark_precision.svg"), height = 4, width = 4)
boxplot(value~level, data=df_subset, xlab="Taxonomic level", ylab="Precision", ylim = c(0, 1))
invisible(dev.off())

df_subset <- subset(df, type == "F1")
df_subset$value <- as.numeric(df_subset$value)
svg(paste0(prefix,"benchmark_F1.svg"), height = 4, width = 4)
boxplot(value~level, data=df_subset, xlab="Taxonomic level", ylab="F1 score", ylim = c(0, 1))
invisible(dev.off())

df_subset <- subset(df, type == "Fbeta")
df_subset$value <- as.numeric(df_subset$value)
svg(paste0(prefix,"benchmark_F",fbeta,".svg"), height = 4, width = 4)
boxplot(value~level, data=df_subset, xlab="Taxonomic level", ylab=paste0("F",fbeta," score"), ylim = c(0, 1))
invisible(dev.off())

