### Karyotype counter v1 -- Andrew R Gross -- 2017-11-09 -- Imports a table of karyotype data to count abnormal lines.
### DO NOT USE -- INCOMPLETE

########################################################################
### Header
ref.id.list <- list()

########################################################################
### Import data

setwd('Z:/Data/Ideograms/')
karyotypes <- read.csv('Karyotypes_list--nov_9.csv')
#joined.dataframe <- cbind(data.frame("Line" = '')[1,],snp.less[1,])

########################################################################
### Format
snp.less <- snp.data[1:20]

### Filter: Select SNPs by coverage

quantile(snp.data$total.coverage, probs = seq(0, 1, 0.1))
#top.tenth.percentile <- quantile(snp.data$total.coverage, probs = seq(0, 1, 0.1))[10]
#top.20th.percentile <- quantile(snp.data$total.coverage, probs = seq(0, 1, 0.1))[9]
top.40th.percentile <- quantile(snp.data$total.coverage, probs = seq(0, 1, 0.1))[5]

pass.cov <- snp.data$total.coverage >=top.40th.percentile
snp.cov <- snp.less[pass.cov==TRUE,]
print(nrow(snp.cov))

### Select SNPs with missense mutations
snp.missense <- snp.cov[snp.cov$ExonicFunc == 'missense SNV',]
nrow(snp.missense)

### Add SNP IDs to list
ref.id.list[[file.selection]] <- as.vector(snp.missense$ID)
print(ref.id.list)

### List unique SNPs
rev(sort(ids.table <- table(unlist(ref.id.list))))
snp.list <- unique(unlist(ref.id.list))

### Generate data frame containing the SNPs of interest
rows <- match(snp.list, snp.less$ID)
rows <- rows[!is.na(rows)]
new.rows <- snp.less[rows,]
Line <- rep(current.file,nrow(new.rows))
new.rows <- cbind(Line,new.rows)
#new.rows[1] <- as.character(new.rows[1])

########################################################################
### Join all relevant SNPs into one df

### Define joined.dataframe the first time
#joined.dataframe <- new.rows

### Add to joined data frame
joined.dataframe <- rbind(joined.dataframe,new.rows)
nrow(joined.dataframe)

########################################################################
### Join all relevant SNPs into one df

write.csv(joined.dataframe,'../Relevant SNPs_40th_percentile.csv')

########################################################################
### Join all relevant SNPs into one df
### List unique SNPs
rev(sort(ids.table <- table(unlist(ref.id.list))))
common.snps <- names(rev(sort(ids.table <- table(unlist(ref.id.list))))[1:13])

### Generate data frame containing the SNPs of interest
common.snps.rows <- which(joined.dataframe$ID%in%  common.snps )
common.snps.rows <- common.snps.rows[!is.na(common.snps.rows)]
common.snps.df <- joined.dataframe[common.snps.rows,]

### Format and output
common.snps.df <- common.snps.df[c(1,4,5,2,3,6,20,8,7,9,10,11,13,21)]
common.snps.df <- common.snps.df[order(common.snps.df$GeneName),]
common.snps.df <- common.snps.df[order(common.snps.df$ID),]

print(common.snps.df[1:12])

write.csv(common.snps.df,'../Shared_SNPs.csv', row.names = FALSE)
