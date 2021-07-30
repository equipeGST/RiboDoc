# Autors : GST team
# This script finds CDS and UTR coverage data when available and makes the corresponding graph

library(stringr)

CDS_total=c()
five_prime_total=c()
three_prime_total=c()

args = commandArgs(trailingOnly=TRUE)
project_name <- args[1]

# Reading the data file
params <- scan(file = paste0("/data/RESULTS/", project_name, ".Analysis_Report.txt"),
               what = "character"
)

CDS <- as.numeric(gsub(" ", "", params[which(params=="CDS")+2], fixed = TRUE))
five_prime <- as.numeric(gsub(" ", "", params[which(params=="5prime")+2], fixed = TRUE))
three_prime <- as.numeric(gsub(" ", "", params[which(params=="3prime")+2], fixed = TRUE))

CDS_mean=mean(CDS)
five_prime_mean=mean(five_prime)
three_prime_mean=mean(three_prime)

tiff(filename = "/data/RESULTS/CDS_vs_UTR.tiff")
barplot(c(CDS_mean,five_prime_mean,three_prime_mean),
        names.arg = c("CDS","5primeUTR","3primeUTR"),
        xlab = "mRNA regions",
        ylab = "RPKM mean"
)
dev.off()