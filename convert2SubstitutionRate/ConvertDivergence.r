#!/usr/bin/env Rscript
library(dplyr)

Ne = 10000
gen = 20
args = commandArgs(trailingOnly=TRUE)

manifestFile <- args[1]
outputFile <- args[2]

compute_mu_div <- function(k, Ne=Ne, t=t, gen=gen){
  mu = k / (2*(t/gen) + 4*Ne)
  return(mu)
}
compute_mu_poly <- function(k, Ne=Ne, gen=gen){
  mu = k / (4*Ne)
  return(mu)
}

df_manifest <- read.delim(manifestFile, sep='')

l_df <- lapply(df_manifest[,2], function(X) read.delim(as.character(X), header=F, sep=''))

for (i in 1:dim(df_manifest)[1]){
  l_df[[i]]$type <- df_manifest[i,1]
  if(df_manifest[i,3] == 0){
    l_df[[i]]$mu <- compute_mu_poly(l_df[[i]]$V8, Ne)
  }else{
    l_df[[i]]$mu <- compute_mu_div(l_df[[i]]$V8, Ne, df_manifest[i,3] , gen)
  }
}

df <- bind_rows(l_df)

df_forGlennis <- df[,c(1,2,3,4,7,8,9,10)]
names(df_forGlennis) <- c("contig","start","end","comparison","fraction_GC","SubstitutionsPerSite","type","mu1")

write.table(df_forGlennis, outputFile, quote = F, row.names = F, sep="\t")
