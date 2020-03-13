#!/usr/bin/env Rscript

library(readxl)
library(data.table)
library(plyr)

args = commandArgs(trailingOnly=TRUE)
file_name=args[1]
print(file_name)
file_ext=gsub(".xlsx","",file_name)
data=read_excel(file_name)
write.csv(data, paste0(file_ext,"_tmp.csv"))
data_tmp=fread(paste0(file_ext,"_tmp.csv"), sep=",", skip="Repeat Sequence")
data_final=data_tmp[,c(2,5,6)]
data_final=rename(data_final, c("Repeat Sequence"="Sequence"))
data_final=subset(data_final, Locus!="Amelogenin")
data_final$SampleID=data[2,2]
data_final$Project=data[3,2]
data_final$Analysis=data[4,2]
write.csv(data_final, paste0(file_ext,".csv"), row.names=F, quote=F)
