#!/usr/bin/env Rscript
if (!require("dplyr")) {
    install.packages("dplyr")
}

library(dplyr)
library(tools)

mod_tss_peaks <- function(input, strand_s, cov_min = 5, merge_w = 5){
    input %>% dplyr::rename(chr = V1, start_peak = V2, end_peak = V3, 
                  prominence = V5, strand_peak = V6, width = V10,
                  start_cov = V12, end_cov = V13, cov = V14, width_cov = V15, mapped_reads = V16) %>%
    dplyr::select(-V4, -V7, -V8, -V11) %>%
    group_by(start_peak, end_peak) %>%
    dplyr::mutate(full_cov_peak = sum(cov))%>%
    dplyr::filter(cov == max(cov)) %>%
    dplyr::mutate(decision_v = ifelse(strand_s == "+", 
                                      min(end_cov), max(end_cov))) %>%
    dplyr::filter(end_cov == decision_v) %>%
    ungroup() %>%
    arrange(end_cov) %>%
    mutate(index = lag(end_cov, default = 1) + as.integer(merge_w),
           index1 = cumsum(ifelse(index >= start_cov, 0, 1))+1) %>%
    dplyr::group_by(index1) %>%
    dplyr::mutate(full_cov_clust = sum(full_cov_peak))%>%
    dplyr::filter(cov == max(cov),
                  cov >= cov_min)%>%
    dplyr::mutate(RPM = 1000000*full_cov_clust/mapped_reads)
}

cluster_peaks <- function(inputDirectory,clusterw){
    directory = paste(inputDirectory,"/*/*.normalized",sep="")
    fileNames <- Sys.glob(directory)
    for (fileName in fileNames) {
        if (grepl(".plus.", fileName, fixed=TRUE)){
            strand <- "+"
        } else if (grepl(".minus.", fileName, fixed=TRUE)){
            strand <- "-"
        }else{
            print("Can't decide between pos and neg")
        }
        counts <- read.table(fileName)
        peaks <- mod_tss_peaks(counts,strand,merge_w=clusterw)
        write.csv(peaks, paste(file_path_sans_ext(fileName),".clustered.csv", sep = ""), row.names = FALSE)
    }
}

directory=commandArgs(TRUE)[1]
clusterw=as.numeric(commandArgs(TRUE)[2])
cluster_peaks(directory,clusterw)
