#!/usr/bin/env Rscript
if (!require("dplyr")) {
    install.packages("dplyr",repos = "http://cran.us.r-project.org")
}

library(dplyr)
library(tools)

mod_tss_peaks <- function(input, strand_s, cov_min = 1, merge_w = 5){
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

cluster_peaks <- function(inputFile,strand,clusterw,covm,out){
        counts <- read.table(inputFile)
        peaks <- mod_tss_peaks(counts,strand,merge_w=clusterw,cov_min=covm)
        write.csv(peaks, out)
}

inputFile=commandArgs(TRUE)[1]
strand=commandArgs(TRUE)[2]
clusterw=as.numeric(commandArgs(TRUE)[3])
covm=as.numeric(commandArgs(TRUE)[4])
out=commandArgs(TRUE)[5]
cluster_peaks(inputFile,strand,clusterw,covm,out)
