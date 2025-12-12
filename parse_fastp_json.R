#! /usr/bin/Rscript

if (!require("jsonlite")) install.packages("jsonlite")
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
json_paths <- args

json_parse <- function(json_paths){  
  data <- fromJSON(json_paths)

  df <- data.frame(
    file = basename(json_paths),

    before_total_reads = as.numeric(data$summary$before_filtering$total_reads),
    before_total_bases = as.numeric(data$summary$before_filtering$total_bases),
    before_q30_bases = as.numeric(data$summary$before_filtering$q30_bases),
    before_q30_rate = as.numeric(data$summary$before_filtering$q30_rate),
    before_gc_content = as.numeric(data$summary$before_filtering$gc_content),

    after_total_reads = as.numeric(data$summary$after_filtering$total_reads),
    after_total_bases = as.numeric(data$summary$after_filtering$total_bases),
    after_q30_bases = as.numeric(data$summary$after_filtering$q30_bases),
    after_q30_rate = as.numeric(data$summary$after_filtering$q30_rate),
    after_read1_mean_length = as.numeric(data$summary$after_filtering$read1_mean_length),

    after_read2_mean_length = ifelse(length(data$summary$after_filtering$read2_mean_length)!=0,
                                     as.numeric(data$summary$after_filtering$read2_mean_length), NA),
    after_gc_content = as.numeric(data$summary$after_filtering$gc_content),

    passed_filter_reads = as.numeric(data$filtering_result$passed_filter_reads),
    low_quality_reads = as.numeric(data$filtering_result$low_quality_reads),
    too_short_reads = as.numeric(data$filtering_result$too_short_reads),

    duplication_rate = as.numeric(data$duplication$rate),
    stringsAsFactors = TRUE
  )
}

tables <- lapply(json_paths, json_parse)
combined.df <- do.call(rbind , tables)
write.csv(combined.df, "fastp_trimmed_summary.csv", row.names = FALSE)
