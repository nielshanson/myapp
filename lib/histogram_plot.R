#!/usr/bin/Rscript

# import libraries
library(ggplot2)
library(grid)

args <- commandArgs(trailingOnly = TRUE)

hist_file_in <- args[1]
output <- args[2]

con=file(hist_file_in,open="r")
lines=readLines(con)

my_parsing_function <- function(x, split="\t") {
    number_list <- strsplit(x, split=split)
    return(number_list[[1]])
}

x <- as.numeric(unlist(lapply(lines, my_parsing_function)))

p1 <- qplot(x, geom="histogram") + theme_bw()

png(output)
p1
dev.off()
