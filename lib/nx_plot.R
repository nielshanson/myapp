#!/usr/bin/Rscript

# import libraries
library(ggplot2)
library(grid)

args <- commandArgs(trailingOnly = TRUE)

nx_file_in <- args[1]
output <- args[2]

nx_df <- read.table(nx_file_in, head=T, sep="\t")

p1 <- ggplot(nx_df, aes(x=N_x, y=value, group=Assembler))
p1 <- p1 + geom_line(size=2, color="blue")
p1 <- p1 + geom_point(size=2)

png(output)
p1
dev.off()
