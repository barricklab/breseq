#!/usr/bin/env Rscript
##*****************************************************************************
##
##AUTHORS
##
##  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
##  David B. Knoester
##
##LICENSE AND COPYRIGHT
##
##  Copyright (c) 2008-2010 Michigan State University
##  Copyright (c) 2011-2012 The University of Texas at Austin
##
##  breseq is free software; you can redistribute it and/or modify it under the  
##  terms the GNU General Public License as published by the Free Software 
##  Foundation; either version 1, or (at your option) any later version.
##
##*****************************************************************************
options <- commandArgs(trailingOnly = T)

table_path <- options[1]
output_path <- options[2]
cutoff <- 3

data <- read.table(table_path, sep = '\t', header = T)


x <- c(data$score)
y <- c(data$TP / (data$TP + data$FP))

png(output_path, height = 600, width = 800, bg = "white")

plot(x, y, col = "blue", ann = F)
lines(x, y)
abline(v = cutoff, col = "red", lty = 22)

title(main = "JC Precision versus Score")
title(xlab = "Score")
title(ylab = "Precision")



