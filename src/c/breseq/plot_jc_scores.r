#!/usr/bin/env Rscript
##
##
## AUTHORS
##
## Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
##
## LICENSE AND COPYRIGHT
##
## Copyright (c) 2008-2010 Michigan State University
## Copyright (c) 2011-2017 The University of Texas at Austin
##
## breseq is free software; you can redistribute it and/or modify it under the
## terms the GNU General Public License as published by the Free Software
## Foundation; either version 1, or (at your option) any later version.
##
##
options <- commandArgs(trailingOnly = T)
par(family="sans")

table_path <- options[1]
prefix <- options[2]
precision_path <- paste(prefix, ".precision.png", sep ="")
sensitivity_path <- paste(prefix, ".sensitivity.png", sep ="")
cutoff <- as.numeric(options[3])
cv_exe <- options[4]


data <- read.table(table_path, sep = '\t', header = T)

##Precision vs Score plot:
x <- c(data$score)
y <- c(data$TP / (data$TP + data$FP))

png(precision_path, height = 600, width = 800, bg = "white")

if (cv_exe == "tophat") {
	plot(x, y, col = "blue", ann = F, xlim = rev(range(x)), ylim = c(0,1))
} else {
	plot(x, y, col = "blue", ann = F, ylim = c(0,1))
}
lines(x, y)
title(main = "JC Precision versus Score")
title(xlab = "Score")
title(ylab = "Precision")

if (cutoff != 0 ) {
	abline(v = cutoff, col = "red", lty = 22)
}


##Sensitivity vs Score plot:
x <- c(data$score)
y <- c(data$TP / (data$TP + data$FN))

png(sensitivity_path, height = 600, width = 800, bg = "white")

if (cv_exe == "tophat") {
	plot(x, y, col = "blue", ann = F, xlim = rev(range(x)), ylim = c(0,1))
} else {
	plot(x, y, col = "blue", ann = F, ylim = c(0,1))
}
lines(x, y)
title(main = "JC Sensitivity versus Score")
title(xlab = "Score")
title(ylab = "Sensitivity")

if (cutoff != 0) {
	abline(v = cutoff, col = "red", lty = 22)
}

