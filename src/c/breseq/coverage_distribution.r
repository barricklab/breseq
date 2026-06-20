##
##
## AUTHORS
##
## Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
##
## LICENSE AND COPYRIGHT
##
## Copyright (c) 2008-2010 Michigan State University
## Copyright (c) 2011-2022 The University of Texas at Austin
##
## breseq is free software; you can redistribute it and/or modify it under the
## terms the GNU General Public License as published by the Free Software
## Foundation; either version 1, or (at your option) any later version.
##
##

## Draws the coverage distribution diagnostic plot. All of the statistics
## (negative binomial fit, deletion propagation cutoff, etc.) are computed in
## C++ -- see CoverageDistribution::fit in coverage_distribution.cpp -- and
## passed in as arguments below, rather than computed here.
##
## Arguments:
##   distribution_file=/path/to/input
##   plot_file=/path/to/output
##   nb_fit_mu=float     (0 if the fit failed)
##   nb_fit_size=float   (0 if the fit failed)
##   nb_fit_scale=float  (scales a raw dnbinom(...) curve to match the histogram's counts)
##   censor_start=int    (left edge of the window fit/censored around the peak)
##   censor_end=int      (right edge of the window fit/censored around the peak)

pdf_output = 1;

for (e in commandArgs(TRUE)) {
  ta = strsplit(e,"=",fixed=TRUE)[[1]]
  if(length(ta)>1) {
    temp = ta[2]
    assign(ta[1],temp)
  } else {
    assign(ta[[1]][1],TRUE)
  }
}

nb_fit_mu = as.numeric(nb_fit_mu);
nb_fit_size = as.numeric(nb_fit_size);
nb_fit_scale = as.numeric(nb_fit_scale);
censor_start = as.numeric(censor_start);
censor_end = as.numeric(censor_end);

#load data
X<-read.table(distribution_file, header=T)

#Nothing to plot if there's no data or no fit window was found (degenerate distribution).
if ((nrow(X) == 0) || (censor_start >= censor_end))
{
  q()
}

## graphics settings
my_pch = 21
my_col = "black";
my_col_censored = "red";

if (pdf_output == 0) {

  ## bitmap() requires ghostscript to be installed.
  ## taa=4, gaa=2 options NOT compatible with earlier R versions!
  ## units = "px" NOT compatible with even earlier R versions!

  if(!capabilities(what = "png"))
  {
    ## fallback to ghostscript
    bitmap(plot_file, height=6, width=7, type = "png16m", res = 72, pointsize=18)
  } else {
    ## use X11 function, which gives better resolution
    png(plot_file, height=6, width=7, units ="in", res = 72, pointsize=18)
    par(family="sans")
  }
} else {
  pdf(plot_file, height=6, width=7)
  par(family="sans")
}

par(mar=c(5.5,7.5,3,1.5));

## don't graph very high values with very little coverage
max_i = which.max(X$n)
max_n = X$n[max_i]
i<-max_i
while (i <= length(X$n) && X$n[i]>0.01*max_n)
{
	i <- i+1;
}
graph_end_i <-i

## This leaves enough room to the right of the peak for the legend
graph_end_i = max(floor(2.2 * max_i), graph_end_i);

fit_nb = c()
if (nb_fit_mu > 0) {
  fit_nb = dnbinom(0:max(X$coverage), mu = nb_fit_mu, size=nb_fit_size) * nb_fit_scale;
}

max_y = max(X$n, fit_nb)

plot(0:10, 0:10, type="n", lty="solid", ylim=c(0, max_y)*1.05, xlim=c(0, graph_end_i), lwd=1, xaxs="i", yaxs="i", axes=F, las=1, main="Coverage Distribution at Unique-Only Positions", xlab="Coverage depth (reads)", ylab="", cex.lab=1.2, cex.axis=1.2)

mtext(side = 2, text = "Number of reference positions", line = 5.5, cex=1.2)

sciNotation <- function(x, digits = 1) {
    if (length(x) > 1) {
        return(append(sciNotation(x[1]), sciNotation(x[-1])))
	}
    if (!x) return(0)

	exponent <- floor(log10(x))
    base <- round(x / 10^exponent, digits)
	as.expression(substitute(base %*% 10^exponent, list(base = base, exponent = exponent)))
}

axis(2, cex.lab=1.2, las=1, cex.axis=1.2, at = axTicks(2), labels = sciNotation(axTicks(2), 1))
axis(1, cex.lab=1.2, cex.axis=1.2, labels=T)
box()

#graph the coverage as points
fit_data <- subset(X, (coverage>=censor_start) & (coverage<=censor_end) );
points(fit_data$coverage, fit_data$n, pch=my_pch, col=my_col, bg="white", cex=1.2)

#graph the censored coverage as red points
censored_data <- subset(X, (coverage<censor_start) | (coverage>censor_end) );
points(censored_data$coverage, censored_data$n, pch=my_pch, col=my_col_censored, bg="white", cex=1.2)

#graph the negative binomial fit
if (nb_fit_mu > 0) {
  lines(0:max(X$coverage), fit_nb, lwd=3, col="black");
}

legend("topright", c("Coverage distribution", "Censored data", "Negative binomial"), lty=c("blank","blank","solid"), lwd=c(1,1,2), pch=c(my_pch, my_pch, -1), col=c("black", "red", "black"), bty="n")

dev.off()
