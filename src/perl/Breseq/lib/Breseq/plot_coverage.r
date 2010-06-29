## Args should be in_file=/path/to/input out_file=/path/to/output

pdf_output = 0;
total_only = 0;

for (e in commandArgs()) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
 #   temp = as.numeric(temp) #Im only inputting numbers so I added this to recognize scientific notation
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    assign(ta[[1]][1],temp)
    cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
  } else {
    assign(ta[[1]][1],TRUE)
    cat("assigned ",ta[[1]][1]," the value of TRUE\n")
  }
}


X<-read.table(in_file, sep="\t", header=T)
X$unique_tot_cov = X$unique_bot_cov + X$unique_top_cov;
X$redundant_tot_cov = X$redundant_bot_cov + X$redundant_top_cov;
maxy=max(X$unique_tot_cov, X$redundant_tot_cov) + 5;
start_pos = X$position[1];
end_pos = X$position[length(X$position)];

if (pdf_output == 0) {
	bitmap(out_file, height=450, width=900, type = "png16m", units = "px", res = 72, pointsize=18, taa=4, gaa=2)
} else {
	pdf(out_file, height=5, width=8)
}

par(mar=c(5.1,4.1,1.2,2));
plot(0:10, 0:10, type="n", lty="solid", ylim=c(0, maxy), xlim=c(start_pos, end_pos), lwd=2, xaxs="i", yaxs="i", xlab="Coordinate in Reference Genome", ylab="Read Coverage Depth")

#### Need to add back the option to gray out the ends!
#rect(pos[start], 0, pos[del_start], maxy, col="grey85", lty=0)
#rect(pos[del_end]+1, 0, pos[end], maxy, col="grey85", lty=0)

lines(X$position, X$redundant_tot_cov, type="s", col="red", lty="solid", lwd=1.5 )
if (total_only == 0)
{
	lines(X$position, X$redundant_top_cov, type="s", col="yellow", lty="solid", lwd=0.7)
	lines(X$position, X$redundant_bot_cov, type="s", col="orange", lty="solid", lwd=0.7)
}

lines(X$position, X$unique_tot_cov, type="s", col="blue", lty="solid", lwd=1.5 )
if (total_only == 0)
{
	lines(X$position, X$unique_top_cov, type="s", col="cyan", lty="solid", lwd=0.7)
	lines(X$position, X$unique_bot_cov, type="s", col="purple", lty="solid", lwd=0.7)
}
dev.off()
