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

## Args should be in_file=/path/to/input out_file=/path/to/output total_length=<total_length_of_sequences>  qual_file=

##error_count_file=/path/to/error_count
par(family="sans")


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


total_length = as.numeric(total_length);


X<-read.table(in_file, sep="\t", header=T)
#print(X)

##allocate output data frame
Y<-data.frame(
  ks_quality_p_value = 1:length(X$major_frequency)
);

qual_dist<-read.table(qual_file, sep="\t", header=F)
print(qual_dist)
qual_dist<-as.vector(qual_dist$V1);

qual_dist = ceiling(qual_dist / 1000)

print(qual_dist)


qual_dist_list = c();
for (i in 1:length(qual_dist))
{
	qual_dist_list = c(qual_dist_list, rep(i, qual_dist[i]));
}

qual_dist = qual_dist / sum(qual_dist)
qual_cdf = c();

pr_sum = 0;
for (i in 1:length(qual_dist))
{
	pr_sum = pr_sum + qual_dist[i];
	qual_cdf[i] = pr_sum;		
}


print(qual_cdf)

qual_pr <- function(qual_cdf = qual_cdf)
{
	if (x > length(qual_cdf))
	{
		return(1)
	}
	return (qual_cdf[trunc(x)]);
}

#print(Y)

if (length(X$major_quals) > 0)
{
	for (i in 1:length(X$major_quals))
	{
		#print (i);
    major_quals_list = c()
    major_quals = c()
    #cat(i, " length: ", length(X$major_quals[i]), "\n")
    if (length(X$major_quals[i]) > 0) {
      major_quals_list <- strsplit(as.character(X$major_quals[i]), ",");
      major_quals <- as.numeric( major_quals_list[[1]] )
    }
    
    minor_quals_list = c()
    minor_quals = c()
    
    #cat(i, " length: ", length(X$minor_quals[i]), "\n")
    if (length(X$minor_quals[i]) > 0) {
      minor_quals_list <- strsplit(as.character(X$minor_quals[i]), ",");
      minor_quals <- as.numeric( minor_quals_list[[1]] )
    }

    ##cat("\n")
    ##print(minor_quals)
        
    ## This code estimates the actual strand and quality score distribution as the total observed.
#		max_qual = max(major_quals, minor_quals)
#		NQ = tabulate(major_quals, nbins=max_qual)
#		RQ = tabulate(minor_quals, nbins=max_qual)
#		TQ = NQ+RQ
#
#		log10_qual_likelihood = log10(dmultinom(NQ, prob=TQ)) + log10(dmultinom(RQ, prob=TQ)) - log10(dmultinom(TQ, prob=TQ))
#		Y$log10_qual_likelihood_position_model[i] = -log10_qual_likelihood
		
#		RS = c(X$major_top_strand[i], X$major_bot_strand[i])
#		NS = c(X$minor_top_strand[i], X$minor_bot_strand[i])
#		TS = RS+NS
		
#		log10_strand_likelihood =  log10(dmultinom(RS, prob=TS)) + log10(dmultinom(NS, prob=TS)) - log10(dmultinom(TS, prob=TS))
#		Y$log10_strand_likelihood_position_model[i] = -log10_strand_likelihood

		#likelihoods are written as NULL model (one base) versus POLYMORPHISM model (mixed bases)
		#log10_base_likelihood should be negative
		#log10_qual_likelihood should be positive
		#log10_strand_likelihood should be positive

		#convert to natural logarithm and back to log10
    #score_combined_log =  -2* ( Y$log10_base_likelihood[i]  + Y$log10_qual_likelihood_position_model[i] + Y$log10_strand_likelihood_position_model[i]) * log(10)
    #Y$quality_position_model[i] = -pchisq(score_combined_log, 1, lower.tail=T, log=T) / log(10)
		
	## This section estimates the actual strand distribution as the reference observations and the quality score distribution from the
	## count across the entire genome for all bases added together (doesn't take into account that there may be more low G's than A's, for example)	

#		max_qual = length(qual_dist)
#		NQ = tabulate(best_quals, nbins=max_qual)
#		RQ = tabulate(second_best_quals, nbins=max_qual)
#		TQ = NQ+RQ
#		EQ = RQ+1
						
#		log10_qual_likelihood_genome_mode = log10(dmultinom(NQ, prob=EQ)) + log10(dmultinom(RQ, prob=EQ)) - log10(dmultinom(TQ, prob=TQ))
#		Y$log10_qual_likelihood_genome_model[i] = -log10_qual_likelihood_genome_mode

#		RS = c(X$ref_top_strand[i], X$ref_bot_strand[i])
#		NS = c(X$new_top_strand[i], X$new_bot_strand[i])
#		TS = RS + NS
		
#		log10_strand_likelihood_genome_mode =  log10(dmultinom(TS, prob=TS+1)) - log10(dmultinom(RS, prob=RS+1)) + log10(dmultinom(NS, prob=RS+1))
#		Y$log10_strand_likelihood_genome_model[i] = -log10_strand_likelihood_genome_mode

#		#convert to natural logarithm and back to 
#		score_combined_log_genome_mode =  -2 * ( Y$log10_base_likelihood[i] + Y$log10_strand_likelihood_genome_model[i]) * log(10)

#		Y$quality_genome_model[i] = -pchisq(score_combined_log_genome_mode, 1, lower.tail=T, log=T) / log(10)
#		score_combined_log_genome_mode =  -2 * ( Y$log10_base_likelihood[i]) * log(10)
		

		## KS test for unusual qualities -- Rewrite to use cumulative distribution! = faster and more accurate
		all_quals = c(minor_quals,major_quals);

		options(warn=-1);

		if (X$major_frequency[i] < 0.5)
		{
			poly_quals = major_quals;
		}
		else
		{
			poly_quals = minor_quals;
		}

## This one tests the quality scores predicting a polymorphism agains the overall distribution
## at all positions in the genome	
#	ks_test_p_value_unusual <- ks.test(qual_dist_list, poly_quals, alternative = "less")
#	ks_test_p_value_unusual <- ks_test_p_value_unusual$p.value
#	Y$ks_quality_p_value_unusual_poly[i] <- ks_test_p_value_unusual


#	ks_test_p_value_unusual <- ks.test(qual_dist_list, second_best_quals, alternative = "less")
#	ks_test_p_value_unusual <- ks_test_p_value_unusual$p.value
#	Y$ks_quality_p_value_unusual_ref[i] <- ks_test_p_value_unusual
	
#	ks_test_p_value_unusual <- ks.test(qual_dist_list, best_quals, alternative = "less")
#	ks_test_p_value_unusual <- ks_test_p_value_unusual$p.value
#	Y$ks_quality_p_value_unusual_new[i] <- ks_test_p_value_unusual

#	ks_test_p_value_unusual <- ks.test(qual_dist_list, all_quals, alternative = "less")
#	ks_test_p_value_unusual <- ks_test_p_value_unusual$p.value
#	Y$ks_quality_p_value_unusual_all[i] <- ks_test_p_value_unusual

	## Oldest code that calculates bias p-values
	
    options(warn=-1);
    
    ks_test_p_value = 1
    if ((length(minor_quals) > 0) && (length(major_quals) > 0)) {
      ks_test_result <- ks.test(minor_quals, major_quals, alternative = "less")
      ks_test_p_value <- ks_test_result$p.value
    }
    options(warn=0);
    
		Y$ks_quality_p_value[i] <- ks_test_p_value

		contingency_table <- matrix(data=c(X$minor_top_strand[i], X$minor_bot_strand[i], X$major_top_strand[i], X$major_bot_strand[i]), nrow=2, ncol=2)
		#print(contingency_table)
	
		fisher_test_p_value <- fisher.test( contingency_table, alternative="two.sided")
		fisher_test_p_value <- fisher_test_p_value$p.value
		#print(fisher_test_p_value)
	
		Y$fisher_strand_p_value[i] <- fisher_test_p_value
	
		# Fisher's method for combining p-values
		combined_log = - 2* ( log(ks_test_p_value) + log(fisher_test_p_value) )
		Y$bias_p_value[i] = pchisq(combined_log, 2*2, lower.tail=F)
	
		Y$bias_e_value[i] = Y$bias_p_value[i] * total_length
	}
}

Y = signif(Y, digits = 6)
write.table(Y, out_file, sep="\t", row.names=F, quote=F)
