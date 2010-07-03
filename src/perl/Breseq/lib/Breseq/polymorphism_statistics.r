## Args should be in_file=/path/to/input out_file=/path/to/output total_length=<total_length_of_sequences>  qual_file=

##error_count_file=/path/to/error_count


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
	log10_base_likelihood = 1:length(X$new_quals)
);

qual_dist<-read.table(qual_file, sep="\t", header=F)
print(qual_dist)
qual_dist<-as.vector(qual_dist$V1);
print(qual_dist)

#print(Y)

if (length(X$new_quals) > 0)
{
	for (i in 1:length(X$new_quals))
	{
		Y$log10_base_likelihood[i] = X$log10_base_likelihood[i]
	
		#print (i);
		new_quals_list <- strsplit(as.vector(X$new_quals[i]), ",");
		new_quals <- as.numeric( new_quals_list[[1]] )
		ref_quals_list <- strsplit(as.vector(X$ref_quals[i]), ",");
		ref_quals <- as.numeric( ref_quals_list[[1]] )


	## This code estimates the actual strand and quality score distribution as the total observed.
	
		max_qual = max(new_quals, ref_quals)
		NQ = tabulate(new_quals, nbins=max_qual)
		RQ = tabulate(ref_quals, nbins=max_qual)
		TQ = NQ+RQ
		
		log10_qual_likelihood = log10(dmultinom(NQ, prob=TQ)) + log10(dmultinom(RQ, prob=TQ)) - log10(dmultinom(TQ, prob=TQ)) 
		Y$log10_qual_likelihood_position_model[i] = -log10_qual_likelihood	
		
		RS = c(X$ref_top_strand[i], X$ref_bot_strand[i])
		NS = c(X$new_top_strand[i], X$new_bot_strand[i])
		TS = RS+NS
		
		log10_strand_likelihood =  log10(dmultinom(RS, prob=TS)) + log10(dmultinom(NS, prob=TS)) - log10(dmultinom(TS, prob=TS))
		Y$log10_strand_likelihood_position_model[i] = -log10_strand_likelihood

		#likelihoods are written as NULL model (one base) versus POLYMORPHISM model (mixed bases)
		#log10_base_likelihood should be negative
		#log10_qual_likelihood should be positive
		#log10_strand_likelihood should be positive

		#convert to natural logarithm and back to 
		score_combined_log =  -2* ( Y$log10_base_likelihood[i]  + Y$log10_qual_likelihood_position_model[i] + Y$log10_strand_likelihood_position_model[i]) * log(10)
		Y$quality_position_model[i] = -pchisq(score_combined_log, 1, lower.tail=T, log=T) / log(10)
		
	## This section estimates the actual strand distribution as the reference observations and the quality score distribution from the
	## count across the entire genome for all bases added together (doesn't take into account that there may be more low G's than A's, for example)	

		max_qual = length(qual_dist)
		NQ = tabulate(new_quals, nbins=max_qual)
		RQ = tabulate(ref_quals, nbins=max_qual)
		TQ = NQ+RQ
		EQ = RQ+1
						
		log10_qual_likelihood_genome_mode = log10(dmultinom(NQ, prob=EQ)) + log10(dmultinom(RQ, prob=EQ)) - log10(dmultinom(TQ, prob=TQ))
		Y$log10_qual_likelihood_genome_model[i] = -log10_qual_likelihood_genome_mode

		RS = c(X$ref_top_strand[i], X$ref_bot_strand[i])
		NS = c(X$new_top_strand[i], X$new_bot_strand[i])
		ES = RS+c(1,1)
		TS = RS + NS
		
		log10_strand_likelihood_genome_mode =  log10(dmultinom(RS, prob=ES)) + log10(dmultinom(NS, prob=ES)) - log10(dmultinom(TS, prob=TS))
		Y$log10_strand_likelihood_genome_model[i] = -log10_strand_likelihood_genome_mode

		#convert to natural logarithm and back to 
		score_combined_log_genome_mode =  -2* ( Y$log10_base_likelihood[i]  + Y$log10_qual_likelihood_genome_model[i] + Y$log10_strand_likelihood_genome_model[i]) * log(10)
		Y$quality_genome_model[i] = -pchisq(score_combined_log_genome_mode, 1, lower.tail=T, log=T) / log(10)


	## Oldest code that calculates bias p-values
	
		options(warn=-1);
		ks_test_p_value <- ks.test(ref_quals, new_quals, alternative = "less")
		options(warn=0);
		ks_test_p_value <- ks_test_p_value$p.value
		Y$ks_quality_p_value[i] <- ks_test_p_value

		contingency_table <- matrix(data=c(X$new_top_strand[i], X$new_bot_strand[i], X$ref_top_strand[i], X$ref_bot_strand[i]), nrow=2, ncol=2)
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

write.table(Y, out_file, sep="\t", row.names=F, quote=F)