## Args should be in_file=/path/to/input out_file=/path/to/output total_length=<total_length_of_sequences>


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
	ks_quality_p_value = 1:length(X$new_quals), 
	fisher_strand_p_value = 1:length(X$new_quals), 
	bias_p_value = 1:length(X$new_quals),
	bias_e_value = 1:length(X$new_quals)
);

#print(Y)

for (i in 1:length(X$new_quals))
{
	#print (i);
	new_quals_list <- strsplit(as.vector(X$new_quals[i]), ",");
	new_quals <- as.numeric( new_quals_list[[1]] )
	ref_quals_list <- strsplit(as.vector(X$ref_quals[i]), ",");
	ref_quals <- as.numeric( ref_quals_list[[1]] )
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

write.table(Y, out_file, sep="\t", row.names=F, quote=F)