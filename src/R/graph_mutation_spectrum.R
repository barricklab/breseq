#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(tidyr)
library(dplyr)

option_list = list(
  make_option(c("-f", "--count"), type="character", default="count.csv", 
              help="count.csv file output by gdtools [default= %default]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="mutation_spectrum.pdf", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## load data
## X is the full table

X = read.csv(opt$count)



#plot defaults

plot.width = 7
plot.height = 5
theme_set(theme_bw(base_size = 12))
line_thickness = 0.8
theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))

##First graph just the ones at each point
unrolled = X %>% 
  select(population, base_substitution.synonymous, base_substitution.nonsynonymous, base_substitution.nonsense, base_substitution.noncoding, base_substitution.pseudogene, base_substitution.intergenic, small_indel, large_deletion, large_insertion, large_amplification, large_substitution, mobile_element_insertion, gene_conversion, inversion) %>%
  gather(mutation.type, count, base_substitution.synonymous:inversion)

unrolled$mutation.type = factor(unrolled$mutation.type, levels=rev(c("base_substitution.synonymous", "base_substitution.nonsynonymous", "base_substitution.nonsense", "base_substitution.noncoding", "base_substitution.pseudogene", "base_substitution.intergenic", "small_indel", "large_insertion", "large_amplification", "large_substitution", "large_deletion", "mobile_element_insertion", "gene_conversion", "inversion")))

unrolled=unrolled %>% subset(count!=0 | mutation.type == "base_substitution.nonsynonymous")



ggplot(unrolled, aes(x=population, y=count, fill=factor(mutation.type))) + geom_bar(stat="identity", position="stack") + scale_y_continuous(limits=c(0,10), expand=c(0,0), breaks = c(1:10))



ggsave(opt$output, height=plot.height, width=plot.width)
