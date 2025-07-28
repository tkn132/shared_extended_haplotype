#!/usr/bin/env Rscript

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
		  stop("Usage: Rscript run_extended_haplotype_sharing_pipeline.R <input_file> <phased_data> <output_file")
}
input_file <- args[1]
phased_data <- args[2]
output_file <- args[3]

# Load libraries
library(vcfR)
library(ggplot2)
library(tidyverse)


# Function to count mismatches
count_mismatches <- function(seq1, seq2) {
		if (nchar(seq1) != nchar(seq2)) {
					  stop("Sequences must be the same length")
	}
	a <- strsplit(seq1, "")[[1]]
		b <- strsplit(seq2, "")[[1]]
		sum(a != b)
}

# Function to extract haplotypes for each sample
get_haplos <- function(sample_id) {
		  sample_genos <- haplo_mat[, sample_id]  # phased genotypes across region
  hap1 <- sapply(strsplit(sample_genos, "\\|"), `[`, 1)  # left haplotype
      hap2 <- sapply(strsplit(sample_genos, "\\|"), `[`, 2)  # right haplotype
            haps = c(paste(hap1, collapse = ""), paste(hap2, collapse = ""))
	            names(haps) = paste0(sample_id, "-", "hap", 1:2)
	            return(haps)
}

# Function to Extract the alt-carrying haplotype for each carrier:
get_alt_haplo <- function(sample_id) {
		  gt <- mutation_genotypes[[sample_id]]
  sample_genos <- haplo_mat[, sample_id]  # phased genotypes across region
      alleles <- strsplit(gt, "\\|")[[1]]
      # Pick the haplotype that contains the alternate allele (1)
      if (alleles[1] == "1") {
	      	        hap <- sapply(strsplit(sample_genos, "\\|"), `[`, 1)  # left haplotype
            } else if (alleles[2] == "1") {
		    	          hap <- sapply(strsplit(sample_genos, "\\|"), `[`, 2)  # right haplotype
	            } else {
			    		    return(NA)  # skip non-carriers
		    	  }
            return(paste(hap, collapse = ""))
            }

# Read the file snp list
list= read.table(input_file, header=T)
head(list)


# Loop over each of the snp and perform extended shared haplotype analysis

snps = c()
genes = c()
windows = c()
nshares = c()
totalcarriers = c()
totalhomos = c()
chrss = c()
poss = c()
nonhapmin = c()


for (line in 1:nrow(list)){
		print(list[line,])
    gene=as.character(list[line,4])
        mutation_pos=as.integer(list[line,3])
        SNP=as.character(list[line,1])

	      chrss = c(chrss, as.integer(list[line,2]))
	      poss = c(poss, mutation_pos)
	              snps = c(snps, SNP)
	              genes = c(genes, gene)
		      	  # Load VCF
		      	  vcf <- read.vcfR(paste0(phased_data,"/",gene, ".",SNP,".5Mbp.phased.vcf.gz"))

		      	  # Extract genotype dosage (0/1/2) for the variant
		      	  geno_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
			  	    # Find the row index of your ABCC8 variant
			  	    variant_pos <- which(vcf@fix[, "POS"] == mutation_pos)
			  	    # Get carrier status per sample
			  	    mutation_dosage <- geno_matrix[variant_pos, ]
				    	      table(mutation_dosage)
				    	      total_carriers = sum(table(mutation_dosage)[-1])
					      	        totalcarriers = c(totalcarriers , total_carriers)
					      	        print(paste0("total carriers for ", gene, " is ",total_carriers))
									  totalhomos = c(totalhomos, sum(mutation_dosage == "1|1"))

									  # Extract phased genotypes
									  haplos <- extract.gt(vcf, element = "GT")
									  		    # Focus on carriers (of the haploblocks)
									  		    carriers <- names(which(mutation_dosage == "1|0" | mutation_dosage == "0|1" | mutation_dosage == "1|1"))
									  		    # Identify the  variant row
									  		    mutation_genotypes <- haplos[variant_pos, carriers]


											    		      # Define window size
											    		      out = c()
											    		      window = c(1, seq(50,1000,50))
													      		        n=0
													      		        # Loop over all the window sizes to extract the maximum number of shared haplotypes per window size
													      		        for (s in window){
																					    print(s)
																			    pos <- as.numeric(vcf@fix[, "POS"])
																			    			        variant_pos <- which(pos == mutation_pos)
																			    			        window_size <- s*1000
																											  # Get region indices
																											    region_idx <- which(pos >= (mutation_pos - window_size) & pos <= (mutation_pos + window_size))
																											    # Subset phased genotypes
																											    haplo_mat <- extract.gt(vcf[region_idx, ], element = "GT")  # SNPs x samples
																											    				        # Apply to all carriers
																											    				        alt_haplotypes <- sapply(carriers, get_alt_haplo)
																											    				        tab = as.data.frame(table(alt_haplotypes))
																																					    tab
																																					    # Find the sequence with the highest freq
																																					    tab$alt_haplotypes = as.character(tab$alt_haplotypes)
																																					    					        maxhap = tab$alt_haplotypes[tab$Freq == max(tab$Freq)][1]
																																					    					        maxhap_idx = which(tab$alt_haplotypes == maxhap)
																																																	    # count mismatches
																																																	    count_mis=c()
																																																	    idx=c()
																																																	    						        for (i in 1:nrow(tab)){
																																																																	      if (tab$alt_haplotypes[i] == maxhap) next
																																																	    						          idx=c(idx, i)
																																																								  							        count_mis = c(count_mis, count_mismatches(maxhap, tab$alt_haplotypes[i]))
																																																								  							      }
																																																	    						        count_mis = cbind(tab[idx,], count_mis)
																																																															    share=max(tab$Freq) + sum(count_mis$Freq[count_mis$count_mis <= n])
																																																															    print(share)
																																																															    							        out = c(out,share)
																																																															    							        # number of mismatches allowed would increase in every loop
																																																															    							        n=n+1
																																																																															  }

																			  out = as.data.frame(cbind(window, out))
																			  names(out) = c("window", "haplo_share")
																			  			    out


																			  			  ##### Looking for the optimal window
																			  			  ## The window has to be at least 250 to make a 500kbp sequence
																			  			  ## The number of people with shared haplotypes should also be over p% of the carriers. Let's choose p as 50
																			  			  ## When looking at the 250 window, if there are at least p of carriers with shared haplotypes, if the next window also have the same freq then keep going to the next windows until the freq falls sharply (delta >= 3)

																			  			    out
																						    			      p=0.5
																						    			      s=NA
																									      			        if (out$haplo_share[out$window == 250] >= p*total_carriers) {
																																			    tmp = out[out$window >= 250 & out$haplo_share >= out$haplo_share[out$window == 250] - 2,]
																									      			          s=max(tmp$window)
																													  				      print(paste0("Shared founder haplotype found with a length of ", s*2))
																													  				      } else {print(paste0("No founder haplotype found for ",   gene, " - ", SNP))}

																									      			        windows = c(windows, s)
																																	  nshares = c(nshares, out$haplo_share[out$window == s][1])

																																	  if (is.na(s)) {nonhapmin = c(nonhapmin, NA)
																																	  				  next}

																																	  				    print(paste0("doing plotting for ", gene) ) 



																																	  				    ## Haplotype sharing vs. window size
																																	  				    png(paste0(gene, ".", SNP, "_Haplotype sharing vs. Window size_allowedmismatches.png"), width=2000, height=1500, res=300)
																																					    				      print(ggplot(out, aes(x=window, y=haplo_share)) + geom_point(col="blue") + geom_line() +
																																																      labs(x="Window size (Mbp)", y="# shared haplotypes",   title=paste0(gene, ": Haplotype sharing vs. Window size")) + geom_hline(yintercept = total_carriers, col='red', linetype='dashed') + geom_vline(xintercept = s, col='blue', linetype='dashed'))
																																					    				      dev.off()

																																									      				        window_size <- s*1000
																																									      				        region_idx <- which(pos >= (mutation_pos - window_size) & pos <= (mutation_pos + window_size))
																																																			  haplo_mat <- extract.gt(vcf[region_idx, ], element = "GT")  # SNPs x samples
																																																			  alt_haplotypes <- sapply(carriers, get_alt_haplo)
																																																			  					    tab = as.data.frame(table(alt_haplotypes))
																																																			  					    tab
																																																								    					      tab$alt_haplotypes = as.character(tab$alt_haplotypes)
																																																								    					      maxhap = tab$alt_haplotypes[tab$Freq == max(tab$Freq)][1]
																																																													      					        maxhap_idx = which(tab$alt_haplotypes == maxhap)
																																																													      					        count_mis=c()
																																																																									  idx=c()
																																																																									  for (i in 1:nrow(tab)){
																																																																										  							      if (tab$alt_haplotypes[i] == maxhap) next
																																																																									  						      idx=c(idx, i)
																																																																															      						          count_mis = c(count_mis, count_mismatches(maxhap, tab$alt_haplotypes[i]))
																																																																															      						        }
																																																																									  						    count_mis = cbind(tab[idx,], count_mis)
																																																																									  						    share=max(tab$Freq) + sum(count_mis$Freq[count_mis$count_mis <= s/50])
																																																																															    						      print(share)

																																																																															    						    ## Prepare data for plotting
																																																																															    						    # write samples that contain the shared extended haploblock
																																																																															    						      exthap = c(maxhap, count_mis$alt_haplotypes[count_mis$count_mis <= s/50])
																																																																																					      						        length(exthap)
																																																																																					      						        exthap_carriers = names(alt_haplotypes)[alt_haplotypes %in% exthap]
																																																																																																			  length(exthap_carriers) == share
																																																																																																			  write.table(exthap_carriers, paste0(gene, ".", SNP, "_extended_haplo_carriers.txt"), row.names = F, quote=F, sep='\t')

																																																																																																			  							  # alt_haplotypes is a named character vector (names = sample IDs, values = 0/1 strings)
																																																																																																			  							  # Convert to a matrix: each row is a sample, each column is a SNP
																																																																																																			  							    haplo_mat2 <- do.call(rbind, strsplit(alt_haplotypes, split = ""))
																																																																																																			  							    rownames(haplo_mat2) <- names(alt_haplotypes)

																																																																																																										    							    # Convert to long format for ggplot
																																																																																																										    							      haplo_df <- as.data.frame(haplo_mat2)
																																																																																																										    							      haplo_df$Sample <- rownames(haplo_df)
																																																																																																																	      							        haplo_long <- haplo_df %>%
																																																																																																																																			    pivot_longer(cols = -Sample, names_to = "SNP", values_to = "Allele") %>%
																																																																																																																																			    									        mutate(SNP = as.integer(gsub("V", "", SNP)),
																																																																																																																																														   										                Allele = as.factor(Allele))


																																																																																																																	      							      # Order samples by haplotypes
																																																																																																																	      							        count_mis = rbind(count_mis, c(unlist(tab[tab$alt_haplotypes == maxhap,]),0))
																																																																																																																																	  count_mis = count_mis[order(as.integer(count_mis$count_mis), decreasing = F),]
																																																																																																																																	  samples_order = c()
																																																																																																																																	  								    for (i in 1:nrow(count_mis)){
																																																																																																																																										    									        samples_order = c(samples_order, names(alt_haplotypes[alt_haplotypes == count_mis$alt_haplotypes[i]]))
																																																																																																																																	  								    }
																																																																																																																																	  								    haplo_long$Sample <- factor(haplo_long$Sample, levels = unique(samples_order))

																																																																																																																																									    								      ## plot
																																																																																																																																									    								      png(paste0(gene, ".",SNP,"_Alt_carrying_Haploblock_allowedmismatches.png"), width=3000, height=2000, res=300)
																																																																																																																																									    								      print(ggplot(haplo_long, aes(x = SNP, y = Sample, color = Allele)) +
																																																																																																																																																												        geom_point(shape = 15, size = 1.8) +
																																																																																																																																																																							    scale_color_manual(values = c("0" = "grey90", "1" = "firebrick")) +
																																																																																																																																																																							    										        #geom_vline(xintercept = mutation_index, color = "blue", linetype = "dashed", linewidth = 0.5) +
																																																																																																																																																																							    										      
																																																																																																																																																																							    										        # Shared haplotype rectangle
																																																																																																																																																																							    										        annotate("rect", 
																																																																																																																																																																																			     												            xmin = -100, xmax = -20,
																																																																																																																																																																																																    													               ymin = 1 - 0.5, ymax = share + 0.5,
																																																																																																																																																																																																    													               fill = "#FEC8C3", alpha = 0.6, color = "#D7263D") +
																																																																																																																																																												        # Other haplotypes rectangle
																																																																																																																																																												        annotate("rect", 
																																																																																																																																																														 											            xmin = -100, xmax = -20,
																																																																																																																																																																										    												               ymin = share + 0.5, ymax = total_carriers + 0.5,
																																																																																																																																																																										    												               fill = "#C7E9F1", alpha = 0.6, color = "#2A9D8F") +
																																																																																																																																																												        coord_cartesian(clip = "off") +
																																																																																																																																																																							    labs(x = "SNP index", y = "Sample", title = "") +
																																																																																																																																																																							    										        theme_minimal() +
																																																																																																																																																																																													    theme(axis.text.y = element_blank(),
																																																																																																																																																																																														  												          axis.ticks = element_blank(),
																																																																																																																																																																																																											  													          panel.grid = element_blank(),
																																																																																																																																																																																																											  													          plot.margin = margin(10, 20, 10, 40),
																																																																																																																																																																																																																									  														          legend.position = "none"))
																																																																																																																																																	      								        dev.off()


																																																																																																																																																	      								      # Check if the identified haplotype is not shared by other samples who are not carrying the 
																																																																																																																																																	      								      # Extract all the haplotypes in the samples that are non-mutation carriers
																																																																																																																																																	      								      non_hap_carriers = colnames(vcf@gt)[-1][!(colnames(vcf@gt)[-1] %in% carriers)]
																																																																																																																																																									      								        
																																																																																																																																																									      								      non_hap_carrier_haplos = c()
																																																																																																																																																									      								      for (i in 1:length(non_hap_carriers)){
																																																																																																																																																																		      									        non_hap_carrier_haplos = c(non_hap_carrier_haplos, get_haplos(non_hap_carriers[i]))
																																																																																																																																																																	      								      }

																																																																																																																																																																	      								      head(non_hap_carrier_haplos)
																																																																																																																																																																	      								      non_hap_carrier_haplos_table = as.data.frame(table(non_hap_carrier_haplos))
																																																																																																																																																																									      								      names(non_hap_carrier_haplos_table )

																																																																																																																																																																									      								      identified_haplotype = maxhap

																																																																																																																																																																																	      								      # Find the distances / mismatches between each of the non-hap_carriers haplotype and the identified haplotype
																																																																																																																																																																																	      								      counts = c()
																																																																																																																																																																																	      								      for (i in 1:nrow(non_hap_carrier_haplos_table)) {
																																																																																																																																																																																										      									        counts = c(counts, count_mismatches(as.character(non_hap_carrier_haplos_table$non_hap_carrier_haplos[i]), identified_haplotype))
																																																																																																																																																																																									      								      }
																																																																																																																																																																																									      								      counts
																																																																																																																																																																																									      								      non_hap_carrier_haplos_table = cbind(non_hap_carrier_haplos_table, counts)
																																																																																																																																																																																																	      								      non_hap_carrier_haplos_table = non_hap_carrier_haplos_table[order(non_hap_carrier_haplos_table$counts),]

																																																																																																																																																																																																	      								      write.table(non_hap_carrier_haplos_table, paste0(gene, ".", SNP, "non_hap_carrier_haplos_table.txt"), row.names=F, quote=F)

																																																																																																																																																																																																									      								      nonhapmin = c(nonhapmin, min(counts))
																																																																																																																																																																																																									      								      min(counts) < s / 50 # does number of mismatches less than what allows for haplotype sharing?

}

res = data.frame(snps, genes, chrss, poss, totalcarriers, totalhomos, windows, nshares, nonhapmin)
res

write.table(res, output_file, quote=F, row.names=F)



