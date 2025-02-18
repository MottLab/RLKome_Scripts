###chunk takes all by all blast and determines primary OG groupings starting with single copy genes as base


###WARNING: this was made to try and catch all cases and combinations of single copy and multicopy genes. There may be instances I missed, so be sure to check output. 

#Split OGs function. Takes as input a blast all by all for all the genes/proteins in the OG we're trying to split, and the output folder. 
SplitOGs = function(blast_filein, OG_name, outfolder){
  
  ####for troubleshooting
  #blast_filein = paste(blastfolder, bfile, sep="")
  #OG_name = "N0.HOG0000000"
 #outfolder = "D:/Analysis_Output/Orthogroup_Assignment/Sep2024/Kinases/Orthofinder_out/Split_HOGs/To_Split/Blast_Refinement/N0.HOG0001157_allbyall"
 #
  #OG_name = gsub("\\.blast","",bfile)
  #print(OGname)
  #outfolder = paste("D:/Analysis_Output/Orthogroup_Assignment/Sep2024/Kinases/Orthofinder_out/Results_Dec9/Refinement/Split/", OGname, sep="")

  require(tidyverse)
  
  #If someone doesn't type in the correct path aka misses a /, then just fix path
  if (str_sub(outfolder, start=-1) != "/"){
    outfolder = paste(outfolder, "/", sep="")
  }
  
#create output directory
dir.create(outfolder, showWarnings = FALSE)

#get our blast file in
blast_file_in = read.table(blast_filein)
#rename columns with balst outfmt 6 default names
colnames(blast_file_in) = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')


#get the total non-redundant list of genes
OG_genes = unique(blast_file_in$qseqid)

counter=1   #set counter for OG numbers
flag = 0    #set our while loop flag

#while there are still multiple groups, split the OGs
while (flag == 0 ){

  #first thing we'll get a list of ecotypes assuming ecotype/species name is separated from gene name by a |
ecotypes = gsub("\\|.*","",OG_genes, perl=TRUE)
#get single copy genes - i.e. ones that are not duplicated in this round
single_copy = ecotypes[!(duplicated(ecotypes) | duplicated(ecotypes, fromLast=TRUE))]
#get ecotypes that have more than one gene in this cluster
multi_copy = unique(ecotypes[duplicated(ecotypes)])




#if we've reached a point where we're out multicopy genes, don't continue. Otherwise, continue 
if (length(multi_copy) < 1 | length(single_copy) == 0) {
  
  flag = 1
  
  #nested insdie here also catch cases to catch the remainders and output them. This happens when single copy 
  if (length(single_copy) >=2){
  singlecopy_blast = blast_file_in %>% filter(grepl(paste("^",single_copy, collapse="|", sep=""), qseqid, perl=TRUE))
  
  write.table(unique(singlecopy_blast$qseqid), paste(outfolder, OG_name, "_", counter, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")
  
  
  #If there are still multicopy remaining but no single copy, then output them. 
  } else if (length(single_copy) == 0 & length(multi_copy) > 0 ) {
      
      multicopy_blast = blast_file_in %>% filter(grepl(paste("^", multi_copy, collapse="|", sep=""), qseqid, perl=TRUE) & (sseqid %in% singlecopy_blast$qseqid))
      
    write.table(unique(multicopy_blast$qseqid), paste(outfolder, OG_name, "_", counter, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")
      
  
  }
  
} else {

  #now so the actual stuff
  
  # get the blast all by all info for the single copy genes in thsi cluseter
singlecopy_blast = blast_file_in %>% filter(grepl(paste("^",single_copy, collapse="|", sep=""), qseqid, perl=TRUE))

#get the duplicated genes from this cluster all by all blast info
multicopy_blast = blast_file_in %>% filter(grepl(paste("^", multi_copy, collapse="|", sep=""), qseqid, perl=TRUE) & (sseqid %in% singlecopy_blast$qseqid))

#get the duplicated genes
duplicated_genes = unique(multicopy_blast$qseqid)

#we're gonna try two ways to do this. First is to find which of the pairs is, on average, more similar to the single copies. 

#the second is we're going to create x groups and figure out which one, on average, has the highest similarity across them. 

genes_to_add = vector(mode = "character")    #set holding vector

#for each ecotype in the multicopy, figure out which one has the highest average similarity to the single copy genes
for (i in 1:length(multi_copy)){
  
    #print("line 187") 
  #find which duplicates to compare
  duplicates_to_compare = duplicated_genes[grepl(multi_copy[i],duplicated_genes)]
  
  #set the average score holding matrix and fill it with the duplicates to compare number
  average_score = matrix(nrow=length(duplicates_to_compare), ncol=2)
  average_score[,1] = duplicates_to_compare
  
  #for each row, calculate the average percent identify and put that in the average score
  for (j in 1:nrow(average_score)){
    
    single_gene_temp_blast = multicopy_blast %>% filter(qseqid == duplicates_to_compare[j])
    
    average_score[j,2] = mean(single_gene_temp_blast$pident)
    
    
    
  }
  #add the highest percent identity gene to a list of genes to add to our cluster for output
  genes_to_add = c(genes_to_add, average_score[order(as.numeric(average_score[,2]), decreasing=TRUE),1][1])
  
}
  #add genes to add to our single copy genes and then sort it A-Z
  holding_data_frame = data.frame(c(unique(singlecopy_blast$qseqid), genes_to_add))
  holding_data_frame[,1] = holding_data_frame[order(holding_data_frame[,1]),]
  
  #write output
  write.table(unique(holding_data_frame), paste(outfolder, OG_name, "_", counter, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")
  
  #find difference between the cluster we just output and the remaining duplicates
  OG_genes = setdiff(duplicated_genes, holding_data_frame[,1])
  
  blast_file_in = blast_file_in %>% filter(qseqid %in% OG_genes)
  
  if (length(OG_genes) <=1){
    flag = 1
    
  
  }
  
  counter = counter + 1

  }
}


} #function end