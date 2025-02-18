#library(tidyverse)
library(dplyr)

#get input dataset from orthofinder. This will be the actual values of counts per OG split into hierarchical OGs
#NOTE: My formatting was weird for my output. First column = HOG name, Second column = Original OG name, Third column = Gene tree OG was derived from, Fourth to end = the species/ecotype data
OG_Data= as.data.frame(read.csv("YOUR_CSV_FILE_WITH_THE_OG_DATA"))          
colnames(OG_Data)[1] = "HOG"

families = unique(OG_Data$Family)
#OG_Data = OG_Data %>% filter(Family == "MLRR")

#prepare new potential OGs

df_to_return = data.frame()
removal_OGs = vector()

#Loop over every row in the OG Gene count data
for (l in 1:nrow(OG_Data)){
  
  #First, we have to split the OGs to ensure each OG ONLY has 1 gene from each ecotype.
  if (max(OG_Data[l,4:ncol(OG_Data)]) > 1){
    
    removal_OGs = c(removal_OGs, OG_Data[l,1])
  
  separation_num= max(OG_Data[l,4:ncol(OG_Data)]) - 1
  
 separation_df = rbind(OG_Data[l,], OG_Data[rep(l,separation_num),])
  
  #separation_matrix = as.data.frame(rep(OG_Data[l,], times=separation_num))
  #separation_matrix = matrix( nrow = separation_num, ncol=ncol(OG_Data))
  #colnames(separation_matrix) = colnames(OG_Data)
  
  
  
  
  
  #Remove duplicate names 
  for (k in 1:nrow(separation_df)){
    
    separation_df[k,4:ncol(separation_df)] = separation_df[k,4:ncol(separation_df)] - (k-1)
    
    separation_df[k,1] = paste(separation_df[k,1], "_", k, sep="")
    
  }
 
  

  # if the value is above 0, change it to 1. If it's below 0, change it to 0
  separation_df[,4:ncol(separation_df)] = ifelse(separation_df[,4:ncol(separation_df)] >0, 1, 0)
  
  df_to_return = rbind(df_to_return, separation_df[1:nrow(separation_df),])
  #df_to_return = filter(df_to_return, )
  
  }
  
}

#Update the OG data to be in the separated format. Now we can start counting 
OG_Data = OG_Data %>% filter(!(HOG %in% removal_OGs))
OG_Data = rbind(OG_Data, df_to_return)


#OG_Data = OG_Data[sum(OG_Data[,4:ncol(OG_Data)]) >=1 ,]
#OG_Data = OG_Data[sum(sapply(OG_Data[,4:ncol(OG_Data)], as.numeric)) >=2 ,]
#OG_Data = OG_Data[rowSums(OG_Data[,4:ncol(OG_Data)]) >=2,]
OG_Data = OG_Data[rowSums(OG_Data[,4:ncol(OG_Data)]) >=1,]    #get total number of genes in each OG 

#set the number of input species/ecotypes and the number of iterations
species_number = 146
num_iterations = 1000

#holding_list = vector("list", length=num_iterations)

#initialize the output matrices
core_genome_matrix = matrix(nrow=num_iterations,ncol=species_number)
rlkome_genome_matrix = matrix(nrow=num_iterations, ncol=species_number)
#rownames(holding_matrix) = c("RLKome_size", "Core_genome_size")

#this for loop is for the number of iterations wanted to the rarefaction. The more numbers, the smoother the average curve. 1000 is default bceause it gives a nice, smooth curve. More the better obviously. 
for (j in 1:num_iterations){
  
print(c("Iteration", j), quote=FALSE)
  
#get a random subsamble starting at the fourth column (first data column) until the last and do it for 146 samples (ecotype number)
OGs_randomized = cbind(OG_Data[,1:3],sample(OG_Data[,4:ncol(OG_Data)],species_number, replace=FALSE))


#set holding variable to hold new OGs
OGs_collected = vector()

#start with vector of all HOGs. We will use this and intersect to get core as we go along
core_ogs = OG_Data$HOG

#loop over columns (num species/ecotypes/observations)
for (i in 1:(ncol(OGs_randomized)-3)){
  
  #filter by those with presence of a gene in the OG
  filtered_OGs = OGs_randomized %>% filter (OGs_randomized[,i+3] == 1)
  
  #find intersection of the HOGs we have and the new ones added. If they intersect, that's the 'core genome' 
  core_ogs = intersect(core_ogs, filtered_OGs$HOG)
  
  
  #OG_list = filtered_OGs$HOG
  # find difference between the new OGs and the one's we've already found to see if there are new ones 
  new_ogs = setdiff(filtered_OGs$HOG,OGs_collected)
  
  
  
  
  #add new ogs to growing collection
  OGs_collected = c(OGs_collected,new_ogs)
  
  #add the NOG count to matrix with rows = iterations, i = genome
  rlkome_genome_matrix[j, i] = length(OGs_collected)
  core_genome_matrix[j,i] = length(core_ogs)
  
  
}


  
}

#Now figure out the number of new genes found per iteration and add it to the output matrix 
new_genes_per_iteration = matrix(nrow=nrow(rlkome_genome_matrix), ncol=ncol(rlkome_genome_matrix))
new_genes_per_iteration[,1] = rlkome_genome_matrix[,1]

for (newgene in 2:ncol(new_genes_per_iteration)){
  
  new_genes_per_iteration[,newgene] = rlkome_genome_matrix[,newgene] - rlkome_genome_matrix[,(newgene-1)]
  
}