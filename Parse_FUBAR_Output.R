#R script to parse fubar JSON output. Requires information about the kinase domain and where to split it. This is a fairly singular function script 

#Load libraries
library(rjson)
library(dplyr)
library(stringr)
library(rhmmer)
library(phylotools)


selection_dir = "INSERT_DIRECTORY_CONTAINING_FUBAR_JSON_FILES"
hmmer_dir = "INSERT_DIRECTORY_CONTAINING_HMMSEARCH_OUTPUT"
output_dir = "INSERT_OUTPUT_DIRECTORY_NAME"

input_files = list.files(selection_dir, pattern="*.json$", full.names=FALSE)
input_pkinase_files = list.files(hmmer_dir, pattern="*Pkinase.hmm.out.txt$", full.names=FALSE)
input_pkinase_tyr_files = list.files(hmmer_dir, pattern="*Pkinase.Tyr.hmm.out.txt$", full.names=FALSE)


#initialize the output matrices 
output_matrix = matrix(nrow=length(input_files), ncol=10)
colnames(output_matrix) = c("Gene", "Num_Pos_Sites", "Num_Neg_Sites", "Ratio", "Pos_ECD", "Neg_ECD", "ECD_Ratio", "Pos_Kinase", "Neg_Kinase", "Kinase_Ratio")

output_matrix_normalized = output_matrix

#initialize output matrix for substition rate matrices gene wide
substitution_rate_matrix = matrix(nrow=length(input_files), ncol=6)
colnames(substitution_rate_matrix) = c("Syn_rate", "NonSyn_Rate", "Syn_ECD", "NonSyn_ECD", "Syn_Kinase", "NonSyn_Kinase")
rownames(substitution_rate_matrix) = input_files

substitution_rate_matrix_normalized = substitution_rate_matrix

filter_cutoff = 0.90    ####this is a decimal posterior probability. max is 1. 


#loop over all files in the directory. This code would be better suited as a function but...
for (i in 1:length(input_files)){
  
   print(paste("Selection:", input_files[i], " hmmerfile", input_pkinase_files[i], sep=""))
  
  #try to read the kinase table into memory. If it doesn't work, it means most likely the file is empty. Catch this as an error
  hmmer_pkinase = try(read.table(paste(hmmer_dir,input_pkinase_files[i], sep="")),silent=TRUE)
  hmmer_kinase_tyr = try(read.table(paste(hmmer_dir,input_pkinase_tyr_files[i], sep="")), silent=TRUE)
  # 
  # 
  # if (class(hmmer_pkinase) != "try-error" && class(hmmer_kinase_tyr) != "try-error") {
  # colnames(hmmer_pkinase) = c("target_name", "accession_species", "tlen", "domain_name", "accession_hmm", "qlen", "E-value",  "sequence_score",  "sequence_bias", "sequence_num", "seq_of", "c-Evalue",  "i-Evalue",  "domain_score",  "domain_bias",  "hmm_from",    "hmm_to",  "ali_from",  "ali_to",  "env_from", "env_to",  "acc", "description_of_target")
  # 
  #  
  # colnames(hmmer_kinase_tyr) = c("target_name", "accession_species", "tlen", "domain_name", "accession_hmm", "qlen", "E-value",  "sequence_score",  "sequence_bias", "sequence_num", "seq_of", "c-Evalue",  "i-Evalue",  "domain_score",  "domain_bias",  "hmm_from",    "hmm_to",  "ali_from",    "ali_to",  "env_from",     "env_to",  "acc", "description_of_target")
  # 
  # hmmer_pkinase = as.data.frame(hmmer_pkinase)
  # hmmer_kinase_tyr = as.data.frame(hmmer_kinase_tyr)
  # 
  # #hmmer_pkinase = read_domtblout(paste(hmmer_dir,input_pkinase_files[i], sep=""))
  # 
  # #hmmer_kinase_tyr = read_domtblout(paste(hmmer_dir,input_pkinase_tyr_files[i],sep=""))
  # 
  # extra_tyr = hmmer_kinase_tyr %>% filter(!(target_name %in% hmmer_pkinase$target_name))
  # 
  # kinase_positions = rbind(hmmer_pkinase, extra_tyr)
  # kinase_positions$domain_name=str_remove(kinase_positions$domain_name,"\\..*")
  # 
  # hmmer_pkinase = try(read.table(paste(hmmer_dir,input_pkinase_files[i], sep="")),silent=TRUE)
  # hmmer_kinase_tyr = try(read.table(paste(hmmer_dir,input_pkinase_tyr_files[i], sep="")), silent=TRUE)
  
  #if loading the HMMSearch output returns a try-error, then we know it's empty or corrupt. 
  if (class(hmmer_pkinase) == "try-error" && class(hmmer_kinase_tyr) == "try-error") {
    
      warning("Input HMM files are empty")
    
    kinase_positions = "Error"
    
	#If pkinase is empty but tyrosine is not, then get positions from there
  } else if (class(hmmer_pkinase) == "try-error" && class(hmmer_kinase_tyr) != "try-error"){
    
     colnames(hmmer_kinase_tyr) = c("target_name", "accession_species", "tlen", "domain_name", "accession_hmm", "qlen", "E-value",  "sequence_score",  "sequence_bias", "sequence_num", "seq_of", "c-Evalue",  "i-Evalue",  "domain_score",  "domain_bias",  "hmm_from",    "hmm_to",  "ali_from",    "ali_to",  "env_from",     "env_to",  "acc", "description_of_target")
     
    kinase_positions = as.data.frame(hmmer_kinase_tyr)
    kinase_positions$domain_name=str_remove(kinase_positions$domain_name,"\\..*")
    
#	If tyrosine kinase is empty but pkianse is not, get positions only from pkinase
  } else if (class(hmmer_pkinase) != "try-error" && class(hmmer_kinase_tyr) == "try-error") {
    
    colnames(hmmer_pkinase) = c("target_name", "accession_species", "tlen", "domain_name", "accession_hmm", "qlen", "E-value",  "sequence_score",  "sequence_bias", "sequence_num", "seq_of", "c-Evalue",  "i-Evalue",  "domain_score",  "domain_bias",  "hmm_from",    "hmm_to",  "ali_from",  "ali_to",  "env_from", "env_to",  "acc", "description_of_target")
    
    kinase_positions = as.data.frame(hmmer_pkinase)
    kinase_positions$domain_name=str_remove(kinase_positions$domain_name,"\\..*")
    
#	If both are not empty, then find whatever is present in tyrosine kinase but not pkinase and merge them. 
  } else if (class(hmmer_pkinase) != "try-error" && class(hmmer_kinase_tyr) != "try-error") {
    
  colnames(hmmer_pkinase) = c("target_name", "accession_species", "tlen", "domain_name", "accession_hmm", "qlen", "E-value",  "sequence_score",  "sequence_bias", "sequence_num", "seq_of", "c-Evalue",  "i-Evalue",  "domain_score",  "domain_bias",  "hmm_from",    "hmm_to",  "ali_from",  "ali_to",  "env_from", "env_to",  "acc", "description_of_target")
  
   
  colnames(hmmer_kinase_tyr) = c("target_name", "accession_species", "tlen", "domain_name", "accession_hmm", "qlen", "E-value",  "sequence_score",  "sequence_bias", "sequence_num", "seq_of", "c-Evalue",  "i-Evalue",  "domain_score",  "domain_bias",  "hmm_from",    "hmm_to",  "ali_from",    "ali_to",  "env_from",     "env_to",  "acc", "description_of_target")
  
  hmmer_pkinase = as.data.frame(hmmer_pkinase)
  hmmer_kinase_tyr = as.data.frame(hmmer_kinase_tyr)
  
  #hmmer_pkinase = read_domtblout(paste(hmmer_dir,input_pkinase_files[i], sep=""))

  #hmmer_kinase_tyr = read_domtblout(paste(hmmer_dir,input_pkinase_tyr_files[i],sep=""))
  
  extra_tyr = hmmer_kinase_tyr %>% filter(!(target_name %in% hmmer_pkinase$target_name))
  
  kinase_positions = rbind(hmmer_pkinase, extra_tyr)
  kinase_positions$domain_name=str_remove(kinase_positions$domain_name,"\\..*")
  
  }
  

#If kinase positions, which we got from the previous if stack, is not a try error then go through with this 
  if (kinase_positions != "Error"){
  
  #I had a thing to make the names nicer, but eh w/e. It works weird with different inputs. 
  #gene_name = str_replace(basename(input_files[i]), "_.*","")
  gene_name = input_files[i]
  
  #kinase_start_pos = as.vector(kinase_positions[grep(gene_name,kinase_positions$domain_name),'env_from'])
  kinase_start_pos = as.vector(min(kinase_positions$env_from))
  
  #get json of selection 
  selection_json = rjson::fromJSON(file=paste(selection_dir,input_files[i],sep=""))
  
  #get number of sites being tested in each section
  gene_length = selection_json$input$`number of sites`
  ecd_length = kinase_start_pos
  kinase_length = (gene_length - ecd_length) 
  
#Get seelction info from the json list then get all names and stuff setup nicely
  selection_df = as.data.frame(t(as.data.frame(selection_json$MLE$content$`0`)))
  rownames(selection_df)=c(1:nrow(selection_df))
  selection_df = cbind(selection_df, rownames(selection_df)) 
  colnames(selection_df) = c("alpha", "beta", "beta-alpha, mean posterior beta-alpha", "negative", "positive", "bayes_factor", "null1", "null2", "Position")

  #get total number of positive and negative selected sites gene wide. Uses posterior probability whatever it is set to 
    
    pos_sites = selection_df %>% filter(positive > filter_cutoff)
    
    neg_sites = selection_df %>% filter(negative > filter_cutoff)
    
    #ratio
    selection_ratio = nrow(pos_sites)/nrow(neg_sites)
    
    
    ###kinase filtering section#####
    
    #because the TM domain is hard to accurately determine, the "ECD" is actually everything N-terminal of the kinase domain. Includes the transmembrane domain 
    ecd_sites = selection_df[1:kinase_start_pos[[1]],]
    kinase_sites = selection_df[kinase_start_pos[[1]]:nrow(selection_df),]  
    
    ##ecd pos and neg
    pos_ecd = ecd_sites %>% filter(positive > filter_cutoff)
  
    neg_ecd = ecd_sites %>% filter(negative > filter_cutoff) 
    
    #kinase pos and neg
    pos_kinase = kinase_sites %>% filter(positive > filter_cutoff)
   
    neg_kinase = kinase_sites %>% filter(negative > filter_cutoff) 
    
    #get simple ratio of pos to neg selection 
    ecd_pos_neg_ratio = nrow(pos_ecd)/nrow(neg_ecd)
    kinase_pos_neg_ratio = nrow(pos_kinase)/nrow(neg_kinase)
  
  if (selection_json$MLE$content$`0`[1]!="NULL"){
    
    
    
    
    
    
   

    
    
    
    ####output section assuming the input json file is not empty
    
    output_matrix[i,1:10] = c(gene_name, (nrow(pos_sites)), (nrow(neg_sites)), selection_ratio, (nrow(pos_ecd)), (nrow(neg_ecd)), ecd_pos_neg_ratio, (nrow(pos_kinase)), (nrow(neg_kinase)), kinase_pos_neg_ratio)
    
    #add the sites found under selection to the output matrix. Here we divide the value of each by the feature length to get a "rate of positive/negatively selected sites per site" 
    output_matrix_normalized[i,1:10] = c(gene_name, (nrow(pos_sites)/gene_length), (nrow(neg_sites)/gene_length), selection_ratio, (nrow(pos_ecd)/ecd_length), (nrow(neg_ecd)/ecd_length), ecd_pos_neg_ratio, (nrow(pos_kinase)/kinase_length), (nrow(neg_kinase)/kinase_length), kinase_pos_neg_ratio)
    
    
    
    
    
    #if it is null, that means there are no sites under selection. Add an empty line to the growing output matrix. 
  } else {
    
    pos_sites = matrix(nrow=0, ncol=9)
    colnames(pos_sites) = c("alpha", "beta", "beta-alpha, mean posterior beta-alpha", "negative", "positive", "bayes_factor", "null1", "null2", "Position")
    
    neg_sites = matrix(nrow=0, ncol=9)
    colnames(neg_sites) = c("alpha", "beta", "beta-alpha, mean posterior beta-alpha", "negative", "positive", "bayes_factor", "null1", "null2", "Position")
    
    #get the average nonsyn rate and syn rate across the genes being compared
    
    syn_sites = mean(selection_df$alpha)
    
    nonsyn_sites = mean(selection_df$beta)
    
    #split the syn rates up by feature
    #ecd
    syn_ecd = mean(ecd_sites$alpha)
    nonsyn_ecd = mean(ecd_sites$beta)
    
    #kinase
    syn_kinase = mean(kinase_sites$alpha)
    nonsyn_kinase = mean(kinase_sites$beta)
    
    #put gene-wide substitution rate into matrix
    substitution_rate_matrix[i, 1:6] = c(syn_sites, nonsyn_sites, syn_ecd, nonsyn_ecd, syn_kinase, nonsyn_kinase)
    
    #put normalized gene-wide substitution rate into matrix
    substitution_rate_matrix_normalized[i, 1:6] = c(syn_sites/gene_length, nonsyn_sites/gene_length, syn_ecd/ecd_length, nonsyn_ecd/ecd_length, syn_kinase/kinase_length, nonsyn_kinase/kinase_length)
    
    
    
    output_matrix[i,1:10] = c(gene_name, 0, 0, 0, 0, 0, 0, 0, 0, 0)
     
     
    }
  
  }
  
    
  
   
  
  
  #output positively selected sites
    write.csv(pos_sites, paste(output_dir, input_files[i], "_posSelection.csv", sep=""))
  
  #output negatively selected sites
   write.csv(neg_sites, paste(output_dir, input_files[i], "_negSelection.csv", sep=""))
  
}

#output the normalized and not-normalized positive selection and subsitutions rates 
write.csv(output_matrix, paste(output_dir, "RLKs_ALL.csv", sep=""), row.names=FALSE)
write.csv(output_matrix_normalized, paste(output_dir, "RLKs_ALL_Normalized.csv", sep=""), row.names=FALSE)

#output the final substitution rate matrix 
write.csv(substitution_rate_matrix, paste(output_dir, "RLKs_ALL_Substitution_Rate.csv", sep=""), row.names=TRUE)
write.csv(substitution_rate_matrix_normalized, paste(output_dir, "RLKs_ALL_Substitution_Rate_Normalized.csv", sep=""), row.names=TRUE)

