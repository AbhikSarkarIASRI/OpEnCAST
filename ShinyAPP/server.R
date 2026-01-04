##Supress Warnings
options(warn=-1)
options(timeout = 600)



library(seqinr)
library(ranger)
library(xgboost)
library(gbm)
library(e1071)
library(ftrCOOL)



##Define Function
my_packages <- c("seqinr","e1071", "ranger", "gbm", "xgboost", "ftrCOOL")

lapply(my_packages, require, character.only = TRUE) 

#rtree<-readRDS("models/Rtree.rds")
#gb<- readRDS("models/GB.rds")
#ann<-readRDS("models/ANN.rds")
Monocot<-function(fasta_file_path){
set.seed(123)
libraries <- c("seqinr",  "e1071",  "ranger", "gbm",  "ftrCOOL")
sapply(libraries, library, character.only = TRUE)
fasta_file <- fasta_file_path
################################Training##############################



######################### Tabular format to Fasta format###############################

#this is a function to convert tabular fasta into plain fasta file
#first column should be squence names
#second column should be sequence

#######################Fasta to Tabular format##################################

FastaToTabular <- function (filename){
  
  #read fasta file
  
  file1 <- readLines(filename)
  
  #find the genename location by grepping >
  
  location <- which((substr(file1,1,1))==">")
  
  #start an empty vector to collect name and sequence
  
  name=c()
  sequence =c()
  
  
  
  #number of genes= number of loops
  #extract name first
  for ( i in 1:length(location)){
    name_line = location[i]
    name1 = file1[name_line]
    name=c(name,name1)
    #extract sequence between the names
    #the last sequence will be missed using this strategy
    #so, we are using if condition to extract last sequence
    start= location[i]+1
    end = location[i+1]-1
    if ( i < length (location)){
      
      end=end
      
    } else {
      
      end=length(file1)
    }
    
    lines = start:end
    sequence1= as.character(paste(file1[lines],collapse = ""))
    sequence =c(sequence,sequence1)
  }
  
  #now create table using name and sequence vector
  
  data <- data.frame(name = name, sequence = sequence, stringsAsFactors = FALSE)
  
  
  
  
  #finally export the file
  #before that remove preexisting file
  unlink(c("dna_table.csv"),force=TRUE)
  as.matrix(data,"dna_table.csv")
  
  #function ends
}
#########################alphabetcheck###########################
alphabetCheck<-function (sequences, alphabet = "aa", label = c())
{
  if (length(sequences) == 0) {
    stop("ERROR: sequence parameter is empty")
  }
  if (length(label) != 0 && length(label) != length(sequences)) {
    stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
  }
  if (alphabet == "rna") {
    alphabet <- c("A", "C", "G", "U")
  }
  else if (alphabet == "dna") {
    alphabet <- c("A", "C", "G", "T")
  }
  else if (alphabet == "aa") {
    alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I",
                  "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                  "W", "Y")
  }
  else {
    stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
  }
  alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                             split = "")[[1]] %in% alphabet))
  flag = 0
  if (length(label) == length(sequences)) {
    flag = 1
    label = label[alphabetCheck]
  }
  else if (length(label) > 0 && length(label) != length(sequences)) {
    stop("ERROR: The number of labels is not equal to the number of sequences!")
  }
  if (is.null(names(sequences))) {
    names(sequences) <- as.character(1:length(sequences))
  }
  nonstanSeq <- names(sequences)[!alphabetCheck]
  if (length(nonstanSeq) != 0) {
    nonstanSeq <- toString(nonstanSeq)
    warMessage <- paste("The sequences (", nonstanSeq, ") were deleted. They contained non-standard alphabets")
    message(warMessage)
  }
  sequences = sequences[alphabetCheck]
  if (length(sequences) == 0) {
    stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
  }
  if (flag == 1) {
    names(label) = names(sequences)
  }
  seq_lab <- list(sequences = sequences, Lab = label)
  return(seq_lab)
}
#################################NCP_DNA############################

ncp_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                   label = c())
{
  if (length(seqs) == 1 && file.exists(seqs)) {
    seqs <- fa.read(seqs, alphabet = "dna")
    seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
    seqs <- seqs_Lab[[1]]
    label <- seqs_Lab[[2]]
  }
  else if (is.vector(seqs)) {
    seqs <- sapply(seqs, toupper)
    seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
    seqs <- seqs_Lab[[1]]
    label <- seqs_Lab[[2]]
  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }
  lenSeqs <- sapply(seqs, nchar)
  nucs <- list(A = c(1, 1, 1), C = c(0, 0, 1), G = c(1, 0,
                                                     0), T = c(0, 1, 0), U = c(0, 1, 0))
  numSeqs <- length(seqs)
  if (outFormat == "mat") {
    if (length(unique(lenSeqs)) > 1) {
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }
    if (binaryType == "strBin") {
      nucs <- c(A = "111", C = "001", G = "100", T = "010",
                U = "010")
      featureMatrix <- sapply(seqs, function(x) {
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        return(cods)
      })
      featureMatrix <- t(featureMatrix)
      colnames(featureMatrix) <- paste("ncp_pos", 1:lenSeqs[1],
                                       sep = "")
      row.names(featureMatrix) <- names(seqs)
    }
    else if (binaryType == "logicBin") {
      nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                  TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                 FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
      featureMatrix <- sapply(seqs, function(x) {
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        cods <- unlist(cods)
        return(cods)
      })
      featureMatrix <- t(featureMatrix)
      temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
      temp2 <- rep(1:lenSeqs[1], each = 3)
      colnames(featureMatrix) <- paste("ncp_pos", temp2, "-",
                                       temp1, sep = "")
      row.names(featureMatrix) <- names(seqs)
    }
    else if (binaryType == "numBin") {
      featureMatrix <- sapply(seqs, function(x) {
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        cods <- unlist(cods)
        return(cods)
      })
      featureMatrix <- t(featureMatrix)
      temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
      temp2 <- rep(1:lenSeqs[1], each = 3)
      colnames(featureMatrix) <- paste("ncp_pos", temp2, "-",
                                       temp1, sep = "")
      row.names(featureMatrix) <- names(seqs)
    }
    else {
      stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
    }
    return(featureMatrix)
  }
  else if (outFormat == "txt") {
    nucs <- c(A = "111", C = "001", G = "100", T = "010",
              U = "010")
    counter <- 0
    namesSeqs <- names(seqs)
    codes <- lapply(seqs, function(x) {
      counter <- counter + 1
      charList <- unlist(strsplit(x, split = ""))
      cods <- nucs[charList]
      namecods <- namesSeqs[counter]
      cods <- unlist(cods)
      cods <- c(namecods, cods)
      temp <- paste(cods, collapse = "\t")
      write(temp, outputFileDist, append = TRUE)
    })
  }
  else {
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }
}



###########################GC_Content##########################
GC.content <- function(fasta_file){
  x <- read.fasta(file=fasta_file)
  tt<-function(x){
    res<-GC(x)
    val=round(res,4)
    return(val)
  }
  
  f_res<-lapply(x,tt)
  s=data.frame(f_res)
  
  rownames(s) <- c("GC-content")
  
  w=t(s)
  return(w)
}
################################ONF#################################
oligo.freq <- function(fasta_file, f){
  
  x <- seqinr::read.fasta(fasta_file, as.string = TRUE, seqtype = "DNA")
  x <- sapply(x, toupper)
  
  bases <- c("A","C","G","T")
  kmers <- expand.grid(rep(list(bases), f))
  kmers <- apply(kmers, 1, paste0, collapse = "")
  
  calc <- function(s){
    tab <- table(factor(sapply(1:(nchar(s)-f+1),
                               function(i) substr(s, i, i+f-1)), levels = kmers))
    as.numeric(tab)
  }
  
  y <- t(sapply(x, calc))
  z <- as.data.frame(y)
  colnames(z) <- kmers
  rownames(z) <- names(x)
  
  return(z)
}

#########################################AMIP################################


#######################mononucleotide_binary_encoding##################################

FastaToTabular <- function (filename){
  
  #read fasta file
  
  file1 <- readLines(filename)
  
  #find the genename location by grepping >
  
  location <- which((substr(file1,1,1))==">")
  
  #start an empty vector to collect name and sequence
  
  name=c()
  sequence =c()
  
  
  
  #number of genes= number of loops
  #extract name first
  for ( i in 1:length(location)){
    name_line = location[i]
    name1 = file1[name_line]
    name=c(name,name1)
    #extract sequence between the names
    #the last sequence will be missed using this strategy
    #so, we are using if condition to extract last sequence
    start= location[i]+1
    end = location[i+1]-1
    if ( i < length (location)){
      
      end=end
      
    } else {
      
      end=length(file1)
    }
    
    lines = start:end
    sequence1= as.character(paste(file1[lines],collapse = ""))
    sequence =c(sequence,sequence1)
  }
  
  #now create table using name and sequence vector
  
  data <- data.frame(name = name, sequence = sequence, stringsAsFactors = FALSE)
  
  
  
  #finally export the file
  #before that remove preexisting file
  unlink(c("dna_table.csv"),force=TRUE)
  as.matrix(data,"dna_table.csv")
  
  #function ends
}
#########################alphabetcheck###########################
alphabetCheck<-function (sequences, alphabet = "aa", label = c())
{
  if (length(sequences) == 0) {
    stop("ERROR: sequence parameter is empty")
  }
  if (length(label) != 0 && length(label) != length(sequences)) {
    stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
  }
  if (alphabet == "rna") {
    alphabet <- c("A", "C", "G", "U")
  }
  else if (alphabet == "dna") {
    alphabet <- c("A", "C", "G", "T")
  }
  else if (alphabet == "aa") {
    alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I",
                  "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                  "W", "Y")
  }
  else {
    stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
  }
  alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                             split = "")[[1]] %in% alphabet))
  flag = 0
  if (length(label) == length(sequences)) {
    flag = 1
    label = label[alphabetCheck]
  }
  else if (length(label) > 0 && length(label) != length(sequences)) {
    stop("ERROR: The number of labels is not equal to the number of sequences!")
  }
  if (is.null(names(sequences))) {
    names(sequences) <- as.character(1:length(sequences))
  }
  nonstanSeq <- names(sequences)[!alphabetCheck]
  if (length(nonstanSeq) != 0) {
    nonstanSeq <- toString(nonstanSeq)
    warMessage <- paste("The sequences (", nonstanSeq, ") were deleted. They contained non-standard alphabets")
    message(warMessage)
  }
  sequences = sequences[alphabetCheck]
  if (length(sequences) == 0) {
    stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
  }
  if (flag == 1) {
    names(label) = names(sequences)
  }
  seq_lab <- list(sequences = sequences, Lab = label)
  return(seq_lab)
}
#################################MBE_DNA############################

mbe_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                   label = c())
{
  if (length(seqs) == 1 && file.exists(seqs)) {
    seqs <- fa.read(seqs, alphabet = "dna")
    seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
    seqs <- seqs_Lab[[1]]
    label <- seqs_Lab[[2]]
  }
  else if (is.vector(seqs)) {
    seqs <- sapply(seqs, toupper)
    seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
    seqs <- seqs_Lab[[1]]
    label <- seqs_Lab[[2]]
  }
  else {
    stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
  }
  lenSeqs <- sapply(seqs, nchar)
  nucs <- list(A = c(1, 0, 0, 0), C = c(0, 1, 0, 0), G = c(0, 0, 1, 0), T = c(0, 0, 0, 1), U = c(0, 0, 0, 1))
  numSeqs <- length(seqs)
  if (outFormat == "mat") {
    if (length(unique(lenSeqs)) > 1) {
      stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
    }
    if (binaryType == "strBin") {
      nucs <- c(A = "1000", C = "0100", G = "0010", T = "0001",
                U = "0001")
      featureMatrix <- sapply(seqs, function(x) {
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        return(cods)
      })
      featureMatrix <- t(featureMatrix)
      colnames(featureMatrix) <- paste("pos_mbe", 1:lenSeqs[1],
                                       sep = "")
      row.names(featureMatrix) <- names(seqs)
    }
    else if (binaryType == "logicBin") {
      nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                  TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                 FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
      featureMatrix <- sapply(seqs, function(x) {
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        cods <- unlist(cods)
        return(cods)
      })
      featureMatrix <- t(featureMatrix)
      temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
      temp2 <- rep(1:lenSeqs[1], each = 3)
      colnames(featureMatrix) <- paste("pos_mbe", temp2, "-",
                                       temp1, sep = "")
      row.names(featureMatrix) <- names(seqs)
    }
    else if (binaryType == "numBin") {
      featureMatrix <- sapply(seqs, function(x) {
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        cods <- unlist(cods)
        return(cods)
      })
      featureMatrix <- t(featureMatrix)
      temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
      temp2 <- rep(1:lenSeqs[1], each = 3)
      colnames(featureMatrix) <- paste("pos_mbe", temp2, "-",
                                       temp1, sep = "")
      row.names(featureMatrix) <- names(seqs)
    }
    else {
      stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
    }
    return(featureMatrix)
  }
  else if (outFormat == "txt") {
    nucs <- c(A = "1000", C = "0100", G = "0010", T = "0001",
              U = "0001")
    counter <- 0
    namesSeqs <- names(seqs)
    codes <- lapply(seqs, function(x) {
      counter <- counter + 1
      charList <- unlist(strsplit(x, split = ""))
      cods <- nucs[charList]
      namecods <- namesSeqs[counter]
      cods <- unlist(cods)
      cods <- c(namecods, cods)
      temp <- paste(cods, collapse = "\t")
      write(temp, outputFileDist, append = TRUE)
    })
  }
  else {
    stop("ERROR: outFormat should be 'mat' or 'txt' ")
  }
}
#########EIIP############


EIIP_final<-EIIP(fasta_file)



# Define a substring or keyword to search for
# find patterns in names
posSeqs <- ftrCOOL::fa.read(url("https://zenodo.org/records/17934726/files/Rice_Pos.fasta?download=1"))
negSeqs <- ftrCOOL::fa.read(url("https://zenodo.org/records/17934726/files/Rice_Neg.fasta?download=1"))
seqs <- fa.read(file = fasta_file, alphabet = "dna")
PSTNPss_Final <- PSTNPss_DNA(seqs = seqs, pos = posSeqs, neg = negSeqs)

rf_model<-readRDS(url("https://zenodo.org/records/17934726/files/rf_rice.rds?download=1"))
gb_model<-readRDS(url("https://zenodo.org/records/17934726/files/gb_rice.rds?download=1"))
svm_model<-readRDS(url("https://zenodo.org/records/17934726/files/svm_rice.rds?download=1"))



#################################data_preparation######################
res<-FastaToTabular(fasta_file)
data<-as.vector(res[,2])
mat<-as.matrix(ncp_dna(seqs = data,binaryType="strBin",outFormat="mat"))
sequence<-rownames(mat)
seq_id<-res[,1]
ncp<-cbind(seq_id,sequence,mat)
rownames(ncp)<-seq_id
ncp_temp<-data.frame(ncp[,-1], stringsAsFactors = FALSE)
ncp_final<-as.data.frame(apply(ncp_temp[,-1], 2, as.numeric))
log_gc_temp<-log((GC.content(fasta_file))*100, base = exp(1))
log_gc<-as.data.frame(as.numeric((ifelse(log_gc_temp>0,log_gc_temp,'0'))))
onf<-oligo.freq(fasta_file, 2)
res_temp_mbe<-FastaToTabular(fasta_file)
data_mbe<-as.vector(res_temp_mbe[,2])
res_mbe<-as.matrix(mbe_dna(seqs = data_mbe,binaryType="strBin",outFormat="mat"))
mbe_temp<-data.frame(res_mbe[,-1], stringsAsFactors = FALSE)
mbe_final<-as.data.frame(apply(mbe_temp[,-1], 2, as.numeric))
# AMIP_final<- NULL
# for (i in 1:6) {
#   for(j in 2:7){
#     if(i<j){
#       AMIP1<-AMIP(fasta_file,i,j)
# 
#       colnames(AMIP1) <- paste('Mean',i,j, sep="_")
#       AMIP_final<-cbind(AMIP_final, AMIP1)
#     }
#   }
# }

temp1<- cbind(onf, gcc =log_gc[,1], ncp_final, mbe_final,  EIIP_final)
my_data_temp<- temp1
inputData <-as.data.frame(my_data_temp)

selected_columns <- c("pos21",	"ncp_pos21",	"pos_mbe21",	"ncp_pos1",	"pos1",
                      "ncp_pos2",	"pos2",	"ncp_pos27",	"pos_mbe27",	"pos27",
                      "ncp_pos25",	"pos_mbe25",	"pos25",	"ncp_pos24",	"pos_mbe24",
                      "pos24",	"AA",	"ncp_pos22",	"pos_mbe22", "pos_mbe28",	"TG",	"pos23",	"ncp_pos23")

data_temp <- inputData[,selected_columns]
test_data<-data_temp
test_data<-cbind(Sequence=ncp_temp$sequence, data_temp)




############################Prediction########################
predicted_prob_svm <- predict(svm_model, newdata = test_data, type = "prob")
#predicted_prob_svm <- attr(predicted_prob_svm_temp, "probabilities")
predicted_value_svm <- predict(svm_model, newdata = test_data)
predicted_prob_rf <- predict(rf_model, newdata = test_data, type = "prob") 
predicted_value_rf <- predict(rf_model, newdata = test_data)
predicted_prob_gb <- predict(gb_model, newdata = test_data, type = "prob")  
predicted_value_gb <- predict(gb_model, newdata = test_data)
test_res_en_prob <-  cbind(SVM = predicted_prob_svm, RF = predicted_prob_rf, GB = predicted_prob_gb)


##################Ensemble###################
# Define weights
weights <- c(SVM = 0.2644348, RF = 0.2532280, GB = 0.4823372)

# Apply weights
weighted_data <- as.data.frame(test_res_en_prob)
weighted_data[1:3] <- test_res_en_prob[1:3] * weights["SVM"]
weighted_data[4:6] <- test_res_en_prob[4:6] * weights["RF"]
weighted_data[7:9] <- test_res_en_prob[7:9] * weights["GB"]

weighted_data$Class_1 <- rowSums(weighted_data[c("SVM.X1", "RF.X1", "GB.X1")])
weighted_data$Class_2 <- rowSums(weighted_data[c("SVM.X2", "RF.X2", "GB.X2")])
weighted_data$Class_3 <- rowSums(weighted_data[c("SVM.X3", "RF.X3", "GB.X3")])



weighted_data$Predicted_Class <- apply(weighted_data[, c("Class_1", "Class_2", "Class_3")], 1, which.max)


final_results <- weighted_data[, c("Class_1", "Class_2", "Class_3", "Predicted_Class")]

final_results$Predicted_Class_Label <- ifelse(final_results$Predicted_Class == 2, "DNA 6mA Methylation",
                                              ifelse(final_results$Predicted_Class == 3, "DNA 5mC Methylation",
                                                     ifelse(final_results$Predicted_Class == 1, "No Methylation", NA)))
Ids <- rownames(ncp_temp)
Ids <- sub("^>", "", Ids)

final_pred <- data.frame(Ids = Ids,
                         Sequence = ncp_temp[,1],
                         Modification = final_results$Predicted_Class_Label,
                         Probability = round(apply(final_results[, c("Class_1", "Class_2", "Class_3")], 1, max), 2))
rownames(final_pred) <- NULL
return(final_pred)
}
Dicot<-function(fasta_file_path){
  set.seed(123)
  libraries <- c("seqinr", "e1071",  "ranger", "xgboost", "ftrCOOL")
  sapply(libraries, library, character.only = TRUE)
  fasta_file <- fasta_file_path
  ################################Training##############################
  
  
  ######################### Tabular format to Fasta format###############################
  
  #this is a function to convert tabular fasta into plain fasta file
  #first column should be squence names
  #second column should be sequence
  
  #######################Fasta to Tabular format##################################
  
  FastaToTabular <- function (filename){
    
    #read fasta file
    
    file1 <- readLines(filename)
    
    #find the genename location by grepping >
    
    location <- which((substr(file1,1,1))==">")
    
    #start an empty vector to collect name and sequence
    
    name=c()
    sequence =c()
    
    
    
    #number of genes= number of loops
    #extract name first
    for ( i in 1:length(location)){
      name_line = location[i]
      name1 = file1[name_line]
      name=c(name,name1)
      #extract sequence between the names
      #the last sequence will be missed using this strategy
      #so, we are using if condition to extract last sequence
      start= location[i]+1
      end = location[i+1]-1
      if ( i < length (location)){
        
        end=end
        
      } else {
        
        end=length(file1)
      }
      
      lines = start:end
      sequence1= as.character(paste(file1[lines],collapse = ""))
      sequence =c(sequence,sequence1)
    }
    
    #now create table using name and sequence vector
    
    data <- data.frame(name = name, sequence = sequence, stringsAsFactors = FALSE)
    
    
    
    
    #finally export the file
    #before that remove preexisting file
    unlink(c("dna_table.csv"),force=TRUE)
    as.matrix(data,"dna_table.csv")
    
    #function ends
  }
  #########################alphabetcheck###########################
  alphabetCheck<-function (sequences, alphabet = "aa", label = c())
  {
    if (length(sequences) == 0) {
      stop("ERROR: sequence parameter is empty")
    }
    if (length(label) != 0 && length(label) != length(sequences)) {
      stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
    }
    if (alphabet == "rna") {
      alphabet <- c("A", "C", "G", "U")
    }
    else if (alphabet == "dna") {
      alphabet <- c("A", "C", "G", "T")
    }
    else if (alphabet == "aa") {
      alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I",
                    "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                    "W", "Y")
    }
    else {
      stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
    }
    alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                               split = "")[[1]] %in% alphabet))
    flag = 0
    if (length(label) == length(sequences)) {
      flag = 1
      label = label[alphabetCheck]
    }
    else if (length(label) > 0 && length(label) != length(sequences)) {
      stop("ERROR: The number of labels is not equal to the number of sequences!")
    }
    if (is.null(names(sequences))) {
      names(sequences) <- as.character(1:length(sequences))
    }
    nonstanSeq <- names(sequences)[!alphabetCheck]
    if (length(nonstanSeq) != 0) {
      nonstanSeq <- toString(nonstanSeq)
      warMessage <- paste("The sequences (", nonstanSeq, ") were deleted. They contained non-standard alphabets")
      message(warMessage)
    }
    sequences = sequences[alphabetCheck]
    if (length(sequences) == 0) {
      stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
    }
    if (flag == 1) {
      names(label) = names(sequences)
    }
    seq_lab <- list(sequences = sequences, Lab = label)
    return(seq_lab)
  }
  #################################NCP_DNA############################
  
  ncp_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                     label = c())
  {
    if (length(seqs) == 1 && file.exists(seqs)) {
      seqs <- fa.read(seqs, alphabet = "dna")
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else if (is.vector(seqs)) {
      seqs <- sapply(seqs, toupper)
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else {
      stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
    }
    lenSeqs <- sapply(seqs, nchar)
    nucs <- list(A = c(1, 1, 1), C = c(0, 0, 1), G = c(1, 0,
                                                       0), T = c(0, 1, 0), U = c(0, 1, 0))
    numSeqs <- length(seqs)
    if (outFormat == "mat") {
      if (length(unique(lenSeqs)) > 1) {
        stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
      }
      if (binaryType == "strBin") {
        nucs <- c(A = "111", C = "001", G = "100", T = "010",
                  U = "010")
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        colnames(featureMatrix) <- paste("ncp_pos", 1:lenSeqs[1],
                                         sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "logicBin") {
        nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                    TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                   FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("ncp_pos", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "numBin") {
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("ncp_pos", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else {
        stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
      }
      return(featureMatrix)
    }
    else if (outFormat == "txt") {
      nucs <- c(A = "111", C = "001", G = "100", T = "010",
                U = "010")
      counter <- 0
      namesSeqs <- names(seqs)
      codes <- lapply(seqs, function(x) {
        counter <- counter + 1
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        namecods <- namesSeqs[counter]
        cods <- unlist(cods)
        cods <- c(namecods, cods)
        temp <- paste(cods, collapse = "\t")
        write(temp, outputFileDist, append = TRUE)
      })
    }
    else {
      stop("ERROR: outFormat should be 'mat' or 'txt' ")
    }
  }
  
  
  
  ###########################GC_Content##########################
  GC.content <- function(fasta_file){
    x <- read.fasta(file=fasta_file)
    tt<-function(x){
      res<-GC(x)
      val=round(res,4)
      return(val)
    }
    
    f_res<-lapply(x,tt)
    s=data.frame(f_res)
    
    rownames(s) <- c("GC-content")
    
    w=t(s)
    return(w)
  }
  ################################ONF#################################
  oligo.freq <- function(fasta_file, f){
    
    x <- seqinr::read.fasta(fasta_file, as.string = TRUE, seqtype = "DNA")
    x <- sapply(x, toupper)
    
    bases <- c("A","C","G","T")
    kmers <- expand.grid(rep(list(bases), f))
    kmers <- apply(kmers, 1, paste0, collapse = "")
    
    calc <- function(s){
      tab <- table(factor(sapply(1:(nchar(s)-f+1),
                                 function(i) substr(s, i, i+f-1)), levels = kmers))
      as.numeric(tab)
    }
    
    y <- t(sapply(x, calc))
    z <- as.data.frame(y)
    colnames(z) <- kmers
    rownames(z) <- names(x)
    
    return(z)
  }
  
  
  #########################################AMIP################################
  
  
  #######################mononucleotide_binary_encoding##################################
  
  FastaToTabular <- function (filename){
    
    #read fasta file
    
    file1 <- readLines(filename)
    
    #find the genename location by grepping >
    
    location <- which((substr(file1,1,1))==">")
    
    #start an empty vector to collect name and sequence
    
    name=c()
    sequence =c()
    
    
    
    #number of genes= number of loops
    #extract name first
    for ( i in 1:length(location)){
      name_line = location[i]
      name1 = file1[name_line]
      name=c(name,name1)
      #extract sequence between the names
      #the last sequence will be missed using this strategy
      #so, we are using if condition to extract last sequence
      start= location[i]+1
      end = location[i+1]-1
      if ( i < length (location)){
        
        end=end
        
      } else {
        
        end=length(file1)
      }
      
      lines = start:end
      sequence1= as.character(paste(file1[lines],collapse = ""))
      sequence =c(sequence,sequence1)
    }
    
    #now create table using name and sequence vector
    
    data <- data.frame(name = name, sequence = sequence, stringsAsFactors = FALSE)
    
    
    
    
    #finally export the file
    #before that remove preexisting file
    unlink(c("dna_table.csv"),force=TRUE)
    as.matrix(data,"dna_table.csv")
    
    #function ends
  }
  #########################alphabetcheck###########################
  alphabetCheck<-function (sequences, alphabet = "aa", label = c())
  {
    if (length(sequences) == 0) {
      stop("ERROR: sequence parameter is empty")
    }
    if (length(label) != 0 && length(label) != length(sequences)) {
      stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
    }
    if (alphabet == "rna") {
      alphabet <- c("A", "C", "G", "U")
    }
    else if (alphabet == "dna") {
      alphabet <- c("A", "C", "G", "T")
    }
    else if (alphabet == "aa") {
      alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I",
                    "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                    "W", "Y")
    }
    else {
      stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
    }
    alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                               split = "")[[1]] %in% alphabet))
    flag = 0
    if (length(label) == length(sequences)) {
      flag = 1
      label = label[alphabetCheck]
    }
    else if (length(label) > 0 && length(label) != length(sequences)) {
      stop("ERROR: The number of labels is not equal to the number of sequences!")
    }
    if (is.null(names(sequences))) {
      names(sequences) <- as.character(1:length(sequences))
    }
    nonstanSeq <- names(sequences)[!alphabetCheck]
    if (length(nonstanSeq) != 0) {
      nonstanSeq <- toString(nonstanSeq)
      warMessage <- paste("The sequences (", nonstanSeq, ") were deleted. They contained non-standard alphabets")
      message(warMessage)
    }
    sequences = sequences[alphabetCheck]
    if (length(sequences) == 0) {
      stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
    }
    if (flag == 1) {
      names(label) = names(sequences)
    }
    seq_lab <- list(sequences = sequences, Lab = label)
    return(seq_lab)
  }
  #################################MBE_DNA############################
  
  mbe_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                     label = c())
  {
    if (length(seqs) == 1 && file.exists(seqs)) {
      seqs <- fa.read(seqs, alphabet = "dna")
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else if (is.vector(seqs)) {
      seqs <- sapply(seqs, toupper)
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else {
      stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
    }
    lenSeqs <- sapply(seqs, nchar)
    nucs <- list(A = c(1, 0, 0, 0), C = c(0, 1, 0, 0), G = c(0, 0, 1, 0), T = c(0, 0, 0, 1), U = c(0, 0, 0, 1))
    numSeqs <- length(seqs)
    if (outFormat == "mat") {
      if (length(unique(lenSeqs)) > 1) {
        stop("ERROR: All sequences should have the same length in 'mat' mode. For sequences with different lengths, please use 'txt' for outFormat parameter")
      }
      if (binaryType == "strBin") {
        nucs <- c(A = "1000", C = "0100", G = "0010", T = "0001",
                  U = "0001")
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        colnames(featureMatrix) <- paste("pos_mbe", 1:lenSeqs[1],
                                         sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "logicBin") {
        nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                    TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                   FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("pos_mbe", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "numBin") {
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("pos_mbe", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else {
        stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
      }
      return(featureMatrix)
    }
    else if (outFormat == "txt") {
      nucs <- c(A = "1000", C = "0100", G = "0010", T = "0001",
                U = "0001")
      counter <- 0
      namesSeqs <- names(seqs)
      codes <- lapply(seqs, function(x) {
        counter <- counter + 1
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        namecods <- namesSeqs[counter]
        cods <- unlist(cods)
        cods <- c(namecods, cods)
        temp <- paste(cods, collapse = "\t")
        write(temp, outputFileDist, append = TRUE)
      })
    }
    else {
      stop("ERROR: outFormat should be 'mat' or 'txt' ")
    }
  }
  #########EIIP############
  
  
  EIIP_final<-EIIP(fasta_file)
  
  ############PSTNP###############
  
  
  
  # Define a substring or keyword to search for
  # find patterns in names
  posSeqs <- suppressWarnings(
    ftrCOOL::fa.read(url("https://zenodo.org/records/17934726/files/AT__Pos.fasta?download=1"))
  )
  
  negSeqs <- suppressWarnings(
    ftrCOOL::fa.read(url("https://zenodo.org/records/17934726/files/AT_Neg.fasta?download=1"))
  )
  seqs <- fa.read(file = fasta_file, alphabet = "dna")
  PSTNPss_Final <- PSTNPss_DNA(seqs = seqs, pos = posSeqs, neg = negSeqs)
  
  
  rf_model <- readRDS(url("https://zenodo.org/records/17934726/files/rf_AT.rds?download=1"))
  xgb_model <- readRDS(url("https://zenodo.org/records/17934726/files/xgb_AT.rds?download=1"))
  svm_model <- readRDS(url("https://zenodo.org/records/17934726/files/svm_AT.rds?download=1"))
  
  
  #################################data_preparation######################
  res<-FastaToTabular(fasta_file)
  data<-as.vector(res[,2])
  mat<-as.matrix(ncp_dna(seqs = data,binaryType="strBin",outFormat="mat"))
  sequence<-rownames(mat)
  seq_id<-res[,1]
  ncp<-cbind(seq_id,sequence,mat)
  rownames(ncp)<-seq_id
  ncp_temp<-data.frame(ncp[,-1], stringsAsFactors = FALSE)
  ncp_final<-as.data.frame(apply(ncp_temp[,-1], 2, as.numeric))
  log_gc_temp<-log((GC.content(fasta_file))*100, base = exp(1))
  log_gc<-as.data.frame(as.numeric((ifelse(log_gc_temp>0,log_gc_temp,'0'))))
  onf<-oligo.freq(fasta_file, 2)
  res_temp_mbe<-FastaToTabular(fasta_file)
  data_mbe<-as.vector(res_temp_mbe[,2])
  res_mbe<-as.matrix(mbe_dna(seqs = data_mbe,binaryType="strBin",outFormat="mat"))
  mbe_temp<-data.frame(res_mbe[,-1], stringsAsFactors = FALSE)
  mbe_final<-as.data.frame(apply(mbe_temp[,-1], 2, as.numeric))
  # AMIP_final<- NULL
  # for (i in 1:6) {
  #   for(j in 2:7){
  #     if(i<j){
  #       AMIP1<-AMIP(fasta_file,i,j)
  # 
  #       colnames(AMIP1) <- paste('Mean',i,j, sep="_")
  #       AMIP_final<-cbind(AMIP_final, AMIP1)
  #     }
  #   }
  # }
  
  temp1<- cbind(onf, gcc =log_gc[,1], ncp_final, mbe_final,  EIIP_final, PSTNPss_Final)
  my_data_temp<- temp1
  inputData <-as.data.frame(my_data_temp)
  
  selected_columns <- c("21",	"20",	"19",	"ncp_pos21",	"pos_mbe21",
                        "pos21",	"22",	"23",	"24",	"25",	"26",	"27",	
                        "ncp_pos24",	"pos24",	"18",	"pos27",	"ncp_pos27",	
                        "pos_mbe27",	"3",	"2",	"pos25",	"ncp_pos25",	"pos_mbe25",
                        "17",	"pos_mbe22",	"pos22",	"ncp_pos22",	"16",	"pos20",
                        "ncp_pos20",	"pos_mbe20",	"gcc",	"15",	"TC",	
                        "TA",	"TT",	"ncp_pos23",	"pos_mbe23")
  
  data_temp <- inputData[,selected_columns]
  test_data<-data_temp
  #test_data<-cbind(Sequence=ncp_temp$sequence, data_temp)
  
  
  
  
  ############################Prediction########################
  predicted_prob_svm <- predict(svm_model, newdata = test_data, type = "prob")
  #predicted_prob_svm <- attr(predicted_prob_svm_temp, "probabilities")
  predicted_value_svm <- predict(svm_model, newdata = test_data)
  predicted_prob_rf <- predict(rf_model, newdata = test_data, type = "prob") 
  predicted_value_rf <- predict(rf_model, newdata = test_data)
  predicted_prob_xgb <- predict(xgb_model, newdata = test_data, type = "prob")  
  predicted_value_xgb <- predict(xgb_model, newdata = test_data)
  test_res_en_prob <-  cbind(SVM = predicted_prob_svm, RF = predicted_prob_rf, XGB = predicted_prob_xgb)
  
  
  ##################Ensemble###################
  # Define weights
  weights <- c(SVM = 0.19070627, RF = 0.78879226, XGB = 0.02050146)
  
  # Apply weights
  weighted_data <- as.data.frame(test_res_en_prob)
  weighted_data[1:3] <- test_res_en_prob[1:3] * weights["SVM"]
  weighted_data[4:6] <- test_res_en_prob[4:6] * weights["RF"]
  weighted_data[7:9] <- test_res_en_prob[7:9] * weights["XGB"]
  
  weighted_data$Class_1 <- rowSums(weighted_data[c("SVM.X1", "RF.X1", "XGB.1")])
  weighted_data$Class_2 <- rowSums(weighted_data[c("SVM.X2", "RF.X2", "XGB.2")])
  weighted_data$Class_3 <- rowSums(weighted_data[c("SVM.X3", "RF.X3", "XGB.3")])
  
  
  
  weighted_data$Predicted_Class <- apply(weighted_data[, c("Class_1", "Class_2", "Class_3")], 1, which.max)
  
  
  final_results <- weighted_data[, c("Class_1", "Class_2", "Class_3", "Predicted_Class")]
  
  final_results$Predicted_Class_Label <- ifelse(final_results$Predicted_Class == 2, "DNA 6mA Methylation",
                                                ifelse(final_results$Predicted_Class == 3, "DNA 4mC Methylation",
                                                       ifelse(final_results$Predicted_Class == 1, "No Methylation", NA)))
  Ids <- rownames(ncp_temp)
  Ids <- sub("^>", "", Ids)
  
  final_pred <- data.frame(Ids= Ids,
                           Sequence = ncp_temp[,1],
                           Modification = final_results$Predicted_Class_Label,
                           Probability = round(apply(final_results[, c("Class_1", "Class_2", "Class_3")], 1, max), 2))
  rownames(final_pred) <- NULL
  return(final_pred)
}

live_users  <- 0
total_users <- 0


# Define server function for the Shiny app
server <- function(input, output, session) {
  
  
  # ---- USER COUNTERS ----
  # ---- USER COUNTERS ----
  live_users  <<- live_users + 1
  total_users <<- total_users + 1
  
  session$onSessionEnded(function() {
    live_users <<- live_users - 1
  })
  
  
  output$downloadButton1 <- downloadHandler(
    filename = function() {
      "Example_Data.zip"  # File name for the downloaded .zip file
    },
    content = function(file_example) {
      # Load the content of the .zip file
      zipContent <- readBin("data/Example_Data.zip", "raw", file.info("data/Example_Data.zip")$size)
      
      # Provide the .zip file for download
      file.copy("data/Example_Data.zip", file_example)
     
    }
  )

  #observe({
    #print(input$folder)}
#)
  # Shiny app logic
  
    
  datasetInput <- reactive({
    req(input$fasta_file)
    
    dat <- input$fasta_file
    data<-dat$datapath
    
    res_monocot <- NULL
    res_dicot <- NULL
    
    if (input$Reference_type %in% c("Rice", "both")) {
      res_monocot <- Monocot(fasta_file_path=data)
    }
    
    if (input$Reference_type %in% c("Arabidopsis", "both")) {
      res_dicot <- Dicot(fasta_file_path=data)
    }
    
    list(PredMonocot=res_monocot, PredDicot=res_dicot)
  })
  
  output$msmono <- renderTable({
    req(input$submitbutton)
    res <- datasetInput()
    res$PredMonocot
  })
  
  output$msdi <- renderTable({
    req(input$submitbutton)
    res <- datasetInput()
    res$PredDicot
  })
  
  output$downloadmonocot <- downloadHandler(
    filename = function() {
      
      # Get uploaded file name (without extension)
      input_name <- tools::file_path_sans_ext(input$fasta_file$name)
      
      paste0(input_name, "_DNA_Methylation_Status_Monocot_", Sys.Date(), ".csv")
    },
    
    content = function(file) {
      res <- datasetInput()
      write.csv(res$PredMonocot, file, row.names = FALSE)
    }
  )
  
  
  output$downloaddicot <- downloadHandler(
    filename = function() {
      
      # Get uploaded file name (without extension)
      input_name <- tools::file_path_sans_ext(input$fasta_file$name)
      
      paste0(input_name, "_DNA_Methylation_Status_Dicot_", Sys.Date(), ".csv")
    },
    
    content = function(file) {
      res <- datasetInput()
      write.csv(res$PredDicot, file, row.names = FALSE)
    }
  )
  
  # Progress status
  observeEvent(input$submitbutton, {
    progress <- Progress$new(session, min = 0, max = 100)
    progress$set(message = 'Predicting...')
    
    # Run the prediction process in a separate reactive conductor
    predictionDone <- reactiveVal(FALSE)
    observeEvent(predictionDone(), {
      progress$close()
    })
    
    observeEvent(input$submitbutton, {
      withProgress(message = "Predicting...", value = 0, {
        
        # break your task into steps
        n <- 100
        for (i in 1:n) {
          Sys.sleep(0.05)  # simulate computation
          incProgress(1/n, detail = paste0("Step ", i, " of ", n))
        }
        
        # now run the actual prediction
        res <- datasetInput()
        
        # save results to a reactive value if needed
        # results(res)
      })
    })
      })
  
  # ---- force live UI refresh for counters ----
  autoInvalidate <- reactiveTimer(1000)   # every 1 second
  
  output$liveCountText <- renderText({
    autoInvalidate()
    paste("👥 Live Users:", live_users)
  })
  
  output$totalCountText <- renderText({
    autoInvalidate()
    paste("🌐 Total Users:", total_users)
  })
  
  observe({
    session$sendCustomMessage("counts", list(
      live  = paste("👥 Live Users:", live_users),
      total = paste("🌐 Total Users:", total_users)
    ))
  })
  
  
  
}

