
# -- Run bumphunter on Illumina 450k datasets w.r.t. age data

setwd("/rafalab/qbet/dna_met")

require(tidyverse)
require(minfi)
require(doParallel)
registerDoParallel(cores = 4)

# -- Methylation data matrix file (normalized) and age data matrix
data_filepath = "D:/Datasets/Horvath 2013 Data/Collated Data/Test_Unnormalized_450k.csv"
age_filepath = "D:/Datasets/Horvath 2013 Data/Collated Data/Test_Age_450k.csv"

# -- Loading in collated methylation data processed by processing/normalization pipeline
data_tibble <- read_csv(data_filepath)
genome_tibble <- select(data_tibble, Probe, Pos, Chr, Probe_type)
beta_tibble <- select(data_tibble, -Probe, -Pos, -Chr, -Probe_type)
NA_vector <- apply(beta_tibble, MARGIN = 1, FUN = function(x) sum(is.na(x)))
print(paste("The data matrix has a total of", sum(NA_vector), "missing values."))
beta_matrix <- as.matrix(beta_tibble)
sampleList <- colnames(beta_matrix)
rownames(beta_matrix) <- data_tibble$Probe
setG <- makeGenomicRatioSetFromMatrix(as.matrix(beta_matrix))
rm(data_tibble, beta_tibble, beta_matrix, NA_vector)


# -- Loading in collated age data and arranging to match order of samples in beta value matrix
age_data = read_csv(age_filepath)
age_data <- age_data %>%
    filter(Sample %in% sampleList) %>%
    arrange(match(Sample, sampleList))

# -- Creating age design matrix for bump hunting. Order of ages 
ages <- age_data$Age
designMatrix <- model.matrix(~ ages)
rownames(designMatrix) <- age_data$Sample

# TODO: Make Bumphunter(blockFinder) actually work
age_dmrs <- minfi::bumphunter(setG, design = designMatrix, cutoff = 0.01, type = "Beta")

clusters <- clusterMaker(genome_tibble$Chr, genome_tibble$Pos, maxGap = 1000)
block_setG <- cpgCollapse(setG, returnBlockInfo = FALSE)
block_age_dmrs <- blockFinder(setG, design = designMatrix, cluster = clusters, cutoff = 0.1, B=0, what = "Beta")

