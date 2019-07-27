
# -- This file exists to automate the process of normalizing Illumina 450k methylation data 
#    Also performs filtering of data based on acceptable # of NA methylation values across sets

setwd("/rafalab/qbet/dna_met")

require("tidyverse")
require("minfi")
require("impute")
require("magrittr")


# -- List of 450k array datasets from Horvath 2013 Publication (With additional datasets to improve robustness)
#    NOTE: GSE42865 did NOT incude (parsable) age data and are samples from patients with Werner Syndrome
meta_data <- read_csv("./cell_tissue_datasets.csv")
dataset_list <- meta_data$Availability[meta_data$Platform == "450K" & meta_data$Availability != "GSE42865"]


# -- Illumina 450k Reference Datasheet (for referencing genome positions)
illumina_sheet <- read_csv("./Platform Information/GPL13534_HumanMethylation450_15017482_v.1.1.csv")


# Terminate if output file already exists
data_filepath = "./Outputs/Normalized_450k.csv"
nonorm_filepath ="./Outputs/Unnormalized_450k.csv"
age_filepath = "./Outputs/Age_450k.csv"


if (file.exists(data_filepath) | file.exists(age_filepath)) {
    print("One or more output files already exist.")
    quit(status = 0)
}

# -- This loop exists solely to join all datasets (and their phenotypic data) together.
for (data_id in dataset_list) {

    # -- Locating and loading in altered data matrix & converting to GenomicRatioSet
    if (startsWith(data_id, "TCGA")) {
        datapath = paste("./Healthy/", data_id, "_450K.csv", sep = "")
        pheno_data = paste("./TCGA Age Data/", data_id, "-Clinical.tsv", sep = "")
        raw_read   <- read_csv(datapath)
        raw_matrix <- as.matrix(raw_read[,2:ncol(raw_read)])
        rownames(raw_matrix) <- as_vector(raw_read[,1])
        setG <- makeGenomicRatioSetFromMatrix(as.matrix(raw_matrix))
        rm(raw_read, raw_matrix)
        
    } else if (startsWith(data_id, "GSE")) {        # NOTE: Most, but not all datasets are uploaded as beta values.
                                                    # The relationship between M and Beta is logistic and computed automatically by minfi
        if (data_id == "GSE34639") { 
            setG <- getGenomicRatioSetFromGEO(data_id, what = "M") 
        } else {
            setG <- getGenomicRatioSetFromGEO(data_id, what = "Beta")
        }
        
    } else {
        print("Dataset must be GSE or modified TCGA format")
        quit(status = 1)
    }
    
    print(paste("Dataset:", data_id, "loaded."))
        
    
    # -- Getting sample and probe names
    sampleNames <- sampleNames(setG)
    probeNames  <- featureNames(setG)
    print("Sample and probe names loaded.")
    
    
    
    
    # -- Phenotype data (including age, sex, etc.)
    if (startsWith(data_id, "TCGA")) {
        pheno <- read.delim(pheno_data, stringsAsFactors = FALSE)
        pheno <- pheno[3:nrow(pheno),]
        
    } else if (startsWith(data_id, "GSE")) {
        pheno <- pData(setG)
    }
    print("Phenotype data loaded")
    
    
    # -- Getting age data for each sample
    if (startsWith(data_id, "TCGA")) {
        tmp <- pheno %>% 
            as_tibble() %>%
            select(bcr_patient_barcode, starts_with("age_at")) %>%
            setNames(c("Sample", "Age")) %>%
            filter(Sample %in% sampleNames) %>%
            arrange(match(Sample, sampleNames)) %>%
            mutate(Age = as.numeric(Age))
        
    } else if (startsWith(data_id, "GSE")) {
        tmp <- pheno %>% 
            as_tibble() %>%
            mutate(Sample = geo_accession[title %in% sampleNames]) %>%
            select(Sample, Age = contains("age")) %>%
            setNames(c("Sample", "Age")) %>%
            mutate(Age = as.numeric(Age)) 
    }
    
    print("Sample ages retrieved")
    sampleList <- tmp$Sample
    
    
    # -- Standardizing GSE dataset age data to years
    if (startsWith(data_id, "GSE")) {
        age_type = colnames(pheno[contains("age", var = colnames(pheno))])
        
        if (grepl("day", age_type)) {
            tmp$Age <- tmp$Age / 365
            
        } else if (grepl("week", age_type)) {
            tmp$Age <- tmp$Age * (7 / 365)
            
        } else if (grepl("month", age_type)) {
            tmp$Age <- tmp$Age / 12
        } 
    } 
    
    
    # -- Saving age data in years for all samples
    if (exists("age_data")) {
        age_data <- age_data %>%
            full_join(tmp)
        
    } else {
        age_data <- tmp
    }
    print("Age data joined to age data matrix") 
    
    
    
    
    # -- Genome position for each probe
    #    This function returns the chromosome & position for each probe.
    gr   <- granges(setG)
    
    
    # -- Obtaining beta values
    #    This are the methylation proportion values
    beta <- getBeta(setG)

    
    # -- Adding probe names to the beta matrix
    colnames(beta) <- sampleList
    tmp <- beta %>%
        as_tibble() %>%
        mutate(Probe = probeNames)
    
    print("Beta matrix labeled")
    
    
    # -- Saving methylation data for all samples organized by Probe
    if (exists("met_dat")) {
        met_dat <- met_dat %>%
            left_join(tmp, by = "Probe")
        
    } else {
        met_dat<- tmp
    }
    print("Beta values joined to methylation matrix.")
    
}



# -- Freeing up data from loop
rm(setG, beta, pheno, pheno_data, age_type, data_id, sampleNames, sampleList)


# -- Removing age and methylation data for samples with missing age values 
missing_age_samples <- age_data$Sample[is.na(age_data$Age)]

age_data <- age_data %>%
    filter(!is.na(Age))

met_dat <- met_dat %>%
    select(-missing_age_samples)
    

# -- Preparing genomic data for Illumina Probe data
tmp <- gr %>%
    as_tibble() %>%
    select(seqnames, start) %>%
    mutate(Probe = probeNames) %>%
    setNames(c("Chr", "Pos", "Probe"))


# -- Getting Illumina probe information and arranging by genomic position
probeInfo <- illumina_sheet %>% 
    select(IlmnID, Infinium_Design_Type) %>%
    setNames(c("Probe", "Probe_type")) %>%
    left_join(tmp, by = "Probe") %>%
    filter(Probe %in% met_dat$Probe) %>%
    arrange(Chr, Pos)

probeInfo$Chr <- str_replace(probeInfo$Chr, "chr", "")


# -- Arranging methylation data by genomic position
met_dat <- met_dat %>%
    left_join(probeInfo, by = "Probe") %>%
    arrange(Chr, Pos)


# -- Documenting probes with less than or equal to the allowed number of NA values in methylation data
#    Probes not meeting this condition will be removed after normalization
NA_ALLOWED = 20
NA_vector <- apply(met_dat, MARGIN = 1, FUN = function(x) sum(is.na(x)))


# -- Imputing NAs before BMIQ normalization
#    Imputation is done using k nearest neighbors designed for gene expression values
met_dat <- met_dat %>%
    select(-Probe, -Probe_type, -Chr, -Pos) %>%
    as.matrix() %>%
    impute.knn() %>%
    use_series(data) %>%
    as_tibble()


# -- Saving age and unnormalized data in the event of normalization failure
write_csv(age_data, path = age_filepath)
probeInfo %>%
    bind_cols(met_dat) %>%
    arrange(Chr, Pos) %>%
    filter(NA_vector <= NA_ALLOWED) %>%
    write_csv(path = nonorm_filepath)

    
# -- Checking some details about the data for the sake of bug-fixing
residual_NAs <- sum(apply(met_dat, MARGIN = 1, FUN = function(x) sum(is.na(x))))
residual_Ms <- sum(apply(met_dat, MARGIN = 1, FUN = function(x) sum(x > 1 | x < 0, na.rm = TRUE)))
probe_NAs <- sum(is.na(probeInfo$Probe_type))
sample_NAs <- sum(is.na(colnames(met_dat)))
print(paste("Residual NAs:", residual_NAs, "Probe Type NAs:", probe_NAs, "Sample Name NAs:", sample_NAs))
print(paste("Residual M (Not beta) values:", residual_Ms))


# -- BMIQ v1.3 normalization function developed by Andrew Teschendorff
#    Corrects for Illumina Type II probe bias using a Beta Mixture Model. Also used in Horvath 2013.
source("../scripts/BMIQ_1.3.R")
met_dat <- apply(met_dat, MARGIN = 2, function(x) use_series(BMIQ(x, probeInfo$Probe_type, plots = FALSE), nbeta))
met_dat <- as_tibble(met_dat)


# -- Joining normalized methylation values to probe information and removing probes with too many NAs
met_dat <- probeInfo %>%
    bind_cols(met_dat) %>%
    arrange(Chr, Pos) %>%
    filter(NA_vector <= NA_ALLOWED)


# -- Saving files, writing dimensions to output for documentation sake
write_csv(met_dat, path = data_filepath)
print(paste("Saving NORMALIZED", nrow(met_dat),"probes from n =", ncol(met_dat), "different samples with less than", NA_ALLOWED, "missing values across all samples."))


