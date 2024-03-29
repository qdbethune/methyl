
# -- Script for subsetting sample methylation data centered at a single CpG of interest

setwd("/rafalab/qbet/dna_met")

require("tidyverse")
require("minfi")

# -- List of 450k array datasets
horvath_meta <- read_csv("./healthy_tissue_datasets.csv")
horvath_450k_sets <- horvath_meta$Availability[horvath_meta$Platform == "450K" & horvath_meta$`Data Use` != "Other" & horvath_meta$Availability != "GSE42865"]
clock_data = "./Horvath_Clock_CpGs.csv"

# -- Horvath CpGs (Vector of identifiers and Clock coefficients)
horvath_clock <- read_csv(clock_data)
horvath_ids <- horvath_clock$CpGmarker[2:length(horvath_clock$CpGmarker)]

# -- Illumina 450k Reference Datasheet (for referencing genome positions)
illumina_sheet <- read_csv("./Platform Information/GPL13534_HumanMethylation450_15017482_v.1.1.csv")

# -- Getting the current CpG of interest from sys.arg (SLURM batch job)
args <- commandArgs()
horvath_cpg = horvath_ids[as.integer(args[length(args)])]

# If resuming a prior operation, skip any plots which have already been produced 
filepath = paste("./Outputs/", horvath_cpg, "_wind100k.csv", sep = "")
if (file.exists(filepath)) {
    print("Data file already exists")
    quit(status = 0)
}

# -- This loop filters data from each dataset before joining them together
for (data_id in horvath_450k_sets) {

    # -- Locating and loading in altered data matrix & converting to GenomicRatioSet
    if (startsWith(data_id, "TCGA")) {
        datapath = paste("./Healthy/", data_id, "_450K.csv", sep = "")
        pheno_data = paste("./TCGA Age Data/", data_id, "-Clinical.tsv", sep = "")
        raw_read   <- read_csv(datapath)
        raw_matrix <- as.matrix(raw_read[,2:ncol(raw_read)])
        rownames(raw_matrix) <- as_vector(raw_read[,1])
        setG <- makeGenomicRatioSetFromMatrix(as.matrix(raw_matrix))
        rm(raw_read, raw_matrix)
        
    } else if (startsWith(data_id, "GSE")) {
        setG <- getGenomicRatioSetFromGEO(data_id)
        
    } else {
        print("Dataset must be GSE or modified TCGA format")
        quit(status = 1)
    }
        
    # -- Name of probes and samples
    probeNames  <- featureNames(setG)
    sampleNames <- sampleNames(setG)
    probeTypes <- illumina_sheet %>% 
                select(IlmnID, Infinium_Design_Type) %>%
                setNames(c("Probe", "Probe_type"))
    
    # -- Phenotype data (including age, sex, etc.)
    #    This function allows you access all phenotypic information available
    if (startsWith(data_id, "TCGA")) {
        pheno <- read.delim(pheno_data, stringsAsFactors = FALSE)
        pheno <- pheno[3:nrow(pheno),]
        
    } else if (startsWith(data_id, "GSE")) {
        pheno <- pData(setG)
    } 
    
    # -- Genome position for each probe
    #     This function returns the position for each probe.
    #     It may be the case that multiple probes were used for a 
    #     a single CpG. Thus, this data structure may be smaller
    gr   <- granges(setG)
    
    # -- Obtaining beta values
    #     This are the methylation proportion values
    beta <- getBeta(setG)
    
    # -- Adding probe names and probe types to the beta matrix
    dat <- beta %>%
        as_tibble() %>%
        mutate(Probe = probeNames) %>%
        right_join(probeTypes, by = "Probe")
    
    # -- Probe genome location
    tmp <- gr %>%
        as_tibble() %>%
        select(seqnames, start) %>%
        mutate(Probe = probeNames) %>%
        setNames(c("Chr", "Pos", "Probe"))
    
    tmp$Chr <- str_replace(tmp$Chr, "chr", "")
    
    # -- Joining dat and tmp
    dat <- dat %>%
        left_join(tmp, by = "Probe")
    
    # -- Filtering based on chromosome and distance from CpG of interest
    horvath_Pos <- illumina_sheet$MAPINFO[illumina_sheet$IlmnID == horvath_cpg]
    window_size = 100000
    dat <- dat %>%
        filter(Chr == horvath_clock$Chr[horvath_clock$CpGmarker == horvath_cpg] & abs(Pos - horvath_Pos) <= window_size) %>%
        mutate(Rel_Pos = Pos - horvath_Pos)
    
    # -- Getting age data for each sample
    if (startsWith(data_id, "TCGA")) {
        tmp <- pheno %>% 
            as_tibble() %>%
            select(bcr_patient_barcode, starts_with("age_at")) %>%
            setNames(c("Sample", "Age")) %>%
            mutate(Age = as.numeric(Age))
        
    } else if (startsWith(data_id, "GSE")) {
        tmp <- pheno %>% 
            as_tibble() %>%
            mutate(Sample = rownames(pheno)) %>%
            select(Sample, Age = contains("age")) %>%
            mutate(Age = as.numeric(Age)) %>%
            setNames(c("Sample", "Age"))
    }
    
    # -- Standardizing GSE age data to years
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
    
    # -- Joining the genomic data with the phenotype data
    dat <- dat %>%
        gather(Sample, Met, -c(Probe, Pos, Rel_Pos, Chr)) %>%
        left_join(tmp, by = "Sample")
    
    # -- Creating/extending tibble containing filtered data for analysis
    if (exists("full_dat")) {
        full_dat <- full_dat %>%
            full_join(dat)
        
    } else {
        full_dat <- dat
    }
    
    # -- Removing objects to prepare for next dataset (only full_dat is needed)
    rm(setG, probeNames, sampleNames, beta, dat, tmp, gr, pheno)
    
}
    
write_csv(full_dat, path = filepath)





