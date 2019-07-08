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
if (file.exists(paste("./Outputs/", horvath_cpg, "_correlation_.png", sep = ""))) {
    print("Plot file already exists")
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
        
    } else if (startsWith(data_id, "GSE")) {
        setG <- getGenomicRatioSetFromGEO(data_id)
        
    } else {
        print("Dataset must be GSE or modified TCGA format")
        quit(status = 1)
    }
        
    # -- Name of probes and samples
    probeNames  <- featureNames(setG)
    sampleNames <- sampleNames(setG)
    
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
    
    # -- Adding probe names to the beta matrix
    dat <- beta %>%
        as_tibble() %>%
        mutate(Probe = probeNames)
    
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
    window_size = 10000
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
    rm(raw_read, raw_matrix, setG, probeNames, sampleNames, beta, dat, tmp, gr, pheno)
    
}
    
# -- Computing correlation between Age and Methylation
tmp <- full_dat %>%
    filter(!is.na(Met)) %>%
    filter(!is.na(Age)) %>%
    group_by(Probe, Rel_Pos) %>%
    mutate(Cor = cor(Met, Age)) %>%
    ungroup()

# -- Only taking relative position on chromosome and correlation value
myd <- tmp %>%
    select(Rel_Pos, Cor) %>%
    unique() 

# -- Plotting the results and saving the plot
horvath_Cor = myd$Cor[myd$Rel_Pos == 0]

plot <- ggplot(data = myd, aes(Rel_Pos, Cor)) +
    geom_point(alpha=0.30) +
    geom_point(aes(0, horvath_Cor), color = "red", size = 3) +
    geom_hline(yintercept = 0, color="red", lty=2) +
    ggtitle(paste("Correlation of CpGs Local to", horvath_cpg)) +
    xlab(paste("Position on Chromosome", horvath_clock$Chr[horvath_clock$CpGmarker == horvath_cpg],
               "Relative to Position", horvath_Pos)) +
    ylab("Pearson Correlation") +
    scale_y_continuous(limits = c(-1,1),
                       breaks = seq(-1,1,by=0.20)) +
    scale_x_continuous(limits = c(-window_size, window_size),
                       breaks = seq(-window_size, window_size, by = window_size / 10)) +
    theme_minimal() +
        theme(plot.margin = unit(c(0.75, 1, 0.75, 0.75), "cm"),
        plot.title=element_text(hjust=0.5, vjust=2, size = 24),
        axis.title.x=element_text(vjust=-2, size = 16),
        axis.title.y=element_text(angle = 90, vjust=3, size = 16))

ggsave(paste(horvath_cpg, "_correlation_.png", sep = ""), 
       width = 16, height = 9, dpi = 150, type = "cairo", path = "./Outputs")






