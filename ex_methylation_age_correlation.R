#NOTE: Much of this is hard-coded right now. Changes will be required upon scaling.

setwd("/rafalab/qbet/dna_met")

# -- Install dependencies on current node
install.packages("tidyverse", repos='http://cran.us.r-project.org')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

require("tidyverse")
require("minfi")

# -- List of 450k array datasets
horvath_meta_path = "./healthy_tissue_datasets.csv"
horvath_meta <- read_csv(horvath_meta_path)
horvath_450k_sets <- horvath_meta$Availability[horvath_meta$Platform == "450K" & horvath_meta$`Data Use` != "Other"]
clock_data = "./Horvath_Clock_CpGs.csv"

# TODO: LOOP THROUGH ALL DATASETS
data_id = horvath_450k_sets[1]

# -- Locating and loading in altered data matrix & converting to GenomicRatioSet
if (startsWith(data_id, "TCGA")) {
    datapath = paste("./Healthy/", data_id, "_450K.csv", sep = "")
    pheno_data = paste("./TCGA Age Data/", data_id, "-Clinical.tsv", sep = "")
    raw_read   <- read_csv(datapath)
    raw_matrix <- as.matrix(raw_read[,2:ncol(raw_read)])
    rownames(raw_matrix) <- as_vector(raw_read[,1])
    setG <- makeGenomicRatioSetFromMatrix(as.matrix(raw_matrix))
} else if (startsWith(data_id, "GSE")) {
    setG <- getGenomicRatioSetFromGEO(datapath)
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

# -- Horvath CpGs (Vector of identifiers and Clock coefficients)
horvath_clock <- read_csv(clock_data)

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
horvath_Probe <- as.character(horvath_clock[2, 1])
horvath_Pos <- dat$Pos[dat$Probe == horvath_Probe]
window_size = 50000
dat <- dat %>%
    filter(Chr == horvath_clock$Chr[2] & abs(Pos - horvath_Pos) <= window_size) %>%
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

# -- Computing correlation between Age and Methylation
tmp <- dat %>%
    filter(!is.na(Met)) %>%
    filter(!is.na(Age)) %>%
    group_by(Probe, Rel_Pos) %>%
    mutate(Cor = cor(Met, Age)) %>%
    ungroup()

# -- Only taking relative position on chromosome and correlation value
myd <- tmp %>%
    select(Rel_Pos, Cor) %>%
    unique() 

# -- Plotting the results
horvath_Cor = myd$Cor[myd$Rel_Pos == 0]

plot <- ggplot(data = myd, aes(Rel_Pos, Cor)) +
    geom_point(alpha=0.30) +
    geom_point(aes(0, horvath_Cor), color = "red", size = 3) +
    geom_hline(yintercept = 0, color="red", lty=2) +
    ggtitle(paste("Correlation of CpGs Local to", horvath_clock$CpGmarker[2])) +
    xlab(paste("Position on Chromosome", horvath_clock$Chr[2])) +
    ylab("Pearson Correlation") +
    scale_y_continuous(limits = c(-1,1),
                       breaks = seq(-1,1,by=0.20)) +
    scale_x_continuous(limits = c(-window_size, window_size),
                       breaks = seq(-window_size, window_size, by = 10000)) +
    theme_minimal() +
        theme(plot.margin = unit(c(0.75, 1, 0.75, 0.75), "cm"),
        plot.title=element_text(hjust=0.5, vjust=2, size = 24),
        axis.title.x=element_text(vjust=-2, size = 16),
        axis.title.y=element_text(angle = 90, vjust=3, size = 16))

ggsave("plot.png", width = 16, height = 9, dpi = 150, type = "cairo", path = "./Outputs")
