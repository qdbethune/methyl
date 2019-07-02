#NOTE: Much of this is hard-coded right now. Changes will be required upon scaling.

require(tidyverse)
require(minfi)

# -- Data location (The "path" for GSE files should just be the identifier [ex: GSE41169])
datapath = "GSE44667"
pheno_data = "D:\\Datasets\\Horvath 2013 Data\\TCGA Age Data\\TCGA-BRCA-Clinical.tsv"
clock_data = "D:\\Datasets\\Horvath 2013 Data\\Horvath_Clock_CpGs.csv"

# -- Loading in altered data matrix & converting to GenomicRatioSet
if (startsWith(basename(datapath), "TCGA")) {
    raw_read   <- read_csv(datapath)
    raw_matrix <- as.matrix(raw_read[,2:ncol(raw_read)])
    rownames(raw_matrix) <- as_vector(raw_read[,1])
    setG <- makeGenomicRatioSetFromMatrix(as.matrix(raw_matrix))
} else if (startsWith(basename(datapath), "GSE")) {
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
if (startsWith(basename(datapath), "TCGA")) {
    pheno <- read.delim(pheno_data, stringsAsFactors = FALSE)
    pheno <- pheno[3:nrow(pheno),]
} else if (startsWith(basename(datapath), "GSE")) {
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

dat <- dat %>%
    filter(Chr == horvath_clock$Chr[2] & abs(Pos - horvath_Pos) <= 50000) %>%
    mutate(Rel_Pos = Pos - horvath_Pos)

# -- Getting age data for each sample
if (startsWith(basename(datapath), "TCGA")) {
    tmp <- pheno %>% 
        as_tibble() %>%
        select(bcr_patient_barcode, age_at_diagnosis) %>%
        mutate(age_at_diagnosis = as.numeric(age_at_diagnosis)) %>%
        setNames(c("Sample", "Age"))
} else if (startsWith(basename(datapath), "GSE")) {
    tmp <- pheno %>% 
        as_tibble() %>%
        mutate(Sample = rownames(pheno)) %>%
        select(Sample, Age = contains("age")) %>%
        mutate(Age = as.numeric(Age)) %>%
        setNames(c("Sample", "Age"))
}

# -- Standardizing GSE age data to years
if (startsWith(basename(datapath), "GSE")) {
    age_type = colnames(pheno[contains("age", var = colnames(pheno))])
    
    if (grepl("day", age_type)) {
        tmp$Age <- tmp$Age / 365
    } else if (grepl("week", age_type)) {
        tmp$Age <- tmp$Age * (7 / 365)
    } else if (grepl("month", age_type)) {
        tmp$Age <- tmp$Age / 12
    } else {}
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

ggplot(data = myd, aes(Rel_Pos, Cor)) +
    geom_point(alpha=0.30) +
    geom_point(aes(0, horvath_Cor), color = "red", size = 3) +
    geom_hline(yintercept = 0, color="red", lty=2) +
    ggtitle(paste("Correlation of CpGs Local to", horvath_clock$CpGmarker[2])) +
    xlab(paste("Position on Chromosome", horvath_clock$Chr[2])) +
    ylab("Pearson Correlation") +
    scale_y_continuous(limits = c(-1,1),
                       breaks = seq(-1,1,by=0.20)) +
    scale_x_continuous(limits = c(-50000, 50000),
                       breaks = seq(-50000, 50000, by = 10000)) +
    theme_minimal() +
    theme(plot.title=element_text( hjust=0.5, vjust=0.5, face='bold'))

