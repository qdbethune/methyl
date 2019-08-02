
# -- Run bumphunter/blockfinder on grouped Illumina 450k datasets w.r.t. age data (use group_450k_data.R)
#    TODO: Use data in other ways. Currently, this script plots the regression coefficients of regions 

setwd("/rafalab/qbet/dna_met")

require(tidyverse)
require(minfi)
require(doParallel)

# -- bumphunter makes use of parallel processing if available
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

# -- Methylation data matrix file (normalized) and age data matrix
data_filepath = "./Collated Data/Unnormalized_450k.csv"
age_filepath = "./Collated Data/Age_450k.csv"

# -- Local function to generate plots of regional regression coefficients
dmr_plot <- function(i, type) {    
    
    # -- Info to be used for filtering and plotting
    dmr_chr = dmr_table$chr[i]
    dmr_start = dmr_table$start[i]
    dmr_end = dmr_table$end[i]
    dmr_number = as.integer(rownames(dmr_table)[i])
    
    # -- Filter collapsed cpgs to only include those in a given "block"
    curr_gr <- gr %>%
        filter(seqnames == dmr_chr, start >= dmr_start, start <= dmr_end)
    
    # Generate and save plot
    plot <- ggplot(data = curr_gr, aes(start, coef)) +
        stat_smooth(se = FALSE, method = "lm") + 
        geom_point(alpha=0.4, size = 5) +
        geom_point(alpha=0.4, size = 5, pch=1) +
        geom_hline(yintercept = 0, color="red", lty=2, size=0.80) +
        ggtitle(paste("Regression coefficients of collapsed cpgs on", dmr_chr)) +
        xlab(paste("Position on", dmr_chr)) +
        ylab("Bumphunter Regression Coefficient") +
        scale_y_continuous(limits = c(-0.01,0.01)) +
        scale_x_continuous(limits = c(dmr_start, dmr_end)) +
        
        theme_minimal() +
        theme(plot.margin = unit(c(0.75, 1, 0.75, 0.75), "cm"),
              plot.title=element_text(hjust=0.5, vjust=2, size = 24),
              axis.title.x=element_text(vjust=-2, size = 16),
              axis.title.y=element_text(angle = 90, vjust=3, size = 16),
              axis.text = element_text(size = 12, color = "black"))
    
    ggsave(paste(type, "_region_", dmr_number, ".png", sep = ""), 
           width = 16, height = 10, dpi = 320, type = "cairo", path = "./Outputs")
}


# -- Loading in and formatting collated methylation data processed by pipeline
data_tibble <- read_csv(data_filepath)
genome_tibble <- select(data_tibble, Probe, Pos, Chr, Probe_type)
beta_tibble <- select(data_tibble, -Probe, -Pos, -Chr, -Probe_type)
NA_vector <- apply(beta_tibble, MARGIN = 1, FUN = function(x) sum(is.na(x)))
print(paste("The data matrix has a total of", sum(NA_vector), "missing values."))
sampleList <- colnames(beta_tibble)


# -- Loading in collated age data and arranging to match order of samples in beta value matrix
age_data = read_csv(age_filepath)
age_data <- age_data %>%
    filter(Sample %in% sampleList) %>%
    arrange(match(Sample, sampleList))

beta_tibble <- beta_tibble %>%
    select(age_data$Sample)

# -- Constructing minfi GenomicRatioSet object using methylation matrix
beta_matrix <- as.matrix(beta_tibble)
rownames(beta_matrix) <- data_tibble$Probe
setG <- makeGenomicRatioSetFromMatrix(as.matrix(beta_matrix))
rm(data_tibble, beta_tibble, beta_matrix, NA_vector)


# -- Creating age design matrix for bump hunting.
ages <- age_data$Age
designMatrix <- model.matrix(~ ages)
rownames(designMatrix) <- age_data$Sample


# -- Blockfinder uses the bumphunter engine to analyze much wider methylation regions.
#    NOTE: Original cpgCollapse function was bugged for objects with no CN data (like GenomicRatioSets)
source("../scripts/cpgcollapse_fix.R")
collapse <- new_cpgCollapse(setG, what = "Beta", maxGap = 500, maxClusterWidth = 5000)
age_dmrs <- blockFinder(collapse$object, design = designMatrix, cutoff = 0.001, B = 1000, what = "Beta")

# -- Filtering data by FWER, currently using a static cutoff
FWER_CUTOFF = 0.05
dmr_table <- age_dmrs$table %>%
    filter(fwer <= FWER_CUTOFF)

print(paste("There are", nrow(dmr_table), "blockfinder bumps with FWER <=", FWER_CUTOFF))

dmr_table$range <- dmr_table$end - dmr_table$start
collapsed_data <- collapse$object

# -- Mapping coefficients produced by bumphunter to collapsed cpgs
gr <- collapsed_data %>%
    granges() %>%
    as_tibble %>%
    bind_cols(coef = age_dmrs$coef) %>%
    mutate(range = end - start + 1)
    

# -- Generate plots utilizing parallel processing
foreach(i = 1:nrow(dmr_table), .packages = "tidyverse") %dopar% dmr_plot(i, type = "blockfinder")




