setwd("/rafalab/qbet/dna_met")

require("tidyverse")

# -- List of 450k array datasets
horvath_meta <- read_csv("./healthy_tissue_datasets.csv")
horvath_450k_sets <- horvath_meta$Availability[horvath_meta$Platform == "450K" & horvath_meta$`Data Use` != "Other" & horvath_meta$Availability != "GSE42865"]
clock_data = "./Horvath_Clock_CpGs.csv"

# -- Horvath CpGs (Vector of identifiers and Clock coefficients)
horvath_clock <- read_csv(clock_data)
horvath_ids <- horvath_clock$CpGmarker[2:length(horvath_clock$CpGmarker)]
horvath_model <- horvath_clock$CoefficientTraining[2:nrow(horvath_clock)]

# -- Illumina 450k Reference Datasheet (for referencing genome positions)
illumina_sheet <- read_csv("./Platform Information/GPL13534_HumanMethylation450_15017482_v.1.1.csv")

data_list <- list.files(path = "./450k_EDA")

for (file in data_list) {
    dat <- read_csv(file = paste("./450k_EDA/", file, sep = ""))
    cpg <- regmatches(file, regexpr("cg\\d{8}", file))
    
    dat <- dat %>%
        filter(!is.na(Met)) %>%
        filter(!is.na(Age)) %>%
        group_by(Probe) %>%
        filter(Probe == cpg) %>%
        select(Probe, Cor, Met, Age) %>%
        ungroup()
    
    if (exists("full_dat")) {
        
        full_dat <- full_dat %>%
            full_join(dat)
        
    } else {
        full_dat <- dat
    }
}

full_dat <- full_dat %>% 
    arrange(Age)

plot <- ggplot(data = full_dat, aes(x = as.factor(Age), y = reorder(Probe, Cor), fill = Met)) +
    geom_tile(color="black") +
    ggtitle("Methylation of Samples by Age") +
    ylab("Illumina Probe") +
    xlab("Age (years)") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.40) +

    theme_minimal() +
    theme(plot.margin = unit(c(0.75, 1, 0.75, 0.75), "cm"),
          plot.title=element_text(hjust=0.5, vjust=2, size = 24),
          axis.title.x=element_text(vjust=-2, size = 16),
          axis.title.y=element_text(vjust = 1, size = 16),
          axis.text.x=element_blank(),
          axis.text.y = element_blank())

ggsave("Methylation_Spearman_Heatmap.png", width = 16, height = 10, dpi = 150, type = "cairo", path = "./Outputs")
