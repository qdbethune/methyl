
# -- Generate plots based on data subsets centered at a single cpg
#    NOTE: DATAPATHS HERE ARE HEAVILY HARD-CODED

setwd("/rafalab/qbet/dna_met")

require("tidyverse")

# -- Horvath CpGs (Vector of identifiers and Clock coefficients)
clock_data = "./Horvath_Clock_CpGs.csv"
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

# TODO: 
data_dir <- list.dirs(path = "./450k_EDA")
full_dat <- read_csv(paste("./450k_EDA/", horvath_cpg, "_wind100k.csv", sep = ""))
    
# -- Computing spearman correlation between Age and Methylation
tmp <- full_dat %>%
    filter(!is.na(Met)) %>%
    filter(!is.na(Age)) %>%
    group_by(Probe, Rel_Pos) %>%
    mutate(Cor = cor(Met, Age, method = "spearman")) %>%
    ungroup()

# -- Only taking relative position on chromosome and correlation value
window_size = 100

myd <- tmp %>%
    select(Rel_Pos, Cor) %>%
    mutate(Rel_Pos = Rel_Pos / 1000) %>%
    filter(abs(Rel_Pos) <= window_size) %>%
    unique() 

# -- Plotting the results and saving the plot
horvath_Pos <- illumina_sheet$MAPINFO[illumina_sheet$IlmnID == horvath_cpg]
horvath_Cor = myd$Cor[myd$Rel_Pos == 0]

plot <- ggplot(data = myd, aes(Rel_Pos, Cor)) +
    stat_smooth(se = FALSE, span = 1) + 
    geom_point(alpha=0.4, size = 5) +
    geom_point(alpha=0.4, size = 5, pch=1) +
    geom_point(aes(0, horvath_Cor), color = "red", size = 8) +
    geom_point(aes(0, horvath_Cor), color = "black", size = 8, pch=1) +
    geom_hline(yintercept = 0, color="red", lty=2, size=0.80) +
    ggtitle(paste("Correlation of CpGs Local to", horvath_cpg)) +
    xlab(paste("Relative Position on Chromosome", horvath_clock$Chr[horvath_clock$CpGmarker == horvath_cpg],
               "(kb)")) +
    ylab("Spearman Correlation Coefficient") +
    scale_y_continuous(limits = c(-0.8,0.8),
                       breaks = seq(-0.8,0.8,by=0.20)) +
    scale_x_continuous(limits = c(-window_size, window_size),
                       breaks = seq(-window_size, window_size, by = window_size / 10)) +
    theme_minimal() +
        theme(plot.margin = unit(c(0.75, 1, 0.75, 0.75), "cm"),
        plot.title=element_text(hjust=0.5, vjust=2, size = 24),
        axis.title.x=element_text(vjust=-2, size = 16),
        axis.title.y=element_text(angle = 90, vjust=3, size = 16),
        axis.text = element_text(size = 12, color = "black"))

ggsave(paste(horvath_cpg, "_correlation_.png", sep = ""), 
       width = 16, height = 10, dpi = 320, type = "cairo", path = "./Outputs")






