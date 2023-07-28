library(ggplot2)
library(dplyr)
library(magrittr)

# Read in the files
prs <- read.table( "prs.best", header = T)
phen <- read.table( "phenotype.phen", header = T
                    
# Rename the phen
phen$Phenotype <- as.factor(phen$Phenotype)
levels(phen$Phenotype) <- c("healthy", "diseased")

# Merge the files
dat <- merge(prs, phen,by = "IID")
med_salary_df <- dat %>%
  group_by(Phenotype) %>%
  summarize(median=median(PRS))

# Start plotting
ggplot(dat, aes( x = PRS, color = Phenotype)) +
  geom_density() +
  theme_classic() +
  labs( x = "Polygenic Score", y="Density") +
  geom_vline( data = med_salary_df, aes( xintercept = median, color = Phenotype), size = 0.5)

# Save the plot
ggsave( "PRSice.png", height = 7, width = 7)
