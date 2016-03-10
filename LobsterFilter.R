# Clean workspace -------
rm(list=ls())

# load the required libraries -------
library(reshape2)
library(ggplot2)
library(stringr)
library(stringi)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(randomForestSRC)
library(doParallel)
library(stackr)
library(purrr)
library(adegenet)


#Populations
Pops <- c("BON","BOO","BRA","BRO","BUZ","CAN","CAR","GAS","LOB",
          "MAG","MAL","MAR","OFF","RHO","SEA","SID","TRI")

#Populations ordered by latitutde (N -> S)
Pops_ordered <- c("TRI","BON","GAS","CAR","MAG","SID","MAL","BRA","CAN",
                  "SEA","LOB","BRO","OFF","BOO","MAR","BUZ","RHO")

## Missing genotypes - apply the 30% filter.
blacklisted.id <- missing_genotypes(haplotypes.file = "data/batch_1.haplotypes.tsv", 
                                    #whitelist.loci = "new.whitelist.txt", 
                                    pop.id.start = 1, pop.id.end = 3,
                                    pop.levels=Pops_ordered, 
                                    missing.geno.threshold = 30)

## impute the data
      genind.lobster <- haplo2genind(data="data/batch_1.haplotypes.tsv", 
                                      #whitelist.loci = "my.whitelist.txt", 
                                      blacklist.id = "blacklisted.id.30.txt", 
                                      pop.levels = Pops_ordered,
                                      pop.id.start = 1, 
                                      pop.id.end = 3, 
                                      imputation.method = "rf", 
                                      imputations.group = "populations",
                                      impute = "genotypes",
                                      num.tree = 100, 
                                      split.number = 100, 
                                      iteration.rf = 10, 
                                      verbose = FALSE)
      
      ##**## Error here 
      
      #Error in UseMethod("select_") : 
      #  no applicable method for 'select_' applied to an object of class "NULL"
      
      ## Create dataframes for imputed and non-imputed data. 
      LobsterData <- genind.lobster$no.imputation
      LobsterData.imputed <- geneind.lobster$imputation


#### Progress further using the vcf data but skipping the data imputation steps 
## tidy vcf
lobster.tidy.id <- read_stacks_vcf(vcf.file = "data/batch_1.vcf",
                                blacklist.id = "blacklisted.id.30.txt",
                                pop.id.start = 1, 
                                pop.id.end = 3, 
                                pop.levels = Pops_ordered,
                                filename= "data/lobster_tidy.vcf")

## Summary of tidy vcf
lobster.tidy.summary <- summary_stats_vcf_tidy(lobster.tidy.id, 
                                               filename= "data/lobster_tidy_summary.vcf")


# Covariance imbalance
cov.imbalance.figure <- plot_coverage_imbalance_diagnostic(tidy.vcf.file = lobster.tidy.id, 
                                                           pop.levels = Pops_ordered, 
                                                           read.depth.threshold = 16, 
                                                           aes.colour = aes(y = ..count..),
                                                           adjust.bin = 1) 
#adjust theme
cov.imbalance.figure <- cov.imbalance.figure + theme_bw()
                                                                                                                      adjust.bin = 1) 
#add facet
p1 <- cov.imbalance.figure + facet_grid(GROUP_GL ~ GROUP_COVERAGE)+theme_bw()

#save plots
ggsave("Figures/Lobster.coverage.imbalance.figure.pdf",
       cov.imbalance.figure+theme_bw(), width = 20, height = 20, dpi = 600, 
       units = "cm", useDingbats = F)
ggsave("Figures/Lobster.coverage.imbalance.figure_facet.pdf",
       p1, width = 20, height = 20, dpi = 600, units = "cm", useDingbats = F)

# You can also explore the imbalance on your own by using the table provided by this function:
low.coverage.summary <- table_low_coverage_summary(tidy.vcf.file = lobster.tidy, 
                                                   read.depth.threshold = 15, 
                                                   filename.low.coverage = "data/low.coverage.summary.table.tsv", 
                                                   filename.low.coverage.imbalance = "data/low.coverage.imbalance.summary.table.tsv")

##**## the One of the output tables is populated with NAs

## Filter data ---------------

## STEP 1 filter genotype likelihood. ----------

# this filter removes loci with an overall poor genotype likelihood (gl). Poor GL is different criteria, 
# the min and max depth threshold, the mean, min and diff (max-min) gl of a loci and the number or 
# percentage of population is required to remove the loci, based on the criteria.

lobster.tidy.id.gl <- filter_genotype_likelihood(lobster.tidy.id, 
                                                 approach="SNP",
                                             allele.min.depth.threshold = 1, 
                                             read.depth.max.threshold = 100, 
                                             gl.mean.threshold = 10, 
                                             gl.min.threshold = 5, 
                                             gl.diff.threshold = 100, 
                                             pop.threshold = 50, percent = T)

## STEP 2 filter loci ~ individual coverage ----------

# filter_individual: this filter will remove loci present in less than x% of the individuals.
lobster.tidy.id.gl.ind <- filter_individual(tidy.vcf=lobster.tidy.id.gl,
                                        pop.id.start = 1, 
                                        pop.id.end = 3, 
                                        pop.levels = Pops_ordered, 
                                        ind.threshold = 65, 
                                        percent = TRUE)

## STEP 3 filter loci ~ population ----------

# filter_population: this filter removes loci based on population representation
lobster.tidy.id.gl.ind.pop <- filter_population(data = lobster.tidy.id.gl.ind, 
                                                pop.threshold = 5,
                                                percent = F)

# STEP 4 update summary vcf object ----------
lobster.tidy.id.gl.ind.pop.sum <- summary_stats_vcf_tidy(data = lobster.tidy.id.gl.ind.pop)


## Step 5 filter Minor Allele frequency ----------

#investigate maf among populations
plot.maf <- plot_density_distribution_maf(data = lobster.tidy.id.gl.ind.pop.sum,
                                          maf.group = aes(x = FREQ_ALT, na.rm = F),
                                          aes.colour = aes(y = ..scaled.., color = POP_ID),
                                          adjust.bin = 0.2, #smooths the curve
                                          x.title = "Local MAF distribution")

plot.maf <- plot.maf + coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 1)) +
            facet_grid(~POP_ID)+theme_bw()+theme(legend.position="none");plot.maf

#Save plot
ggsave("Figures/Lobster_maf_plot.pdf",
       plot.maf, width = 30, height = 10, dpi = 600, units = "cm", useDingbats = F)


## Filter Minor Allele Frequency
lobster.tidy.id.gl.ind.pop.sum.maf <-    filter_maf(lobster.tidy.id.gl.ind.pop.sum,
                               approach="SNP",
                               local.maf.threshold = 2,
                               global.maf.threshold = 5,
                               pop.threshold = 1, 
                               filename = "data/lobster_tidy_id_gl_ind_pop_sum_maf.vcf")


## STEP 6 observed heterozygosity filter ----------

lobster.tidy.id.gl.ind.pop.sum.maf.het <- filter_het(data = lobster.tidy.id.gl.ind.pop.sum.maf,
                                                 approach = "SNP",
                                                 het.threshold = 0.5, 
                                                 het.diff.threshold = 0.5, #can turn it off by using 1
                                                 pop.threshold = 5, 
                                                 percent = F,
                                                 filename = "data/lobster.tidy.summary.maf.het.vcf")


## STEP 7 filter SNP number ----------
lobster.tidy.id.gl.ind.pop.sum.maf.het.snp <- filter_snp_number(data = lobster.tidy.id.gl.ind.pop.sum.maf.het, 
                                                                max.snp.number = 6, 
                                                                pop.threshold = 100,
                                                                percent = TRUE)
# , filename = "vcf.tidy.id.gl.ind.pop.sum.maf.het.fis.tsv")

# The inbreeding coef (Fis) ## do I need this...
#lobster.tidy.id.gl.ind.pop.sum.maf.het.fis <- filter_fis(data = lobster.tidy.id.gl.ind.pop.sum.maf.het, 
#                                                     approach = "haplotype",
#                                                     fis.min.threshold = -0.3, 
#                                                     fis.max.threshold = 0.3, 
#                                                     fis.diff.threshold = 0.5, 
#                                                     pop.threshold = 5, 
#                                                     percent = F)

## save the workspace
save.image("data/March10_Workspace.RData")

