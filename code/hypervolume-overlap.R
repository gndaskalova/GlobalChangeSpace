library(factoextra)
library(tidyverse)
#library(ggbiplot)
library(hypervolume)

# There are five global change spaces (and thus 5 PCAs, the files for each are attached too). 
# I think what would work well is for each global change space combination 
# (e.g., global plus marine LPD, global plus terr LPD) 
# to progressively sample bigger and bigger random samples and calculate the resulting 
# overlap until the actual sample size is reached. 
# And it could be done in percentages (e.g., first 5% of the data, then 10%, etc.
# up until 100% of the data) so that then I can plot “global change accumulation curves” for 
# all the databases together (so that the scale is the same, as otherwise the BioTIME marine sample size is much larger than all the others). 
# Or I could plot the curves separately in which case the sample size can stay in their normal units 
# and not be transformed to percentages.

#set up all possible tasks

#all files
allFiles <- c("combined_bio_marine_simple.RData","combined_bio_terr_simple.RData",
              "combined_pop_marine_simple.RData","combined_pop_terr_simple.RData",
              "combined_predicts_terr_simple.RData")

#all proportions
allProps <- c(0.01,seq(0.05,1,by=0.05))#start with 0.01 instead of 0 proportion

#choose combination for this task:

currentTask <- expand_grid(allFiles,allProps) %>%
                add_column(n = 1:nrow(.)) %>% #105 tasks
                filter(n == as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))

currentFile <- currentTask$allFiles
currentProp <- currentTask$allProps

#load file for this task
myFolder <- "/data/idiv_ess/Gergana"
assign('combined_pop', get(load(paste(myFolder,currentFile, sep="/"))))

#do pca
pop_pca <- prcomp(combined_pop[, -6], scale = TRUE, center = TRUE)
df <- as.data.frame(get_pca_ind(pop_pca)$coord[,1:2])
df$sampling2 <- combined_pop$sampling2

#get overlap - write as a function

overlapFun <- function(df, currentProp){

#get random sample
group1 <- df %>%
            filter(grepl("Random",df$sampling2)) %>%
            select(Dim.1, Dim.2) %>%
            hypervolume_box()

#get data sample
group2 <- df %>%
            filter(!grepl("Random",df$sampling2)) %>%
            select(Dim.1, Dim.2) %>%
            slice_sample(prop = currentProp) %>% 
            hypervolume_box()

#get overlap
mySet <- hypervolume_set(group1,group2, check.memory = FALSE)
hypervolume_overlap_statistics(mySet)

}

#run function 100 times for each random sample proportion
output <- replicate(100,
                    overlapFun(df, currentProp)) %>%
                    t() %>%
                    as_tibble() %>%
                    add_column(run = 1:nrow(.))

#save output
saveRDS(output, file=paste('hypervolumn_overlap', 
                           gsub('.RData','',currentFile), #jut get rid of file extension
                           currentProp,
                           ".rds", 
                           sep="_"))
