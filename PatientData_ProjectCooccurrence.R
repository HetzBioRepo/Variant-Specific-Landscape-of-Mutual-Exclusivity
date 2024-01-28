library(dplyr)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(ggtext)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(foreach)
library(doParallel)
library(svglite)
library(officer)
library(rsvg)

t1 <- Sys.time()

### All genes to analyze
genesToDo <- c("KRAS", "BRAF", "EGFR")

################################################################################
# Preprocessing of patient data
################################################################################
# Patient data contains duplicates, three redundant columns and addition of term
# ' (driver)' to certain mutations. Solve all this!

### Load patient data
# All files has to be stored in one folder to be automatically read.
rawdatafolder <- choose.dir(default = getwd(), caption = "Select raw data folder")
inputData <- read.table(file = paste0(rawdatafolder, "\\2022_12_alterations_across_samples.tsv"),
                        header = TRUE,
                        sep = "\t")
# create result folders
resultfolder <- paste0(dirname(rawdatafolder),
                       "/results_samples")
if(!dir.exists(resultfolder)) dir.create(resultfolder)
fileresultfolder <- paste0(resultfolder,
                       "/files")
if(!dir.exists(fileresultfolder)) dir.create(fileresultfolder)

### Remove redundant columns
pp_inputData <- inputData
for(i in genesToDo){
  pp_inputData <- pp_inputData[, -which(colnames(pp_inputData) == paste0(i, "..MUT"))]
}

### Remove 'Altered' column, since it counts the unspecific term 'MUTATED' as altered!
pp_inputData <- subset(pp_inputData, select = -c(Altered))

### remove duplicates
pp_inputData <- pp_inputData[!duplicated(pp_inputData$Sample.ID),]
pp_inputData <- pp_inputData[!duplicated(pp_inputData$Patient.ID),]

### Replace synomynous mutations (e.g. L210=) with "no alteration"
for(i in genesToDo){
  pp_inputData[,i] <- gsub("[A-Z][0-9]{1,}=", "no alteration", pp_inputData[,i])
  # If synonymous mutation co-occurred with other mutations, remove 'no alteration'
  pp_inputData[,i] <- gsub("no alteration, ", "", pp_inputData[,i])
  pp_inputData[,i] <- gsub(", no alteration", "", pp_inputData[,i])
}

### Remove ' (driver)' term
for(i in genesToDo){
  pp_inputData[,i] <- gsub(" (driver)", "", pp_inputData[,i], fixed = TRUE)
}

### Which patient data is valid for each gene? (Unspecific mutational information 'MUTATED'
# discards sample as valid for this mutation!)
validList <- list()
for(i in genesToDo){
  validList[[i]] <- (pp_inputData[, i] != "MUTATED" & pp_inputData[, i] != "not profiled")
}

### Save Alteration status of genes
alterationList <- list()
for(i in genesToDo){
  alterationList[[i]] <- (pp_inputData[, i] != "no alteration" & validList[[i]] == TRUE)
}


################################################################################
# Preprocessing of EGFR mutations
################################################################################
# For classificaton of EGFR mutations, certain classes are defined by an Exon 19
# deletion (Ex19del). However, for future single mutation analysis, this Ex19del
# should be replaced by all mutations occuring in the dataset defined as Ex19del.

### Load KRAS/BRAF/EGFR class definitions
# All files has to be stored in one folder to be automatically read.
classes <- list()
for(i in genesToDo){
  classes[[i]] <- read.table(file = paste0(rawdatafolder, "\\", i, "_Mutation_Classes.txt"),
                             header = TRUE,
                             sep = "\t")
}
classesEGFR_temp <- classes$EGFR

### Regular expression for all Exon 19 deletions (amino acid 729-761)
regexEGFR19del <- "[A-Z](729|7[3-5][0-9]|76[0-1])(_[A-Z]){0,1}(729|7[3-5][0-9]|76[0-1]){0,1}del"

### Identify all single deletions in Exon 19
# First get all single mutations in EGFR column
single_mutations_egfr <- unique(unlist(strsplit(x = pp_inputData$EGFR, split = ", ")))
# Then grab all exon 19 deletions
sdel_ex19 <- grep(pattern = regexEGFR19del,
                  x = single_mutations_egfr,
                  fixed = FALSE,
                  value = TRUE)
sdel_ex19 <- sort(sdel_ex19)

### Replace every Ex19del in EGFR classification by single mutation
# Where to find term 'Ex19del'?
ids_to_replace <- grep(pattern = "Ex19del",
                       x = classesEGFR_temp$Mutation,
                       fixed = TRUE)
# save non-interesting rows as classesEGFR
classes$EGFR <- classesEGFR_temp[-ids_to_replace,]
# get rows containing 'Ex19del'
subTemp <- classesEGFR_temp[ids_to_replace,]
# replace every 'Ex19del' by single mutations in dataset
for(i in seq(nrow(subTemp))){
  actReplaced <- sapply(X = sdel_ex19,
                        FUN = function(x) return(gsub("Ex19del",
                                                      x,
                                                      subTemp$Mutation[i],
                                                      fixed = TRUE)))
  actReplaced <- actReplaced[!(actReplaced %in% classes$EGFR$Mutation)]
  
  if(length(actReplaced) > 0){
    temp <- data.frame(Mutation = actReplaced,
                       Class = rep(subTemp$Class[i], length(actReplaced)))
    classes$EGFR <- rbind(classes$EGFR, temp)
  }
}
# order new data frame alphabetically
classes$EGFR <- classes$EGFR[order(classes$EGFR$Class,
                                   classes$EGFR$Mutation,
                                 decreasing = FALSE),]
# remove duplicated 'E746_A750del, A647T' in classification...
classes$EGFR <- classes$EGFR[!duplicated(classes$EGFR$Mutation),]

# save EGFR classes
write.table(x = classes$EGFR,
            file = paste0(rawdatafolder, "\\EGFR_Mutation_Classes_v2_complete.txt"),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

################################################################################
# Create dependency lists
################################################################################
# Especially for EGFR mutations, where class assignment is in some cases defined
# by cooccurrence of multiple mutations, certain definitions infer with each other.
# E.g. T790M (T790M_like_3S) and G724S, T790M (T790M_like_3R). In this case,
# if a more specific definition (defined as relying on cooccurrence of more mutations)
# is available, this definition is assigned instead of the more general definition.
# In this example, G724S, T790M is prioritized against T790M alone.
# This relationships between mutations are now to be found.
depList <- list()
for(i in genesToDo){
  depList[[i]] <- list()
}

### Helper function to identify presence of mutations
# Are all mutations of string a present in string b?
compareMutations <- function(a, b, psplit = ", "){

  # decomposite a
  decA <- unlist(strsplit(a, split = psplit))
  
  # decomposite b
  decB <- unlist(strsplit(b, split = psplit))
  
  boolEqual <- all(decA %in% decB)
  
  return(boolEqual)
}

### Now create dependency lists
# Iterate through every gene
for(actGene in genesToDo){
  for(i in seq(nrow(classes[[actGene]]))){
    # Iterate through all remaining mutations
    for(j in seq(nrow(classes[[actGene]]))[-i]){
      # If all mutations are matched, add dependency
      if(compareMutations(classes[[actGene]]$Mutation[i],
                          classes[[actGene]]$Mutation[j]) == TRUE){
        depList[[actGene]][[classes[[actGene]]$Mutation[i]]] <- c(depList[[actGene]][[classes[[actGene]]$Mutation[i]]],
                                                                  classes[[actGene]]$Mutation[j])
      }
    }
  }
}


################################################################################
# Class assignment
################################################################################

### Class assignment
# Iterate through every gene
for(actGene in genesToDo){
  
  print(actGene)
  
  # Vector to store class assignment
  tempVec <- rep("", nrow(pp_inputData))
  
  # Iterate through every mutated patient
  for(i in which(alterationList[[actGene]] == TRUE)){
    
    # For progress
    #print(i)
    
    actClass <- c()
    
    # Iterate through every class definition for actual gene
    for(j in seq(nrow(classes[[actGene]]))){
      # Is current class definition matched in current patient?
      if(compareMutations(classes[[actGene]]$Mutation[j],
                          pp_inputData[i, actGene]) == TRUE){
        addClass <- TRUE
        
        # Check if class assignment is not prohibited by a prioritized class definition
        if(classes[[actGene]]$Mutation[j] %in% names(depList[[actGene]])){
          # If prioritized class assignment is present, reject current class assignment
          for(k in depList[[actGene]][[classes[[actGene]]$Mutation[j]]]){
            if(compareMutations(k, pp_inputData[i, actGene]) == TRUE) addClass <- FALSE
          }
        }
        
        # If addClass, add class to sample
        if(addClass){
          actClass <- c(actClass, paste0(classes[[actGene]]$Class[j],
                                         ": ",
                                         classes[[actGene]]$Mutation[j]))
        }
      }
    }
    
    actClass <- sort(actClass)
    tempVec[i] <- paste(actClass, collapse = "; ")
  }
  
  pp_inputData[, paste0("classes", actGene)] <- tempVec
}

### Check how often in EGFR mutations the same mutation contributes to two different class assignments
# Priorization via dependency list prohibits e.g. assignment via T790M (T790M_like_3S)
# if G724S, T790M (T790M_like_3R) as combination is present. However, how often
# do cases occur with a shared mutation but in an equal priority? (e.g. E709K, L858R (Classical_like)
# and L858R, T790M (T790M_like_3S) in the same patient).
print(pp_inputData[sapply(X = pp_inputData$classesEGFR,
                          FUN = function(x){
                            if(x == "") return(FALSE)
                            
                            actMuts <- unlist(strsplit(x = x, split = "; ", fixed = TRUE))
                            actMuts <- gsub("^.*: ", "", actMuts, fixed = FALSE)
                            actMuts <- unlist(strsplit(x = actMuts, split = ", ", fixed = TRUE))
                            return(any(table(actMuts) > 1))
                          }), c("Study.ID", "Sample.ID", "Patient.ID", "EGFR", "classesEGFR")])

# Result: four cases in which the same single mutation contributes to two class assignments
# Study.ID         Sample.ID Patient.ID                       EGFR                                               classesEGFR
# 5331 luad_oncosg_2020              A349       A349        S768I, G719C, V769L                    PACC: G719C, S768I; PACC: S768I, V769L
# 6140  msk_impact_2017 P-0001328-T01-IM3  P-0001328   T790M, L858R, E709K, G5A Classical_like: E709K, L858R; T790M_like_3S: L858R, T790M
# 6280  msk_impact_2017 P-0002085-T02-IM5  P-0002085        T790M, L858R, T854S           PACC: L858R, T854S; T790M_like_3S: L858R, T790M
# 6700  msk_impact_2017 P-0003923-T01-IM5  P-0003923 S768I, T790M, G719C, C775Y           PACC: G719C, S768I; T790M_like_3S: S768I, T790M


################################################################################
# Class statistics
################################################################################
# Get counts and cooccurrence of classes and genes

### Create class counts list
classesCounts <- list()
for(i in genesToDo){
  temp <- unique(classes[[i]]$Class)
  classesCounts[[i]] <- data.frame(Index = seq(temp),
                                   Gene = rep(i, length(temp)),
                                   Class = temp)
  classesCounts[[i]] <- rbind(classesCounts[[i]],
                              data.frame(Index = nrow(classesCounts[[i]])+1,
                                         Gene = i,
                                         Class = "Unclassified"))
}

### Create class matrices
# Rows represent each patient, Column represent a Class of one gene.
# 0 indicates no mutation of according class in this patient, 1 indicates the occurrence of this class.
classesMatrices <- list()
# Iterate over every gene
for(i in genesToDo){
  classesMatrices[[i]] <- matrix(rep(0, nrow(classesCounts[[i]])*nrow(pp_inputData)),
                                 ncol = nrow(classesCounts[[i]]))
  # Iterate over every class
  for(j in 1:(ncol(classesMatrices[[i]])-1)){
    classesMatrices[[i]][,j] <- sapply(pp_inputData[,paste0("classes", i)],
                                       FUN = function(x){
                                         if(x == "") return(0)
                                         
                                         actClasses <- unlist(strsplit(x = x, split = "; ", fixed = TRUE))
                                         actClasses <- gsub(":.*$", "", actClasses, fixed = FALSE)
                                         actClasses <- unique(actClasses)
                                         
                                         # Is current class present in current patient?
                                         if(classesCounts[[i]]$Class[j] %in% actClasses) return(1) else return(0)
                                       })
    classesMatrices[[i]][,ncol(classesMatrices[[i]])] <- as.integer(alterationList[[i]] == TRUE &
                                                                      pp_inputData[,paste0("classes", i)] == "")
  }
}

# Save class matrices
for(i in genesToDo){
  write.table(x = classesMatrices[[i]],
              file = paste0(fileresultfolder, "\\Class_Matrix_", i, ".txt"),
              row.names = TRUE,
              col.names = TRUE,
              sep = "\t",
              quote = FALSE)
}

### Get counts and frequency of every class
for(i in genesToDo){
  classesCounts[[i]]$Counts <- sapply(seq(nrow(classesCounts[[i]])),
                                      FUN = function(x){
                                        return(sum(classesMatrices[[i]][,x]))
                                      })
  classesCounts[[i]]$Freq <- classesCounts[[i]]$Counts/sum(classesCounts[[i]]$Counts)
}

### Draw pie chart for all classes of each gene
# First a summarizing dataframe with all class counts is needed
summaryClassCounts <- data.frame(Gene = c(), Class = c(), Counts = c(), Freq = c())
for(i in genesToDo){
  summaryClassCounts <- rbind(summaryClassCounts,
                              classesCounts[[i]][, c("Gene", "Class", "Counts", "Freq")])
}

write.table(x = summaryClassCounts,
            file = paste0(fileresultfolder, "\\Class counts.txt"),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

# Add labels and position in pie chart
summaryClassCounts <- summaryClassCounts %>% 
                      group_by(Gene)%>%
                      # arrange(perc) %>%
                      mutate(Labels = paste0(Class, "\n", Counts, "/", sum(`Counts`), " (", scales::percent(Freq), ")")) %>%
                      mutate(cumulativeSum = rev(cumsum(rev(Counts))), 
                      Position = Counts/2 + lead(cumulativeSum, 1),
                      Position = if_else(is.na(Position), Counts/2, Position))
# Load Color code
colorCode <- read.table(file = paste0(rawdatafolder, "\\Color_Code.txt"),
                        header = TRUE,
                        sep = "\t",
                        comment.char = "@")
# Draw Pie function
drawPie <- function(pGene, nudgeX = 1, nudgeY = 0, title = "", angle = 0, direction = 1){
  
  # Color assignment and preprocessing
  tDataframe <- summaryClassCounts %>% filter(Gene == pGene)
  unqValuesOfClass <- sort(unique(tDataframe$Class))
  unqColorsOfClass <- sapply(unqValuesOfClass, FUN = function(x){
    return(colorCode$Color[colorCode$Gene == pGene & colorCode$Class == x])
  })
  tDataframe$Labels <- gsub("_", "-", tDataframe$Labels, fixed = TRUE)
  
  tDataframe$Position <- tDataframe$Position/max(tDataframe$Counts)
  tDataframe$Counts <- tDataframe$Counts/max(tDataframe$Counts)
  
  
  # Generate ggplot object
  gg <- ggplot(tDataframe, aes(x = "" , y = Counts, fill = fct_inorder(Class))) +
    ggtitle(title) +
    geom_col(color = "black") +
    coord_polar(theta = "y", start = angle, direction = direction) +
    geom_label_repel(aes(y = Position, label = Labels, segment.size = 1),
                     fontface = "bold",
                     size = 3,
                     nudge_x = nudgeX,
                     nudge_y = nudgeY, show.legend = FALSE) +
    scale_fill_manual(breaks = unqValuesOfClass, 
                      values= unqColorsOfClass) +
    # geom_text_repel(aes(y = pos, label = labels),
    #                  size = 4, nudge_x = nudgeX,
    #                  nudge_y = nudgeY, show.legend = FALSE) +
    theme_void() +
    theme(legend.position="none",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  
  return(gg)
}
# Generate result folder
pieresultfolder <- paste0(resultfolder,
                          "/pie")
if(!dir.exists(pieresultfolder)) dir.create(pieresultfolder)
# Draw it!
svg(filename = paste0(pieresultfolder, "/KRAS_classes.svg"), width = 8/2.54, height = 8/2.54)
drawPie(pGene = "KRAS", nudgeX = 0.75, title = "KRAS", direction = 1, angle = 0)
dev.off()

svg(filename = paste0(pieresultfolder, "/BRAF_classes.svg"), width = 8/2.54, height = 8/2.54)
drawPie(pGene = "BRAF", nudgeX = 0.75, title = "BRAF", direction = 1, angle = 0)
dev.off()

svg(filename = paste0(pieresultfolder, "/EGFR_classes.svg"), width = 8/2.54, height = 8/2.54)
drawPie(pGene = "EGFR", nudgeX = 0.75, title = "EGFR", direction = 1, angle = pi*3/4)
dev.off()

# Save in PowerPoint
my_doc <- read_pptx()

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(pieresultfolder, "/KRAS_classes.svg"), width = 8, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(pieresultfolder, "/BRAF_classes.svg"), width = 8, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(pieresultfolder, "/EGFR_classes.svg"), width = 8, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)


print(my_doc, target = paste0(pieresultfolder, "/piecharts.pptx"))

### Now calculate cooccurrence and mutual exclusivity statistics between classes and genes
# Function to calculate OddsRatio and standard error
calcORandSE <- function(pDataframe){
  retdf <- pDataframe
  
  # To calculate Odds Ratio, account for prohibited 0 values by adding 0.5 if a 0 occurs (Haldane correction).
  AaB <- retdf$AandB
  nAaB <- retdf$notAandB
  AanB <- retdf$AandnotB
  nAanB <- retdf$notAandnotB
  
  # v1.2: To avoid bias, add 0.5 to all data!
  # if(AaB == 0 | nAaB == 0 | AanB == 0 | nAanB == 0){
  # }
  
  AaB <- AaB+0.5
  nAaB <- nAaB+0.5
  AanB <- AanB+0.5
  nAanB <- nAanB+0.5
  
  # Calculate Odds Ratio
  retdf$oddsRatio <- (AaB*nAanB)/(nAaB*AanB)
  retdf$standardError <- sqrt(1/AaB + 1/nAaB + 1/AanB + 1/nAanB)
  
  return(retdf)
}
# Function to calculate statistics between two boolean vectors representing
# occurrence of event A and B, respectively.
calcStats <- function(boolEventA, boolEventB, FUN_additional_Stats = calcORandSE){
  
  # Calculate four-field table
  temp <- data.frame(CountsA = sum(boolEventA),
                     CountsB = sum(boolEventB),
                     AandB = sum((boolEventA == TRUE) & (boolEventB == TRUE)),
                     notAandB = sum((boolEventA == FALSE) & (boolEventB == TRUE)),
                     AandnotB = sum((boolEventA == TRUE) & (boolEventB == FALSE)),
                     notAandnotB = sum((boolEventA == FALSE) & (boolEventB == FALSE)))
  # Calculate Fisher's Exact Test
  temp$fisherTwoSided <- fisher.test(matrix(as.numeric(temp[,3:6]), ncol = 2, byrow = TRUE),
                                     conf.level = 0.95,
                                     alternative = "two.sided")$p.value
  # For following purposes
  temp$significantTwoSided <- ""
  
  # Calculate additional stats
  # Default: Odds ratio and standard error
  if(!is.null(FUN_additional_Stats)){
    temp <- FUN_additional_Stats(temp)
  }
  
  return(temp)
}
# Get Dataframe with all combinations of classes and genes
temp <- summaryClassCounts[summaryClassCounts$Class != "Unclassified", c("Gene", "Class")]
coocStats <- NULL
for(i in genesToDo){
  tempdf <- temp[temp$Gene != i,]
  colnames(tempdf) <- c("GeneA", "ClassA")
  tempdf$GeneB <- rep(i, nrow(tempdf))
  
  coocStats <- rbind(coocStats,
                     tempdf)
}
# Statsv1: Now calculate statistics of cooccurrence of ClassA and GeneB in all patients with mutated GeneA
tempdf <- NULL
for(i in seq(nrow(coocStats))){
  actA <- coocStats$GeneA[i]
  actB <- coocStats$GeneB[i]
  classIndex <- classesCounts[[actA]]$Index[classesCounts[[actA]]$Class == coocStats$ClassA[i]]
  
  # all patients with mutation of classA if GeneA is mutated
  boolEventA <- classesMatrices[[actA]][, classIndex][alterationList[[actA]] == TRUE & validList[[actB]] == TRUE]
  
  # all patients with GeneB mutated AND GeneA mutated
  boolEventB <- alterationList[[actB]][alterationList[[actA]] == TRUE & validList[[actB]] == TRUE]
  
  tempdf <- rbind(tempdf,
                  calcStats(boolEventA = boolEventA,
                            boolEventB = boolEventB))
}

# Bonferroni adjustment for each analysis seperately
# Significance declaration
actSignComb <- unique(coocStats[, c("GeneA", "GeneB")])
for(i in 1:nrow(actSignComb)){
  ids <- which(coocStats$GeneA == actSignComb$GeneA[i] & coocStats$GeneB == actSignComb$GeneB[i])
  tempdf$significantTwoSided <- tempdf$fisherTwoSided < 0.05/length(ids)
}

coocStatsv1 <- cbind(coocStats, tempdf)

# Save statistics
write.table(x = coocStatsv1,
            file = paste0(fileresultfolder, "\\Cooccurrence Statistics Subset GeneA mutated.txt"),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

################################################################################
# Single mutation or tuple counts
################################################################################
# Get counts of all single mutations ot tuples of mutations (in case of e.g. EGFR)

### Create mutation counts list
mutationCounts <- list()
for(i in genesToDo){
  
  # Remove all single mutations of tuples used for classification and add tuple information
  tempMutations <- sapply(seq(nrow(pp_inputData))[alterationList[[i]] == TRUE],
                          FUN = function(x){
                                actMutations <- pp_inputData[x,i]
                                actTuples <- pp_inputData[x,paste0("classes", i)]
                                if(actMutations == "") return("")
                                
                                # Get all single mutations of patient x
                                actMutations <- unique(unlist(strsplit(x = actMutations, split = ", ")))
                                # Get all single mutations/tuples used for classification of patient x
                                actTuples <- unique(unlist(strsplit(x = actTuples, split = "; ")))
                                actTuples <- gsub("^.*: ", "", actTuples, fixed = FALSE)
                                # Get all single mutations of tuples used for classification of patient x
                                actMutTuples <- unique(unlist(strsplit(x = actTuples, split = ", ")))
                                
                                # Remove single mutations of tuples from all single mutations of patient x
                                actMutations <- actMutations[!(actMutations %in% actMutTuples)]
                                
                                # Combine non-classified mutations and tuples of mutations used for classification
                                return(paste(sort(c(actMutations, actTuples)), collapse = "; "))
                              })
  pp_inputData[,paste0("tuples", i)] <- rep("", nrow(pp_inputData))
  pp_inputData[alterationList[[i]] == TRUE, paste0("tuples", i)] <- tempMutations
  
  # Get all single mutations and tuples
  temp <- sort(unique(unlist(strsplit(pp_inputData[,paste0("tuples", i)], split = "; "))))
  mutationCounts[[i]] <- data.frame(Index = seq(temp),
                                    Gene = rep(i, length(temp)),
                                    Mutation = temp)
}

### Create mutation matrices
# Rows represent each patient, Column represent a mutation/tuple of one gene.
# 0 indicates no occurence of this mutation in this patient, 1 indicates the occurrence of this mutation.
mutationMatrices <- list()
# Iterate over every gene
for(i in genesToDo){
  mutationMatrices[[i]] <- matrix(rep(0, nrow(mutationCounts[[i]])*nrow(pp_inputData)),
                                 ncol = nrow(mutationCounts[[i]]))
  # Iterate over every class
  for(j in 1:ncol(mutationMatrices[[i]])){
    mutationMatrices[[i]][,j] <- sapply(pp_inputData[,paste0("tuples", i)],
                                       FUN = function(x){
                                         if(x == "") return(0)
                                         
                                         actTuples <- unlist(strsplit(x = x, split = "; ", fixed = TRUE))
                                         actTuples <- unique(actTuples)
                                         
                                         # Is current mutation/tuple present in current patient?
                                         if(mutationCounts[[i]]$Mutation[j] %in% actTuples) return(1) else return(0)
                                       })
  }
}

# Save mutation matrices
for(i in genesToDo){
  write.table(x = mutationMatrices[[i]],
              file = paste0(fileresultfolder, "\\Mutation_Matrix_", i, ".txt"),
              row.names = TRUE,
              col.names = TRUE,
              sep = "\t",
              quote = FALSE)
}


### Get counts and frequency of every mutation
for(i in genesToDo){
  mutationCounts[[i]]$Class <- sapply(mutationCounts[[i]]$Mutation, FUN = function(x){
    for(j in seq(nrow(classes[[i]]))){
      if(x == classes[[i]]$Mutation[j]) return(classes[[i]]$Class[j])
    }
    
    return("Unclassified")
  })
  mutationCounts[[i]]$Counts <- sapply(seq(nrow(mutationCounts[[i]])),
                                      FUN = function(x){
                                        return(sum(mutationMatrices[[i]][,x]))
                                      })
  mutationCounts[[i]]$Freq <- mutationCounts[[i]]$Counts/sum(mutationCounts[[i]]$Counts)
}

### Ranked plot: Counts of single mutations
# Function for ranked plot
printRankedPlot <- function(pGene, title = "", jitter = 0){
  
  # Color assignment and preprocessing
  tDataframe <- mutationCounts[[pGene]]
  tDataframe <- tDataframe[order(tDataframe$Counts, decreasing = TRUE),]
  tDataframe$rank <- seq(nrow(tDataframe))
  unqValuesOfClass <- sort(unique(tDataframe$Class))
  unqColorsOfClass <- sapply(unqValuesOfClass, FUN = function(x){
    return(colorCode$Color[colorCode$Gene == pGene & colorCode$Class == x])
  })
  
  tDataframe$Class <- gsub("_", "-", tDataframe$Class, fixed = TRUE)
  unqValuesOfClass <- gsub("_", "-", unqValuesOfClass, fixed = TRUE)
  names(unqColorsOfClass) <- gsub("_", "-", names(unqColorsOfClass), fixed = TRUE)
  
  highlighted <- tDataframe %>% filter(rank <= 10 & Class == "Unclassified")
  
  tDataframe <- tDataframe %>% filter(!(rank <= 10 & Class == "Unclassified"))
  
  ggplot(tDataframe %>% arrange(fct_rev(Class)), aes(x=rank, y=Counts, color = Class)) +
    ggtitle(title) +
    geom_point(size=1)+
    geom_jitter(width = jitter, height = jitter) +
    geom_point(data = highlighted %>% arrange(fct_rev(Class)), aes(x=rank, y=Counts, color = Class))+
    scale_color_manual(breaks = unqValuesOfClass, 
                       values = unqColorsOfClass) +
    xlab("Rank") + ylab("Count") +
    theme_classic() +
    geom_label_repel(
      data = subset(highlighted, rank <= 10 & Class == "Unclassified"),
      nudge_x = max(tDataframe$rank)/2-subset(highlighted, rank <= 10 & Class == "Unclassified")$rank,
      hjust = 0,
      aes(label = as.character(Mutation), color = Class, segment.size = 0.5),
      size = 5,
      show.legend = FALSE
    ) +
    scale_fill_manual(breaks = unqValuesOfClass, 
                      values = unqColorsOfClass) +
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text = element_text (size = 9, face = "plain"),
          axis.title = element_text (size = 9, face = "bold"),
          legend.title = element_text (size = 9, face = "bold"),
          legend.text = element_text (size = 9, face = "plain"))
}
# Generate result folder
singleMutresultfolder <- paste0(resultfolder,
                          "/Single mutations")
if(!dir.exists(singleMutresultfolder)) dir.create(singleMutresultfolder)
# Draw it!
svg(filename = paste0(singleMutresultfolder, "/BRAF mutation counts.svg"), width = 10/2.54, height = 8/2.54)
printRankedPlot(pGene = "BRAF", title = "BRAF", jitter = 0.01)
dev.off()

svg(filename = paste0(singleMutresultfolder, "/KRAS mutation counts.svg"), width = 10/2.54, height = 8/2.54)
printRankedPlot(pGene = "KRAS", title = "KRAS", jitter = 0.01)
dev.off()

svg(filename = paste0(singleMutresultfolder, "/EGFR mutation counts.svg"), width = 10/2.54, height = 8/2.54)
printRankedPlot(pGene = "EGFR", title = "EGFR", jitter = 0.01)
dev.off()

# Save in PowerPoint
my_doc <- read_pptx()

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/KRAS mutation counts.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/BRAF mutation counts.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/EGFR mutation counts.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)


print(my_doc, target = paste0(singleMutresultfolder, "/SingleMutationsOccurrence.pptx"))

################################################################################
# Single Mutation vs Gene statistics
################################################################################
# Now the cooccurrence of every single mutation/tuple with each gene should be tested.

### Now calculate cooccurrence and mutual exclusivity statistics between single mutations and genes
# First a summarizing dataframe with all class counts is needed
summaryMutationCounts <- data.frame(Gene = c(), Mutation = c(), Counts = c(), Freq = c())
for(i in genesToDo){
  summaryMutationCounts <- rbind(summaryMutationCounts,
                              mutationCounts[[i]][, c("Gene", "Mutation", "Class", "Counts", "Freq")])
}

write.table(x = summaryMutationCounts,
            file = paste0(fileresultfolder, "\\Single Mutation counts.txt"),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

# Get Dataframe with all combinations of mutations and genes
temp <- summaryMutationCounts[, c("Gene", "Mutation", "Class")]
coocStatsMutations <- NULL
for(i in genesToDo){
  tempdf <- temp[temp$Gene != i,]
  colnames(tempdf) <- c("GeneA", "MutationA", "ClassA")
  tempdf$GeneB <- rep(i, nrow(tempdf))
  
  coocStatsMutations <- rbind(coocStatsMutations,
                     tempdf)
}
# Statsv1: Now calculate statistics of cooccurrence of MutationA and GeneB in all patients with mutated GeneA
tempdf <- NULL
for(i in seq(nrow(coocStatsMutations))){
  actA <- coocStatsMutations$GeneA[i]
  actB <- coocStatsMutations$GeneB[i]
  mutationIndex <- mutationCounts[[actA]]$Index[mutationCounts[[actA]]$Mutation == coocStatsMutations$MutationA[i]]
  
  # all patients with mutation of classA if GeneA is mutated
  boolEventA <- mutationMatrices[[actA]][, mutationIndex][alterationList[[actA]] == TRUE & validList[[actB]] == TRUE]
  
  # all patients with GeneB mutated AND GeneA mutated
  boolEventB <- alterationList[[actB]][alterationList[[actA]] == TRUE & validList[[actB]] == TRUE]
  
  tempdf <- rbind(tempdf,
                  calcStats(boolEventA = boolEventA,
                            boolEventB = boolEventB))
}
# Bonferroni adjustment for each analysis seperately
# Significance declaration
actSignComb <- unique(coocStatsMutations[, c("GeneA", "GeneB")])
for(i in 1:nrow(actSignComb)){
  ids <- which(coocStatsMutations$GeneA == actSignComb$GeneA[i] & coocStatsMutations$GeneB == actSignComb$GeneB[i])
  tempdf$significantTwoSided <- tempdf$fisherTwoSided < 0.05/length(ids)
}

coocStatsMutationsv1 <- cbind(coocStatsMutations, tempdf)

# Save statistics
write.table(x = coocStatsMutationsv1,
            file = paste0(fileresultfolder, "\\Cooccurrence Statistics_Single Mutations vs Genes_Subset GeneA mutated.txt"),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)


### Plot cooccurrence of signle mutations vs. one entire gene
# Function for plot
printQuadrantPlot <- function(pGene, pDataframe, title = "", jitter = 0){
  
  # Calculate log values
  tDataframe <- pDataframe
  tDataframe$logvalFisher <- -log10(tDataframe$fisherTwoSidedAdjusted)
  tDataframe$logvalOdd <- log10(tDataframe$oddsRatio)
  
  # Color assignment and preprocessing
  unqValuesOfClass <- sort(unique(tDataframe$ClassA))
  unqColorsOfClass <- sapply(unqValuesOfClass, FUN = function(x){
    return(colorCode$Color[colorCode$Gene == pGene & colorCode$Class == x])
  })
  
  tDataframe$ClassA <- gsub("_", "-", tDataframe$ClassA, fixed = TRUE)
  unqValuesOfClass <- gsub("_", "-", unqValuesOfClass, fixed = TRUE)
  names(unqColorsOfClass) <- gsub("_", "-", names(unqColorsOfClass), fixed = TRUE)
  
  # ggplot(tDataframe %>% arrange(fct_rev(ClassA)), aes(x=logvalOdd, y=logvalFisher, color = fct_inorder(ClassA))) +
  #   ggtitle(title) +
  #   geom_hline(yintercept= -log10(0.05/nrow(tDataframe)), color = "black") +
  #   geom_vline(xintercept= log10(1), color = "black") +
  #   geom_jitter(width = jitter, height = jitter) +
  #   #geom_point() +
  #   xlab("Log odds ratio (Comutated/Non-comutated)") + ylab("-log10(p-value)") +
  #   scale_color_manual(breaks = unqValuesOfClass, 
  #                      values= unqColorsOfClass, name = "Class") +
  #   theme_classic() +
  #   geom_text_repel(
  #     data = subset(tDataframe, fisherTwoSided < 0.05/nrow(tDataframe)),
  #     aes(label = MutationA, color = ClassA),
  #     size = 5,
  #     max.overlaps = 100,
  #     show.legend = FALSE
  #   ) +
  #   scale_fill_manual(breaks = unqValuesOfClass, 
  #                     values= unqColorsOfClass) +
  #   scale_x_continuous(breaks = breaks_pretty(n = 6)(tDataframe$logvalOdd),
  #                      labels = paste0("10^", breaks_pretty(n = 6)(tDataframe$logvalOdd))) +
  #   theme(axis.text.x = element_markdown(),
  #         axis.text.y = element_markdown()) +
  #   theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
  #         axis.text = element_text (size = 8, face = "plain"),
  #         axis.title = element_text (size = 8, face = "bold"),
  #         legend.title = element_text (size = 8, face = "bold"),
  #         legend.text = element_text (size = 8, face = "plain"))
    #theme(legend.position = "none")
  
  ggplot(tDataframe %>% arrange(fct_rev(ClassA)), aes(x=oddsRatio, y=logvalFisher, color = fct_inorder(ClassA))) +
    ggtitle(title) +
    geom_hline(yintercept= -log10(0.05/nrow(tDataframe)), color = "black") +
    geom_vline(xintercept= 1, color = "black") +
    geom_jitter(width = jitter, height = jitter) +
    #geom_point() +
    xlab("Odds ratio (Comutated/Non-comutated)") + ylab("-log10(p-value)") +
    scale_color_manual(breaks = unqValuesOfClass, 
                       values= unqColorsOfClass, name = "Class") +
    theme_classic() +
    geom_label_repel(
      data = subset(tDataframe, significantTwoSided == TRUE),
      aes(label = MutationA, color = ClassA, segment.size = 0.5),
      size = 4,
      show.legend = FALSE
    ) +
    scale_fill_manual(breaks = unqValuesOfClass, 
                      values= unqColorsOfClass) +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    # theme(axis.text.x = element_markdown(),
    #       axis.text.y = element_markdown()) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text = element_text (size = 9, face = "plain"),
          axis.title = element_text (size = 9, face = "bold"),
          legend.title = element_text (size = 9, face = "bold"),
          legend.text = element_text (size = 9, face = "plain"))
}
# Draw quadrant plots with statistics of subset for one altered gene
svg(filename = paste0(singleMutresultfolder, "/QuadrantPlot_BRAF_Comutation EGFR_Only BRAF mutated.svg"), width = 10/2.54, height = 8/2.54)
printQuadrantPlot(pGene = "BRAF", pDataframe = coocStatsMutationsv1 %>% filter(GeneA == "BRAF" & GeneB == "EGFR"),
                  title = "BRAF mutations\nEGFR co-mutated", jitter = 0.01)
dev.off()

svg(filename = paste0(singleMutresultfolder, "/QuadrantPlot_BRAF_Comutation KRAS_Only BRAF mutated.svg"), width = 10/2.54, height = 8/2.54)
printQuadrantPlot(pGene = "BRAF", pDataframe = coocStatsMutationsv1 %>% filter(GeneA == "BRAF" & GeneB == "KRAS"),
                  title = "BRAF mutations\nKRAS co-mutated", jitter = 0.01)
dev.off()

svg(filename = paste0(singleMutresultfolder, "/QuadrantPlot_KRAS_Comutation EGFR_Only KRAS mutated.svg"), width = 10/2.54, height = 8/2.54)
printQuadrantPlot(pGene = "KRAS", pDataframe = coocStatsMutationsv1 %>% filter(GeneA == "KRAS" & GeneB == "EGFR"),
                  title = "KRAS mutations\nEGFR co-mutated", jitter = 0.01)
dev.off()

svg(filename = paste0(singleMutresultfolder, "/QuadrantPlot_KRAS_Comutation BRAF_Only KRAS mutated.svg"), width = 10/2.54, height = 8/2.54)
printQuadrantPlot(pGene = "KRAS", pDataframe = coocStatsMutationsv1 %>% filter(GeneA == "KRAS" & GeneB == "BRAF"),
                  title = "KRAS mutations\nBRAF co-mutated", jitter = 0.01)
dev.off()

svg(filename = paste0(singleMutresultfolder, "/QuadrantPlot_EGFR_Comutation KRAS_Only EGFR mutated.svg"), width = 10/2.54, height = 8/2.54)
printQuadrantPlot(pGene = "EGFR", pDataframe = coocStatsMutationsv1 %>% filter(GeneA == "EGFR" & GeneB == "KRAS"),
                  title = "EGFR mutations\nKRAS co-mutated", jitter = 0.01)
dev.off()

svg(filename = paste0(singleMutresultfolder, "/QuadrantPlot_EGFR_Comutation BRAF_Only EGFR mutated.svg"), width = 10/2.54, height = 8/2.54)
printQuadrantPlot(pGene = "EGFR", pDataframe = coocStatsMutationsv1 %>% filter(GeneA == "EGFR" & GeneB == "BRAF"),
                  title = "EGFR mutations\nBRAF co-mutated", jitter = 0.01)
dev.off()

# Save in PowerPoint
my_doc <- read_pptx()

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/QuadrantPlot_BRAF_Comutation EGFR_Only BRAF mutated.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/QuadrantPlot_BRAF_Comutation KRAS_Only BRAF mutated.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/QuadrantPlot_KRAS_Comutation EGFR_Only KRAS mutated.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/QuadrantPlot_KRAS_Comutation BRAF_Only KRAS mutated.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/QuadrantPlot_EGFR_Comutation KRAS_Only EGFR mutated.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)

my_doc <- add_slide(my_doc, layout = "Title and Content", master = "Office Theme")
my_doc <- ph_with(my_doc, value = external_img(src = paste0(singleMutresultfolder, "/QuadrantPlot_EGFR_Comutation BRAF_Only EGFR mutated.svg"), width = 10, height = 8, unit = "cm", guess_size = FALSE), location = ph_location_type(type = "body"), use_loc_size = FALSE)


print(my_doc, target = paste0(singleMutresultfolder, "/Quadrant plots with subset of alterations in primary gene.pptx"))


################################################################################
# Single Mutation vs Single Mutation statistics
################################################################################
# Now the cooccurrence of every single mutation/tuple of one gene withevery
# single mutation/tuple of another gene should be tested.

# Minimal count threshold for a mutation to be part of the analysis
minOccurrence <- 5

### Now calculate cooccurrence and mutual exclusivity statistics between single mutations

# Get Dataframe with all combinations of mutations
temp <- t(combinat::combn(seq(genesToDo), 2))
coocStatsMutationvsMutation <- NULL
for(i in seq(nrow(temp))){
  combns <- expand.grid(seq(nrow(mutationCounts[[temp[i,1]]])),
                        seq(nrow(mutationCounts[[temp[i,2]]])))
  tempdf <- cbind(mutationCounts[[temp[i,1]]][combns[,1], c("Gene", "Mutation", "Class")],
                  mutationCounts[[temp[i,2]]][combns[,2], c("Gene", "Mutation", "Class")])
  colnames(tempdf) <- c("GeneA", "MutationA", "ClassA", "GeneB", "MutationB", "ClassB")
  
  coocStatsMutationvsMutation <- rbind(coocStatsMutationvsMutation,
                                       tempdf)
}

## Stats: Now calculate statistics of cooccurrence of GeneA and GeneB in all patients with mutated GeneA or mutated GeneB
tempdf <- NULL

# Parallel computation
n_cores <- 1
# To enable parallel computation run the following line
n_cores <- detectCores()-2

#go
actCluster <- makeCluster(n_cores, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(actCluster)

tempdf <- foreach(i = seq(nrow(coocStatsMutationvsMutation)),
                  .combine = c) %dopar% {
  if(i%%10000 == 0) cat(paste0(i/nrow(coocStatsMutationvsMutation)*100), "%\n") #Sadly only works on Linux
  actA <- coocStatsMutationvsMutation$GeneA[i]
  actB <- coocStatsMutationvsMutation$GeneB[i]
  mutationIndexA <- mutationCounts[[actA]]$Index[mutationCounts[[actA]]$Mutation == coocStatsMutationvsMutation$MutationA[i]]
  mutationIndexB <- mutationCounts[[actB]]$Index[mutationCounts[[actB]]$Mutation == coocStatsMutationvsMutation$MutationB[i]]
  
  # all patients with mutation of MutationA AND GeneA or GeneB mutated
  boolEventA <- mutationMatrices[[actA]][, mutationIndexA][alterationList[[actA]] == TRUE | alterationList[[actB]] == TRUE]
  
  # all patients with mutation of MutationB AND GeneA or GeneB mutated
  boolEventB <- mutationMatrices[[actB]][, mutationIndexB][alterationList[[actA]] == TRUE | alterationList[[actB]] == TRUE]
  
  acttemp <- calcStats(boolEventA = boolEventA,
                       boolEventB = boolEventB)
  acttemp$significantTwoSided <- 0
  
  return(unlist(acttemp))
}

# Stop parallel computation
stopCluster(actCluster)

# Convert output to data.frame
dummyStats <- calcStats(boolEventA = rep(1,10), boolEventB = rep(1,10))
tempdf <- matrix(tempdf, ncol = ncol(dummyStats), byrow = TRUE)
tempdf <- as.data.frame(tempdf)
colnames(tempdf) <- colnames(dummyStats)

# Only consider mutations with counts at least of minimum threshold
coocStatsMutationvsMutation <- coocStatsMutationvsMutation[tempdf$CountsA >= minOccurrence &
                                                             tempdf$CountsB >= minOccurrence,]
tempdf <- tempdf[tempdf$CountsA >= minOccurrence &
                 tempdf$CountsB >= minOccurrence,]

# Bonferroni adjustment for each analysis seperately
# Significance declaration
actSignComb <- unique(coocStatsMutationvsMutation[, c("GeneA", "GeneB")])
for(i in 1:nrow(actSignComb)){
  ids <- which(coocStatsMutationvsMutation$GeneA == actSignComb$GeneA[i] & coocStatsMutationvsMutation$GeneB == actSignComb$GeneB[i])
  tempdf$significantTwoSided <- tempdf$fisherTwoSided < 0.05/length(ids)
}

# Calculate log10 of odds ratio
tempdf$log10OddsRatio <- log10(tempdf$oddsRatio)

# Calculate observed cooccurrence and expected cooccurrence
sumOfFourSquareTable <- tempdf$AandB+tempdf$notAandB+tempdf$AandnotB+tempdf$notAandnotB
tempdf$observedFreq <- (tempdf$AandB + 0.5)/sumOfFourSquareTable # v1.2: Do Haldane correction as well as did for odds ratio
tempdf$expectedFreq <- (tempdf$CountsA/sumOfFourSquareTable)*(tempdf$CountsB/sumOfFourSquareTable)

# Calculate ratio of observed cooccurrence vs. expected cooccurrence
tempdf$ratioObsExp <- tempdf$observedFreq/tempdf$expectedFreq

# Calculate log10 of observed cooccurrence vs. expected cooccurrence
tempdf$log10ratioObsExp <- log10(tempdf$ratioObsExp)

# Combine all information
coocStatsMutationvsMutationv2 <-  cbind(coocStatsMutationvsMutation, tempdf)

# Save statistics
write.table(x = coocStatsMutationvsMutationv2,
            file = paste0(fileresultfolder, "\\Cooccurrence Statistics_Single Mutations vs Single Mutations_Subset GeneA or GeneB mutated.txt"),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

# Function to plot heatmap of cooccurrence
plot_heatmap_cooc <- function(geneA, geneB,
                              pDataframe = coocStatsMutationvsMutationv2,
                              columnAsData = "log10OddsRatio",
                              ylab = "log10(Odds ratio)",
                              pminOccurrence = minOccurrence,
                              filename){
  
  # Create legend titles for classes
  legend_column_title = paste0(geneA, " classes")
  legend_row_title =  paste0(geneB, " classes")
  
  # Create axis titles 
  column_title = paste0(geneA, " variants")
  row_title =  paste0(geneB, " variants")
  
  # Get data for heatmap
  subdata <- pDataframe %>% filter(GeneA == geneA & GeneB == geneB & CountsA >= pminOccurrence & CountsB >= pminOccurrence)
  
  # Convert data to matrix
  heatmap_data <- matrix(data = subdata[,columnAsData],
                         ncol = length(unique(subdata$MutationA)),
                         byrow = TRUE)
  colnames(heatmap_data) <- unique(subdata$MutationA)
  rownames(heatmap_data) <- unique(subdata$MutationB)
  
  # Get boolean matrix whether p-value of cell is significant
  signif_heatmap_data <- matrix(data = subdata$significantTwoSided == TRUE,
                                ncol = length(unique(subdata$MutationA)),
                                byrow = TRUE)
  colnames(signif_heatmap_data) <- unique(subdata$MutationA)
  rownames(signif_heatmap_data) <- unique(subdata$MutationB)
  
  # Get color values for classes
  col_vec_column <- colorCode$Color[colorCode$Class %in% unique(subdata$ClassA) &
                                    colorCode$Gene == geneA]
  names(col_vec_column) <- sapply(col_vec_column, FUN = function(x){
    return(colorCode$Class[colorCode$Color == x & colorCode$Gene == geneA])
  })
  col_vec_row <- colorCode$Color[colorCode$Class %in% unique(subdata$ClassB) &
                                 colorCode$Gene == geneB]
  names(col_vec_row) <- sapply(col_vec_row, FUN = function(x){
    return(colorCode$Class[colorCode$Color == x & colorCode$Gene == geneB])
  })
  
  lgd_column <- Legend(at = names(col_vec_column), legend_gp = gpar(fill = col_vec_column), title = legend_column_title)
  lgd_row <- Legend(at = names(col_vec_row), legend_gp = gpar(fill = col_vec_row), title = legend_row_title)
  
  breaks <- 10^(0:(floor(log10(max(c(subdata$CountsA, subdata$CountsB))))))
  
  # Calculate Heatmap
  ret <- Heatmap(heatmap_data,
                 column_title = column_title,
                 row_title = row_title,
                 name = ylab,
                 col = colorRamp2(c(min(heatmap_data), 0, max(heatmap_data)), c("blue", "white", "red")),
                 column_dend_height = unit(2, "cm"), 
                 row_dend_width = unit(2, "cm"),
                 column_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                 row_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                 column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                 row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                 right_annotation = rowAnnotation(" " = anno_simple(x = subdata$ClassB[!duplicated(subdata$MutationB)],
                                                                    col = col_vec_row,
                                                                    border = TRUE, width = unit(x = 0.35, units = "cm")),
                                                  "  " = anno_barplot(x = log10(subdata$AandB[!duplicated(subdata$MutationB)] + subdata$notAandB[!duplicated(subdata$MutationB)]),
                                                                      gp = gpar(fill = "black"),
                                                                      axis_param = list(at = log10(breaks), labels =breaks),
                                                                      width = unit(x = 1, units = "cm"))),
                 bottom_annotation = columnAnnotation(" " = anno_simple(x = subdata$ClassA[!duplicated(subdata$MutationA)],
                                                                        col = col_vec_column,
                                                                        border = TRUE, height = unit(x = 0.35, units = "cm")),
                                                      "  " = anno_barplot(x = log10(subdata$AandB[!duplicated(subdata$MutationA)] + subdata$AandnotB[!duplicated(subdata$MutationA)]),
                                                                          gp = gpar(fill = "black"),
                                                                          axis_param = list(at = log10(breaks), labels =breaks),
                                                                          height = unit(x = 1, units = "cm"))),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.rect(x = x, y = y, width = width, height = height, 
                             gp = gpar(col = "black", fill = NA))
                   if(signif_heatmap_data[i, j] == TRUE)
                     grid.text(".", x = x, y = y, 
                               gp = gpar(fontsize = 12, fontface = "bold"))
                 })
  
  #svg(filename = filename, width = (ncol(heatmap_data)*10+600)*0.035/2.54, height = (nrow(heatmap_data)*10+350)*0.035/2.54)
  svg(filename = filename, width = (ncol(heatmap_data)*0.3+8.5)/2.54, height = (nrow(heatmap_data)*0.3+4.5)/2.54)
  draw(ret)
  dev.off()
}

# Generate result folder
heatmapresultfolder <- paste0(resultfolder,
                          "/heatmap")
if(!dir.exists(heatmapresultfolder)) dir.create(heatmapresultfolder)

# Plot heatmaps with mutations occurring minimum 5 times and log10 of odds ratio
plot_heatmap_cooc(geneA = "KRAS", geneB = "BRAF",
                  filename = paste0(heatmapresultfolder,
                                    "/KRAS_and_BRAF_log odds ratio.svg"))

plot_heatmap_cooc(geneA = "KRAS", geneB = "EGFR",
                  filename = paste0(heatmapresultfolder,
                                    "/KRAS_and_EGFR_log odds ratio.svg"))

plot_heatmap_cooc(geneA = "BRAF", geneB = "EGFR",
                  filename = paste0(heatmapresultfolder,
                                    "/BRAF_and_EGFR_log odds ratio.svg"))

################################################################################
# Save sample data
################################################################################
write.table(x = pp_inputData,
            file = paste0(fileresultfolder, "\\Sample Data Class assignment and tuples.txt"),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

print(Sys.time()-t1)
