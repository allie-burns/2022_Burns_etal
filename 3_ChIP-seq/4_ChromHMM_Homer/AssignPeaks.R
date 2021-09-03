library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(nationalparkcolors)

################################################################################
## Load data
################################################################################
## Load Homer
homer  <-  readRDS(
    "../annotation_files/Homer_peakassignments/annotate_peaks/homer_peaks_all.rds"
)

## load ChromHMM
files <- list.files("../annotation_files/ChromHMM_out/chromHMM_beds/",full.names = TRUE)
chrom <- lapply(seq(1:length(files)), function(i) {
    ## Setup files and naming
    file <- files[[i]]
    names <- strsplit(basename(files[[i]]),"-")[[1]][[2]]
    names <- substr(names,1,nchar(names) -4)
    
    ## read table
    tab <- read.table(file)
    colnames(tab)  <- c("chr","start","end","strand","score","unknown")
    tab$enhancer_type  <- names
    tab
})

chrom <- do.call(rbind,chrom)
chrom$chr <- gsub("^","chr",chrom$chr)

## Convert to GRanges
homer  <- GRanges(seqnames    = homer$chr,
                  ranges      = IRanges(start = homer$start,
                                        end = homer$end),
                  strand      = "*",
                  id          = homer$id,
                  annotation  = homer$annotation,
                  detailed_annotation = homer$detailed_annotation,
                  dist_tss    = homer$dist_tss,
                  ensembl     = homer$ensembl,
                  gene        = homer$gene,
                  alias       = homer$alias,
                  description = homer$description,
                  biotype     = homer$biotype)
homer <- sortSeqlevels(homer)
homer <- sort(homer)

chrom <- GRanges(seqnames    = chrom$chr,
                 ranges      = IRanges(start = chrom$start + 2,
                                       end   = chrom$end),
                 strand      = "*",
                 chrom_state = chrom$enhancer_type)
chrom <- sortSeqlevels(chrom)
chrom <- sort(chrom)

################################################################################
## Find overlaps
################################################################################
options(dplyr.summarise.inform = FALSE)
assignRegions <- function(i) {
    ## Find overlaps > 100bp
    ols <- findOverlaps(homer[i],chrom, minoverlap = 100)
    ## Find %overlapping
    overlaps <- pintersect(chrom[subjectHits(ols)], homer[i][queryHits(ols)])
    ## Merge same states
    states <- data.frame(chromHMM_state = chrom[subjectHits(ols)]$chrom_state,
                         perc = width(overlaps) / width(homer[i][queryHits(ols)]))
    states  <- states %>%
      group_by(chromHMM_state) %>%
      summarise(chromHMM_perc = round(sum(perc)*100,3))
    ## Which state does peak have highest overlap with?
    states[states$chromHMM_perc == max(states$chromHMM_perc),]
}

states <- lapply(seq(1:length(homer)), assignRegions)

##save.states <- states
states <- save.states

## Randomly assign states with equal overlap
## There are only 300 and all have less than 50% so for now, this is the easiest
set.seed(775)
states[unlist(lapply(states,nrow)) > 1] <- 
    lapply(states[unlist(lapply(states,nrow)) > 1], function(x) {
        x[sample(x = 1:nrow(x),1),]
    })

states <- do.call(rbind, states)  
peak_anns <- cbind(homer,states) ## add chromHMM to homer information

## save homer
## saveRDS(peak_anns,
##         "../annotation_files/Homer_peakassignments/homer_chromHMM_peakassignments.rds")

################################################################################
## Homer/ChromHMM characterization (and sanity checks)
################################################################################
cols <- c("#F9F5EA","#E3EEF4","#6BBAE5", "#81974C", "#454B68")
peak_anns <- readRDS(
    "../annotation_files/Homer_peakassignments/homer_chromHMM_peakassignments.rds"
)

peak_anns <- data.frame(peak_anns)

## merge promoter groups
peak_anns$chromHMM_state <- recode(peak_anns$chromHMM_state,
                                   "weak_promoter" = "promoter",
                                   "strong_promoter" = "promoter")

peak_anns$chromHMM_state <- factor(peak_anns$chromHMM_state,
                                   levels = c("control",
                                              "repressed_regions",
                                              "promoter",
                                              "active_enhancer",
                                              "primed_enhancer"))


## Remove Control
peak_anns <- peak_anns[!peak_anns$chromHMM_state == "control",]
peak_anns$chromHMM_state <- droplevels(peak_anns$chromHMM_state)

## Percentage of overlaps based on chromatin state
ggplot(peak_anns, aes(x = chromHMM_state, y = chromHMM_perc, fill = chromHMM_state)) +
    geom_boxplot() +
    ylim(0,100) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    ylab("% peak overlap with ChromHMM") + 
    scale_fill_manual(values = cols) ##park_palette("Yosemite"))

peak_anns <- peak_anns[peak_anns$chromHMM_perc > 50,]

## compare chomHMM state lengths
state_length <- data.frame(state = chrom$chrom_state,
                           length = width(chrom))
state_length$state <- recode(state_length$state,
                             "weak_promoter" = "promoter",
                             "strong_promoter" = "promoter")
state_length$state <- factor(state_length$state,
                             levels = c("control",
                                        "repressed_regions",
                                        "promoter",
                                        "active_enhancer",
                                        "primed_enhancer"))

ggplot(state_length, aes(x = state, y = length, fill = state)) +
    geom_boxplot() +
    ylim(0,1000000) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    scale_fill_manual(values = cols)##park_palette("GreatBasin"))


## Percentages of ChromHMM types
## % by label
ggplot(data.frame(table(peak_anns$chromHMM_state)), aes(x="", y = Freq, fill = Var1)) + 
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    theme_void() +
    theme(axis.text.x=element_blank()) +
    scale_fill_manual(values = cols)##park_palette("GreatBasin"))

## % by coverage
perc_cov  <-  peak_anns %>%
    select(chromHMM_state, chromHMM_perc) %>%
    group_by(chromHMM_state) %>%
    summarise(p_freq = sum(chromHMM_perc)) %>%
    ungroup() %>%
    mutate(perc = p_freq/sum(p_freq) * 100) 

ggplot(perc_cov, aes(x="", y = perc, fill = chromHMM_state)) + 
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    theme_void() +
    theme(axis.text.x=element_blank()) +
    scale_fill_manual(values = cols)##park_palette("GreatBasin"))

## ratio of total state vs peak states
state_cov <-  state_length %>% 
    group_by(state) %>%
    summarise(tot_cov = sum(length))
state_cov <- state_cov[!state_cov$state == "control",]

peak_cov  <- peak_anns %>%
    select(seqnames, start, end, strand, chromHMM_state) %>%
    mutate(width = end - start) %>%
    group_by(chromHMM_state) %>%
    summarise(tot_cov = sum(width))

peak_cov$perc  <- peak_cov$tot_cov / state_cov$tot_cov * 100

ggplot(peak_cov, aes(x="", y = perc, fill = chromHMM_state)) + 
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(values = cols)##park_palette("GreatBasin"))

ggplot(peak_cov, aes(x="", y = perc, fill = chromHMM_state)) + 
    geom_bar(width = 1, stat = "identity", position = position_dodge(), color = "black") +
    coord_polar("y", start=0) +
    theme_minimal() +
    scale_y_continuous(limits = c(0,100)) + ##name, breaks, labels, limits, trans)
    theme(axis.text.x=element_blank()) +
    scale_fill_manual(values = cols)

