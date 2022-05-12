<a href="https://zenodo.org/badge/latestdoi/402792893"><img src="https://zenodo.org/badge/402792893.svg" alt="DOI"></a>

# The HDAC inhibitor CI-994 acts as a molecular memory aid by facilitating synaptic and intra-cellular communication after learning
Here is the code used to perform the analyses presented in the manuscript, "The HDAC inhibitor CI-1 994 acts as a molecular memory aid by facilitating synaptic and intra-cellular communication after learning" [(Burns et al., 2021, BioRxiv)](https://www.biorxiv.org/content/10.1101/2021.09.21.460970v1). 

The raw data (fastq files) for this analysis is in the process of being uploaded to the GEO database and will be linked once that is accepted.

### Abstract
Long-term memory formation relies on synaptic plasticity, activity-dependent transcription and epigenetic modifications. Multiple studies have shown that HDAC inhibitor (HDACi) treatments can enhance individual aspects of these processes, and thereby act as putative cognitive enhancers. However, their mode of action is not fully understood. In particular, it is unclear how systemic application of HDACis, which are devoid of substrate specificity, can target pathways that promote memory formation. In this study, we explore the electrophysiological, transcriptional and epigenetic responses that are induced by CI-994, a class I HDAC inhibitor, combined with contextual fear conditioning (CFC) in mice. We show that CI-994-mediated improvement of memory formation is accompanied by enhanced long-term potentiation in the hippocampus, a brain region recruited by CFC, but not in the striatum, a brain region not primarily implicated in contextual memory formation. Furthermore, using a combination of bulk and single cell RNA sequencing, we find that synaptic plasticity-promoting gene expression cascades are more strongly engaged in the hippocampus than in the striatum, but only when HDACi treatment co-occurred with CFC, and not by either treatment alone. Lastly, using ChIP-sequencing, we show that the combined action of HDACi application and conditioning is required to elicit enhancer histone acetylation in pathways that may underlie improved memory performance. Together, our results indicate that systemic HDACi administration amplifies brain-region specific processes that are naturally induced by learning. These findings shed light onto the mode of action of HDACis as cognitive enhancers.

### RNA-seq analysis
We used bulk RNA-sequencing in the hippocampus and the striatum to show that a systemic HDACi treatment differentially alters the transcriptional programs in two different brain regions after HDACi treatment and contextual fear conditioning. Code should be run in the following order, however, the data provided in the files directory should be sufficient to start analysis from `3_DEanalysis.R`: 

- `1_runAlignment.sh`: align .fastq files to mouse genome (mm10, Ensembl 93) with STAR. This script is specifically used to run alignments on a cluster using a slurm workload manager.

- `2_CountGenes.R`: count number of reads aligning to exonic regions. 

- `3_DEanalysis.R`: run the Differential Expression (DE) analysis usinge DEseq2.

- `4_TrajectoryAnalysis.R`: cluster genes based on their log2FC values in each group compared to the baseline group (VEH-Context).


### snRNA-seq analysis
We used snRNA-sequencing in the hippocampus to show that systemic HDACi treatment alters differential transcriptional programs in the different cell types after treatment paired with contextual fear conditioning. Code should be run in the following order:

- `1_runAlignment.sh`: align fastq files to mouse genome (mm10-premrna; created with `cellranger mkfastq`) using CellRanger.

- `2_SeuratAnalysis.R`: remove doublets, filter nuclei and normalize nuclear gene counts for downstream UMAP clustering with Seurat. After clustering, assign clusters to cell types and locations within the nucleus. 

- `3_DEanalysis.R`: differential expression analysis (logistic framework) between HDACi and Vehicle treatments paired with contextual fear conditioning for each cell type. 

- `4_Augur_cellComposition.R`: determine cell type perturbation prioritization (not biased by nuclear numbers) and cell number changes after HDACi treatment.

- `5_DEG_removal_test.R`: analysis associated with Figures 3E-F and Supplemental Figure 8 in the manuscript. Remove up and down-regulated genes (DEGs) in each cell type to determine whether or not DEGs are driving the HDACi mediated split of Excitatory Neurons of the DG.


### ChIP-seq analysis
Based on the results from the snRNA-sequencing analysis, we performed ChIP-sequencing in Neun+ cells (neurons) of the dentate gyrus (DG). Animals underwent similar treatment and behavior as in the RNA-sequencign experiment (drug: Veh/HDACi; behavior: context/CFC) and an hour after behavior the DG of each animal was collected. Nuclei were extracted and Neun+ nuclei were sorted before being combined and IPed fo H3K27ac. Code for the analysis of the ChIP libraries should be run in the following order:

- `1_appendMIDs.R`: Append MID information from Read_1 to Read_2 for removal of duplicates and build adapter  sequence files for adapter trimming (2_trimAdapt.sh).

- `2_trimAdapt.sh`: Trim adapters (from library preparation) and output updated fastq files.

- `3_runAlignment.sh`:  Align trimmed fastq files to mm10 genome with bowtie2 then convert sam to bam and sort bam files.

- `4_filterBamFiles.sh`:  Call the perl script `4_rmDupByMids.pl` (provided by Active Motif) to remove duplicated reads based on MID sequences. Remove low mapping alignments (mapq <= 40) and re-index bam files.

- `5_peakCalls.sh`: Run MACS2 for peak calling in broad peak mode. 

- `6_DiffBind.R`: Differential Expression Analysis using DiffBind. Creates a diffbind object that accounts for peak locations and reads within those peaks.

- `7_ChromHMM_Homer.sh`: Run chromHMM on bam files for multiple histone post-translational modifications (data coming from Halder et al. 2016 collected 1 hour after CFC). Run Homer `annotatePeaks.pl` to assign peaks to gene regions based on proximity. 

- `8_AssignPeaks.R`: Assign the peaks that were defined in `6_DiffBind.R` to the chromatin states and genes defined in `7_ChromHmm_Homer.sh` and check the assignment statistics.
