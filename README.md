# The HDAC inhibitor CI-994 acts as a molecular memory aid by facilitating synaptic and intra-cellular communication after learning
Here is the code used to perform the analyses presented in the manuscript, "The HDAC inhibitor CI-1 994 acts as a molecular memory aid by facilitating synaptic and intra-cellular communication after learning" [(Burns et al., 2021, BioRxiv)](https://www.biorxiv.org/content/10.1101/2021.09.21.460970v1). 

The raw data (fastq files) for this analysis is in the process of being uploaded to the GEO database and will be linked once that is accepted.

### Abstract
Long-term memory formation relies on synaptic plasticity, activity-dependent transcription and epigenetic modifications. Multiple studies have shown that HDAC inhibitor (HDACi) treatments can enhance individual aspects of these processes, and thereby act as putative cognitive enhancers. However, their mode of action is not fully understood. In particular, it is unclear how systemic application of HDACis, which are devoid of substrate specificity, can target pathways that promote memory formation. In this study, we explore the electrophysiological, transcriptional and epigenetic responses that are induced by CI-994, a class I HDAC inhibitor, combined with contextual fear conditioning (CFC) in mice. We show that CI-994-mediated improvement of memory formation is accompanied by enhanced long-term potentiation in the hippocampus, a brain region recruited by CFC, but not in the striatum, a brain region not primarily implicated in contextual memory formation. Furthermore, using a combination of bulk and single cell RNA sequencing, we find that synaptic plasticity-promoting gene expression cascades are more strongly engaged in the hippocampus than in the striatum, but only when HDACi treatment co-occurred with CFC, and not by either treatment alone. Lastly, using ChIP-sequencing, we show that the combined action of HDACi application and conditioning is required to elicit enhancer histone acetylation in pathways that may underlie improved memory performance. Together, our results indicate that systemic HDACi administration amplifies brain-region specific processes that are naturally induced by learning. These findings shed light onto the mode of action of HDACis as cognitive enhancers.

### RNA-seq analysis
We used bulk RNA-sequencing in the hippocampus and the striatum to show that a systemic HDACi treatment differentially alters the transcriptional programs in two different brain regions after HDACi treatment and contextual fear conditioning. Code should be run in the following order, however, the data provided in the files directory should be sufficient to start analysis from `3_DEanalysis.R`: 

- `1_runAlignment.sh`: align .fastq files to mouse genome (mm10, Ensembl 93) with STAR. This script is specifically used to run alignments on a cluster using a slurm workload manager.

- `2_CountGenes.R`: Count number of reads aligning to exonic regions. 

- `3_DEanalysis.R`: Run the Differential Expression (DE) analysis usinge DEseq2.

- `4_TrajectoryAnalysis.R`: Cluster genes based on their log2FC values in each group compared to the baseline group (VEH-Context).


### snRNA-seq analysis

### ChIP-seq analysis
