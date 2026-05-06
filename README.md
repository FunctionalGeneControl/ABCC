# Activity-By-Captured-Contact (ABCC)

ABCC (also known as CHiC-ABC) is a computational model for predicting the effects of enhancers on gene expression based on chromatin activity information and chromosomal contact data from high-resolution Promoter Capture Hi-C (PCHi-C). It is an modification of the Activity-By-Contact (ABC) model developed by Jesse Engreitz's lab [1] [2] that uses Hi-C data as a source of chromosomal contact information.

Promoter Capture Hi-C (PCHi-C) is a method of enriching Hi-C libraries for contacts involving (at least on one end) gene promoters using hybridisation probes. Please see Freire-Pritchett et al. [1] for information about PCHi-C and its analysis  tools.

ABCC leverages the statistical modelling of PCHi-C data by CHiCAGO software [1][2] for normalisation and imputation, and uses the results as a source of contact information, as described in more detail below.

The ABC code used for this implementation was downloaded from the repo https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) on 22 September 2022 (latest commit: a74aa73e66329968bc2d4068464ee9e8cacf86a1). 

## Requirements
For each cell-type, the inputs to the ABC model are:

Inputs:

 	* bam file for Dnase-Seq or ATAC-Seq (indexed and sorted)
 	* bam file for H3K27ac ChIP-Seq (indexed and sorted) (optional)
 	* PCHi-C data (see the PCHi-C imputation section below)
 	* A measure of gene expression (gene-level TPM)

In addition the following (non-cell-type specific) genome annotation files are required:

 * bed file containing gene annotations (may change across cell types if using cell-type specific TSS's)
 * bed file containing chromosome annotations

## The principle of the ABC(C) Model

The Activity by Contact (ABC) model [1][2] is designed to represent a mechanistic model in which enhancers activate gene transcription upon enhancer-promoter contact. In a simple conception of such a model, the quantitative effect of an enhancer depends on the frequency with which it contacts a promoter multiplied by the strength of the enhancer (i.e., the ability of the enhancer to activate transcription upon contacting a promoter). Moreover, the contribution of a specific enhancer to a gene’s expression should depend on the surrounding context (ie, the strength and contact frequency of other enhancers for the gene). 

To convert this conceptual framework into a practical score (which can be applied genome-wide), we formulated the ABC score:

ABC score for effect of element E on gene G = Activity of E × Contact frequency between E and G /  Sum of (Activity × Contact Frequency) over all candidate elements within 5 Mb.

Activity (A) is defined as the geometric mean of the read counts of DNase-seq and H3K27ac ChIP-seq at an element E, and Contact (C) frequency between E and the promoter of gene G. Elements are defined as ~500bp regions centered on DHS/ATAC peaks. In the original ABC, the KR normalized Hi-C contact frequency is used as Contact. In ABCC, CHiCAGO-normalised and imputed contact frequencies are used instead.

## Running the PCHiC-ABC Model
Running the PCHiC-ABC model consists of the following steps:

 1. Define candidate enhancer and gene regions (the same step as in the original ABC model)
 2. Running PCHiC imputation (a new step for ABCC)
 3. Compute ABC Scores (using a predict.py script modified from the original ABC codebase)

### Step 1. Defining Candidate Enhancers and Genes.

Either the original ABC code used in this implementation or the newly developed ABC snakemake pipeline can be used for this purpose, as described in https://abc-enhancer-gene-prediction.readthedocs.io/.

* If running the origial ABC version for this step using ``makeCandidateRegions.py`` and ``run.neighborhoods.py``, you may find this notebook useful: https://hoellin.github.io/eg/notes_ABC/generic_notebooks/turnkey_notebook_to_run_ABC_with_example_over_GM12878.html.

#### Information on defining candidate elements from the original ABC method is copied below:
'Candidate elements' are the set of putative enhancers; ABC scores will be computed for all 'Candidate elements' within 5Mb of each gene. In computing the ABC score, the product of DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq reads will be counted in each candidate element. Thus the candidate elements should be regions of open (nucleasome depleted) chromatin of sufficient length to capture H3K27ac marks on flanking nucleosomes. In [1], we defined candidate regions to be 500bp (150bp of the DHS peak extended 175bp in each direction). Given that the ABC score uses absolute counts of Dnase-seq reads in each region, ```makeCandidateRegions.py ``` selects the strongest peaks as measured by absolute read counts (not by pvalue). In order to do this, we first call peaks using a lenient significance threshold (.1 in the above example) and then consider the peaks with the most read counts. This procedure implicitly assumes that the active karyotype of the cell type is constant. We recommend removing elements overlapping regions of the genome that have been observed to accumulate anomalous number of reads in epigenetic sequencing experiments (‘block-listed regions’). For convenience, we provide the list of block-listed regions available from <https://sites.google.com/site/anshulkundaje/projects/blacklists>. We also force the candidate enhancer regions to include gene promoters, even if the promoter is not among the candidate elements with the strongest signals genome-wide in a cell type, by specifying `--region_includelist`. 

### Step 2. Running PCHiC imputation

To adapt ABC for PCHi-C data, we took advantage of the CHiCAGO normalisation algorithm and developed an imputation procedure in the normalised counts space based on the inferred decay of interaction read counts with distance. 

As we do not expect the frequency of enhancer-promoter contacts to fall below levels expected due to Brownian collision, for a given pair of fragments involving a baited promoter we select the maximum between the CHiCAGO-normalised observed read counts (Nobs) and expected read counts Nexp estimated as Nexp = Bmean/(si*sj), where Bmean is the CHiCAGO-estimated Brownian noise level and si and sj and the bait- and other end-specific scaling factors. For promoters that could not be baited in the Capture Hi-C design and those, the reads for which were filtered out due to QC fail, we estimate the expected normalised read count directly from the interaction distance d, using the distance function f(d), estimated by CHiCAGO. Please refer to Additional File 1 in the publication presenting the CHiCAGO pipeline [2] for the formal definition of these parameters and their estimation procedures.

To esitmate the N imputed for candidate regions, as presented in Candidate Enhancer list, please use the following procedure:

```
k562_roadmap=/path/to/ABC_data/input/K562_roadmap_DpnII_frag/
pchic_rds=${k562_roadmap}K562_fres_merged_reweighted.Rds
design=/path/to/Design_for_PCHIC/Human_DpnII_binsize1500_maxL75K/
enhancerdir=${k562_roadmap}
split_pchic=${k562_roadmap}K562_frag_dpnii_
distout_rds=${k562_roadmap}K562_frag_dpnii_dist.rds
imputed_pchic_prefix=${k562_roadmap}imputed_contact/
mkdir $imputed_pchic_prefix

Rscript imputation_script_whole_genome.r ${pchic_rds} ${design} ${enhancerdir} ${split_pchic} ${distout_rds} ${imputed_pchic_prefix}

```
The main output from this will be a directory, with a subfolder for each chromosome, where will contain imputed PCHiC contact frequency file in bedpe format (eg. chr1.bedpe.gz), as well as a separate file with CHiCAGO scores for contacts. This directory structure can be used down the pipeline in the next step.

Imputation is a computationally heavy script, for speed we devised such that it is done by chromosome, so each candidate lists need to be split into individual chromosomes. We attach the script to separate the files chromosome by chromosome. Chromosome splitting and saving distance parameters are also needed to be performed before running the imputations. The following command using our scripts will separate the files into individual CHiCAGO object chromosome files as well as will save a distance parameters file. These steps are run in the background, the supporting files are present in the src folder. 

### Step 3. Computing the ABC Score

Compute ABC scores by combining Activity (as calculated by ```run.neighborhoods.py```) and imputed PCHi-C. We modified predict.py to work with output from PCHiC imputation, so please refer to our version of predictor.py when running the ABC Scores using our version of imputed PCHiC. Please include ```--hic_type bedpe``` flag.

Sample Command:

```
python src_mod/predict.py \
--enhancers ${k562_roadmap}/EnhancerList.txt \
--genes ${k562_roadmap}/GeneList.txt \
--HiCdir ${k562_roadmap}imputed_contact/ \
--chrom_sizes reference/chrom.sizes.og \
--threshold .02 \
--cellType K562 \
--hic_type bedpe \
--outdir ABC_scores \
--make_all_putative
```

## Human Genome assembly conversion.
For conversion of EnahncerList and GeneList files, we have added the script enhancerliftover that performs liftover the files between Hg19 and Hg38. The use is shown below. The .chain file can be found separately. 

```
Rscript enhancerliftoverwrapper.R EnhancerList.txt GeneList.txt
```

## Tips
* We found that an optimal ABCC score cutoff is highly data-dependent. One way to estimate a cutoff is plot a correlation between ABCC numerators (i.e., unnormalised ABCC scores) per gene and gene expression levels across a range of ABCC score cutoffs and choose the lowest cutoff, at which an appreciable correlation is achieved.
* Please note the gene expression file should be tab separated with no column names.

## Contributors
* ABCC development was led by Pavel Artemov supervised by Mikhail Spivakov, with valuable contributions by Joseph Ellaway, Valeriya Malysheva and Helen Ray-Jones.

## Acknowledgements
* We thank the Engreitz lab for developing the original ABC model and making the codebase accessible.
* We thank Jonathan Cairns for leading the development of CHiCAGO in Mikhail Spivakov's lab and helpful discussions about CHiCAGO-based imputation.
* We thank Stephen Rong for testing the ABCC software and helpful suggestions on improving its usability. 

## Contact
Please submit a github issue with any questions or if you experience any issues/bugs. 

### Dependencies
The codebase relies on the following dependencies (tested version provided in 
parentheses):

```
Python (3.6.4)
CHiCAGO
samtools (0.1.19)
bedtools (2.26.0)
Tabix (0.2.5) - Partial dependancy
MACS2 (2.1.1.20160309) - Partial dependancy
Java (1.8) - Partial dependancy
Juicer Tools (1.7.5) - Partial dependancy

Python packages:
pyranges (0.0.55)
numpy (1.15.2)
pandas (0.23.4)
scipy (0.18.1)
pysam (0.15.1)
pyBigWig (0.3.2) - Partial dependancy
```

### References 
[1] Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0

[2] Nasser J, Bergman DT, Fulco CP, Guckelberger P, Doughty BR, Patwardhan TA, Jones TR, Nguyen TH, Ulirsch JC, Lekschas F, Mualim K, Natri HM, Weeks EM, Munson G, Kane M, Kang HY, Cui A, Ray JP, Eisenhaure TM, Collins RL, Dey K, Pfister H, Price AL, Epstein CB, Kundaje A, Xavier RJ, Daly MJ, Huang H, Finucane HK, Hacohen N, Lander ES, Engreitz JM. Genome-wide enhancer maps link risk variants to disease genes. Nature. 2021 May;593(7858):238-243. doi: 10.1038/s41586-021-03446-x

[3] Freire-Pritchett, P., Ray-Jones, H., Della Rosa, M. et al. Detecting chromosomal interactions in Capture Hi-C data with CHiCAGO and companion tools. Nat Protoc 16, 4144–4176 (2021) https://doi.org/10.1038/s41596-021-00567-5

[4] Cairns, J., Freire-Pritchett, P., Wingett, S.W. et al. CHiCAGO: robust detection of DNA looping interactions in Capture Hi-C data. Genome Biol 17, 127 (2016). https://doi.org/10.1186/s13059-016-0992-2

### Citation
Malysheva et al., High-resolution promoter interaction analysis implicates genes involved in the activation of Type 3 Innate Lymphoid Cells in autoimmune disease risk. bioRxiv doi: https://doi.org/10.1101/2022.10.19.512842.
