![haplomagic2](https://github.com/GoliczGenomeLab/haploMAGIC/assets/134378980/7e5d9334-7760-456e-ab8e-68e86156e7bc)

# README
haploMAGIC is a pipeline for phasing and CO detection in multiparental populations (MPPs) with single-nucleotide polymorphism (SNP) genotypes of individuals derived from inbred founder lines. haploMAGIC was specifically designed for MPPs and outperformed similar tools, duoHMM and LINKPHASE3, in F1 scores obtained detecting recombination events in simulated MPPs (Montero-Tena et al. 2024). Notably, haploMAGIC performs consistently with high genotyping error rates.

## Publication
Montero-Tena, J. A., Abdollahi Sisi, N., Kox, T., Abbadi, A., Snowdon, R. J., & Golicz, A. A. (2024). haploMAGIC: accurate phasing and detection of recombination in multiparental populations despite genotyping errors. G3 (Bethesda, Md.), 14(8), jkae109. https://doi.org/10.1093/g3journal/jkae109

## Quick setup

#### Download repository and execute on the example data:
```
Rscript 0.haplomagic.R example 1 3 imputeTHonly correctFalseHom 10000
```
#### For more information, read the following manual.

## User's manual

#### For running haploMAGIC, install the scripts in the same directory with the input.

#### The 2 inputs required by haploMAGIC for each run consist in:
* One PED file in 12 format, containing the pedigree info and the SNP genotypes of all the individuals of the population. The naming format must be *[population]_[chromosome].ped*. PED files show genotypes by rows with space-separated columns, where the first 6 contain pedigree information and the rest X columns contain the alleles in each SNP (0, missing ; 1, major ; 2, minor), with X being twice the number of rows in the MAP file. Pedigree columns are generation, named *GX* with X from 0 (founder lines) until the last generation (G0, G1, G2...), offspring, father and mother IDs, sex and phenotype. These two last columns are not required. haploMAGIC populations are required to have (1) complete pedigree (both parents known for every individual except founders), (2) complete genotypes, (3) homozygous founder lines. If your population does not meet this criteria, consider subsetting into subfamilies and filtering homozygous loci.
* One MAP file with the SNP marker information. The naming format must be *[population]_[chromosome].map*. MAP files show SNPs by rows with four tab-separated columns: chromosome, SNP marker ID, genetic position (cM) and physical position in the chromosome (bp). Genetic positions are not required and can be filled with 0.

#### The standard haploMAGIC command prompt looks like this:
```
Rscript 0.haplomagic.R <population(s)> <chromosome(s)> <min> <imp> <cor> <thr>
```
#### Arguments explained:
* 1. Population(s). The IDs of one or a list of populations to analyze. The PED files of these populations must be present for the chromosomes provided. If multiple, write IDs between "" and split by space.
* 2. Chromosome(s). The IDs of one or a list of chromosomes to analyze. The PED files of these chromosomes must be present for the populations provided. If multiple, write IDS between "" and split by space.
* 3. Min threshold, *min* (any integer >0). Minimum number of informative alleles per haploblock. Haplotype origins of the alleles within the haploblock <min are imputed, thus they will not be contribute to recombination events. *min* is a filtering method that increases precision, with *min=1* being equivalent to no filtering and higher *min* values increasing stringency.
* 4. Phase imputation method, *imp* (*imputeAll*/*imputeTHonly*/*imputeNot*). At phasing, some loci cannot be resolved if for them all three trio members are heterozygous (TH), any of them has missing data (MD) or follow incorrect Mendelian inheritance patterns, i.e., Mendelian errors (ME). Imputing missing phases increases recall, but might reduce precision. Users can choose between these options:
	* *imputeTHonly*: Only impute the phases of TH loci.
	* *imputeAll*: Impute the phase of all unresolved loci (TH, MD & ME).
	* *imputeNot*: Do not impute any unresolved phase.
* 5. Post-imputation phase correction method, *cor* (*correctAll*/*correctFalseHom*/*reImpute*/*correctNot*):
	* *correctFalseHom*: The phase from triply heterozygous loci that were incorrectly imputed as homozygous remain unresolved. This method increases precision.
	* *reImpute*: The phase of unresolved alleles are imputed if the phase of the homologous alelle is known. This method increases recall, but it is not recommended alone. Instead, use *correctAll*.
	* *correctAll*: *correctFalseHom*+*reImpute*.
	* *correctNot*: no correction applied. Default when *imp=imputeNot*.
* 6. Base pair threshold, *thr* (any integer >0). For classifying recombination events as gene conversions (<) or crossovers (>) based on the length of the flanking haploblocks.

#### When running the standard prompt with lists of files, each PED/MAP pair is analyzed in series. For parallelizing haploMAGIC runs, we recommend the following method:
* 1. Create a list (INPUT) where each line represents a different haploMAGIC input*
```
pop1 chr1 min imp cor thr #Run1
pop2 chr2 min imp cor thr #Run2
pop3 chr3 min imp cor thr #Run3
```
* 2. Use xargs to run haploMAGIC on each line simultaneously.
```
cat INPUT | xargs -L1 -P3 Rscript 0.haplomagic.R
```
-P for adjusting the number of cores available. -L1 should not be changed.

#### Output explained:
All output files are space-separated. Here is a description of the files users might find most important:
| File extension	| Description |
| --------------	| ----------- |
| *.phase*	| Sequence of phased alleles (0, 1, 2) on same chromosome with row names indicating *[individual]*_*[parental phase]*, where P and M stand for paternal or maternal phase respectively. |
| *.f.i.origin*	| Sequence of filtered and imputed haplotype origins per allele (P, M, *, ?, !) assigned after phasing by comparing parental and offspring phases by trios. Symbol meanings explained in Montero-Tena et al. 2024. Row names like in *.phase* files. |
| *.haplo*	| Sequence of labels indicating the founder line annotated to each allele. Generated by sequentially assigning founders by generation using the origins assigned and pedigree information. Non-assigned alleles expressed with *. Row names like in *.phase* files. |
| ***.f.reco***	| Information on the recombination events detected in the population. Recombination events split by rows. Columns for the population, chromosome, generation and individual where the recombination was detected, as well as the father and mother. Also, whether it was detected in maternal or paternal meiosis, the transition type (M-->P/P-->M),  the founder haploblock transition, start and end coordinates in the SNP array as well as with genomic position, and the event classification based on *thr* (CO/GC) |
| *.stats*	| Statistics per individual per meiosis about the recombination events detected by haploMAGIC. The first 6 columns are shared with *.f.reco*. It shows the number of recombination events, crossovers and gene conversions detected before and after filtering and provides the rates of filtering for each type of event, as well as the rate of non-informative alleles, missing data and Mendelian error. |
| files with *.nf.*	| Show non-filtered events or origins. |
| files with *.ni.*	| Show non-imputed events or origins. |

