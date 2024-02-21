![haplomagic2](https://github.com/GoliczGenomeLab/haploMAGIC/assets/134378980/7e5d9334-7760-456e-ab8e-68e86156e7bc)
## Quick setup

#### Download repository and execute on the example data:
```
Rscript 0.haplomagic.R example 1 3 imputeTHonly correctFalseHom 10000
```
#### For more information, read the following manual.

## User's manual

#### For running haploMAGIC, install the scripts in the same directory with the input.

#### The 2 inputs required by haploMAGIC for each run consist in:
* One PED file in 12 format, containing the pedigree info and the SNP genotypes of all the individuals of the population. The name format must be [population]_[chromosome].ped
* One MAP file with the SNP marker information. The name format must be [population]_[chromosome].map

#### The standard haploMAGIC command prompt looks like this:
```
Rscript 0.haplomagic.R <population(s)> <chromosome(s)> <min threshold> <phase imputation method> <post-imputation phase correction method> <CO-GC discrimination threshold>
```
#### Arguments explained:
* 1. Population(s). The IDs of one or a list of populations to analyze. The PED files of these populations must be present for the chromosomes provided. If multiple, write IDs between "" and split by space.
* 2. Chromosome(s). The IDs of one or a list of chromosomes to analyze. The PED files of these chromosomes must be present for the populations provided. If multiple, write IDS between "" and split by space.
* 3. Min threshold. Minimum number of informative alleles that reconstructed haplotype blocks must count with. If this number lies under min, the haplotypes are imputed. Min is a filtering method that increases precision.
* 4. Phase imputation method (imputeAll/imputeTHonly/imputeNot). At phasing, some loci cannot be resolved if for them all three trio members are heterozygous (TH), any of them has missing data (MD) or mendelian transmission incompabilities are observed, i.e., mendelian error (ME). Imputing missing phases increases recall, but might reduce precision. Users can choose between these options:
	* imputeTHonly: Only impute the phases of TH loci.
	* imputeAll: Impute the phase of all unresolved loci (TH, MD & ME).
	* imputeNot: Do not impute any unresolved phase.
* 5. Phase correction method (correctAll/correctFalseHom/reImpute/correctNot).
	* correctFalseHom: The phase from triply heterozygous loci that were incorrectly imputed as homozygous remain unresolved. This method increases precision.
	* reImpute: The phase of unresolved alleles are imputed if the phase of the homologous alelle is known. This method increases recall, but it is not recommended alone. Instead, use 'correctAll'.
	* correctAll: correctFalseHom+reImpute.
	* correctNot: no correction applied. Default when imp=imputeNot.
* 6. CO-GC discrimination threshold. Number for the basepair threshold to discriminate haploblocks as gene conversions (<) or cross-overs (>) based on the haploblock size.

#### When running the standard prompt with lists of files, each PED/MAP pair is analyzed in series. For parallelizing haploMAGIC runs, we recommend the following method:
* 1. Create a list (INPUT) where each line represents a different haploMAGIC input
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
