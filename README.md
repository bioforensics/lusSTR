# lusSTR

lusSTR is a tool written in Python to convert NGS sequence data of forensic STR loci to different sequence representations (sequence bracketed form) and allele designations (CE allele, LUS/LUS+ alleles) for ease in downstream analyses. See the below section ```Converting STR sequences to other sequence representations and allele designations``` for more information. 

Further, lusSTR can perform filtering and stutter identification using the CE allele or the bracketed sequence form for autosomal loci and create files for direct input into two probabilistic genotyping software packages, EuroForMix (EFM) and STRmix (both CE and NGS). 

lusSTR also processes SNP data from the Verogen ForenSeq and Kintelligence panels and create evidence and/or reference files for use in EFM. See the below section ```SNP Data Processing``` for more information.

This Python package has been written for use with either:  
* ForenSeq Signature Prep panel
  * 27 autosomal STR loci
  * 24 Y-chromosome STR loci
  * 7 X-chromosome STR loci
  * 94 identity SNPs
  * 22 phenotype (hair/eye color) SNPs
  * 54 ancestry SNPs
  * 2 phenotype and ancestry SNPs
* ForenSeq Kintelligence panel
  * 10,230 SNPs for forensic genetic genealogy purposes
* Promega PowerSeq panel
  * 22 autosomal STR loci
  * 22 Y-chromosome loci
  
The package accomodates either the Sample Details Report/Phenotype Report/Sample Report from the ForenSeq Universal Analysis Software (UAS) or STRait Razor output. If STRait Razor output is provided, sequences are filtered to the UAS sequence region for conversion.


## Installation

For best results, install from bioconda.
```
conda install -c bioconda lusstr
```
Contributors/developers who want to create a dedicated environment on their machine can do so:
```
conda create --name lusSTR -y python=3.7 pandas
git clone https://www.github.com/bioforensics/lusSTR.git
cd lusSTR
make devenv
```

## Usage

lusSTR accomodates three different input formats:  
(1) UAS Sample Details Report, UAS Sample Report, and UAS Phenotype Report (for SNP processing) in .xlsx format (a single file or directory containing multiple files)  
(2) STRait Razor output with one sample per file (a single file or directory containing multiple files)  
(3) Sample(s) sequences in CSV format; first four columns must be Locus, NumReads, Sequence, SampleID; Optional last two columns can be Project and Analysis IDs.  

*These individual sample files or directory of files must be specified in the config file (see below).*


lusSTR utilizes the ```lusstr``` command to invoke various Snakemake workflows. The ```lusstr strs``` command invokes the STR analysis workflow while the ```lusstr snps``` command invokes the SNP analysis workflow. Please see below for further information on processing SNP data.
___
### Creating the STR config file

Running ```lusstr config``` creates a config file containing the default settings for the lusSTR STR analysis pipeline. The settings can be changed with command line arguments (see below) or by manually editing the config file. The default settings, along with their descriptions, are as follows:

### general settings
uas: ```True``` (True/False); if ran through UAS (invoke ```--straitrazor``` flag if STRait Razor was used)  
sex: ```False``` (True/False); include sex-chromosome STRs (invoke ```--sex``` flag)  
samp_input: ```/path/to/input/directory/or/samples``` input directory or sample; if not provided, will be current working directory (indicate using ```--input path/to/dir``` )  
output: ```lusstr_output``` output file/directory name (indicate using ```--out dir/sampleid e.g. --out test_030923```)

### convert settings  
kit: ```forenseq``` (forenseq/powerseq) (invoke the ```--powerseq``` flag if using PowerSeq data)  
nocombine: ```False``` (True/False); do not combine identical sequences during the ```convert``` step, if using STRait Razor data. (invoke the ```--nocombine``` flag)  

### filter settings  
output_type: ```strmix``` (strmix/efm) (invoke ```--efm``` flag if creating output for EuroForMix)  
profile_type: ```evidence``` (evidence/reference) (invoke ```--reference``` flag if creating a reference output file)  
data_type: ```ngs``` (ce/ngs) (invoke ```--ce``` if using CE allele data)  
info: ```True``` (True/False); create allele information file (invoke ```--noinfo``` flag to not create the allele information file)  
separate: ```False``` (True/False); for EFM only, if True will create individual files for samples; if False, will create one file with all samples (invoke ```--separate``` flag to separate EFM output files)  
nofilters: ```False``` (True/False); skip all filtering steps but still creates EFM/STRmix output files (invoke ```--nofilters``` flag)  
strand: ```uas``` (uas/forward); indicates the strand orientation in which to report the sequence in the final output table for STRmix NGS only (indicate using ```--strand```)

One additional argument can be provided with ```lusstr config```:  
```-w```/```-workdir``` sets the working directory (e.g. ```-w lusstr_files/```) and all created files are stored in that directory.

**Once the config file is created with all the desired settings, the STR workflow can be run. The config file must be located in the working directory.**
___
## Running the lusSTR STR workflow

The lusSTR STR workflow consists of three steps:  
(1) ```format```: formatting input  
(2) ```convert```: converting sequences to other sequence representations and allele designations  
(3) ```filter```: performing several filtering steps and creating appropriately formatted files for use in EuroForMix or STRmix.

Any or all steps can be run. In order to run all three steps, the following command can be used:  
```
lusstr strs all
```

One additional argument can be provided, a working directory.  
**This working directory must contain the config file.**  
If not specified, the working directory is the current directory.
```
lusstr strs all -w lusstr_files/
```

Individual steps can also be run
```
lusstr strs format -w lusstr_files/
```

```
lusstr strs convert -w lusstr_files/
```

**In order to run the ```convert``` step, the appropriately formatted ```.csv``` file containing the sequences normally created in the ```format``` step must be present in the working directory. See the below ```Formatting input for STR loci sequences``` section for specific information about that file (required columns, etc.).** 

----

## Additional information about each step  


### Formatting input for STR loci sequences

If inputting data from either the UAS Sample Details Report or STRait Razor output, the user must first invoke the ```format``` step to extract necessary information and format for the ```convert``` step.

The ```format``` command removes unnecessary rows/columns and outputs a table in CSV format containing the following columns:
*  Locus
*  Number of Reads observed with the specified sequence
*  Sequence
*  Sample ID
*  Project ID (if provided)
*  Analysis ID (if provided)

If including the sex chromosome loci as specified in the config file, the ```format``` command will output a second table for the sex loci with the same columns (```*_sexloci.csv```).

---

### Converting STR sequences to other sequence representations and allele designations

The ```convert``` step produces a tab-delineated table with the following columns:
*  Sample ID
*  Project ID (if provided)
*  Analysis ID (if provided)
*  Locus
*  UAS Output sequence: can be forward or reverse strand
*  Forward strand sequence: will be same as UAS Output sequence for those loci reported on forward strand
*  UAS Output Bracketed notation: Bracketed sequence form for the reported UAS sequence output (will be same for those loci which report the forward strand)
*  Forward Strand Bracketed notation: Bracketed notation for forward strand sequence
*  CE allele: common length-based CE allele designation (also called the repeat unit, or RU, allele)
*  LUS: Longest uninterrupted stretch
*  LUS+: Notation combining multiple allele designations including CE, LUS, secondary motif (if applicable) and tertiary motif (if applicable)
*  Reads: number of reads observed with the specified sequence

If including the sex chromosome loci as specified in the config file, a second table with the above columns for the sex chromosome loci will be outputted as well.


If STRait Razor data is specified, several additional processes occur with the ```convert``` step:
*  The full sequences are filtered to the UAS region before the translation step. The number of bases to remove is determined based on the specified kit.
*  Once the sequences are filtered to the UAS region, any duplicated sequences are removed and their reads are summed in with the remaining sequence ```Reads``` column. NOTE: This step can be skipped with the ```nocombine``` setting in the config file.
  
Further, a second table (labeled as ```*_flanks.txt```) containing information related to the flanking sequences surrounding the UAS sequence region is also produced with the following columns:
*  Sample ID
*  Project ID
*  Analysis ID (same as Project ID)
*  Locus
*  Reads: number of reads observed for the specified sequence
*  Length-based (CE) Allele
*  Full Sequence
*  5' Flanking Sequence Bracketed Notation
*  UAS Region Sequence Bracketed Notation (same as column ```UAS Output Bracketed Notation``` in the main table)
*  3' Flanking Sequence Bracketed Notation
*  Potential Issues (such as: Possible indel or partial sequence)

The ```Potential Issues``` column in this report is to draw attention to potential problem sequences (due to perhaps an indel or partial sequence) and requires the attention of the user to further evaluate the sequence for it's authenticity.
 
---

### Filtering Alleles/Sequences and Creation of Files for Use in ProbGen Software

The ```filter``` step provides the opportunity to filter sequences using thresholds such as:
* Detection threshold (both static and dynamic)
* Analytical threshold (both static and dynamic)
* Same size threshold (dynamic)

Custom static and dynamic thresholds for each locus are stored in the ```filters.json``` file. This file should be updated to utilize validated thresholds for individual labs.

In addition, stutter alleles can be identified using the ```info``` setting in the config file. This creates a separate file containing information about each allele, including an allele classification (```real allele```, ```stutter``` or ```BelowAT```). Stutter alleles are classified as either ```-1 stutter```, ```-2 stutter```, or ```+1 stutter```. For these stutter alleles, the stuttering allele is reported along with the percent stutter (# of reads for that allele/# of reads for stuttering allele). In instances where a stutter allele could be multiple different types of stutter, all potential designations will be reported as such: ```-1 stutter/-2 stutter```, ```-1 stutter/+1 stutter```, or ```-2 stutter/+1 stutter```. No percent stutter is calculated for these alleles. If a sequence is identified as noise, the percent noise is calculated (# of reads for that sequence/total locus reads).

Each locus is checked for containing greater than 2 alleles (indicating a potential mixture) and for intralocus imbalance. If either are identified, a separate file (```Flagged_Loci.csv```) is created, containing the SampleID, Locus and either ```>2Alleles``` or ```IntraLocusImbalance```.

When using STRmix data, the data type can be specified using the ```data-type``` setting as either ```ce``` or ```ngs``` (default is ```ngs```). If ```ngs``` is specified, the same size filter is applied following the stutter filter. Further, the columns and column names in the output file differ based on the data type.

Finally, output files are created for direct use in EuroForMix (EFM) or STRmix. If EFM is specified, a single file is created containing all samples in the input file (however, separate output files for each sample can be created with the ```separate``` setting specified in the config file). If STRmix is specified, a directory containing files for each individual sample is created. The ```profile-type``` setting allows for the creation of either a ```reference``` or ```evidence``` profile. Both EuroForMix and STRmix require different formatting depending on the type of sample. 

___

 ## SNP Data Processing

lusSTR is able to process SNPs derived from the ForenSeq Signature Prep assay and the ForenSeq Kintelligence assay. SNPs from the ForenSeq Signature Prep assay could be analyzed using either the Verogen UAS or STRait Razor. SNPs from the ForenSeq Kintelligence assay must first be analyzed using the UAS.

___
### Creating the SNP config file

Running ```lusstr config --snps``` creates a config file containing the default settings for the lusSTR SNP workflow. The settings can be changed with command line arguments (see below) or by manually editing the config file. The default settings, along with their descriptions, are as follows:  


### general settings  
uas: ```True```  (True/False); if ran through UAS  (invoke ```--straitrazor``` flag if STRait Razor was used)  
samp_input: ```/path/to/input/directory/or/samples``` input directory or sample; if not provided, will be current working directory (indicate using ```--input path/to/dir```)  
output: ```lusstr_output``` output file/directory name; (indicate using ```--out dir/sampleid e.g. --out test_030923```)   
kit: ```sigprep``` (sigprep/kintelligence) (invoke using the ```--kintelligence``` flag if using Kintelligence data)  

### format settings  
types: ```all``` choices are "all", "i" (identity SNPs only), "p" (phenotype only), "a" (ancestry only) or any combination (indicate using the ```--snp-type e.g. --snp-type i, p```)  
nofilter: ```False``` (True/False); if no filtering is desired at the format step; if False, will remove any allele designated as Not Typed (invoke using the ```--nofiltering```)  

### convert settings  
strand: ```forward``` (forward/uas); indicates which orientation to report the alleles for the SigPrep SNPs; uas indicates the orientation as reported by the UAS or the forward strand
references: ```None```; list IDs of the samples to be run as references in EFM; default is no reference samples  
separate: ```False``` (True/False); if want to separate samples into individual files for use in EFM  
thresh: ```0.03```; Analytical threshold value    


One additional argument can be provided with ```lusstr config --snps```:  
```-w```/```-workdir``` sets the working directory (e.g. ```-w lusstr_files/```) and all created files are stored in that directory.

**Once the config file is created with all the desired settings, the SNP workflow can be run. The config file must be located in the working directory.**
___
## Running the lusSTR SNP workflow

The lusSTR SNP workflow consists of three steps:  
(1) ```format```: formatting input and calling alleles if using STRait Razor data  
(2) ```convert```: applying analytical threshold; converting data to correct format for input into EuroForMix;  

Any or all steps can be run. In order to run all three steps, the following command can be used:  
```
lusstr snps all
```

One additional argument can be provided, a working directory.  
**This working directory must contain the config file.**  
The default working directory is the current directory.  
```
lusstr snps all -w lusstr_files/  
```

Individual steps can also be run  
```
lusstr snps format -w lusstr_files/  
```

```
lusstr snps convert -w lusstr_files/  
```

**In order to run the ```convert``` step, the appropriately formatted ```.csv``` file containing the sequences normally created in the ```format``` step must be present in the working directory. See the below ```Usage``` section for specific information about that file (required columns, etc.).**  

----

## Additional information about each step  


### Formatting input for SNP data  

If inputting data from either the UAS Sample Details Report/Phenotype Report/Sample Report or STRait Razor output, the user must first invoke the ```format``` step to extract necessary information and format for the ```convert``` step.  

The ```format``` command removes unnecessary rows/columns and outputs a table in CSV format containing the following columns:  
*  Sample ID  
*  Project ID  
*  Analysis ID  
*  SNP (rsID)  
*  Reads  
*  Forward Strand Allele  
*  UAS orientation Allele  
*  Type (ancestry/identity/phenotype/kintelligence)  
*  Issues  

### Converting to appropriately formatted files for use in EuroForMix  

This step will convert the table generated in the ```format``` step into the correct format for use in EuroForMix. An analytical threshold can be applied (this is especially useful for data analyzed using STRait Razor) in this step.  

If any samples are to be used as references, their IDs can be provided in the config file to create a separate file appropriately formatted for use as reference profiles in EFM. Any samples not specified as references are assumed to be evidence samples and will be formatted as such.  

There is the option to create separate evidence files for each sample (as specified in the config file); this is especially useful for Kintelligence profiles given their larger size.  

This command also changes the alleles to numeric (```A```=```1```, ```C```=```2```, ```G```=```3```, ```T```=```4```)


----


lusSTR is still under development and any suggestions/issues found are welcome!
