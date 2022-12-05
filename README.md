# lusSTR

lusSTR is a tool written in Python to convert NGS sequence data of forensic STR loci to different sequence representations and allele designations for ease in downstream analyses.

This Python package has been written for use with either: (1) the 27 autosomal STR loci, 24 Y-chromosome STR loci and 7 X-chromosome STR loci from the Verogen ForenSeq panel, or (2) the 22 autosomal STR loci and 22 Y-chromosome loci from the Promega PowerSeq panel. The package accomodates either the Sample Details Report from the ForenSeq Universal Analysis Software (UAS) or STRait Razor output. If STRait Razor output is provided, sequences are filtered to the UAS sequence region for translation.

lusSTR can perform filtering and stutter identification using the RU allele for autosomal loci and create files for direct input into two probabilistic genotyping software packages, including EuroForMix (EFM) and STRmix.

lusSTR also processes SNP data from the Verogen ForenSeq panel. ForenSeq consists of 94 identity SNPs, 22 phenotype (hair/eye color) SNPs, 54 ancestry SNPs and 2 phenotype and ancestry SNPs. Identity SNP data is provided in the UAS Sample Details Report; phenotype and ancestry SNP data is provided in the UAS Phenotype Report. All SNP calls are also reported in the STRait Razor output.


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
(1) UAS Sample Details Report and UAS Phenotype Report (for SNP processing) in .xlsx format
(2) STRait Razor output with one sample per file
(3) Sample(s) sequences in CSV format; first four columns must be Locus, NumReads, Sequence, SampleID; Optional last two columns can be Project and Analysis IDs.

### Formatting input for STR loci sequences

If inputting data from either the UAS Sample Details Report or STRait Razor output, the user must first invoke the ```format``` command to extract necessary information and format for the ```annotate``` command.

The ```format``` command removes unnecessary rows/columns and outputs a table in CSV format containing the following columns:
*  Locus
*  Number of Reads observed with the specified sequence
*  Sequence
*  Sample ID
*  Project ID (if provided)
*  Analysis ID (if provided)

If including the sex chromosome loci (using the ```--include-sex``` flag), the ```format``` command will output a second table for the sex loci with the same columns.


#### **UAS Sample Details Report**

If using the UAS Sample Details Report, the user must specify the input file or folder as well an output file and the ```--uas``` flag:
```
lusstr format <input> -o <output> --uas
```

Example:
```
lusstr format UAS_Sample_Details_Report.xlsx -o UAS_test_file.csv --uas
```
Example using a folder of UAS Sample Details Reports:
```
lusstr format Run01/ -o Run01_compiled_file.csv --uas
```

Including the sex chromosome loci using the ```--include-sex```:

```
lusstr format UAS_Sample_Details_Report.xlsx -o UAS_test_file.csv --uas --include-sex
```

The above command will output two tables which are used in the ```annotate``` command: ```UAS_test_file.csv``` and ```UAS_test_file_sexloci.csv```.  

#### **STRait Razor**

If using lusSTR version 0.4 or above, STRait Razor data **must** be produced using the STRait Razor config file released in January 2021 (ForenSeqv1.25.config and PowerSeqv2.1.config). These config files are available here: https://github.com/Ahhgust/STRaitRazor/tree/103ef68746f010add8f21266fa8bf8fb9f4a076e/.

If using the output from STRait Razor, the files **must** be labeled as ```SampleID.txt``` (example: ```Sample0001.txt```) and can either be specified as a single file or as a folder of multiple STRait Razor output files (folder labeled with the project ID). The user must specify the file or folder name for the ```format``` command as well as an output filename (all sample files will be compiled into one file):
```
lusstr format <input> -o <output>
```

Examples:

```
lusstr format STRaitRazorOutputFolder/ -o STRaitRazor_test_file.csv
```
```
lusstr format A001.txt -o A001.csv
```

Again, sex loci can be included using the ```--include-sex``` flag.
```
lusstr format STRaitRazorOutputFolder/ -o STRaitRazor_test_file.csv --include-sex
```
With this, two tables will be produced: ```STRaitRazor_test_file.csv``` and ```STRaitRazor_test_file_sex_loci.csv```.


### Translation of STR loci sequences

The ```annotate``` command produces a tab-delineated table with the following columns:
*  Sample ID
*  Project ID (if provided)
*  Analysis ID (if provided)
*  Locus
*  UAS Output sequence: can be forward or reverse strand
*  Forward strand sequence: will be same as UAS Output sequence for those loci reported on forward strand
*  RU allele: common length-based repeat unit (RU) allele designation
*  Forward Strand Bracketed notation: Bracketed notation for forward strand sequence
*  UAS Output Bracketed notation: Bracketed annotation for the reported UAS sequence output (will be same for those loci which report the forward strand)
*  LUS: Longest uninterrupted stretch
*  LUS+: Notation combining multiple allele designations including RU, LUS, secondary motif (if applicable) and tertiary motif (if applicable)
*  Reads: number of reads observed with the specified sequence

If the ```--include-sex``` flag is included, a second table with the above columns for the sex chromosome loci will be outputted as well.

**NOTE** on including the sex chromosome STR loci: in the ```annotate``` step, lusSTR requires two files for input: (1) the properly formatted file of autosomal STR loci produced from the ```format``` command (or a file with the appropriate format) with a label such as ```lusSTRinput.csv```, and (2) the properly formatted file of X- and Y-STR loci produced from the ```format``` command with the ```--include-sex``` flag (or a file with the appropriate format) labeled as ```lusSTRinput_sexloci.csv```. The file containing the X- and Y-STR loci *must* have the identical file name to the file containing the autosomal STRs but with ```_sexloci.csv``` (see above for precise examples). These two files are automatically created (and named appropriately) when using the ```--include-sex``` flag with the ```format``` command.

For the ```annotate``` command, the following must be specified:
*  Input filename
*  Output filename
*  Kit (forenseq or powerseq; default is forenseq)

```--uas``` flag indicates the sequences are only of the UAS region; otherwise, lusSTR assumes full length sequences.  
```--include-sex``` flag indicates to include the sex chromosome loci.

```
lusstr annotate <input> -o <output> --kit forenseq --uas --include-sex
```
Example:
```
lusstr annotate UAS_test_file.csv -o UAS_final_table.txt --kit forenseq --uas
```

If no ```--uas``` flag is provided, several additional processes occur with the ```annotate``` command:
*  The full sequences are filtered to the UAS region before the translation step. The number of bases to remove is determined based on the specified kit.
*  Once the sequences are filtered to the UAS region, any duplicated sequences are removed and their reads are summed in with the remaining sequence ```Reads``` column. NOTE: This step can be skipped with the ```--nocombine``` flag.
  
Further, a second table (labeled as ```*_flanks_anno.txt```) containing information related to the flanking sequences surrounding the UAS sequence region is also produced with the following columns:
*  Sample ID
*  Project ID
*  Analysis ID (same as Project ID)
*  Locus
*  Reads: number of reads observed for the specified sequence
*  Length-based Allele
*  Full Sequence
*  5' Flanking Sequence Bracketed Notation
*  UAS Region Sequence Bracketed Notation (same as column ```UAS Output Bracketed Notation``` in the main table)
*  3' Flanking Sequence Bracketed Notation
*  Potential Issues (such as: Possible indel or partial sequence)

The ```Potential Issues``` column in this report is to draw attention to potential problem sequences (due to perhaps an indel or partial sequence) and requires the attention of the user to further evaluate the sequence for it's authenticity.
 
Example:
```
lusstr annotate STRaitRazor_test_file.csv -o STRaitRazor_powerseq_final.txt --kit powerseq
```
The above example would produce two files: (1) ```STRaitRazor_powerseq_final.txt``` and (2) ```STRaitRazor_powerseq_final_flanks_anno.txt```. 

If the ```--include-sex``` flag is included, as below:
```
lusstr annotate STRaitRazor_test_file.csv -o STRaitRazor_powerseq_final.txt --kit powerseq --include-sex
```
 Two additional tables will be produced: (1) ```STRaitRazor_powerseq_final_sexloci.txt``` and (2) ```STRaitRazor_powerseq_final_sexloci_flanks_anno.txt``` for translation of the sex chromosome loci and their flanking regions.

 ## SNP Data Processing

 The ```snp``` command produces tab-delineated table with the following columns:
 * Sample ID
 * Project ID
 * Analysis ID (same as Project ID)
 * SNP (rsID)
 * Reads: number of reads observed for the specified allele
 * Foward Strand Allele: allele call on the forward strand
 * UAS Allele: allele call as reported from the UAS
 * Type: SNP type (identity/phenotype/ancestry)
 * Issues: Indicates if called allele is one of two expected alleles for SNP

If STRait Razor data is used as input, the number of reads for identical alleles within a SNP are combined in the above table. Further, if STRait Razor data is used as input, a second table (```*_full_output.txt```) is produced containing information for each sequence (not combined) with the following columns:
 * Sample ID
 * Project ID
 * Analysis ID
 * SNP
 * Sequence: sequence containing the SNP of interest
 * Reads
 * Forward Strand Allele
 * UAS Allele
 * Type
 * Potential issues: flags sequences which may contains errors, such as an unexpected allele call or short than expected sequence length.

 ### Usage

 ```
 lusstr snps <input_directory> -o <output file name> --type <all, i, p> --uas
 ```

The ```snp``` command requires a folder of either UAS Reports (Sample Details Report(s) and/or Phenotype Report(s)) or STRait Razor output file(s).
The ```-o``` flag specifies the name of the output file (should end in ```.txt```)
The ```--type``` flag specifies the type of SNPs to include in the output file(s). The options are: ```all``` (all SNPs), ```i``` (identity SNPs only), or ```p``` (ancestry and phenotype SNPs only).  The default is ```i```.
Similar to the processing of STR loci sequences, the ```--uas``` flag indicates the input files are Reports from the UAS. Absence of this flag indicates the provided files are STRait Razor output files.

**Examples**:
```
lusstr snps UAS_files/ -o uas_output_all.txt --type all --uas
```
```
lusstr snps STRait_Razor_output/ -o strait_razor_p.txt --type p
```

## Filtering RU alleles and Creation of Files for Use in ProbGen Software

**Currently, lusSTR is only set up to filter and identify stutter based on RU alleles for autosomal loci; future work will expand into the use of the LUS+ allele as well as for sex chromosome loci and SNPs.*

The ```filter``` command provides the opportunity to filter sequences using thresholds such as:
* Detection threshold (both static and dynamic)
* Analytical threshold (both static and dynamic)

Custom static and dynamic thresholds for each locus are stored in the ```filters.json``` file. This file should be updated to utilize validated thresholds for individual labs.

In addition, stutter alleles can be identified using the ```--info``` flag. This creates a separate file containing information about each allele, including an allele classification (```real allele```, ```stutter``` or ```noise```). Stutter alleles are classified as either ```-1 stutter```, ```-2 stutter```, or ```+1 stutter```. For these stutter alleles, the stuttering allele is reported along with the percent stutter (# of reads for that allele/# of reads for stuttering allele). In instances where a stutter allele could be multiple different types of stutter, all potential designations will be reported as such: ```-1 stutter/-2 stutter```, ```-1 stutter/+1 stutter```, or ```-2 stutter/+1 stutter```. No percent stutter is calculated for these alleles. If a sequence is identified as noise, the percent noise is calculated (# of reads for that sequence/total locus reads).

Each locus is checked for containing greater than 2 alleles (indicating a potential mixture) and for intralocus imbalance. If either are identified, a separate file (```Flagged_Loci.csv```) is created, containing the SampleID, Locus and either ```>2Alleles``` or ```IntraLocusImbalance```.

Finally, output files are created for direct use in EuroForMix (EFM) or STRmix. If EFM is specified, a single file is created containing all samples in the input file (however, separate output files for each sample can be created with the ```--separate``` flag). If STRmix is specified, a directory containing files for each individual sample is created. The ```--profile-type``` flag allows for the creation of either a ```reference``` or ```evidence``` profile. Both EuroForMix and STRmix require different formatting depending on the type of sample. 

### Usage
```
lusstr filter <input file> -o <output file/directory> --output-type <efm or strmix> --profile-type <evidence or reference> --info --no-filters --separate
```
The ```filter``` command requires the input of a ```.txt``` file produced by the ```lusstr annot``` command.
The ```-o/--out``` flag specifies the name of the output file (for EFM) or output directory (for STRmix)
```--output-type``` specifies the type of output file created, either ```efm``` or ```strmix```. ```efm``` is the default.
```--profile-type``` specifies the sample type, either ```evidence``` or ```reference```. ```evidence``` is the default.
```--info``` creates the allele information file, containing allele designations (e.g. stutter, noise or real allele) as well as stutter/noise percentages.
The ```--no-filters``` flag will not apply any filters and therefore all alleles present in the input file will be in the created output file(s).
The ```--separate``` flag will indicate to separate samples into individual output files for EFM. STRmix creates separate files by default.

**Examples**:

```
lusstr filter experiment01.txt -o experiment01_efm.csv --output-type efm --info
```

```
lusstr filter experiment01.txt -o STRmix_files/ --output-type strmix --profile-type reference --info
```

----

lusSTR is still under development and any suggestions/issues found are welcome!
