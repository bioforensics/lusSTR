# lusSTR

lusSTR is a tool written in Python to convert NGS sequence data of forensic STR loci to different annotation types for ease in downstream analyses.

This Python package has been written for use with either the 27 autosomal STR loci from the Verogen ForenSeq panel or the 22 autosomal STR loci from the Promega PowerSeq panel. The
package accomodates either the Sample Details Report from the ForenSeq Universay Analysis
Software (UAS) or STRait Razor output. If STRait Razor output is provided, sequences are
filtered to the UAS sequence region for annotation.

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
(1) UAS Sample Details Report in .xlsx format
(2) STRait Razor output with one sample per file
(3) Sample(s) sequences in CSV format; first four columns must be Locus, NumReads, Sequence, SampleID; Optional last two columns can be Project and Analysis IDs.

### Formatting input

If inputting data from either the UAS Sample Details Report or STRait Razor output, the user must first invoke the ```format``` command to extract necessary information and format for the ```annotate``` command.

The ```format``` command removes unnecessary rows/columns and outputs a table in CSV format containing the following columns:
*  Locus
*  Number of Reads observed with the specified sequence
*  Sequence
*  Sample ID
*  Project ID (if provided)
*  Analysis ID (if provided)


#### UAS Sample Details Report

If using the UAS Sample Details Report, the user must specify the input file as well an output file and the ```--uas``` flag:
```
lusstr format <input> -o <output> --uas
```

Example:
```
lusstr format UAS_Sample_Details_Report.xlsx -o UAS_test_file.csv --uas
```

#### STRait Razor

If using the output from STRait Razor, the files **must** be labeled as ```SampleID_STRaitRazor.txt``` (example: ```Sample0001_STRaitRazor.txt```) and **must** be compiled in a separate folder (labeled with the project ID). The user must specify the folder name for the ```format``` command as well as an output filename (all sample files will be compiled into one file):
```
lusstr format <input> -o <output>
```

Example:

```
lusstr format STRaitRazorOutputFolder/ -o STRaitRazor_test_file.csv
```


### Annotation

The ```annotate``` command produces a tab-delineated table with the following columns:
*  Sample ID
*  Project ID (if provided)
*  Analysis ID (if provided)
*  Locus
*  UAS Output sequence: can be forward or reverse strand
*  Forward strand sequence: will be same as UAS Output sequence for those loci reported on forward strand
*  Traditional STR allele: common repeat unit based annotation
*  Forward Strand Bracketed annotation: Bracketed annotation for forward strand sequence
*  UAS Output Bracketed annotation: Bracketed annotation for the reported UAS sequence output (will be same for those loci which report the forward strand)
*  LUS: Longest uninterrupted stretch
*  LUS+: annotation combining multiple annotations including traditional STR allele designation, LUS, secondary motif (if applicable) and tertiary motif (if applicable)
*  Reads: number of reads observed with the specified sequence


For the ```annotate``` command, the following must be specified:
*  Input filename
*  Output filename
*  Kit (forenseq or powerseq; default is forenseq)

```--uas``` flag indicates the sequences are only of the UAS region; otherwise, lusSTR assumes full length sequences.

```
lusstr annotate <input> -o <output> --kit forenseq --uas
```
Example:
```
lusstr annotate UAS_test_file.csv -o UAS_final_table.txt --kit forenseq --uas
```

If no ```--uas``` flag is provided, several additional processes occur with the ```annotate``` command:
*  The full sequences are filtered to the UAS region before the annotation step. The number of bases to remove is determined based on the specified kit.
*  Once the sequences are filtered to the UAS region, any duplicated sequences are removed and their reads are summed in with the remaining sequence ```Reads``` column. NOTE: This step can be skipped with the ```--nocombine``` flag.
  
Further, a second table (labeled as ```*_flanks_anno.txt```) containing information related to the flanking sequences surrounding the UAS sequence region is also produced with the following columns:
*  Sample ID
*  Project ID
*  Analysis ID (same as Project ID)
*  Locus
*  Reads: number of reads observed for the specified sequence
*  Length-based Allele
*  Full Sequence
*  5' Flanking Sequence Annotation
*  UAS Region Sequence Annotation (same as column ```UAS Output Bracketed annotation``` in the main table)
*  3' Flanking Sequence Annotation
*  Potential Issues (such as: Possible indel or partial sequence)

The ```Potential Issues``` column in this report is to draw attention to potential problem sequences (due to perhaps an indel or partial sequence) and requires the attention of the user to further evaluate the sequence for it's authenticity.
 
Example:
```
lusstr annotate STRaitRazor_test_file.csv -o STRaitRazor_powerseq_final.txt --kit powerseq
```
The above example would produce two files: (1) ```STRaitRazor_powerseq_final.txt``` and (2) ```STRaitRazor_powerseq_final_flanks_anno.txt```. 

----

lusSTR is still under development and any suggestions/issues found are welcome!
