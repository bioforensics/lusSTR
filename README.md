# lusSTR

lusSTR is a tool written in Python to convert NGS sequence data of forensic STR loci to different annotation types for ease in downstream analyses.

This Python package has been written for use with the 27 autosomal STR loci from the ForenSeq panel and the sequence range output by the ForenSeq Universal Analysis Software (UAS).

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

lusSTR accomodates two different input formats:
(1) UAS Sample Details Report in .xlsx format
(2) Sample(s) sequences in CSV format; first four columns must be Locus, NumReads, Sequence, SampleID; Optional last two columns can be Project and Analysis IDs.

If inputting the UAS Sample Details Report, the user must first invoke the ```format``` command to extract necessary information from the UAS Sample Details Report and format for the ```annotate``` command. The user must specify the input file as well an output file:
```
lusstr format <input> -o <output>
```
The ```format``` command removes unnecessary rows/columns and outputs a table in CSV format containing the following columns:
*  Locus
*  Number of Reads observed with the specified sequence
*  Sequence
*  Sample ID
*  Project ID (if provided)
*  Analysis ID (if provided)

This CSV file format is required for the ```annotate``` command:
```
lusstr annotate <input> -o <output>
```
The ```annotate``` command produces a table with the following columns:
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


lusSTR is still under development and any suggestions/issues found are welcome!
