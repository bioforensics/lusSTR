lusSTR is a tool written to convert NGS sequence data to different annotation types for forensic STR loci.

These python scripts have been written for use with the 27 autosomal STR loci from the ForenSeq panel and the sequence range output by the ForenSeq Universal Analysis Software (UAS).

The format_UAS_output.R script requires the Sample Details Report output directly from the UAS software. The R script removes unnecessary rows/columns and outputs a table in CSV format containing the following columns:
*  Locus
*  Number of Reads observed with the specified sequence
*  Sequence
*  Sample ID
*  Project ID (if provided)
*  Analysis ID (if provided)

The STR_annotation.py script is run on the output from the above R script (or any provided .csv files in the correct format) and currently outputs a table with the following columns:
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
 

The provided shell script (submit_for_annotation.sh) is to be run in the directory containing either .xlsx files (direct from the UAS software) or the provided .csv files.

User should copy the scripts to a directory in their `$PATH`.

lusSTR is still under development and any suggestions/issues found are welcome!
