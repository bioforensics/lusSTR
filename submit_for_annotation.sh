#!/usr/bin/env bash

files=*xlsx

for file in $files
do
	file_name=$(echo $file | sed 's/.xlsx//g')
	echo $file_name
	Rscript --vanilla format_UAS_output.R "$file"
	rm "$file_name_tmp.csv"
done

csv_files=*csv

for csv_file in $csv_files
do
	python STR_annotation.py "$csv_file"
done
