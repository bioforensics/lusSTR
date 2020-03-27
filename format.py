#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import os
import pandas as pd
import re
import sys



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o', '--out', metavar='FILE',
        help='file to which output will be written; default is terminal (stdout)'
    )
    parser.add_argument(
        'input', help='UAS Sample Details Report; .xlsx format'
    )
    args = parser.parse_args()
    file = pd.read_excel(io = args.input, sheet_name = 0)
    well_index = file[file["Sample Autosomal STR Report"] == "Coverage Information"].index.tolist()
    results_newdf = file[(well_index[0] + 2):]
    results_newdf.columns = file.iloc[(well_index[0] + 1)]
    results_filt = results_newdf[results_newdf.Locus != "Amelogenin"]
    results_final = results_filt[['Locus', 'Reads', 'Repeat Sequence']]
    results_final['SampleID'] = file.iloc[1,1]
    results_final['Project'] = file.iloc[2,1]
    results_final['Analysis'] = file.iloc[3,1]
    output_file = sys.stdout
    if args.out is not None:
        results_final.to_csv(args.out, index = False)
            #output_file = open(args.out, "w")
            #output_file.write("Locus", "Reads", "Sequence", "SampleID", "Project", "Analysis")
            #output_file.write(results_final)

if __name__ == "__main__":
    main()