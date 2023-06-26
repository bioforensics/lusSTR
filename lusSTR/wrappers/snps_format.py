import pandas as pd
import os


def kintelligence_filtering(input):
    return input


def create_evidence(evid_df, orientation, separate):
    if orientation == "uas":
        allele_col = "UAS_Allele"
    else:
        allele_col = "Forward_Strand_Allele"
    all_samples_df = pd.DataFrame()
    for sample in evid_df["SampleID"].unique():
        sample_df = evid_df[evid_df["SampleID"] == sample]
        compiled_table = (
            sample_df.groupby([sample_df.groupby(["SampleID", "SNP"]).cumcount() + 1, "SNP"])
            .first()[[allele_col, "Reads"]]
            .unstack(0)
            .reset_index()
        )
        compiled_table.columns = ["Marker", "Allele 1", "Allele 2", "Height 1", "Height 2"]
        compiled_table.insert(0, "Sample Name", sample)
        all_samples_df = all_samples_df.append(compiled_table)
        if separate:
            compiled_table.to_csv(
                f"evidence_samples/{sample}_snp_evidence.csv", index=False, sep="\t"
            )
    return all_samples_df


def main(input, output, kit, strand, separate, refs):
    output = str(output)
    input = str(input)
    output_name = os.path.splitext(output)[0]
    input_file = pd.read_csv(input, sep="\t")
    if kit == "kintelligence":
        results = kintelligence_filtering(input)
    else:
        results = input_file
    if refs is None:
        evidence_table = create_evidence(results)
    else:
        ref_samples = results[results["SampleID"].isin([refs])]
        if len(ref_samples) > 0:
            ref_table = create_ref(ref_samples)
            ref_table.to_csv(f"{output_name}_snp_references.csv", index=False, sep="\t")
        evid_samples = results[~results["SampleID"].isin([refs])]
        if len(evid_samples) > 0:
            evid_table = create_evidence(evid_samples, strand, separate)
            evid_table.to_csv(output, index=False, sep="\t")


if __name__ == "__main__":
    main(
        snakemake.input,
        snakemake.output,
        kit=snakemake.params.kit,
        strand=snakemake.params.strand,
        separate=snakemake.params.separate,
        refs=snakemake.params.refs,
    )
