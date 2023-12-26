import lusSTR


configfile: "snp_config.yaml"
output_name = config["output"]
input_name = config["samp_input"]
refs = config["references"]


def format_filename(id, refs):
    if refs is None:
        return f"{id}_snp_evidence"
    else:
        return f"{id}_snp_reference"


rule all:
    input:
        expand("{name}.txt", name=output_name),
        expand("{name}.csv", name=format_filename(output_name, refs))


rule format:
    input:
       expand("{samp_input}", samp_input=input_name)
    output:
        expand("{name}.txt", name=output_name)
    params:
       software=config["analysis_software"],
       kit=config["kit"],
       types=config["types"],
       nofilter=config["nofilter"]
    script:
        lusSTR.wrapper("snps_format")


rule convert:
    input:
        rules.format.output
    output:
        expand("{name}.csv", name=format_filename(output_name, refs))
    params:
        strand=config["strand"],
        separate=config["separate"],
        kit=config["kit"],
        refs=refs,
        outputid=output_name,
        software=config["analysis_software"],
        thresh=config["thresh"]
    script:
        lusSTR.wrapper("snps_convert")
