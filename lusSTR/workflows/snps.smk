import lusSTR


configfile: "snp_config.yaml"
output_name = config["output"]
input_name = config["samp_input"]
refs = config["references"]

def format_filename(id, refs):
    if refs == "":
        return f"{id}_snp_evidence"
    else:
        return f"{id}_snp_reference"


rule all:
    input:
        expand("{name}.txt", name=output_name),
        expand("{name}.csv", name=format_filename(output_name, refs))


rule convert:
    input:
       expand("{samp_input}", samp_input=input_name)
    output:
        expand("{name}.txt", name=output_name)
    params:
       uas=config["uas"],
       kit=config["kit"],
       types=config["types"],
       nofilter=config["nofilter"]
    script:
        lusSTR.wrapper("snps_convert")


rule format:
    input:
        rules.convert.output
    output:
        expand("{name}.csv", name=format_filename(output_name, refs))
    params:
        strand=config["strand"],
        separate=config["separate"],
        kit=config["kit"],
        refs=refs,
        outputid=output_name
    script:
        lusSTR.wrapper("snps_format")
