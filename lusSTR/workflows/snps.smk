import lusSTR


configfile: "snp_config.yaml"
output_name = config["output"]
input_name = config["samp_input"]


rule all:
    input:
        expand("{name}.txt", name=output_name),
        expand("{name}_snp_evidence.csv", name=output_name)


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
        expand("{name}_snp_evidence.csv", name=output_name)
    params:
        strand=config["strand"],
        separate=config["separate"],
        kit=config["kit"],
        refs=config["references"]
    script:
        lusSTR.wrapper("snps_format")
