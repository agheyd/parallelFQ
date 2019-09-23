import pandas as pd
from os import path

configfile: "config.yaml"

def GetAccession(wildcards):
    entry = SAMPLE_SHEET.loc["{wildcard}".format(wildcard=wildcards)]
    accession = entry.accession
    return accession

## Input variables
SAMPLE_SHEET = pd.read_csv(config["locations"]["sample_sheet"], sep="\t").set_index("sample")
SAMPLES = sorted(SAMPLE_SHEET.index)

## Output variables
OUTPUT_DIR = config["locations"]["output_dir"]
LOG_DIR = path.join(OUTPUT_DIR, "logs")
FASTQ_DIR = path.join(OUTPUT_DIR, "fastq")
FASTP_DIR = path.join(OUTPUT_DIR, "fastp")
FASTQC_DIR = path.join(OUTPUT_DIR, "fastqc")
MULTIQC_DIR = path.join(OUTPUT_DIR, "multiqc")

localrules: all, RenameFQ

rule all:
    input:
        path.join(MULTIQC_DIR, "multiqc_report.html")

rule ParallelFQD:
    output:
        touch(path.join(FASTQ_DIR, "{sample}", "download_done.txt"))
    params:
        fastq_dir = path.join(FASTQ_DIR, "{sample}"),
        accession = GetAccession,
        partition = "main",
        nodes = "1",
        ntasks = "6",
        time = "05:00:00",
        mem = "8G",
        mail_user = "didrio87@zedat.fu-berlin.de",
        mail_type = "end"
    conda:
        "envs/pfastq.yml"
    log:
        path.join(LOG_DIR, "fastq", "{sample}.log")
    script:
        "scripts/parallel_fastq_dump.py"

rule RenameFQ:
    input:
        rules.ParallelFQD.output
    output:
        touch(path.join(FASTQ_DIR, "{sample}", "rename_done.txt"))
    params:
        fastq_dir = path.join(FASTQ_DIR, "{sample}")
    shell:
        """
        arr=({params.fastq_dir}/*.gz)

        if [ ${{#arr[@]}} == 1 ];

        then
        	mv ${{arr[0]}} {params.fastq_dir}/{wildcards.sample}.fq.gz

        elif [ ${{#arr[@]}} == 2 ];

        then
        	mv ${{arr[0]}} {params.fastq_dir}/{wildcards.sample}_1.fq.gz
        	mv	${{arr[1]}} {params.fastq_dir}/{wildcards.sample}_2.fq.gz
        fi
        """

rule Fastp:
    input:
        rules.RenameFQ.output
    output:
        path.join(FASTP_DIR, "{sample}", "fastp.json")
    params:
        fastq_dir = path.join(FASTQ_DIR, "{sample}"),
        partition = "main",
        nodes = "1",
        ntasks = "8",
        time = "01:00:00",
        mem = "8G",
        mail_user = "didrio87@zedat.fu-berlin.de",
        mail_type = "end"
    conda:
        "envs/fastp.yml"
    log:
        path.join(LOG_DIR, "fastp", "{sample}.log")
    script:
        "scripts/fastp_before.py"

# rule FastQC:
#     input:
#         rules.Fastp.output
#     output:
#         touch(path.join(FASTQC_DIR, "{sample}", "fastqc_done.txt"))
#     params:
#         fastq_dir = path.join(FASTP_DIR, "{sample}"),
#         output_dir = path.join(FASTQC_DIR, "{sample}"),
#         partition = "main",
#         nodes = "1",
#         ntasks = "8",
#         time = "01:00:00",
#         mem = "8G",
#         mail_user = "didrio87@zedat.fu-berlin.de",
#         mail_type = "end"
#     conda:
#         "envs/fastqc.yml"
#     log:
#         path.join(LOG_DIR, "fastqc", "{sample}.log")
#     shell:
#         """
#         fastqc {params.fastq_dir}/*.gz -o {params.output_dir} -d {params.output_dir}
#         """

rule GetReadDistributions:
    input:
        rules.Fastp.output
    output:
        touch(path.join(FASTQ_DIR, "{sample}", "distribution_done.txt")),
        touch(path.join(FASTP_DIR, "{sample}", "distribution_done.txt"))
    params:
        fastq_dir = path.join(FASTQ_DIR, "{sample}"),
        fastp_dir = path.join(FASTP_DIR, "{sample}"),
        partition = "main",
        nodes = "1",
        ntasks = "4",
        time = "01:00:00",
        mem = "8G",
        mail_user = "didrio87@zedat.fu-berlin.de",
        mail_type = "end"
    shell:
        """
        ## For each .gz file in the fastq directory
        for f in {params.fastq_dir}/*.gz; do
            ## Name of file
            name=$(echo ${{f##*/}} | cut -d "." -f1)

            ## Output read distribution to file
            zcat $f | \
            awk -v name="$name" 'BEGIN{{OFS="\t"}} NR%4 == 2 {{lengths[length($0)]++; counter++}} END {{for (l in lengths){{print name, l, lengths[l]}}}}' > {params.fastq_dir}/${{name}}.read_distribution.tab
        done

        ## For each .gz file in the fastp directory
        for f in {params.fastp_dir}/*.gz; do
            ## Name of file
            name=$(echo ${{f##*/}} | cut -d "." -f1)

            ## Output read distribution to file
            zcat $f | \
            awk -v name="$name" 'BEGIN{{OFS="\t"}} NR%4 == 2 {{lengths[length($0)]++; counter++}} END {{for (l in lengths){{print name, l, lengths[l]}}}}' > {params.fastp_dir}/${{name}}.read_distribution.tab
        done
        """

rule MultiQC:
    input:
        expand(path.join(FASTQ_DIR, "{sample}", "distribution_done.txt"), sample=SAMPLES)
    output:
        path.join(MULTIQC_DIR, "multiqc_report.html")
    params:
        fastp_dir = FASTP_DIR,
        fastqc_dir = FASTQC_DIR,
        output_dir = MULTIQC_DIR,
        partition = "main",
        nodes = "1",
        ntasks = "4",
        time = "01:00:00",
        mem = "8G",
        mail_user = "didrio87@zedat.fu-berlin.de",
        mail_type = "end"
    conda:
        "envs/multiqc.yml"
    log:
        path.join(LOG_DIR, "multiqc", "multiqc.log")
    shell:
        """
        multiqc {params.fastqc_dir} {params.fastp_dir} -o {params.output_dir}
        """
