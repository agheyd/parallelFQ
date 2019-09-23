from os import path
from snakemake.shell import shell
from glob import glob

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fastq_dir = snakemake.params.fastq_dir
fastq_list = sorted(glob(path.join(fastq_dir, "*.gz")))
output_dir = path.dirname(snakemake.output[0])
json_prefix = path.join(output_dir, "fastp.json")
html_prefix = path.join(output_dir, "fastp.html")

extra = ""

## Provide correct input command if single or paired end
if len(fastq_list) == 1:
    read1 = fastq_list[0]
    output1 = path.join(output_dir, "{}.fastq.gz".format(snakemake.wildcards.sample))
    input_cmd = "-i {} -o {}".format(read1, output1)
elif len(fastq_list) == 2:
    read1 = fastq_list[0]
    output1 = path.join(output_dir, "{}_1.fastq.gz".format(snakemake.wildcards.sample))
    read2 = fastq_list[1]
    output2 = path.join(output_dir, "{}_2.fastq.gz".format(snakemake.wildcards.sample))
    input_cmd = "-i {} -o {} -I {} -O {}".format(read1, output1, read2, output2)
    extra += "--detect_adapter_for_pe "
    extra += "--correction " # Base correction in overlapping regions, only for PE


shell("fastp "
""
"-j {json_prefix} "
"-h {html_prefix} "
"-w {snakemake.params.ntasks} "
"--overrepresentation_analysis "
"{extra} "
"{input_cmd} "
"{log}")
