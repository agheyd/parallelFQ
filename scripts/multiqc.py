from os import path
from snakemake.shell import shell
from glob import glob

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fastp_dir = snakemake.params.fastp_dir
output_dir = path.dirname(snakemake.output[0])

extra = ""

shell("multiqc {fastp_dir} "
"-o {output_dir} "
"{log}")
