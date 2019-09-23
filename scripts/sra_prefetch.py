from os import path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

accession = snakemake.params.accession
output_dir = snakemake.output.output_dir
extra = ""

shell("prefetch "
"{accession} "
"-O {output_dir} "
"{log}")
