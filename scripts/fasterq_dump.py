from os import path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

threads = snakemake.params.ntasks
memory = snakemake.params.mem
accession = snakemake.params.accession
output_dir = snakemake.output.output_dir
tmp_dir = path.join(output_dir, "tmp")

extra = ""

shell("fasterq-dump "
"{accession} "
"-O {output_dir} "
"-m {threads} "
"-t {tmp_dir} "
"-e {memory} "
"{log}")
