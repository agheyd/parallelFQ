from os import path
from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

threads = snakemake.params.ntasks
memory = snakemake.params.mem
accession = snakemake.params.accession
fastq_dir = snakemake.params.fastq_dir
tmp_dir = path.join(fastq_dir, "tmp")
makedirs(tmp_dir)

extra = ""

shell("parallel-fastq-dump "
"-s {accession} "
"-t {threads} "
"-O {fastq_dir} "
"--tmpdir {tmp_dir} "
"--split-3 "
"-I "
"--gzip "
"{log}")
