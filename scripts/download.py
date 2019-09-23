import pandas as pd
from subprocess import Popen

def StartDL(ftp):

    cmd = ["wget", ftp]

    try:
        p = Popen(cmd)
        p.wait()
    except (KeyboardInterrupt,  SystemExit):
        p.kill()
        raise

    return None

accession = snakemake.input.accession
request = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run".format(accession)
data = pd.read_csv(request, sep="\t")
ftp_list = data["fastq_ftp"].loc[0].split(";")

for ftp in ftp_list:
    StartDL(ftp)
