#!/usr/bin/env python

import os

# Get config file
configfile: "config.yaml"

# Dependencies
os.environ['PATH'] += ':' + os.path.abspath("./scripts")
os.environ['PATH'] += ':' + os.path.abspath("./bin")

## Make all scripts executable

def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)

for filename in os.listdir("./scripts"):
    f = os.path.join("./scripts", filename)
    # checking if it is a file
    if os.path.isfile(f):
        make_executable(filename)

# Get parameters
hifi_prefix = os.path.basename(config["hifi_reads"]).replace(".bam", "")
hic_R1_prefix = os.path.basename(config["hiC_read1"]).replace(".gz", "").replace(".fastq", "").replace(".fq", "")
hic_R2_prefix = os.path.basename(config["hiC_read2"]).replace(".gz", "").replace(".fastq", "").replace(".fq", "")

# Opts
threads = config["threads"]
out_dir = config["outdir"]

# Wildcards
types = config["types"]
hic_reads_prefixs = [hic_R1_prefix, hic_R2_prefix]
wildcard_constraints:
    type = "|".join(types),
    hic_reads_prefix = "|".join(hic_reads_prefixs)


# Include rule files
include: "rules/00.reads_preprocessing.smk"
include: "rules/01.hifi_assembly.smk"
include: "rules/02.arima_hic_mapping.smk"
include: "rules/03.yahs.smk"

# Main
rule all:
    input:
        # Reads Preprocessing
        expand(os.path.join(out_dir,"qc","fastqc","{hic_reads_prefix}_fastqc.html"),
            hic_reads_prefix=[hic_R1_prefix, hic_R2_prefix]),
        expand(os.path.join(out_dir,"trimmed_reads","{hic_reads_prefix}.clean.fastq.gz"),
            hic_reads_prefix=[hic_R1_prefix, hic_R2_prefix]),
        os.path.join(out_dir, "qc", "fastqc", "{}_fastqc.html".format(hifi_prefix)),
        os.path.join(out_dir, "qc", "fastqc", "{}.filt_fastqc.html".format(hifi_prefix)),
        # HiFi Assembly
        expand(os.path.join(out_dir,"qc","hifi_stat","hifi_assembly.hic.{type}.stat"), type = types),
        # Hi-C Assembly
        expand(os.path.join(out_dir, "assembly", "hic_contact_map", "{type}_scaffolds_final.hic"), type = types)
