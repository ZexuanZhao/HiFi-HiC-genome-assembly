## Hi-C
rule fastqc_hic_before_trimming:
    conda:
        os.path.join(workflow.basedir,"envs/preprocessing.yaml")
    input:
        config["hiC_read1"],
        config["hiC_read2"]
    threads: 2
    params:
        outdir = os.path.join(out_dir, "qc", "fastqc")
    output:
        os.path.join(out_dir, "qc", "fastqc", "{hic_reads_prefix}_fastqc.html")
    shell:
        """
        fastqc --quiet --outdir {params.outdir} --noextract {input} -t {threads}
        """

rule fastp_hic:
    conda:
        os.path.join(workflow.basedir,"envs/preprocessing.yaml")
    input:
        config["hiC_read1"],
        config["hiC_read2"]
    threads: config["threads"]
    output:
        expand(os.path.join(out_dir, "trimmed_reads", "{hic_reads_prefix}.clean.fastq.gz"), hic_reads_prefix = hic_reads_prefixs)
    shell:
        """
        fastp \
            -i {input[0]} -I {input[1]} \
            -o {output[0]} -O {output[1]} \
            --thread {threads} \
            -j /dev/null -h /dev/null
        """

rule fastqc_hic_after_trimming:
    conda:
        os.path.join(workflow.basedir,"envs/preprocessing.yaml")
    input:
        os.path.join(out_dir,"trimmed_reads","{hic_reads_prefix}.clean.fastq.gz")
    threads: 2
    params:
        outdir=os.path.join(out_dir,"qc","fastqc")
    output:
        os.path.join(out_dir, "qc", "fastqc","{hic_reads_prefix}.clean_fastqc.html")
    shell:
        """
        fastqc --quiet --outdir {params.outdir} --noextract {input} -t {threads}
        """


## HiFi
rule fastqc_hifi_before_trimming:
    conda:
        os.path.join(workflow.basedir,"envs/preprocessing.yaml")
    input:
        config["hifi_reads"]
    output:
        os.path.join(out_dir, "qc", "fastqc", "{}_fastqc.html".format(hifi_prefix))
    threads: 2
    params:
        outdir=os.path.join(out_dir,"qc","fastqc")
    shell:
        """
        fastqc --quiet --outdir {params.outdir} --noextract {input} -t {threads}
        """

rule hifiAdapterFilt:
    conda:
        os.path.join(workflow.basedir,"envs/preprocessing.yaml")
    input:
        bam = config["hifi_reads"],
        db = config["hifi_adapter_db"]
    output:
        os.path.join(out_dir,"trimmed_reads", "{}.filt.fastq.gz".format(hifi_prefix))
    threads:
        config["threads"]
    params:
        mem = config["mem"],
        outdir=os.path.join(out_dir,"trimmed_reads")
    shell:
        """
        hifiadapterfilt_modified.sh \
            -b {input.bam} \
            -d {input.db} \
            -t {threads} \
            -o {params.outdir} \
            -M {params.mem}
        """

rule fastqc_hifi_after_trimming:
    conda:
        os.path.join(workflow.basedir,"envs/preprocessing.yaml")
    input:
        os.path.join(out_dir,"trimmed_reads", "{}.filt.fastq.gz".format(hifi_prefix))
    output:
        os.path.join(out_dir, "qc", "fastqc", "{}.filt_fastqc.html".format(hifi_prefix))
    threads: 2
    params:
        outdir=os.path.join(out_dir,"qc","fastqc")
    shell:
        """
        fastqc --quiet --outdir {params.outdir} --noextract {input} -t {threads}
        """
