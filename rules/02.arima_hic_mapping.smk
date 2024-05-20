rule index_assembly:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.fasta")
    output:
        os.path.join(out_dir,"assembly","hifi","hifi_assembly.hic.{type}.fasta.bwt"),
        os.path.join(out_dir,"assembly","hifi","hifi_assembly.hic.{type}.fasta.fai")
    shell:
        """
        bwa index {input}
        samtools faidx {input}
        """

rule hic_mapping:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        bwt = os.path.join(out_dir,"assembly","hifi","hifi_assembly.hic.{type}.fasta.bwt"),
        assembly = os.path.join(out_dir,"assembly","hifi","hifi_assembly.hic.{type}.fasta"),
        reads = os.path.join(out_dir,"trimmed_reads", "{hic_reads_prefix}.clean.fastq.gz")
    output:
        os.path.join(out_dir,"hic_mapping", "{type}_{hic_reads_prefix}.mapped.bam")
    threads:
        threads
    shell:
        """
        bwa mem \
            -t {threads}\
            {input.assembly} \
            {input.reads}| \
            samtools view -@ {threads} -Sb - \
            > {output}
        """

rule filter5end:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping","{type}_{hic_reads_prefix}.mapped.bam")
    output:
        os.path.join(out_dir,"hic_mapping", "{type}_{hic_reads_prefix}.mapped.5endFiltered.bam")
    threads:
        threads
    shell:
        """
        samtools view -h -@ {threads} {input} | \
            perl filter_five_end.pl | \
            samtools view -Sb -@ {threads} - \
            > {output}
        """

rule conbine_and_filter_bams:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        bams = expand(os.path.join(out_dir,"hic_mapping", "{{type}}_{hic_reads_prefix}.mapped.5endFiltered.bam"),
            hic_reads_prefix=hic_reads_prefixs),
        fai = os.path.join(out_dir,"assembly","hifi","hifi_assembly.hic.{type}.fasta.fai")
    output:
         os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.bam")
    threads:
        threads
    params:
        mapq_filter=config["mapq_filter"]
    shell:
        """
        perl two_read_bam_combiner.pl {input.bams} samtools {params.mapq_filter} | \
            samtools view -bS -@ {threads} -t {input.fai} - | \
            samtools sort -@ {threads} -o {output}
        """

rule mark_duplicate:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.bam")
    output:
        metric = os.path.join(out_dir,"qc", "{type}_markDuplicate.metric.txt"),
        bam = os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.purged.bam")
    params:
        mem = config["mem"]
    shell:
        """
        picard MarkDuplicates \
            -Xmx{params.mem} -XX:-UseGCOverheadLimit \
            INPUT={input} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metric} \
            ASSUME_SORTED=TRUE \
            VALIDATION_STRINGENCY=LENIENT\
            REMOVE_DUPLICATES=TRUE
        """

rule sort_by_name_bam:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.purged.bam")
    output:
         os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.purged.sorted.bam")
    threads:
        threads
    shell:
        """
        samtools sort -@ {threads} -o {output} -n {input}
        """
