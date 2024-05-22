rule index_assembly:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.fasta")
    output:
        os.path.join(out_dir,"assembly","hifi","hifi_assembly.hic.{type}.fasta.bwt"),
        os.path.join(out_dir,"assembly","hifi","hifi_assembly.hic.{type}.fasta.fai")
    threads:
        1
    log:
        os.path.join(out_dir, "log", "{type}.index.log")
    shell:
        """
        bwa index {input} 2>{log}
        samtools faidx {input} 2>{log}
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
    log:
        os.path.join(out_dir, "log", "{type}_{hic_reads_prefix}.mapping.log")
    threads:
        threads
    shell:
        """
        bwa mem \
            -t {threads}\
            {input.assembly} \
            {input.reads} 2>{log} | \
            samtools view -@ {threads} -Sb - \
            > {output} 2>{log}
        """

rule filter5end:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping","{type}_{hic_reads_prefix}.mapped.bam")
    output:
        os.path.join(out_dir,"hic_mapping", "{type}_{hic_reads_prefix}.mapped.5endFiltered.bam")
    threads:
        5
    log:
        os.path.join(out_dir, "log", "{type}_{hic_reads_prefix}.filter5end.log")
    shell:
        """
        samtools view -h -@ {threads} {input} | \
            filter_five_end.pl | \
            samtools view -Sb -@ {threads} - \
            > {output} 2>{log}
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
        5
    params:
        mapq_filter=config["mapq_filter"]
    log:
        os.path.join(out_dir, "log", "{type}.conbine_and_filter_bams.log")
    shell:
        """
        two_read_bam_combiner.pl {input.bams} samtools {params.mapq_filter} | \
            samtools view -bS -@ {threads} -t {input.fai} - | \
            samtools sort -@ {threads} -o {output} \
            2>{log}
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
    threads:
        2
    log:
        os.path.join(out_dir, "log", "{type}.mark_duplicate.log")
    shell:
        """
        picard MarkDuplicates \
            -Xmx{params.mem} -XX:-UseGCOverheadLimit \
            INPUT={input} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metric} \
            ASSUME_SORTED=TRUE \
            VALIDATION_STRINGENCY=LENIENT\
            REMOVE_DUPLICATES=TRUE \
            2>{log}
        """

rule sort_by_name_bam:
    conda:
        os.path.join(workflow.basedir,"envs/arima_hic_mapping.yaml")
    input:
        os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.purged.bam")
    output:
         os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.purged.sorted.bam")
    threads:
        min(threads, 16)
    log:
        os.path.join(out_dir, "log", "{type}.sort_by_name_bam.log")
    shell:
        """
        samtools sort -@ {threads} -o {output} -n {input} 2>{log}
        """
