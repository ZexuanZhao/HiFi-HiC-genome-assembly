rule hifiasm:
    conda:
        os.path.join(workflow.basedir,"envs/hifi_assembly.yaml")
    input:
        expand(os.path.join(out_dir, "trimmed_reads", "{hic_reads_prefix}.clean.fastq.gz"),
            hic_reads_prefix = hic_reads_prefixs),
        hifi = os.path.join(out_dir,"trimmed_reads", "{}.filt.fastq.gz".format(hifi_prefix))
    output:
        expand(os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.gfa"),
            type = types)
    threads:
        threads
    log:
        os.path.join(out_dir, "log", "{}.asm.log".format(hifi_prefix))
    params:
        out_dir = os.path.join(out_dir, "assembly", "hifi")
    shell:
        """
        hifiasm \
	        -o {params.out_dir}/hifi_assembly \
	        -t {threads} \
	        --h1 {input[0]} \
	        --h2 {input[1]} \
	        {input.hifi} \
	        2>{log}
        """

rule gfa2fasta:
    input:
        os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.gfa")
    output:
        os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.fasta")
    threads:
        1
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input} > {output}
        """

rule hifi_assembly_stats:
    input:
        os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.fasta")
    output:
        os.path.join(out_dir, "qc", "hifi_stat", "hifi_assembly.hic.{type}.stat")
    threads:
        1
    shell:
        """
        sequence-stats.sh -a {input} > {output}
        """
