rule yahs:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        assembly = os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.fasta"),
        bam = os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.purged.sorted.bam")
    output:
        assembly = os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.fa"),
        bin = os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.bin"),
        agp = os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.agp")
    params:
        prefix = os.path.join(out_dir, "assembly", "hic_hifi", "{type}")
    threads:
        threads
    shell:
        """
        yahs \
            --no-contig-ec \
            --no-mem-check \
            -o {params.prefix} \
            {input.assembly} {input.bam}
        """

rule juicer_pre:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        bin = os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.bin"),
        agp = os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.agp"),
        fai = os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.fasta.fai")
    output:
        os.path.join(out_dir, "assembly", "hic_contact_map", "{type}.alignments_sorted.txt")
    threads: config["threads"]
    shell:
        """
        (juicer pre {input.bin} {input.agp} {input.fai} | \
            sort -k2,2d -k6,6d -T ./ --parallel={threads} -S24G | \
            awk 'NF' \
            > alignments_sorted.txt.part) \
            && (mv alignments_sorted.txt.part {output})
        """

rule scaffolds_final_chrom_sizes:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.fa")
    output:
        os.path.join(out_dir, "assembly", "hic_contact_map", "{type}_scaffolds_final.chrom.sizes")
    shell:
        """
        samtools faidx {input}
        cut -f 1,2 {input}.fai > {output}
        """

rule contact_matrix:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        aln = os.path.join(out_dir, "assembly", "hic_contact_map", "{type}.alignments_sorted.txt"),
        chrom_size = os.path.join(out_dir, "assembly", "hic_contact_map", "{type}_scaffolds_final.chrom.sizes")
    output:
        os.path.join(out_dir, "assembly", "hic_contact_map", "{type}_scaffolds_final.hic")
    params:
        juicer_tools_jar = config["juicer_tools_jar"]
    shell:
        """
        (java -jar -Xmx24G {params.juicer_tools_jar} pre {input.aln} out.hic.part {input.chrom_size}) \
        && (mv out.hic.part {output})
        """
