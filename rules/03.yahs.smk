rule yahs:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        assembly = os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.fasta"),
        bam = os.path.join(out_dir,"hic_mapping", "{type}.combined.filtered.purged.sorted.bam")
    output:
        assembly = os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.fa"),
        bin = os.path.join(out_dir, "assembly", "hic_hifi", "{type}.bin"),
        agp = os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.agp")
    params:
        prefix = os.path.join(out_dir, "assembly", "hic_hifi", "{type}")
    log:
        os.path.join(out_dir, "log", "{type}.yahs.log")
    threads:
        threads
    shell:
        """
        yahs \
            --no-contig-ec \
            --no-mem-check \
            -o {params.prefix} \
            {input.assembly} {input.bam} \
            2>{log}
        """

rule juicer_pre:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        bin = os.path.join(out_dir, "assembly", "hic_hifi", "{type}.bin"),
        agp = os.path.join(out_dir, "assembly", "hic_hifi", "{type}_scaffolds_final.agp"),
        fai = os.path.join(out_dir, "assembly", "hifi", "hifi_assembly.hic.{type}.fasta.fai")
    output:
        txt = os.path.join(out_dir, "assembly", "hic_contact_map", "{type}.txt")
    log:
        os.path.join(out_dir, "log", "{type}.juicer_pre.log")
    params:
        prefix = os.path.join(out_dir, "assembly", "hic_contact_map", "{type}")
    threads:
        threads
    shell:
        """
        juicer pre -a -o {params.prefix} {input.bin} {input.agp} {input.fai} >{log} 2>&1
        """

rule contact_matrix:
    conda:
        os.path.join(workflow.basedir,"envs/yahs.yaml")
    input:
        txt = os.path.join(out_dir, "assembly", "hic_contact_map", "{type}.txt"),
        log = os.path.join(out_dir, "log", "{type}.juicer_pre.log")
    output:
        os.path.join(out_dir, "assembly", "hic_contact_map", "{type}.hic")
    log:
        os.path.join(out_dir, "log", "{type}.contact_matrix.log")
    params:
        temp_file = os.path.join(out_dir, "assembly", "hic_contact_map", "{type}.hic.part"),
        mem = config["mem"],
        juicer_tools_jar = config["juicer_tools_jar"]
    shell:
        """
        (java -jar -Xmx{params.mem} {params.juicer_tools_jar} pre -r 1000 {input.txt} {params.temp_file} 1>{log} 2>{log} <(cat {input.log} | grep PRE_C_SIZE | awk '{{print $2" "$3}}')) && (mv {params.temp_file} {output})
        """
