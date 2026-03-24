#this is a germline variant calling pipeline from whole-exome sequencing data using GATK v4.5.0.0 HaplotypeCaller 

workflow.name = "Germline Variant Calling in Whole-Exome Sequencing"
workflow.description = "Whole exome sequencing variant calling and annotation for Multiple Myelmoma"
workflow.author = "Qualia Hooker"
workflow.version = "1.0"

configfile: "config.yaml"

SAMPLES = glob_wildcards(f"{config['bam_dir']}/{{sample}}.bqsr.bam").sample

#rule all define all outputs
rule all:
    input:
        f"{config['results_dir']}/annotation/gatk.cohort.vqsr.annotated.vcf.gz",
        f"{config['results_dir']}/qc/gatk_stats.txt"

#haplotypecaller
rule haplotypecaller:
    input:
        bam = expand(f"{config['bam_dir']}/{{sample}}.bqsr.bam"),
        ref = config['ref']
    output:
        vcf = f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.vcf.gz",
        tbi = f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.vcf.gz.tbi",
        haplo_bam = f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.bam",
        haplo_bai = f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.bai"
    log:
        f"{config['results_dir']}/logs/{{sample}}.haplotypecaller.log"
    resources:
        mem_gb = 80,
        cpus = 16,
        time = "12:00:00"
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0",
        "bbc2/samtools/samtools-1.21"
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \ 
        HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        -L {config['regions']['gatk_intervals']} \
        -O {output.vcf} \
        --bam-output {output.haplo_bam} \
        --create-output-bam-index true \
        --ERC GVCF \
        --dont-use-soft-clipped-bases true \
        --standard-min-confidence-threshold-for-calling 20 \
        --native-pair-hmm-threads {resources.cpus} \
        --sample-ploidy 2 \
        2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """

#joint genotyping
rule genomicsdb:
    input:
        gvcfs = expand(f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.vcf.gz",
        sample=SAMPLES),
        tbis = expand(f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.vcf.gz.tbi",
        sample=SAMPLES)
    output:
        db = directory(f"{config['results_dir']}/jointcalling/genomicsdb")
    log:
        f"{config['results_dir']}/logs/genomicsdb.log"
    resources:
        mem_gb = 150,
        cpus = 16
    params:
        sample_map = f"{config['results_dir']}/sample_map.txt",
        sample_map_lines = lambda wildcards: "\n".join(
            f"{config['results_dir']}/haplotypecaller/{s}/{s}.hc.vcf.gz\t{s}"
            for s in SAMPLES)
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0"
    shell:
        """
        mkdir -p ./tmp
        mkdir -p $(dirname {params.sample_map})

        # create mapping file for samples
        cat > {params.sample_map} << 'EOF'
{params.sample_map_lines}
EOF

        #genomicsdb better for joint genotyping of large cohort per GATK best practices
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        --sample-name-map {params.sample_map} \
        --db-workspace-path {output.db} \
        -L {config['regions']['gatk_intervals']} \
        2> {log}
        """
rule joint_genotyping:
    input:
        db = rules.genomicsdb.output.db,
        ref = config['ref']
    output:
        vcf = f"{config['results_dir']}/jointcalling/cohort.joint.vcf.gz",
        tbi = f"{config['results_dir']}/jointcalling/cohort.joint.vcf.gz.tbi"
    log:
        f"{config['results_dir']}/logs/joint_genotyping.log"
    resources:
        mem_gb = 96,
        cpus = 4
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0",
        "bbc2/samtools/samtools-1.21"
    shell:
        """
         mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {input.ref} \
        -V gendb://{input.db} \
        -O {output.vcf} \
        2> {log}
       
        tabix -p vcf {output.vcf} 2>> {log}
        """
#vqsr for indels
rule vqsr_indel:
#apply vqsr for indels
#vqsr for gatk
rule vqsr_snp:
    input:
        vcf = f"{config['results_dir']}/jointcalling/cohort.joint.vcf.gz"
    output:
        recal = f"{config['results_dir']}/vqsr/cohort.recal",
        tranches = f"{config['results_dir']}/vqsr/cohort.tranches",
        rscript = f"{config['results_dir']}/vqsr/cohort.Rplots.pdf"
    log:
        f"{config['results_dir']}/logs/vqsr.log"
    resources:
        mem_gb = 96,
        cpus = 4
    params:
        hapmap = config['vqsr']['hapmap'],
        omni = config['vqsr']['omni'],
        thousand_genomes = config['vqsr']['thousand_genomes'],
        dbsnp = config['vqsr']['dbsnp'],
        axiom = config['vqsr']['axiom']
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0"
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        VariantRecalibrator \
        -R {config['ref']} \
        -V {input.vcf} \
        -O {output.recal} \
        --tranches-file {output.tranches} \
        --rscript-file {output.rscript} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.thousand_genomes} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
        -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an MQ \
        -mode SNP \
        2> {log}
        """
#apply snp vqsr
rule apply_vqsr_snp:
    input:
        vcf = rules.joint_genotyping.output.vcf,
        recal = rules.vqsr.output.recal,
        tranches = rules.vqsr.output.tranches
    output:
        vcf = f"{config['results_dir']}/vqsr/cohort.joint.vqsr_filtered.vcf.gz",
        tbi = f"{config['results_dir']}/vqsr/cohort.joint.vqsr_filtered.vcf.gz.tbi"
    log:
        f"{config['results_dir']}/logs/apply_vqsr.log"
    resources:
        mem_gb = 96,
        cpus = 4
    params:
        ts_filter_level = 99.5
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0",
        "bbc2/samtools/samtools-1.21"
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        ApplyVQSR \
        -V {input.vcf} \
        --recal-file {input.recal} \
        --tranches-file {input.tranches} \
        -O {output.vcf} \
        -mode SNP \
        --truth-sensitivity-filter-level {params.ts_filter_level} \
        2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """
#vqsr for indels
rule vqsr_indel
    input:
        vcf= rules.apply_vqsr_snp.output.vcf
    output:
        recal = f"{config['results_dir']}/vqsr/cohort.recal",
        tranches = f"{config['results_dir']}/vqsr/cohort.tranches",
        rscript = f"{config['results_dir']}/vqsr/cohort.Rplots.pdf"
    params:
        mills =
        dbsnp =
        axiom =
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        VariantRecalibrator \
        -R {config['ref']} \
        -V {input.vcf} \
        -O {output.recal} \
        --tranches-file {output.tranches} \
        --rscript-file {output.rscript} \
        --max-gaussians 4 \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 {config['known_sites']['mills_indels']} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {config['vqsr']['dbsnp']} \
        -an QD \
        -an DP \
        -an FS \
        -an SOR \
        -an MQRankSum \
        -an ReadPosRankSum \
        -mode INDEL
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
        -L {config['regions']['gatk_intervals']}
        """
#now apply vqsr
rule apply_vqsr_indel:

#vep annotation
rule annotate_gatk_vcf:
    input:
        vcf = rules.apply_vqsr.output.vcf,
        ref = config['ref']
    output:
        vcf = f"{config['results_dir']}/annotation/gatk.cohort.vqsr.annotated.vcf.gz",
        html = f"{config['results_dir']}/annotation/gatk.cohort.vqsr.annotated.html"
    log:
        f"{config['results_dir']}/annotated_gatk.log"
    resources:
        mem_gb = 64,
        cpus = 8
    envmodules:
        "bbc2/vep/ensembl-vep-singularity-115.0",
        "bbc2/samtools/samtools-1.21"
    shell:
        """
        singularity exec --bind /varidata/research/ $VEP_SIF vep \
        -i {input.vcf} \
        -o {output.vcf} \
        --vcf \
        --assembly GRCh38 \
        --dir_cache $VEP_CACHE \
        --offline \
        --fork {resources.cpus} \
        --everything \
        --plugin REVEL,file=$VEP_REVEL_GRCH38 \
        --stats {output.html} \
        2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """

#summary_stats
rule variant_stats:
    input:
        gatk_vcf = rules.apply_vqsr.output.vcf
    output:
        summary = f"{config['results_dir']}/qc/variant_calling_summary.txt",
        gatk_stats = f"{config['results_dir']}/qc/gatk_stats.txt"
    log:
        f"{config['results_dir']}/logs/variant_stats.log"
    envmodules:
        "bbc2/bcftools/bcftools-1.17"
    shell:
        """
        #gatk stats
        bcftools stats {input.gatk_vcf} > {output.gatk_stats} 2> {log}
        """
