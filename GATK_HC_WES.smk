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
        f"{config['results_dir']}/qc/gatk_stats.txt",
        f"{config['results_dir']}/qc/multiqc_variantcalling.html"

#haplotypecaller
rule haplotypecaller:
    input:
        bam = f"{config['bam_dir']}/{{sample}}.bqsr.bam",
        bai = f"{config['bam_dir']}/{{sample}}.bqsr.bai",
        ref = config['ref']
    output:
        vcf = f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.vcf.gz",
        tbi = f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.vcf.gz.tbi",
        haplo_bam = f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.bam",
        haplo_bai = f"{config['results_dir']}/haplotypecaller/{{sample}}/{{sample}}.hc.bai"
    benchmark:
        f"{config['results_dir']}/benchmarks/haplotypecaller/{{sample}}.tsv"
    log:
        f"{config['results_dir']}/logs/{{sample}}.haplotypecaller.log"
    resources:
        mem_mb = 32768,
        runtime = 17280
    threads: 6
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0"
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx28g -Djava.io.tmpdir=./tmp" \
        HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        -L {config[regions][gatk_intervals]} \
        -O {output.vcf} \
        --bam-output {output.haplo_bam} \
        --create-output-bam-index true \
        --ERC GVCF \
        --dont-use-soft-clipped-bases true \
        --standard-min-confidence-threshold-for-calling 20 \
        --native-pair-hmm-threads {threads} \
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
    benchmark:
        f"{config['results_dir']}/benchmarks/genomicsdb/genomicsdb.tsv"
    log:
        f"{config['results_dir']}/logs/genomicsdb.log"
    resources:
        mem_mb = 8192,
        runtime = 4320
    params:
        sample_map = f"{config['results_dir']}/sample_map.txt",
        sample_map_lines = lambda wildcards: "\n".join(
            f"{s}\t{config['results_dir']}/haplotypecaller/{s}/{s}.hc.vcf.gz"
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
        gatk --java-options "-Xms2g -Xmx6g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        --sample-name-map {params.sample_map} \
        --genomicsdb-workspace-path {output.db} \
        -L {config[regions][gatk_intervals]} \
        --reader-threads 2 \
        2> {log}
        """
rule joint_genotyping:
    input:
        db = rules.genomicsdb.output.db,
        ref = config['ref'],
        dbsnp = config['vqsr']['dbsnp']
    output:
        vcf = f"{config['results_dir']}/jointcalling/cohort.joint.vcf.gz",
        tbi = f"{config['results_dir']}/jointcalling/cohort.joint.vcf.gz.tbi"
    benchmark:
        f"{config['results_dir']}/benchmarks/joint_genotyping/joint_genotyping.tsv"
    log:
        f"{config['results_dir']}/logs/joint_genotyping.log"
    resources:
        mem_mb = 98304,
        runtime = 2880
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0"
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
         mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx92g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {input.ref} \
        -V gendb://{input.db} \
        --dbsnp {input.dbsnp} \
        -O {output.vcf} \
        2> {log}
       
        tabix -p vcf {output.vcf} 2>> {log}
        """
#vqsr for gatk
rule vqsr_snp:
    input:
        vcf = f"{config['results_dir']}/jointcalling/cohort.joint.vcf.gz"
    output:
        recal = f"{config['results_dir']}/vqsr/snp/cohort.recal",
        tranches = f"{config['results_dir']}/vqsr/snp/cohort.tranches"
    benchmark:
        f"{config['results_dir']}/benchmarks/vqsr/vqsr_snp.tsv"
    log:
        f"{config['results_dir']}/logs/vqsr.log"
    resources:
        mem_mb = 98304,
        runtime = 1440
    params:
        hapmap = config['vqsr']['hapmap'],
        omni = config['vqsr']['omni'],
        thousand_genomes = config['vqsr']['thousand_genomes'],
        dbsnp = config['vqsr']['dbsnp']
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0"
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx92g -Djava.io.tmpdir=./tmp" \
        VariantRecalibrator \
        -R {config[ref]} \
        -V {input.vcf} \
        -O {output.recal} \
        --tranches-file {output.tranches} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.thousand_genomes} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
        -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an MQ \
        -mode SNP \
        -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        2> {log}
        """
#apply snp vqsr
rule apply_vqsr_snp:
    input:
        vcf = rules.joint_genotyping.output.vcf,
        recal = rules.vqsr_snp.output.recal,
        tranches = rules.vqsr_snp.output.tranches
    output:
        vcf = f"{config['results_dir']}/vqsr/snp/cohort.joint.snp.vqsr.vcf.gz",
        tbi = f"{config['results_dir']}/vqsr/snp/cohort.joint.snp.vqsr.vcf.gz.tbi"
    benchmark:
        f"{config['results_dir']}/benchmarks/vqsr/apply_vqsr_snp.tsv"
    log:
        f"{config['results_dir']}/logs/apply_vqsr_snp.log"
    resources:
        mem_mb = 98304,
        runtime = 720
    params:
        ts_filter_level = 99.5
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0"
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx92g -Djava.io.tmpdir=./tmp" \
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
rule vqsr_indel:
    input:
        vcf= rules.apply_vqsr_snp.output.vcf
    output:
        recal = f"{config['results_dir']}/vqsr/indel/cohort.recal",
        tranches = f"{config['results_dir']}/vqsr/indel/cohort.tranches"
    benchmark:
        f"{config['results_dir']}/benchmarks/vqsr/vqsr_indel.tsv"
    log:
        f"{config['results_dir']}/logs/vqsr_indel.log"
    resources:
        mem_mb = 65536,
        runtime = 5760 
    params:
        mills = config['vqsr']['mills_indels'],
        dbsnp = config['vqsr']['dbsnp'],
        axiom = config['vqsr']['axiom_indels']
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0"
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx64g -Djava.io.tmpdir=./tmp" \
        VariantRecalibrator \
        -R {config[ref]} \
        -V {input.vcf} \
        -O {output.recal} \
        --tranches-file {output.tranches} \
        --max-gaussians 4 \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
        --resource:axiom,known=false,training=true,truth=false,prior=10.0 {params.axiom} \
        -an QD \
        -an DP \
        -an FS \
        -an SOR \
        -an MQRankSum \
        -an ReadPosRankSum \
        -mode INDEL \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        -L {config[regions][gatk_intervals]} \
        2> {log}
        """
#now apply vqsr
rule apply_vqsr_indel:
    input:
        vcf= rules.apply_vqsr_snp.output.vcf,
        recal= rules.vqsr_indel.output.recal,
        tranches= rules.vqsr_indel.output.tranches
    output:
        vcf= f"{config['results_dir']}/vqsr/indel/cohort.vqsr.final.vcf.gz",
        tbi= f"{config['results_dir']}/vqsr/indel/cohort.vqsr.final.vcf.gz.tbi"
    benchmark:
        f"{config['results_dir']}/benchmarks/vqsr/apply_vqsr_indel.tsv"
    log:
        f"{config['results_dir']}/logs/apply_vqsr_indel.log"
    resources:
        mem_mb = 32768,
        runtime = 720
    params:
        ts_filter_level = 99.0
    envmodules:
        "bbc2/gatk/gatk-4.5.0.0"
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        mkdir -p ./tmp
        gatk --java-options "-Xms8g -Xmx28g -Djava.io.tmpdir=./tmp" \
        ApplyVQSR \
        -V {input.vcf} \
        --recal-file {input.recal} \
        --tranches-file {input.tranches} \
        -O {output.vcf} \
        -mode INDEL \
        --truth-sensitivity-filter-level {params.ts_filter_level} \
        2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """
###this rule will catch sites with a missing filter field i.e. sites that VQSR could not model, so I will use hard filters for them
rule gatk_hard_filter:
    input:
        vcf = rules.apply_vqsr_indel.output.vcf
    output:
        vcf = f"{config['results_dir']}/vqsr/cohort.vqsr.hardfilter.vcf.gz",
        tbi = f"{config['results_dir']}/vqsr/cohort.vqsr.hardfilter.vcf.gz.tbi"
    benchmark:
        f"{config['results_dir']}/benchmarks/hard_filter/hard_filter_fallback.tsv"
    log:
        f"{config['results_dir']}/logs/hard_filter_fallback.log"
    resources:
        mem_mb  = 16384,
        runtime = 240
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        bcftools filter \
            -e 'FILTER="." && TYPE="snp" && (QD < 2.0 || FS > 60.0 || MQ < 40.0)' \
            -s "HardFilter_SNP" \
            -O u {input.vcf} | \
        bcftools filter \
            -e 'FILTER="." && TYPE="indel" && (QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0)' \
            -s "HardFilter_INDEL" \
            -O z -o {output.vcf} \
            2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """
#norm vcfs
rule normalize_vcf:
    input:
        vcf = rules.gatk_hard_filter.output.vcf,
        ref = config['ref']
    output:
        vcf= f"{config['results_dir']}/filtered/cohort.vqsr.norm.vcf.gz",
        tbi= f"{config['results_dir']}/filtered/cohort.vqsr.norm.vcf.gz.tbi"
    benchmark:
        f"{config['results_dir']}/benchmarks/norm_vcf/normalize_vcf.tsv"
    log:
        f"{config['results_dir']}/logs/normalized_vcf.log"
    resources:
        mem_mb = 16384,
        runtime = 240
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        bcftools norm \
        -m + \
        -f {input.ref} \
        -Oz \
        -o {output.vcf} \
        {input.vcf} 2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """
rule fill_tags:
    input:
        vcf = rules.normalize_vcf.output.vcf
    output:
        vcf = f"{config['results_dir']}/filtered/cohort.vqsr.norm.tagged.vcf.gz",
        tbi = f"{config['results_dir']}/filtered/cohort.vqsr.norm.tagged.vcf.gz.tbi"
    benchmark:
        f"{config['results_dir']}/benchmarks/fill_tags/fill_tags.tsv"
    log:
        f"{config['results_dir']}/logs/fill_tags.log"
    resources:
        mem_mb = 16384,
        runtime = 240
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        bcftools +fill-tags \
        {input.vcf} \
        -Oz \
        -o {output.vcf} \
        -- -t VAF,AF,AC,AN,F_MISSING,MAF 2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """
#want to filter based on depth >=10 and VAF >=0.3, if they do not meet these thresholds they will be assigned as missing
rule genotype_filter:
    input:
        vcf = rules.fill_tags.output.vcf
    output:
        vcf = f"{config['results_dir']}/filtered/cohort.vqsr.norm.tagged.gtfiltered.vcf.gz",
        tbi = f"{config['results_dir']}/filtered/cohort.vqsr.norm.tagged.gtfiltered.vcf.gz.tbi"
    benchmark:
        f"{config['results_dir']}/benchmarks/genotype_filter/genotype_filter.tsv"
    log:
        f"{config['results_dir']}/logs/genotype_filter.log"
    resources:
        mem_mb  = 16384,
        runtime = 240
    params:
        min_dp  = 10,
        min_vaf = 0.3
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        #set genotypes to missing if DP < 10 or VAF < 0.3
        bcftools filter \
            -S . \
            -e 'FORMAT/DP < {params.min_dp} | FORMAT/VAF < {params.min_vaf}' \
            -O u {input.vcf} | \

        #drop sites where all samples are now missing
        bcftools view \
            -e 'AC=0' \
            -O z -o {output.vcf} \
            2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """
#vep annotation
rule annotate_gatk_vcf:
    input:
        vcf = rules.genotype_filter.output.vcf,
        ref = config['ref']
    output:
        vcf = f"{config['results_dir']}/annotation/gatk.cohort.vqsr.annotated.vcf.gz",
        html = f"{config['results_dir']}/annotation/gatk.cohort.vqsr.annotated.html"
    benchmark:
        f"{config['results_dir']}/benchmarks/annotation/annotate_gatk_vcf.tsv"
    log:
        f"{config['results_dir']}/annotated_gatk.log"
    resources:
        mem_mb = 65536,
        runtime = 1440
    threads: 8
    envmodules:
        "bbc2/vep/ensembl-vep-singularity-115.0"
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        singularity exec --bind /varidata/research/ $VEP_SIF vep \
        -i {input.vcf} \
        -o {output.vcf} \
        --vcf \
        --assembly GRCh38 \
        --dir_cache $VEP_CACHE \
        --offline \
        --fork {threads} \
        --everything \
        --plugin REVEL,file=$VEP_REVEL_GRCH38 \
        --compress_output bgzip \
        --force_overwrite \
        --stats_file {output.html} \
        2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """

#summary_stats
rule variant_stats:
    input:
        gatk_vcf = rules.genotype_filter.output.vcf
    output:
        gatk_stats = f"{config['results_dir']}/qc/gatk_stats.txt"
    benchmark:
        f"{config['results_dir']}/benchmarks/variant_stats/variant_stats.tsv"
    log:
        f"{config['results_dir']}/logs/variant_stats.log"
    resources:
        mem_mb = 8192,
        runtime = 540
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        #gatk stats
        bcftools stats {input.gatk_vcf} > {output.gatk_stats} 2> {log}
        """
#get per sample stats-- i.e. do certain samples have poor call rates, missing genotypes, etc Ti/Tv should be ~2.8 for exonic variants
rule per_sample_stats:
    input:
        vcf = rules.genotype_filter.output.vcf
    output:
        stats = expand(
            f"{config['results_dir']}/qc/per_sample/{{sample}}.stats.txt",
            sample=SAMPLES)
    benchmark:
        f"{config['results_dir']}/benchmarks/per_sample_stats/per_sample_stats.tsv"
    log:
        f"{config['results_dir']}/logs/per_sample_stats.log"
    resources:
        mem_mb  = 16384,
        runtime = 480
    params:
        outdir = f"{config['results_dir']}/qc/per_sample",
        sample = SAMPLES
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        for sample in {SAMPLES}; do
            bcftools stats \
                -s $sample \
                {input.vcf} \
                > {params.outdir}/${{sample}}.stats.txt \
                2>> {log}
        done
        """
#aggregate stats
rule multiqc_variantcalling:
    input:
        cohort_stats = rules.variant_stats.output.gatk_stats,
        per_sample   = expand(
            f"{config['results_dir']}/qc/per_sample/{{sample}}.stats.txt",
            sample=SAMPLES)
    output:
        f"{config['results_dir']}/qc/multiqc_variantcalling.html"
    benchmark:
        f"{config['results_dir']}/benchmarks/multiqc/multiqc_variantcalling.tsv"
    log:
        f"{config['results_dir']}/logs/multiqc_variantcalling.log"
    resources:
        mem_mb  = 8192,
        runtime = 120
    params:
        outdir      = f"{config['results_dir']}/qc",
        search_dirs = f"{config['results_dir']}/qc"
    envmodules:
        "bbc2/multiqc/multiqc-1.14"
    shell:
        """
        multiqc {params.search_dirs} \
            -o {params.outdir} \
            -n multiqc_variantcalling.html \
            --force \
            2> {log}
        """
