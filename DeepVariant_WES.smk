#this is a germline variant calling pipeline from whole-exome sequencing data using Google DeepVariant

workflow.name = "Germline Variant Calling in Whole-Exome Sequencing"
workflow.description = "Whole exome sequencing variant calling and annotation"
workflow.author = "Qualia Hooker"
workflow.version = "1.0"

configfile: "config.yaml"

SAMPLES = glob_wildcards(f"{config['bam_dir']}/{{sample}}.bqsr.bam").sample

#rule all define all outputs
rule all:
    input:
        f"{config['results_dir']}/annotation/dv.cohort.annotated.vcf.gz",
        f"{config['results_dir']}/qc/deepvariant_stats.txt"

#deepvariant
rule deepvariant_vc:
    input:
        bam = config['bam_dir'],
        ref = config['ref']
    output:
        vcf = f"{config['results_dir']}/deepvariant/{{sample}}/{{sample}}.dv.vcf.gz",
        tbi = f"{config['results_dir']}/deepvariant/{{sample}}/{{sample}}.dv.vcf.gz.tbi",
        gvcf = f"{config['results_dir']}/deepvariant/{{sample}}/{{sample}}.dv.g.vcf.gz"
    log:
        f"{config['results_dir']}/logs/dv/{{sample}}.deepvariant.log"
    params:
        model_type="WES"
    threads: 8
    resources:
        cpus = 8,
        mem_gb = 96,
        time = "12:00:00"
    envmodules:
        "bbc2/DeepVariant/deepvariant-1.5.0"
    shell:
        """
        singularity exec --bind /varidata/research/ $DEEPVAR /opt/deepvariant/bin/run-deepvariant \
        --model_type={params.model_type} \
        --ref={input.ref} \
        --reads={input.bam} \
        --output_vcf={output.vcf} \
        --output_gvcf={output.gvcf} \
        --num_shards={resources.cpus} \
        --make_example_extra_args=split_skip_reads=true \
        2> {log}
        """
#joint calling dv with GLnexus
rule dv_joint_calling:
    input:
        dv_vcfs = expand(
            f"{config['results_dir']}/deepvariant/{{sample}}/{{sample}}.dv.g.vcf.gz",
            sample=SAMPLES
        )
    output:
        vcf = f"{config['results_dir']}/deepvariant/dv.merged.g.vcf.gz",
        tbi = f"{config['results_dir']}/deepvariant/dv.merged.g.vcf.gz.tbi"
    log:
        f"{config['results_dir']}/logs/dv/merged_dv_vcfs.log"
    resources:
        mem_gb = 64,
        cpus = 4
    envmodules:
        "bbc2/GLnexus/GLnexus-1.4.1"
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        mkdir -p $(dirname {output.vcf})

        bcftools merge \
        --threads {resources.cpus} \
        -o {output.vcf} \
        -O z \
        {input.dv_vcfs} \
        2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """
#convert bcf to vcf for vep
rule bcf_to_vcf:
    input:
    output:
        vcf = f"{config['results_dir']}/deepvariant/dv.merged.g.vcf.gz",
        tbi = f"{config['results_dir']}/deepvariant/dv.merged.g.vcf.gz.tbi"
    log:
        f"{config['results_dir']}/logs/dv/bcf_to_vcf.log"
    resources:
        mem_gb = ,
        cpus = 4
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        bcftools view \
        --threads {resources.cpus}
        -O v {input.bcf} > {output.vcf} \
        2>> {log}
        --config DeepVariant WES \
        --list \
        --bed {rule.dv_joint_calling.input.bed}\
        --mem-gbytes \
        --threads \
        2>> {log}
        """ 
#normalize vcf 
rule vcf_norm:

#vep annotation
rule annotate_dv_vcf:
    input:
        vcf = rules.bcf_to_vcf.output.vcf,
        ref = config['ref']
    output:
        vcf = f"{config['results_dir']}/annotation/dv.cohort.annotated.vcf.gz",
        html = f"{config['results_dir']}/annotation/dv.cohort.annotated.html"
    log:
        f"{config['results_dir']}/annotated_deepvariant.log"
    resources:
        mem_gb = 64,
        cpus = 8
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
        dv_vcf = rules.merge_dv.output.vcf,
        dv_annotated = rules.annotate_dv_vcf.output.vcf
    output:
        summary = f"{config['results_dir']}/qc/variant_calling_summary.txt",
        dv_stats = f"{config['results_dir']}/qc/deepvariant_stats.txt"
    log:
        f"{config['results_dir']}/logs/variant_stats.log"
    conda:
        "envs/vcf_utils.yml"
    shell:
        """
        #deepvariant stats
        bcftools stats {input.dv_vcf} > {output.dv_stats} 2>> {log}
        """
