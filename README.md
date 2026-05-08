# Workflows for Germline and Somatic Variant Calling
Automated variant calling pipelines for whole-exome sequencing (WES) data using GATK, DeepVariant, and other callers. Designed for large cohorts with best-practice quality filtering and annotation.

## Overview
Snakemake workflows for reproducible genomic analyses of large cohort data designed to be executed on HPC. It makes use of conda environments for versioning of specific tools.

## Pipelines

### GATK HaplotypeCaller (Germline)
Traditional variant calling. Follows GATK Best Practices.
**File:** `GATK_HC_WES.smk`
### DeepVariant (Germline) *in progress*

Deep learning-based variant calling with GLnexus joint genotyping. Used for high-confidence consensus variants when paired with other callers.
**File:** `DeepVariant_WES.smk`
