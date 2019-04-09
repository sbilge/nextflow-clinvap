# nf-core/clinvap: Output

This document describes the output produced by the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [VEP](#vep) - Variant effect prediction via Ensembl VEP
* [Report Generation](#report) - Aggregate clinical variant annottion report

## VEP
[VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html) is used to annotated raw VCF files. 

**Output directory: `results`**

*  `<VCF>.vcf`
  * Annotated VCF file 

## Report Generation
[Report](#report) is the report generation step precessing the annotated VCF file. 

**Output directory: `results`**

* `<REPORT>.json`
  * Clinical report 
