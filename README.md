# nf-core/clinvap
**Clinical Variant Annotation Pipeline**

[![Build Status](https://travis-ci.org/nf-core/clinvap.svg?branch=master)](https://travis-ci.org/nf-core/clinvap)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/clinvap.svg)](https://hub.docker.com/r/nfcore/clinvap)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->
### Introduction

**ClinVAP** retrieves information from simple somatic mutations (SNVs) given in VCF files and creates structured clinical reports by annotating, prioritizing, and filtering the genomic variants. The report is designed to assist Molecular Tumor Boards (MTB) in making therapeutic decisions by providing them with information on the molecular mechanisms initiating carcinogenesis and on actionable genes.


### Quick Start
To test the pipeline you may run:   
`nextflow run nf-core/clinvap -profile test,docker`

To run the pipeline with your data:

`nextflow run nf-core/clinvap --vcf '/input/folder' -profile docker`

If you already have Ensembl VEP cache files: 

`nextflow run nf-core/clinvap --vep_cache 'PATH' --vcf '/input/folder' -profile docker`

### Documentation
The nf-core/clinvap pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
    * [Reference genomes](docs/configuration/reference_genomes.md)  
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)



### Credits
nf-core/clinvap was originally written by Bilge Sürün, Julian Heinrich, Charlotta Schärfe, Mathew Divine.

### Disclaimer 
The report created by ClinVAP is intended as a hypothesis generating framework and thus for research use only. It is not intended for diagnostic or clinical purposes. Information provided in the report does not replace a physician’s medical judgment and usage is entirely at your own risk. The providers of this resource shall in no event be liable for any direct, indirect, incidental, consequential, or exemplary damages.
