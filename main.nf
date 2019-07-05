#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/clinvap
========================================================================================
 nf-core/clinvap Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/clinvap
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/clinvap v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/clinvap --vcf '/input/folder' -profile docker

    Mandatory arguments:
      --vcf  [Path]                 Path to the input data
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Other options:
      --outdir                      The output directory where the results will be saved
      --vep_cache [Path]            Path to Ensemble VEP cache files
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on

    Skipping steps:
      --skip_vep_cache              Skip downloading Ensemble VEP cache files
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// // TODO nf-core: Add any reference files that are needed
// // Configurable reference genomes
// fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
// if ( params.fasta ){
//     fasta = file(params.fasta)
//     if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
// }

params.vcf = params.vcf ?: { log.error "No input data folder is provided. Make sure you have used the '--vcf' option.": exit 1 }()
params.outdir = params.outdir ?: {log.warn "No ouput directory is provided. Results will be saved into './results'"; return "$baseDir/results"}()
// params.vep_cachedir = params.vep_cachedir ?: {log.warn "No VEP cache directory is provided. Cache files will be downloaded into './vep_cache'"; return "./vep_cache"}()

/*
 * Define the default parameters
 */


// Configurable variables

params.name = false
params.email = false
params.plaintext_email = false

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

/* 
 * Create a channel for input read files
*/
Channel
  .fromPath(params.vcf)
  .ifEmpty { exit 1, "params.vcf was empty - no input files supplied"}
  .set {input_vcf}

/* 
 * Create a channel for ensembl vep cache files when provided by user
*/
params.vep_cache = false
if (params.vep_cache){
  params.skip_vep_cache = true
  Channel
    .frompath("${params.vep_cache}")
    .set(vep_offline_files)
} else {
  params.skip_vep_cache = false
}

// /*
//  * Create a channel for input read files
//  */
//  if(params.readPaths){
//      if(params.singleEnd){
//          Channel
//              .from(params.readPaths)
//              .map { row -> [ row[0], [file(row[1][0])]] }
//              .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
//              .into { read_files_fastqc; read_files_trimming }
//      } else {
//          Channel
//              .from(params.readPaths)
//              .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
//              .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
//              .into { read_files_fastqc; read_files_trimming }
//      }
//  } else {
//      Channel
//          .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
//          .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
//          .into { read_files_fastqc; read_files_trimming }
//  }


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/clinvap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/clinvap'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Variants'] = params.vcf
summary['VEP Cache dir'] = params.vep_cachedir
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-clinvap-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/clinvap Workflow Summary'
    section_href: 'https://github.com/nf-core/clinvap'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


// /*
//  * Parse software version numbers
//  */
// process get_software_versions {

//     output:
//     file 'software_versions_mqc.yaml' into software_versions_yaml

//     script:
//     // TODO nf-core: Get all tools to print their version number here
//     """
//     echo $workflow.manifest.version > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     fastqc --version > v_fastqc.txt
//     multiqc --version > v_multiqc.txt
//     scrape_software_versions.py > software_versions_mqc.yaml
//     """
// }


/*
 * STEP 1 - Vep Cache Files
 */

process ensembl_vep_files {

  storeDir "${params.outdir}/offline"
  output:
  file("*") into vep_offline_files 

  when:
  !params.skip_vep_cache

  script:
  """
  git clone -b release/95 https://github.com/Ensembl/ensembl-vep.git
  cd ensembl-vep
  perl INSTALL.pl --NO_HTSLIB -n --CACHE_VERSION 95 --CACHEDIR './offline_cache' --VERSION 95 -a acf -s homo_sapiens -y GRCh37
  wget 'https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/90/LoFtool_scores.txt'
  """
}

/*
 * STEP 2 - Ensembl VEP Annotation
 */

process vep_on_input_file {

  publishDir "${params.outdir}"
  input:
  file vcf_file from input_vcf
  file('ensembl-vep') from vep_offline_files
  
  output:
  file "${vcf_file.baseName}_out.vcf" into annotated_vcf

  script:
  if (!params.skip_vep_cache)
  """
  vep -i ${vcf_file} -o ${vcf_file.baseName}_out.vcf --dir_cache "${params.outdir}/offline/ensembl-vep/offline_cache" --config $baseDir/assets/vep.ini
  """
  else
  """
  vep -i ${vcf_file} -o ${vcf_file.baseName}_out.vcf --dir_cache ${params.vep_cache} --config $baseDir/assets/vep.ini
  """
 }

/*
 * STEP 3 - Report Generation
 */

process report_generation {

  publishDir "${params.outdir}/reports"

  input:
  file out_vcf from annotated_vcf

  output:
  file "${out_vcf.baseName}.json"

  script:
  """
  Rscript --no-save --no-restore --no-init-file --no-site-file $baseDir/bin/reporting.R -f ${out_vcf} -r "${out_vcf.baseName}.json" -d $baseDir/assets/driver_db_dump.json
  """
}


/*
 * STEP 4 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/clinvap] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/clinvap] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/clinvap] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/clinvap] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/clinvap] Pipeline Complete"
}
