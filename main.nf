#!/usr/bin/env nextflow
/*
========================================================================================
                         KohlbacherLab/nextflow-clinvap
========================================================================================
 KohlbacherLab/nextflow-clinvap Clinical variant annotation pipeline
 #### Homepage / Documentation
 https://github.com/KohlbacherLab/nextflow-clinvap
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --skip_vep false --vcf '/input/folder/*.vcf' -profile docker

    Mandatory arguments:
      --vcf  [Path]                 Path to the input directory.
      --annotated_vcf [Path]        Path to the VEP-annotated vcf.
      --genome [str]                Genome version. Supported are: 'GRCh37' or 'GRCh38'.
      -profile [str]                Configuration profile to use. Can use multiple (comma separated).
                                    Available: conda, docker, singularity, awsbatch, test and more.
    Other options:
      --vep_cache [Path]            Path to Ensemble VEP cache files
      --skip_vep_cache [bool]       Skip downloading Ensemble VEP cache files
      --skip_vep [bool]             Skip variant effect prediction step
      --metadata_json [Path]        Json file with patient information such as diagnosis
      --cnv [Path]                  Path to cnv file
      --docx_template [Path]        DOCX template to render JSON report
      --diagnosis_filter_option [str]       Diagnosis based filter type     
      --outdir [file]                 The output directory where the results will be saved
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// vcf input check  
if (!params.skip_vep) {
    params.vcf = params.vcf ?: { log.error "No input data folder is provided. Make sure you have used the '--vcf' option.": exit 1 }()
}

params.outdir = params.outdir ?: {log.warn "No ouput directory is provided. Results will be saved into './results'"; return "$baseDir/results"}()

//  workflow run name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files 
ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input read files
 */
if (!params.skip_vep) {
    Channel
    .fromPath(params.vcf)
    .ifEmpty { exit 1, "params.vcf was empty - no input files supplied"}
    .set {input_vcf}
} else {
    input_vcf = Channel.empty()
}


/*
* Create a channel for annotated vcf files
*/

if (params.skip_vep){
    Channel
        .fromPath(params.annotated_vcf)
        .ifEmpty {exit 1, "Cannot find any vcf matching: ${params.annotated_vcf}.\nTry enclosing paths in quotes!\nTry adding a * wildcard!"}
        .set {annotated_input_vcf}
} else {
    annotated_input_vcf = Channel.empty()
}



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

/* 
 * Create a channel for ensembl vep cache files when provided by user
*/
if (params.metadata_json){
    Channel
    .fromPath(params.metadata_json)
    .ifEmpty {exit 1, "Cannot find any json matching: ${params.metadata_json}.\nTry enclosing paths in quotes!"}
    .set {ch_metadata}
} else {
    ch_metadata = Channel.empty()
}


/*
* Create a channel for cnv file
*/
if (params.cnv){
    Channel
    .fromPath(params.cnv)
    .ifEmpty {exit 1, "Cannot find any tsv matching: ${params.cnv}.\nTry enclosing paths in quotes!"}
    .set {ch_cnv}
} else {
    ch_cnv = Channel.value("EMPTY")
}


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Variants'] = params.vcf
summary['VEP Cache dir'] = params.vep_cachedir
summary['DOCX Template'] = params.docx_template
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nextflow-clinvap-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'KohlbacherLab/nextflow-clinvap Workflow Summary'
    section_href: 'https://github.com/KohlbacherLab/nextflow-clinvap'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}



/*
 * STEP 1 - Vep Cache Files
 */

process ensembl_vep_files {

  storeDir "${params.outdir}/offline"

  output:
  file("*") into vep_offline_files 

  when:
  !params.skip_vep_cache && !params.skip_vep

  script:
  """
  git clone -b release/95 https://github.com/Ensembl/ensembl-vep.git
  cd ensembl-vep
  perl INSTALL.pl --NO_HTSLIB -n --CACHE_VERSION 95 --CACHEDIR './offline_cache' --VERSION 95 -a acf -s homo_sapiens -y ${params.genome}
  wget 'https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/90/LoFtool_scores.txt'
  """
}

/*
 * STEP 2 - VCF Filter 
 */

process filter_vcf {

    publishDir "${params.outdir}/reports/vcf", mode: 'copy'

    input:
    file vcf_file from input_vcf

    output:
    file "${vcf_file.baseName}.filter.vcf" into vep

    when:
    !params.skip_vep

    script:
    """
    filter_vcf.py ${vcf_file} ${vcf_file.baseName}.filter.vcf
    """
}

/*
 * STEP 3 - Ensembl VEP Annotation
 */

process vep_on_input_file {

  publishDir "${params.outdir}/reports/vcf", mode: 'copy'

  input:
  file vcf from vep
  file('ensembl-vep') from vep_offline_files
  
  output:
  file "${vcf.simpleName}.out.vcf" into ch_annotated_vcf, rep_ch_annotated_vcf

  when:
  !params.skip_vep

  script:
  if (!params.skip_vep_cache)
  """
  vep -i ${vcf} -o ${vcf.simpleName}.out.vcf --dir_cache "${params.outdir}/offline/ensembl-vep/offline_cache" --assembly ${params.genome}  --config $baseDir/assets/vep.ini
  """
  else
  """
  vep -i ${vcf} -o ${vcf.simpleName}.out.vcf --dir_cache ${params.vep_cache} --assembly ${params.genome} --config $baseDir/assets/vep.ini
  """
}

/*
 * STEP 4 - Report Generation
 */
process snv_report_generation {

  publishDir "${params.outdir}/reports/json", mode: 'copy'

  input:
  file vcf from ch_annotated_vcf.mix(annotated_input_vcf)

  output:
  file "${vcf.simpleName}.vcf.out.json" into snv_metadata, snv_report_generate
  
  when:
  !params.cnv

  script:
  """
  snv_reporting.py -i ${vcf} -o ${vcf.simpleName}.vcf.out.json -g ${params.genome} -k $baseDir/assets/cancerDB_final.json
  """
}

/*
 * STEP 5 - Report Generation
 */


process cnv_report_generation {

  publishDir "${params.outdir}/reports/json", mode: 'copy'

  input:
  file vcf from rep_ch_annotated_vcf
  file cnv from ch_cnv

  output:
  file "${cnv.baseName}.cnv.out.json" optional true into cnv_metadata, cnv_metadata_dummy, cnv_report_generate
  
  when:
  params.cnv

  script:
  if (cnv != "EMPTY")
  """
  cnv_reporting.py -i ${vcf} -c ${cnv} -o ${cnv.baseName}.cnv.out.json -g ${params.genome} -k $baseDir/assets/cancerDB_final.json
  """
}

/*
 * STEP 6 - MERGE METADATA - FILTER DIAGNOSIS
 */

process metadata_diagnosis {
    publishDir "${params.outdir}/reports/json", mode: 'copy'

    input:
    file metadata from ch_metadata
    file main_json from snv_metadata
    val dummy from cnv_metadata_dummy.ifEmpty("EMPTY")
    file cnv_json from cnv_metadata.ifEmpty("EMPTY")

    output:
    file "${main_json.simpleName}.json" into ch_snv_diagnosis
    file "${cnv_json.simpleName}.cnv.json" optional true into ch_cnv_diagnosis


    when:
    params.metadata_json

    script:
    if (dummy =="EMPTY")
    """
    process_metadata.py ${main_json} ${metadata} ${main_json.simpleName}.json $baseDir/assets/database_diagnosis_lookup_table.txt $baseDir/assets/icd10_lookup_table.txt ${params.diagnosis_filter_option}
    """
    else
    """
    process_metadata.py ${main_json} ${metadata} ${main_json.simpleName}.json $baseDir/assets/database_diagnosis_lookup_table.txt $baseDir/assets/icd10_lookup_table.txt ${params.diagnosis_filter_option}
    process_metadata.py ${cnv_json} ${metadata} ${cnv_json.simpleName}.cnv.json $baseDir/assets/database_diagnosis_lookup_table.txt $baseDir/assets/icd10_lookup_table.txt ${params.diagnosis_filter_option}
    """
}

/*
 * STEP 7 - DOCX
 */
process render_report_snv {

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    file out_json from snv_report_generate
    file diagnosis_json from ch_snv_diagnosis.ifEmpty("EMPTY").collect()

    output:
    file "${out_json.baseName}.docx"

    script:
    if (!params.metadata_json)
    """
    docx_generate.py ${out_json} ${params.docx_template} ${out_json.baseName}.docx
    """
    else
    """
    docx_generate.py ${diagnosis_json} ${params.docx_template} ${out_json.baseName}.docx
    """
}

/*
 * STEP 8 - DOCX
 */

process render_report_cnv {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    file cnv_json from cnv_report_generate
    file cnv_diagnosis_json from ch_cnv_diagnosis.ifEmpty("EMPTY")

    output:
    file "${cnv_json.baseName}.docx" optional true

    when:
    params.cnv

    script:
    if (!params.metadata_json)
    """
    docx_generate.py ${cnv_json} ${params.docx_template} ${cnv_json.baseName}.docx
    """
    else
    """
    docx_generate.py ${cnv_diagnosis_json} ${params.docx_template} ${cnv_json.baseName}.docx
    """
}



/*
 * STEP 9 - Output Description HTML
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
    def subject = "[KohlbacherLab/nextflow-clinvap] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[KohlbacherLab/nextflow-clinvap] FAILED: $workflow.runName"
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
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[KohlbacherLab/nextflow-clinvap] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[KohlbacherLab/nextflow-clinvap] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

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
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[KohlbacherLab/nextflow-clinvap] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[KohlbacherLab/nextflow-clinvap] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[KohlbacherLab/nextflow-clinvap]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[KohlbacherLab/nextflow-clinvap]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  KohlbacherLab/nextflow-clinvap v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}