#!/usr/bin/env nextflow
/*
========================================================================================
                         bigbio/pride-molecules-indexer
========================================================================================
 bigbio/pride-molecules-indexer Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/bigbio/pride-molecules-indexer
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run bigbio/pride-molecules-indexer --project_accession "PXD029360" -profile conda

    Main arguments:
      --project_accession           Project accession to convert the identifications to json files.
      --outdir                      Output directory containing the information (json) of the project

      Advanced Options:
      --pride_production_folder     (Optional) Local folder that contains the submitted files for PRIDE

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Stage config files
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

// Validate input
if (isCollectionOrArray(params.project_accession)){
  tocheck = params.project_accession[0]
} else {
  tocheck = params.project_accession
}

params.project_accession = params.project_accession ?: { log.error "No project accession provided. Make sure you have used the '--project_accession' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()

process project_get_result_files{

  publishDir "${params.outdir}/result_files", mode: 'copy', pattern: '*.tsv'

  input:

  output:
   file "*.tsv" into result_file_summary

  script:
  """
  java -jar ${projectDir}/bin/pride-molecules-indexer-1.0.0-SNAPSHOT.jar get-result-files --app.project-accession=${params.project_accession} \
       --app.file-output=${params.project_accession}-result_files.tsv
  """
}

result_file_summary.splitCsv(skip: 1, sep: '\t')
  .multiMap{ row -> id = row[0]
                    result_files: tuple(id, !params.root_folder ? row[3].replace("ftp://", "http://") :params.root_folder + "/" + row[0])
           }
  .set{ch_result_files}

//ch_result_files.subscribe { println "value: $it" }

process uncompress_result_files{

   label 'downloading_thread'

   input:
     tuple result_id, result_file_path from ch_result_files.result_files

   output:
     tuple result_id, file("*") into ch_result_uncompress
     tuple result_id, file("*") into ch_result_uncompress_process

   script:
   """
   wget '${result_file_path}'
   gunzip '${result_id}'
   """
}

process project_get_related_spectra{

publishDir "${params.outdir}/result_files", mode: 'copy', pattern: '*.tsv'

  input:
    tuple result_id, file(uncompress_result) from ch_result_uncompress

  output:
    tuple result_id, file("*.tsv") into ch_spectra_summary

  script:
  """
  java -jar ${projectDir}/bin/pride-molecules-indexer-1.0.0-SNAPSHOT.jar get-related-files --app.project-accession=${params.project_accession} \
       --app.file-output=${params.project_accession}-${result_id}-result_spectra.tsv --app.result-file=${uncompress_result}
  """
}

// ch_spectra_summary.subscribe { println "value: $it" }

ch_spectra_summary.map { tuple_summary ->
                         def key = tuple_summary[0]
                         def summary_file = tuple_summary[1]
                         def list_spectra = tuple_summary[1].splitCsv(skip: 1, sep: '\t')
                         .collect{"'" + it[5] + "'"}
                         .flatten()
                         return tuple(key.toString(), list_spectra.findAll{ it != "'null'"})
                        }
                   .groupTuple()
                   .into { ch_spectra_tuple_results; ch_spectra_pride_xml}

// ch_spectra_tuple_results.subscribe { println "value: $it" }

process download_spectra_files{

  input:
  tuple val(result_id), spectra from ch_spectra_tuple_results

  output:
  tuple val(result_id), path("*") into ch_spectra_files

  when:
  !spectra.flatten().isEmpty()

  script:
  """
  wget ${spectra.flatten().join(" ")}
  """
}

ch_spectra_files_process = (!ch_spectra_files)? ch_spectra_files.map { id, files -> files instanceof List ? [ id, files.collect{"'" + it + "'"} ] : [ id, [ "'" + files + "'" ] ] }: ch_spectra_pride_xml

ch_result_uncompress_process.combine(ch_spectra_files_process, by:0).set{ch_final_map}

//ch_final_map.subscribe { println "value: $it" }

process generate_json_index_files{

  label 'process_high'

  publishDir "${params.outdir}/result_files", mode: 'copy', pattern: '**.json'

  input:
    val(result_id) from ch_final_map

  output:
    file("**.json") into final_index_json

  script:
  """
  java -jar ${projectDir}/bin/pride-molecules-indexer-1.0.0-SNAPSHOT.jar generate-index-files --app.result-file=${result_id[1]} --app.folder-output=`pwd` --app.spectra-files=${result_id[2].join(",")} --app.project-accession=${params.project_accession}
  """
}

final_index_json.subscribe { println "value: $it" println " "}

//--------------------------------------------------------------- //
//---------------------- Nextflow specifics --------------------- //
//--------------------------------------------------------------- //


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
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
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-proteomicslfq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'bigbio/pride-molecules-indexer Workflow Summary'
    section_href: 'https://github.com/bigbio/pride-molecules-indexer'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    echo $workflow.manifest.version &> v_msstats_plfq.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[bigbio/pride-molecules-indexer] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[bigbio/pride-molecules-indexer] FAILED: $workflow.runName"
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

    // On success try attach the multiqc report
    def mqc_report = ""
    try {
        if (workflow.success && ch_ptxqc_report.println()) {
            mqc_report = ch_ptxqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[bigbio/pride-molecules-indexer] Found multiple reports from process 'ptxqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
        else {
          mqc_report = ""
        }
    } catch (all) {
        log.warn "[bigbio/pride-molecules-indexer] Could not attach report to summary email"
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
            log.info "[bigbio/pride-molecules-indexer] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc != "" && mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[bigbio/pride-molecules-indexer] Sent summary e-mail to $email_address (mail)"
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
        log.info "-${c_purple}[bigbio/pride-molecules-indexer]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[bigbio/pride-molecules-indexer]${c_red} Pipeline completed with errors${c_reset}-"
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
    ${c_purple}  bigbio/pride-molecules-indexer v${workflow.manifest.version}${c_reset}
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


//--------------------------------------------------------------- //
//---------------------- Utility functions  --------------------- //
//--------------------------------------------------------------- //

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Check class of an Object for "List" type
boolean isCollectionOrArray(object) {
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}